/*************************************************************************
 *  Timestamp free synchronization
 *  Slave node code with coarse delay estimation
 *  DRB Apr 14, 2015
 *  Plays out a modulated sinc pulse on the left channel every 4 clock periods

 *  States:
 *  State -1: warmup (needed since codec inputs seem to glitch on startup)
 *  State 0: searching **
 *  State 1: recording **
 *  State 2: coarse/fine delay estimation
 *  State 3: set up adjusted virtual clock buffer?
 *
 *  TODO: Need a timeout on state 0. We should be entering state 0 only
 *  after a slave->master pulse has been played out. If searching doesn't
 *  enter state 1 within a certain amount of time, we should go back to
 *  state -1.
 *************************************************************************/


#define CHIP_6713 1

// length of searching window in samples
#define M 40

// threshold value for searching window
#define T1 100000

// sinc pulse normalized bandwidth
#define BW 0.0125
#define CBW 0.25 	//2kHz@8k Fs  carrier frequency

// 2*N+1 is the number of samples in the sinc function
#define N 512

// virtual clock counter
#define L 4096

// modulated sinc pulse buffer length (mostly zeros)
#define LL 16384

// number of saves for testing/debugging
#define SAVES 30

// Number of fine delay buffers
#define FER	1024

// sampling period
#define TS 0.000125

// define PI and INVPI
#define PI 3.14159265358979323846
#define INVPI 0.318309886183791

#define STATE_WARMUP -1
#define STATE_SEARCH 0
#define STATE_RECORD 1
#define STATE_ESTIMATE 2
#define STATE_BUFFER 3

// set this to one to pass through the left channel
// (useful for debugging master node functionality)
// or zero to not pass through (should normally be zero)
#define PASSTHROUGH 0

#include <stdio.h>
#include <c6x.h>
#include <csl.h>
#include <csl_mcbsp.h>
#include <csl_irq.h>
#include <math.h>

#include "dsk6713.h"
#include "dsk6713_aic23.h"
#include "dsk6713_led.h"
#include "shared_functions.h"

// ------------------------------------------
// start of variables
// ------------------------------------------
 //The following must be volatile 
 	// max_samp state zc zs buf zz max_recbuf recbuf_start_clock recbufindex bufindex recbuf delay_est_done 
 	//vclock_counter swap_pending buffer_just_swapped current_clockbuf sincpulsebuf_counter
volatile float buf[M];       // search buffer
volatile float mfc[M];		// in-phase correlation buffer
volatile float mfs[M];       // quadrature correlation buffer

float corr_max, corr_max_s, corr_max_c; // correlation variables
float corr_c[2*M];
float corr_s[2*M];
float s[2*M];
short corr_max_lag;
volatile short bufindex = 0;

short i,j,ell,imax,k;
volatile float zc, zs, zz; //zc, zs not needed in main
float z, zi, zq, zmax;
double t,x,y;
double fractionalShift = 0.0;
short wait_count = 0;

short sincpulsebuf[LL];  // sinc pulse buffer (left channel)
short clockbuf[L];  // clock buffer (right channel)

volatile short current_clockbuf = 0;  // selects which clock buffer to play
volatile float recbuf[2*N+2*M]; // recording buffer
//volatile short recbuf[2*N+2*M]; // recording buffer
short clockbuf_shifted[2][L];  // double clock buffer (right channel)

float si[2*N+1];  // in-phase sinc pulse
float sq[2*N+1];  // quadrature sinc pulse

volatile short recbufindex = 0;
volatile short max_recbuf = 0;
short playback_scale = 1;

//volatile short state = -1;  // start in state -1 to let codec settle
volatile int state = -1;
volatile short vclock_counter = 0; // virtual clock counter
volatile int sincpulsebuf_counter = 0;  // sinc pulse buffer counter
volatile short recbuf_start_clock = 0; // virtual clock counter for first sample in recording buffer
char r = 0;

volatile char swap_pending;  // flag for buffer swap pending
volatile char buffer_just_swapped = 0; // flag to tell main code that the buffer was just swapped
volatile short max_samp = 0;

float pcf = 0.0;		// phase correction factor
short cde = 0;			// coarse delay estimate (integer)
float fde = 0.0;		// fine delay estimate

short recbuf_start_clock_save[SAVES];
short imax_save[SAVES];
short cde_save[SAVES];     // save fine delay estimates for testing/debugging
float fde_save[SAVES];     // save fine delay estimates for testing/debugging
float pcf_save[SAVES];     // save fine delay estimates for testing/debugging
float ppm_estimate_save[SAVES];

float clockoffset;
float clockoffset_save[SAVES]; // save clock offset estimates for testing/debugging

double ppm_estimate = 0.0;
short save_index = 0;
short buffer_swap_index = L;  // initialize to L to prevent swapping until new buffer is ready
double one_over_beta;
double adjustedclockoffset;
volatile char delay_est_done = 0;  // flag to tell ISR when delay estimates have been calculated
DSK6713_AIC23_CodecHandle hCodec;							// Codec handle
DSK6713_AIC23_Config config = DSK6713_AIC23_DEFAULTCONFIG;  // Codec configuration with default settings

// lookup table of shifted modulated sinc pulses (xxx temporary (wasteful))
#pragma DATA_SECTION(allMyDelayedWaveforms,".mydata")
far short allMyDelayedWaveforms[FER][L];
// ------------------------------------------
// end of variables
// ------------------------------------------

double sin(double);
double cos(double);
interrupt void serialPortRcvISR(void);

void main()
{
	// initialize the save buffers
	for (i=0;i<SAVES;i++) {
		cde_save[i] = 0;
		pcf_save[i] = 0;
		clockoffset_save[i] = 0.0;
		fde_save[i] = 0.0;
		imax_save[i] = 0;
		recbuf_start_clock_save[i] = 0;
		ppm_estimate_save[i] = 0.0;
	}

	// set up the fractionally shifted buffers
	for(k = 0; k < FER; k++){
		fractionalShift = ((double) k )/((double) FER); // number of samples to shift
		for (i=-N;i<=N;i++){
			x = ((double) i - fractionalShift)*BW;
			if (x==0.0) {
				y = 32767.0;
			}
			else
			{
				t = ((double) i - fractionalShift)*CBW;
				y = 32767.0*cos(2*PI*t)*sin(PI*x)/(PI*x); // double
			}
			j = i;
			if (j<0) {
				j += L; // wrap
			}
			allMyDelayedWaveforms[k][j] = (short) y;
		}
	}

	// set up the cosine and sin matched filters for searching
	// also initialize searching buffer

	calc_filter(mfc, mfs, buf); //Initialize the searching buffer.

	// initialize clock buffers
	for (i=0;i<L;i++) {
		clockbuf[i] = 0;
		clockbuf_shifted[0][i] = 0;
		clockbuf_shifted[1][i] = 0;
	}

	// initialize sinc pulse buffer
	for (i=0;i<LL;i++)
		sincpulsebuf[i] = 0;

	// set up clock buffer and sinc pulse buffer
	// to play modulated sinc centered at zero

	/////////////////////////////////////////////
	for (i=-N;i<=N;i++){
		x = i*BW;
		if (i!=0) {
			t = i*CBW;
			y = 32767.0*cos(2*PI*t)*sin(PI*x)/(PI*x); // double
		}
		else {
			y = 32767.0;
		}
		j = i;
		if (j<0) {
			j += L; // wrap
		}
		clockbuf[j] = (short) y;
		j = i;
	/////////////////////////////////////////////
		if (j<0) {
			j += LL; // wrap
		}
		sincpulsebuf[j] = (short) y;
	}

	// set up inphase and quadrature sinc pulses for coarse and fine delay estimators
	j = 0;
	for (i=-N;i<=N;i++){
		x = i*BW;
		if (i!=0) {
			t = i*CBW;
			si[j] = (float) (cos(2*PI*t)*sin(PI*x)/(PI*x));
			sq[j] = (float) (sin(2*PI*t)*sin(PI*x)/(PI*x));
		}
		else {
			si[j] = 1.0;
			sq[j] = 0.0;
		}
		j++;
	}

	//initialize codec.
	init_codec(config, hCodec);


	while(1)						// main loop
	{

		if ((state==STATE_ESTIMATE)&&(delay_est_done==0)){  // TIME TO COMPUTE DELAY ESTIMATES

			DSK6713_LED_toggle(1);	// toggle LED here for diagnostics

			// compute coarse delay estimate
			zmax = 0.0;				// maximum
			imax = 0;				// index of maximum
			for (i=0;i<(2*M-1);i++){  // lag index
				z = 0;
				for (j=0;j<(2*N+1);j++) {
					z+= si[j]*recbuf[i+j];  // correlation at lag i
				}
				if (abs(z)>zmax) {
					zmax = abs(z);  // store maximum
					imax = i;       // store index of maximum
				}
			}
			cde = recbuf_start_clock + imax + N;  // coarse delay estimate (DRB: +N here because si is already shifted by N)
			// cde  is the number of samples elapsed since we launched the S->M sinc pulse

		    // compute fine delay estimate
			zi = 0.0;  // in phase
			zq = 0.0;  // quadrature
			for (j=0;j<(2*N+1);j++) {
				zi+= si[j]*recbuf[imax+j];  // correlation at lag imax
				zq+= sq[j]*recbuf[imax+j];  // correlation at lag imax
			}
			pcf = atan2(zq,zi)*(2*INVPI); // assumes wc = pi/2
			fde = (float) cde + pcf;

			// compute actual clock offset//////////////////////////////////////////////////////////////////////////////
			// the value computed here is always non-negative and
			// represents the number of samples the master clock is ahead of the slave clock
			// (or the number of samples the slave clock is behind the master clock)
			// S->M sinc pulse was launched at time 0
			// M->S sinc pulse was received at time fde (should be positive)
			// to synchronize, we want slave clock ticks to appear at fde/2 + k*L for k=0,1,....
			clockoffset = fde*0.5;

			while (clockoffset>((float) L))
				clockoffset = clockoffset - (float) L;

			// testing/debugging
			imax_save[save_index] = imax;
			recbuf_start_clock_save[save_index] = recbuf_start_clock;
			cde_save[save_index] = cde;
			pcf_save[save_index] = pcf;
			fde_save[save_index] = fde;
			clockoffset_save[save_index] = clockoffset;

			// ppm estimator (quick and dirty)
			// actual formula is
			// (clockoffset_save[save_index]-clockoffset_save[save_index-1])/(LL)*1E6
			// xxx this will not work if clock offset wraps
			if (save_index>=1) {
				ppm_estimate = (clockoffset_save[save_index]-clockoffset_save[save_index-1])*61.03515625; //
			}

			//*** I have to read from here
			ppm_estimate = 6.198883; // xxx temporary
			ppm_estimate_save[save_index] = ppm_estimate;

			// update save index
			save_index++;
			if (save_index>=SAVES)  // wrap
				save_index = 0;


			// copy appropriate fractionally shifted clock buffer to shifted clock buffer
			one_over_beta = ((double) LL)/((double) LL - ppm_estimate/61.03515625); // approx 1+ppm_estimate/1E6

			//Need to account for the possible 1-tick situation
			adjustedclockoffset = clockoffset + ((double) 2*L)*one_over_beta;  // latency for first pulse is 2*L ********(xxx revisit)   

			//Wrapping is wrong in this case.
			//while (adjustedclockoffset>((double) L))
			//	adjustedclockoffset = adjustedclockoffset - (double) L;

			j = (short) adjustedclockoffset; // casting from float to short just truncates, e.g., 3.99 becomes 3 
			fractionalShift = adjustedclockoffset - (double) j; // fractionalShift >= 0 and < 1 tells us how much of a sample is left
			k = (short) (fractionalShift * (double) FER + 0.5);  // this now rounds and givse results on 0,...,FER

			if (k==FER) {  // we rounded up to the next whole sample
				k = 0;
				j++;
			}
			
			for (i=0;i<L;i++) {
				ell = j+i;
				while (ell>=L) {
					ell = ell-L;
				}
				if (current_clockbuf==0) {
					clockbuf_shifted[1][ell] = allMyDelayedWaveforms[k][i];  // write other buffer
				}
				else {
					clockbuf_shifted[0][ell] = allMyDelayedWaveforms[k][i];  // write other buffer
				}

			}

			// when can we swap buffers?
			buffer_swap_index = j+L/2;  // midpoint of the silent part (I think)
			while (buffer_swap_index>=L)
				buffer_swap_index = buffer_swap_index-L;

			// tell the ISR the calculations are done
			delay_est_done = 1;

			DSK6713_LED_toggle(1);	// toggle LED here for diagnostics

		}
		else if (buffer_just_swapped==1) {
			// this code computes the next buffer and attempts to corect for the frequency offset
			DSK6713_LED_toggle(2);
			buffer_just_swapped = 0;
			// copy appropriate fractionally shifted clock buffer to shifted clock buffer
			adjustedclockoffset = adjustedclockoffset + ((double) L)*one_over_beta;  // adjust latency of next pulse
			while (adjustedclockoffset>((double) L))
				adjustedclockoffset = adjustedclockoffset - (double) L;
			j = (short) adjustedclockoffset; // casting from float to short just truncates, e.g., 3.99 becomes 3
			fractionalShift = adjustedclockoffset - (double) j; // fractionalShift >= 0 and < 1 tells us how much of a sample is left
			k = (short) (fractionalShift * (double) FER);  // this also truncates and should give result on 0,...,FER-1
			for (i=0;i<L;i++) {
				ell = j+i;
				if (ell>=L) {
					ell = ell-L;
				}
				if (current_clockbuf==0) {
					clockbuf_shifted[1][ell] = allMyDelayedWaveforms[k][i];  // write other buffer
				}
				else {
					clockbuf_shifted[0][ell] = allMyDelayedWaveforms[k][i];  // write other buffer
				}
				DSK6713_LED_toggle(2);	
			}
		} // if ((state==2)&&(delay_est_done==0))
	}  // while(1)
}  // void main


interrupt void serialPortRcvISR()
{ 
	union {Uint32 combo; short channel[2];} temp;
/*
	short ii=0;
	short jj=0;
*/
	temp.combo = MCBSP_read(DSK6713_AIC23_DATAHANDLE);
	// Note that right channel is in temp.channel[0]
	// Note that left channel is in temp.channel[1]

	// keep track of largest sample for diagnostics
	if (temp.channel[0]>max_samp)
		max_samp = temp.channel[0];

///////////////////////////////////////////////////////
	if (state==STATE_SEARCH) {  // SEARCHING STATE

//		/*
		/*
		search_for_thresh(buf,
						  recbuf,
						  &bufindex,
						  &(temp.channel[0]),
						  &max_recbuf,
						  &state,
						  &recbufindex,
						  PASSTHROUGH,
						  mfc,
						  mfs);
						  */
	//	search_for_thresh(&(temp.channel[0]),PASSTHROUGH);


//	*/
///*
		buf[bufindex] = (float) temp.channel[0];  // right channel

				int ii, jj;
			    // compute incoherent correlation
			    zc = 0;
			    zs = 0;
				for (ii=0;ii<M;ii++) {
					zc += mfc[ii]*buf[ii];
					zs += mfs[ii]*buf[ii];
				}
				zz = zc*zc+zs*zs;

				if (zz>T1) {  				// threshold exceeded?
					max_recbuf = 0;
					state = 1; 				// enter "recording" state (takes effect in next interrupt)
					DSK6713_LED_toggle(0);	// toggle LED here for diagnostics
					// record time of first sample (DRB fixed 4/19/2015: added +1)
					recbuf_start_clock = sincpulsebuf_counter-M+1; // should not be negative since we started counting when we launched the S->M sinc pulse
					recbufindex = M;		// start recording new samples at position M
					jj = bufindex;			//
					for (ii=0;ii<M;ii++){  	// copy samples from buf to first M elements of recbuf
						jj++;   				// the first time through, this puts us at the oldest sample
						if (jj>=M)
							jj=0;
						recbuf[ii] = buf[jj];
						buf[jj] = 0;  		// clear out searching buffer to avoid false trigger
					}
				}
				else {
					// increment and wrap pointer
				    bufindex++;
				    if (bufindex>=M)
				    	bufindex = 0;
				}
//				*/
	}

	///////////////////////
	else if (state==STATE_RECORD) { // RECORDING STATE

		// put sample in recording buffer
		recbuf[recbufindex] = (float) temp.channel[0];  // right channel
		recbufindex++;
		if (recbufindex>=(2*N+2*M)) {
			state = STATE_ESTIMATE;  		// buffer is full
			delay_est_done = 0; // clear flag
			DSK6713_LED_toggle(0);	// toggle LED here for diagnostics
			recbufindex = 0; 	// shouldn't be necessary
		}
	}
	else if (state==STATE_ESTIMATE) { // CALCULATING DELAY ESTIMATES STATE
		if (delay_est_done==1) {  // are the delay estimates done calculating?
			state = STATE_BUFFER; // next state
		}
	}
	else if (state==STATE_BUFFER) { // WRITE ADJUSTED VIRTUAL CLOCK BUFFER STATE
		delay_est_done = 0; // clear flag
		state = STATE_SEARCH;  // xxx temporary
	}

	if (state==STATE_WARMUP) { // WARMUP STATE (NO OUTPUT)
		if (sincpulsebuf_counter>LL/2) {
			state = 0;
		}
		temp.channel[1] = 0;
		temp.channel[0] = 0;
	}
	else { 
		if (vclock_counter==buffer_swap_index)
			swap_pending = 1;
		if ((swap_pending==1)&&(state!=STATE_ESTIMATE)) // ok to swap buffer
		{
			swap_pending = 0;
			buffer_just_swapped = 1;
			if (current_clockbuf==0)
				current_clockbuf = 1;
			else
				current_clockbuf = 0;
		}
		temp.channel[1] = clockbuf_shifted[current_clockbuf][vclock_counter];  // slave *shifted* clock signal (always played)
//		temp.channel[1] = clockbuf[vclock_counter];  // slave *unshifted* clock signal (for debug)
		temp.channel[0] = sincpulsebuf[sincpulsebuf_counter];  // this initiates the sinc pulse exchange with the master
	}

	MCBSP_write(DSK6713_AIC23_DATAHANDLE, temp.combo); // output L/R channels

	// update virtual clock (cycles from 0 to L-1)
	vclock_counter++;
	if (vclock_counter>=L) {
		vclock_counter = 0; // clock tick occurred, wrap
	}

	// update sinc pulse counter (cycles from 0 to LL-1)
	// this is for sinc pulses from the slave to the master
	// sinc pulses from the slave to the master have their peak at
	// sincpulsebuf_counter = 0 (this makes clock offset calculations simple)
	sincpulsebuf_counter++;
	if (sincpulsebuf_counter>=LL) {
		sincpulsebuf_counter = 0; // wrap
	}

}


