/*
 * shared_functions.c
 *
 *  Created on: May 20, 2015
 *      Author: IEUser
 */
#include <math.h>
#include <stdio.h>
#include <c6x.h>
#include <csl.h>
#include <csl_mcbsp.h>
#include <csl_irq.h>

#include "dsk6713.h"
#include "dsk6713_aic23.h"
#include "dsk6713_led.h"


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

// set this to one to pass through the left channel
// (useful for debugging master node functionality)
// or zero to not pass through (should normally be zero)
#define PASSTHROUGH 0

// define PI and INVPI
#define PI 3.14159265358979323846
#define INVPI 0.318309886183791


#define STATE_WARMUP -1 //Slave
#define STATE_SEARCH 0
#define STATE_RECORD 1
#define STATE_WAIT_CLOCK 2 //Master
#define STATE_ESTIMATE 2 //Slave
#define STATE_WAIT_PLAYBACK 3 //Master
#define STATE_BUFFER 3 //Slave
#define STATE_RESPOND 4 //Master

void init_codec(DSK6713_AIC23_Config config, DSK6713_AIC23_CodecHandle hCodec){

	DSK6713_init();		// Initialize the board support library, must be called first
	DSK6713_LED_init(); // initialize LEDs
	hCodec = DSK6713_AIC23_openCodec(0, &config);	// open codec and get handle

	// Configure buffered serial ports for 32 bit operation
	// This allows transfer of both right and left channels in one read/write
	MCBSP_FSETS(SPCR1, RINTM, FRM);
	MCBSP_FSETS(SPCR1, XINTM, FRM);
	MCBSP_FSETS(RCR1, RWDLEN1, 32BIT);
	MCBSP_FSETS(XCR1, XWDLEN1, 32BIT);

	// set codec sampling frequency
	DSK6713_AIC23_setFreq(hCodec, DSK6713_AIC23_FREQ_8KHZ);

	// interrupt setup
	IRQ_globalDisable();			// Globally disables interrupts
	IRQ_nmiEnable();				// Enables the NMI interrupt
	IRQ_map(IRQ_EVT_RINT1,15);		// Maps an event to a physical interrupt
	IRQ_enable(IRQ_EVT_RINT1);		// Enables the event
	IRQ_globalEnable();				// Globally enables interrupts

	DSK6713_LED_toggle(3); //Debug
	return;
}

void calc_filter(volatile float *mfc, volatile float *mfs, volatile float *buf){
	double t, y;
	int i;

	for (i=0;i<M;i++){
			t = i*CBW;				// time
			y = cos(2*PI*t);		// cosine matched filter (double)
			mfc[i] = (float) y;		// cast and store
			y = sin(2*PI*t);		// sine matched filter (double)
			mfs[i] = (float) y;     // cast and store
			buf[i] = 0;             // clear searching buffer
		}
}

//Inputs are buf, bufindex, temp.channel[0], PASSTHROUGH
extern volatile float *buf, *recbuf, *mfc, *mfs;
extern volatile short bufindex, max_recbuf, recbufindex;
extern volatile int state;

void search_for_thresh(short *channel, short pass){
    // put sample in searching buffer
	float zc, zs, z;
	int i, j;

	buf[bufindex] = (float) (*channel);  // right channel

	if (pass !=1) {
		channel[0] = 0;  // can comment this out for debugging at home
	}

    // compute incoherent correlation
    zc = 0;
    zs = 0;
	for(i=0;i<M;i++) {
		zc += mfc[i]*buf[i];
		zs += mfs[i]*buf[i];
	}
	z = zc*zc+zs*zs;

	if (z>T1) {  				// threshold exceeded?
		max_recbuf = 0;
		state = STATE_RECORD; 	// enter "recording" state (takes effect in next interrupt)
		recbufindex = M;		// start recording new samples at position M
		j = bufindex;			//
		for (i=0;i<M;i++){  	// copy samples from buf to first M elements of recbuf
			j++;   				// the first time through, this puts us at the oldest sample
			if (j>=M)
				j=0;
			recbuf[i] = (short) buf[j];
			buf[j] = 0;  		// clear out searching buffer to avoid false trigger
		}
	}
	else {
		// increment and wrap pointer
	    bufindex++;
	    if (bufindex>=M)
	    	bufindex = 0;
	}
}
