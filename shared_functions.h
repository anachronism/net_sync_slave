/*
 * shared_functions.h
 *
 *  Created on: May 20, 2015
 *      Author: IEUser
 */

#ifndef SHARED_FUNCTIONS_H_
#define SHARED_FUNCTIONS_H_

#include <math.h>
#include <stdio.h>
#include <c6x.h>
#include <csl.h>
#include <csl_mcbsp.h>
#include <csl_irq.h>

#include "dsk6713.h"
#include "dsk6713_aic23.h"
#include "dsk6713_led.h"

void calc_filter(volatile float *mfc, volatile float *mfs, volatile float *buf);
void init_codec(DSK6713_AIC23_Config config, DSK6713_AIC23_CodecHandle hCodec);
//void search_for_thresh(volatile float *buf, volatile short *recbuf, volatile short *bufindex, short *channel, volatile short *max_recbuf, volatile int *state, volatile short *recbufindex, short pass, float *mfc, float *mfs);
//void search_for_thresh(volatile float *buf, volatile float *recbuf, volatile short *bufindex, short *channel, volatile short *max_recbuf, volatile int *state, volatile short *recbufindex, short pass, float *mfc, float *mfs);
void search_for_thresh(short *channel,short pass);

#endif /* SHARED_FUNCTIONS_H_ */
