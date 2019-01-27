// 
//  Project: Experiment 6.6.4 128-point intrinsics FFT - Chapter 6
//  File name: intrinsic_fft.h   
//
//  Description: This is the header file for intrinsic FFT experiment
//
//  For the book "Real Time Digital Signal Processing: 
//                Implementation and Application, 2nd Ed"
//                By Sen M. Kuo, Bob H. Lee, and Wenshun Tian
//                Publisher: John Wiley and Sons, Ltd
//
//  Tools used: CCS v.2.12.07
//              TMS320VC5510 DSK Rev-C
//
#include "types.h"

extern void kuo_fixp_fft_init(int16_complex_t *, unsigned short);
extern void kuo_fixp_fft_bit_rev(int16_complex_t *, unsigned short);
extern void kuo_fixp_fft(int16_complex_t *, unsigned short, int16_complex_t *, unsigned short);

