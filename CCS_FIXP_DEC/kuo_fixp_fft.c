// 
//  Project: Experiment 6.6.4 128-point intrinsics FFT - Chapter 6
//  File name: intrinsic_fft.c   
//
//  Description: Complex radix-2 decimation-in-time intrinsics FFT
//
//  For the book "Real Time Digital Signal Processing: 
//                Implementation and Application, 2nd Ed"
//                By Sen M. Kuo, Bob H. Lee, and Wenshun Tian
//                Publisher: John Wiley and Sons, Ltd
//
//  Perform in place FFT the output overwrite the input array
//
//  Using C55x C intrinsics to perform in place FFT, 
//  the output overwrite the input array
//
//  Tools used: CCS v.2.12.07
//              TMS320VC5510 DSK Rev-C
//

#include "types.h"
#include <math.h>
#include <intrindefs.h> 

#define pi  3.1415926535897
#define SFT16 16

void kuo_fixp_fft_init(int16_complex_t *W, unsigned short EXP)
{
	  unsigned short L,LE,LE1;

	  for (L=1; L<=EXP; L++)       // Create twiddle factor table
	  {
	    LE=1<<L;                   // LE=2^L=points of sub DFT
	    LE1=LE>>1;     	           // Number of butterflies in sub-DFT
	    W[L-1].re = (int16_t)((0x7fff*cos(pi/LE1))+0.5);
	    W[L-1].im = -(int16_t)((0x7fff*sin(pi/LE1))+0.5);
	  }
}

void kuo_fixp_fft_bit_rev(int16_complex_t *X, unsigned short EXP)
{
  unsigned short i,j,k;
  unsigned short N=1<<EXP; // Number of points for FFT
  unsigned short N2=N>>1;
  int16_complex_t temp;            // Temporary storage of int complex variable

  for (j=0,i=1;i<N-1;i++)
  {
    k=N2;
    while(k<=j)
    {
      j-=k;
      k>>=1;
    }
    j+=k;

    if (i<j)
    {
      temp = X[j];
      X[j] = X[i];
      X[i] = temp;
    }
  }
}

void kuo_fixp_fft(int16_complex_t *X, unsigned short EXP, int16_complex_t *W, unsigned short SCALE)
{
  int16_complex_t temp;            // Temporary storage of int complex variable
  int32_complex_t ltemp;          // Temporary storage of long complex variable
  int16_complex_t U;               // Twiddle factor W^k
  unsigned short i,j;
  unsigned short id;       // Index for lower point in butterfly 
  unsigned short N=1<<EXP; // Number of points for FFT 
  unsigned short L;        // FFT stage 
  unsigned short LE;       // Number of points in sub DFT at stage L
                           // and offset to next DFT in stage 
  unsigned short LE1;      // Number of butterflies in one DFT at
                           // stage L.  Also is offset to lower point
                           // in butterfly at stage L 
  short scale;

  kuo_fixp_fft_bit_rev(X, EXP);

  scale = 1;     
  if (SCALE == 0)         
  {
    scale = 0;
  }
  for (L=1; L<=EXP; L++)   // FFT butterfly 
  {
    LE=1<<L;               // LE=2^L=points of sub DFT 
    LE1=LE>>1;             // Number of butterflies in sub DFT 
    U.re = 0x7fff;
    U.im = 0;

    for (j=0; j<LE1;j++)
    {
      for(i=j; i<N; i+=LE) // Do the butterflies 
      {
        id=i+LE1;                 
        ltemp.re = _lsmpy(X[id].re, U.re);
        temp.re = (_smas(ltemp.re, X[id].im, U.im)>>SFT16);  
        temp.re = _sadd(temp.re, 1)>>scale; // Rounding & scale 
        ltemp.im = _lsmpy(X[id].im, U.re);
        temp.im = (_smac(ltemp.im, X[id].re, U.im)>>SFT16);
        temp.im = _sadd(temp.im, 1)>>scale; // Rounding & scale 
        X[id].re = _ssub(X[i].re>>scale, temp.re);
        X[id].im = _ssub(X[i].im>>scale, temp.im);
        X[i].re = _sadd(X[i].re>>scale, temp.re);
        X[i].im = _sadd(X[i].im>>scale, temp.im);
      }
          
      // Recursive compute W^k as W*W^(k-1) 
      ltemp.re = _lsmpy(U.re, W[L-1].re);
      ltemp.re = _smas(ltemp.re, U.im, W[L-1].im);  
      ltemp.im = _lsmpy(U.re, W[L-1].im);
      ltemp.im = _smac(ltemp.im, U.im, W[L-1].re);
      U.re = ltemp.re>>SFT16;            
      U.im = ltemp.im>>SFT16;
    }
  }
}  
