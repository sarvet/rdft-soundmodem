/*
    DFT routines used by Wyman1x digital communications programs
    Copyright (C) 2001  Barry Sanderson

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/*  A mailing address for the author is:
        Barry Sanderson
        1725 N. Bolton Ave.
        Indianapolis, IN 46218  USA
*/
/* Rev 25 Nov  2001	version 1.1.0, for initial public release */
/* Rev 17 Dec 1992  Function to peform foward or inverse Discrete
        Fast Fourier Transform on an integer-power-of-2 number of
        complex floating point data values. */

#include "dcom.h"
#include <math.h>

/* calculate quarter cycle sine table and return frequency step */
double init_dft3(unsigned int dft_siz, float t_step, float *snqt)
{
  unsigned int i;
  double       sstp;

  sstp = (double) (PI2 / (double) (dft_siz));
  for(i = 0; i <= (dft_siz / 4); i++)
    snqt[i] = (float) (sin(sstp * i));
  return ((double) (1.0 / (t_step * dft_siz)));
}

void dfft3f(float *inr, float *ini, float *outr, float *outi, float *sint,
            unsigned int N, unsigned int inv, unsigned int real_in)
{ /* inr		points to real part of input
     ini		points to imaginary part of input
     outr	points to real part of result
     outi	points to imaginary part of result
     sint	points to a quarter period sine table  [0:(N+1)/4]
     N		is maximum index on arrays: inr, ini, outr, and outi [0:N]
     inv		is 1 if inverse DFT is to be computed
     real_in	is 1 if input imaginary part is 0 */

  float         Rr, Ri, R1r, R1i, S1, S2;
  unsigned int  I, J, K, L, M, NT, NH, NQ, IR, IT;
  unsigned long ID, ND;

  /* copy input data to output array, setting imaginary part to zero, if
         input is not complex */
  for(I = 0; I <= N; I++)
    *(outr + I) = *(inr + I);
  if((real_in == 1) && (inv != 1))
  {
    for(I = 0; I <= N; I++)
      *(outi + I) = 0.0;
  }
  else
  {
    for(I = 0; I <= N; I++)
      *(outi + I) = *(ini + I);
  }

  NT = N;
  NQ = NT / 4 + 1;
  NH = NQ + NQ;

  /* reorder initial data, so result will be in desired order */
  J = 0;
  for(I = 0; I <= NT - 2; I++)
  {
    if(J > I)
    {
      Rr          = *(outr + J);
      Ri          = *(outi + J);
      *(outr + J) = *(outr + I);
      *(outi + J) = *(outi + I);
      *(outr + I) = Rr;
      *(outi + I) = Ri;
    }
    K = NH;
    while(J >= K)
    {
      J = J - K;
      K = K / 2;
    }
    J = J + K;
  }

  /* calculate the discrete Fourier transform */
  ID = 1;
  IR = 1;
  M  = 0;
  IT = 0;
tran:
  I  = ID;
  ID = ID + ID;
  for(J = 0; J < I; J++)
  {
    S2 = *(sint + M);
    if(inv != 1)
      S2 = -S2;
    S1 = *(sint + NQ - M);
    if(J >= IR)
    {
      M  = M - IT;
      S1 = -S1;
    }
    else
      M = M + IT;
    for(ND = J; ND <= NT; ND += ID)
    {
      K           = ND;
      L           = K + I;
      /* butterfly computation */
      R1r         = *(outr + L);
      R1i         = *(outi + L);
      Rr          = S1 * R1r - S2 * R1i;
      Ri          = S1 * R1i + S2 * R1r;
      R1r         = *(outr + K);
      R1i         = *(outi + K);
      *(outr + K) = R1r + Rr;
      *(outi + K) = R1i + Ri;
      *(outr + L) = R1r - Rr;
      *(outi + L) = R1i - Ri;
    }
  }
  IR = I;
  IT = NQ / I;
  if(I < NH)
    goto tran;
  if(inv != 1)
    return;
  Rr = N + 1;
  for(I = 0; I <= N; I++)
  {
    *(outr + I) /= Rr;
    *(outi + I) /= Rr;
  }
  return;
}
