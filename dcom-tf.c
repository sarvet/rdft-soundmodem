/*
    filters used by Wyman1x digital communications programs
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
/* Rev 25 Nov 2001	version 1.1.0, for initial public release */

/* Rev 26 Feb 2001	introduced filters for 122.5 baud, 230 Hz spacing */
/* Rev 22 Jan 2001	fixed gain on 2nd stage of interpolation filter */
/* Rev 21 Jan 2001	filters for 230 Hz spacing, 115 baud introduced */
/* Rev 29 Dec 2000	narrowband filters with 11025 Hz sample rate, for
	fine tuning introduced */
/* Rev 09 Oct 2000	frevt 23, Qt 12, introduced: for use with phase
	modulated version */

#include "dcom-t.h"

/**   Sat Dec 16 10:30:02 2000
25X stage of 25X-4X-2X decimator for fine tuning of digital HF signal,
	for 11025 Hz sample rate.  Overall filter:
	passes 0-2 Hz	with a gain of 2,
	stops > 15 Hz.
        Filter Length =  59
********* IMPULSE RESPONSE ***********/
ftype	ds1nb2cf[]={
	1.86502549e-03,	/*   H(  1) =  1.86502549e-03 =  H( 59)*/
	1.58623036e-03,	/*   H(  2) =  1.58623036e-03 =  H( 58)*/
	2.22568284e-03,	/*   H(  3) =  2.22568284e-03 =  H( 57)*/
	3.00472416e-03,	/*   H(  4) =  3.00472416e-03 =  H( 56)*/
	3.93414963e-03,	/*   H(  5) =  3.93414963e-03 =  H( 55)*/
	5.02365502e-03,	/*   H(  6) =  5.02365502e-03 =  H( 54)*/
	6.27943175e-03,	/*   H(  7) =  6.27943175e-03 =  H( 53)*/
	7.70607591e-03,	/*   H(  8) =  7.70607591e-03 =  H( 52)*/
	9.30345710e-03,	/*   H(  9) =  9.30345710e-03 =  H( 51)*/
	1.10695558e-02,	/*   H( 10) =  1.10695558e-02 =  H( 50)*/
	1.29967174e-02,	/*   H( 11) =  1.29967174e-02 =  H( 49)*/
	1.50743620e-02,	/*   H( 12) =  1.50743620e-02 =  H( 48)*/
	1.72882490e-02,	/*   H( 13) =  1.72882490e-02 =  H( 47)*/
	1.96186025e-02,	/*   H( 14) =  1.96186025e-02 =  H( 46)*/
	2.20438521e-02,	/*   H( 15) =  2.20438521e-02 =  H( 45)*/
	2.45377775e-02,	/*   H( 16) =  2.45377775e-02 =  H( 44)*/
	2.70707738e-02,	/*   H( 17) =  2.70707738e-02 =  H( 43)*/
	2.96112988e-02,	/*   H( 18) =  2.96112988e-02 =  H( 42)*/
	3.21257263e-02,	/*   H( 19) =  3.21257263e-02 =  H( 41)*/
	3.45787816e-02,	/*   H( 20) =  3.45787816e-02 =  H( 40)*/
	3.69347185e-02,	/*   H( 21) =  3.69347185e-02 =  H( 39)*/
	3.91585529e-02,	/*   H( 22) =  3.91585529e-02 =  H( 38)*/
	4.12154645e-02,	/*   H( 23) =  4.12154645e-02 =  H( 37)*/
	4.30725142e-02,	/*   H( 24) =  4.30725142e-02 =  H( 36)*/
	4.46998477e-02,	/*   H( 25) =  4.46998477e-02 =  H( 35)*/
	4.60705683e-02,	/*   H( 26) =  4.60705683e-02 =  H( 34)*/
	4.71617989e-02,	/*   H( 27) =  4.71617989e-02 =  H( 33)*/
	4.79548611e-02,	/*   H( 28) =  4.79548611e-02 =  H( 32)*/
	4.84359413e-02,	/*   H( 29) =  4.84359413e-02 =  H( 31)*/
	4.85971384e-02,	/*   H( 30) =  4.85971384e-02 =  H( 30)*/
};
/**   Sat Dec 16 10:35:07 2000
        Filter Length =  13
********* IMPULSE RESPONSE ***********/
ftype	ds2nb2cf[]={
	1.38167304e-03,	/*   H(  1) =  1.38167304e-03 =  H( 13)*/
	1.39118424e-02,	/*   H(  2) =  1.39118424e-02 =  H( 12)*/
	4.66636084e-02,	/*   H(  3) =  4.66636084e-02 =  H( 11)*/
	1.04117669e-01,	/*   H(  4) =  1.04117669e-01 =  H( 10)*/
	1.75465882e-01,	/*   H(  5) =  1.75465882e-01 =  H(  9)*/
	2.36136258e-01,	/*   H(  6) =  2.36136258e-01 =  H(  8)*/
	2.60084450e-01,	/*   H(  7) =  2.60084450e-01 =  H(  7)*/
};
/**   Sat Dec 16 10:37:51 2000
        Filter Length =  23
********* IMPULSE RESPONSE ***********/
ftype	ds3nb2cf[]={
	-4.06523142e-03,	/*   H(  1) = -4.06523142e-03 =  H( 23)*/
	-7.66611099e-03,	/*   H(  2) = -7.66611099e-03 =  H( 22)*/
	-1.16024539e-02,	/*   H(  3) = -1.16024539e-02 =  H( 21)*/
	-1.28093325e-02,	/*   H(  4) = -1.28093325e-02 =  H( 20)*/
	-8.06417502e-03,	/*   H(  5) = -8.06417502e-03 =  H( 19)*/
	5.38879819e-03,	/*   H(  6) =  5.38879819e-03 =  H( 18)*/
	2.85850838e-02,	/*   H(  7) =  2.85850838e-02 =  H( 17)*/
	5.98970354e-02,	/*   H(  8) =  5.98970354e-02 =  H( 16)*/
	9.48844925e-02,	/*   H(  9) =  9.48844925e-02 =  H( 15)*/
	1.27134368e-01,	/*   H( 10) =  1.27134368e-01 =  H( 14)*/
	1.49918631e-01,	/*   H( 11) =  1.49918631e-01 =  H( 13)*/
	1.58144489e-01,	/*   H( 12) =  1.58144489e-01 =  H( 12)*/
};
/**#$   Mon Feb 26 13:06:43 2001
5X stage of 5X-4X-2X decimator for phase modulation version of digital HF signal,
	prior to transmission.  Overall filter:
	passes 0-85.546 Hz	with a gain of 2,
	stops > 137.8125 Hz.
        Filter Length =  17
********* IMPULSE RESPONSE *********$#**/
ftype	ds1cf[]={
	-4.66673682e-03,	/*   H(  1) = -4.66673682e-03 =  H( 17)*/
	-9.77225602e-03,	/*   H(  2) = -9.77225602e-03 =  H( 16)*/
	-9.81219579e-03,	/*   H(  3) = -9.81219579e-03 =  H( 15)*/
	6.17947802e-03,	/*   H(  4) =  6.17947802e-03 =  H( 14)*/
	4.70137782e-02,	/*   H(  5) =  4.70137782e-02 =  H( 13)*/
	1.11866206e-01,	/*   H(  6) =  1.11866206e-01 =  H( 12)*/
	1.86312556e-01,	/*   H(  7) =  1.86312556e-01 =  H( 11)*/
	2.46336490e-01,	/*   H(  8) =  2.46336490e-01 =  H( 10)*/
	2.69410163e-01,	/*   H(  9) =  2.69410163e-01 =  H(  9)*/
};
/**#$   Mon Feb 26 13:14:20 2001
4X stage of 5X-4X-2X decimator for phase modulation version of digital HF signal,
	prior to transmission.  Overall filter:
	passes 0-85.546 Hz	with a gain of 2,
	stops > 137.8125 Hz.
        Filter Length =  19
********* IMPULSE RESPONSE *********$#**/
ftype	ds2cf[]={
	-1.51501581e-05,	/*   H(  1) = -1.51501581e-05 =  H( 19)*/
	-8.60584527e-03,	/*   H(  2) = -8.60584527e-03 =  H( 18)*/
	-2.04488263e-02,	/*   H(  3) = -2.04488263e-02 =  H( 17)*/
	-2.86971144e-02,	/*   H(  4) = -2.86971144e-02 =  H( 16)*/
	-1.73100699e-02,	/*   H(  5) = -1.73100699e-02 =  H( 15)*/
	2.78713070e-02,	/*   H(  6) =  2.78713070e-02 =  H( 14)*/
	1.07877158e-01,	/*   H(  7) =  1.07877158e-01 =  H( 13)*/
	2.04094678e-01,	/*   H(  8) =  2.04094678e-01 =  H( 12)*/
	2.83449382e-01,	/*   H(  9) =  2.83449382e-01 =  H( 11)*/
	3.14247578e-01,	/*   H( 10) =  3.14247578e-01 =  H( 10)*/
};
/**#$   Mon Feb 26 13:15:16 2001
2X stage of 5X-4X-2X decimator for phase modulation version of digital HF signal,
	prior to transmission.  Overall filter:
	passes 0-85.546 Hz	with a gain of 2,
	stops > 137.8125 Hz.
        Filter Length =  29
********* IMPULSE RESPONSE *********$#**/
ftype	ds3cf[]={
	-2.21228367e-03,	/*   H(  1) = -2.21228367e-03 =  H( 29)*/
	-3.76911904e-03,	/*   H(  2) = -3.76911904e-03 =  H( 28)*/
	3.18166800e-03,	/*   H(  3) =  3.18166800e-03 =  H( 27)*/
	8.79500899e-03,	/*   H(  4) =  8.79500899e-03 =  H( 26)*/
	1.97633961e-03,	/*   H(  5) =  1.97633961e-03 =  H( 25)*/
	-1.53443608e-02,	/*   H(  6) = -1.53443608e-02 =  H( 24)*/
	-1.52472183e-02,	/*   H(  7) = -1.52472183e-02 =  H( 23)*/
	1.47586232e-02,	/*   H(  8) =  1.47586232e-02 =  H( 22)*/
	3.78787741e-02,	/*   H(  9) =  3.78787741e-02 =  H( 21)*/
	3.85557208e-03,	/*   H( 10) =  3.85557208e-03 =  H( 20)*/
	-6.46372661e-02,	/*   H( 11) = -6.46372661e-02 =  H( 19)*/
	-6.13453612e-02,	/*   H( 12) = -6.13453612e-02 =  H( 18)*/
	8.66528228e-02,	/*   H( 13) =  8.66528228e-02 =  H( 17)*/
	3.01640928e-01,	/*   H( 14) =  3.01640928e-01 =  H( 16)*/
	4.04814333e-01,	/*   H( 15) =  4.04814333e-01 =  H( 15)*/
};


/**   Fri Jun 15 21:47:24 2001
        Filter Length =  33
********* IMPULSE RESPONSE ***********/
ftype	lp01cf[lp1dly+1]={
	1.23847800e-03,	/*   H(  1) =  1.23847800e-03 =  H( 33)*/
	-6.40280123e-05,	/*   H(  2) = -6.40280123e-05 =  H( 32)*/
	-3.26283998e-03,	/*   H(  3) = -3.26283998e-03 =  H( 31)*/
	-3.12315254e-03,	/*   H(  4) = -3.12315254e-03 =  H( 30)*/
	3.97119531e-03,	/*   H(  5) =  3.97119531e-03 =  H( 29)*/
	9.35720000e-03,	/*   H(  6) =  9.35720000e-03 =  H( 28)*/
	2.12917250e-04,	/*   H(  7) =  2.12917250e-04 =  H( 27)*/
	-1.66874249e-02,	/*   H(  8) = -1.66874249e-02 =  H( 26)*/
	-1.38472766e-02,	/*   H(  9) = -1.38472766e-02 =  H( 25)*/
	1.74594857e-02,	/*   H( 10) =  1.74594857e-02 =  H( 24)*/
	3.75460498e-02,	/*   H( 11) =  3.75460498e-02 =  H( 23)*/
	3.75332660e-04,	/*   H( 12) =  3.75332660e-04 =  H( 22)*/
	-6.63000345e-02,	/*   H( 13) = -6.63000345e-02 =  H( 21)*/
	-5.83123863e-02,	/*   H( 14) = -5.83123863e-02 =  H( 20)*/
	9.02222842e-02,	/*   H( 15) =  9.02222842e-02 =  H( 19)*/
	3.00431788e-01,	/*   H( 16) =  3.00431788e-01 =  H( 18)*/
	4.00438458e-01,	/*   H( 17) =  4.00438458e-01 =  H( 17)*/
};

void	dfilter(ftype *in, ftype *out, unsigned int count,
		unsigned int delay, ftype *coef)
{ftype		acum;
 register int	si, k;
 unsigned int	j;
/*** For valid results in[-delay] through in[-1] must contain valid data.
     For valid results in[(count-1)+1] through in[(count-1)+delay]
	 must contain valid data. ***/

 k=0;
 for (j=0; j< count; j++)
 {acum= in[k] * coef[delay];	/* special case */
  for (si=0; si < delay; si++)
   acum+= (in[si-delay+k] + in[delay-si+k])* coef[si];
  out[j]= acum;			/* save filtered result */
  k++;
 }
}



/* Rev 11 May 1998 converted to C using doubles for data and coefficients
dfilterx:	;Subroutine to calculate and store filter results,
		; with decimation
;	"in" points to element of input data array, whose filtered value is
;		to be calculated first.
;	"out" points to where results are to be stored.
;	"count" is the number of output values to calculate and store.
;	"dcf" is decimation factor
;	An odd number of filter taps is assumed.
*/
void	ddfilter(ftype *in, ftype *out, unsigned int count,
		unsigned int delay, unsigned int dcf, ftype *coef)
{ftype		acum;
 register int	si, k;
 unsigned int	j;
/*** For valid results in[-delay] through in[-1] must contain valid data.
     For valid results in[dcf*(count-1)+1] through in[dcf*(count-1)+delay]
	 must contain valid data. ***/

 k=0;
 for (j=0; j< count; j++)
 {acum= in[k] * coef[delay];	/* special case */
  for (si=0; si < delay; si++)
   acum+= (in[si-delay+k] + in[delay-si+k])* coef[si];
  out[j]= acum;			/* save filtered result */
  k+= dcf;
 }
}


void	ds1t(ftype *in, ftype *out, unsigned int count)
{ddfilter(in, out, count, ds1dly, ds1dcf, ds1cf);}

void	ds2t(ftype *in, ftype *out, unsigned int count)
{ddfilter(in, out, count, ds2dly, ds2dcf, ds2cf);}

void	ds3t(ftype *in, ftype *out, unsigned int count)
{ddfilter(in, out, count, ds3dly, ds3dcf, ds3cf);}

void	ds1nb2(ftype *in, ftype *out, unsigned int count)
{ddfilter(in, out, count, ds1nb2dly, ds1nb2dcf, ds1nb2cf);}

void	ds2nb2(ftype *in, ftype *out, unsigned int count)
{ddfilter(in, out, count, ds2nb2dly, ds2nb2dcf, ds2nb2cf);}

void	ds3nb2(ftype *in, ftype *out, unsigned int count)
{ddfilter(in, out, count, ds3nb2dly, ds3nb2dcf, ds3nb2cf);}

void  lpf01(ftype *in, ftype *out, unsigned int count)
{dfilter(in, out, count, lp1dly, &lp01cf[0]);}

/* Rev 27 Aug 1998	changed for Linux */
/* Rev 11 Jun 1998 */
/* Rev 13 May 1998	conversion from assembly to C */


#include <malloc.h>


/*  2 times interpolation using 24 original samples

   Sun Jan 21 12:19:19 2001
Overall filter:
	Passes: 0-103.125 Hz
	Stops: > 137.8125 Hz
        Filter Length =  47
********* IMPULSE RESPONSE **********/
ftype is1cI[12]={
/* Start of coefficients for  1/ 2 interpolation */
	8.75213591e-04,	/*   H(  1) =  8.75213591e-04 =  H( 47) */
	-2.68876483e-03,	/*   H(  3) = -2.68876483e-03 =  H( 45) */
	6.24762010e-03,	/*   H(  5) =  6.24762010e-03 =  H( 43) */
	-1.15779210e-02,	/*   H(  7) = -1.15779210e-02 =  H( 41) */
	1.73623580e-02,	/*   H(  9) =  1.73623580e-02 =  H( 39) */
	-2.13883221e-02,	/*   H( 11) = -2.13883221e-02 =  H( 37) */
	2.03674845e-02,	/*   H( 13) =  2.03674845e-02 =  H( 35) */
	-1.00583965e-02,	/*   H( 15) = -1.00583965e-02 =  H( 33) */
	-1.50937997e-02,	/*   H( 17) = -1.50937997e-02 =  H( 31) */
	6.48910925e-02,	/*   H( 19) =  6.48910925e-02 =  H( 29) */
	-1.71027765e-01,	/*   H( 21) = -1.71027765e-01 =  H( 27) */
	6.22227669e-01,	/*   H( 23) =  6.22227669e-01 =  H( 25) */
#if (1==2)
	6.22227669e-01,	/*   H( 23) =  6.22227669e-01 =  H( 25) */
	-1.71027765e-01,	/*   H( 21) = -1.71027765e-01 =  H( 27) */
	6.48910925e-02,	/*   H( 19) =  6.48910925e-02 =  H( 29) */
	-1.50937997e-02,	/*   H( 17) = -1.50937997e-02 =  H( 31) */
	-1.00583965e-02,	/*   H( 15) = -1.00583965e-02 =  H( 33) */
	2.03674845e-02,	/*   H( 13) =  2.03674845e-02 =  H( 35) */
	-2.13883221e-02,	/*   H( 11) = -2.13883221e-02 =  H( 37) */
	1.73623580e-02,	/*   H(  9) =  1.73623580e-02 =  H( 39) */
	-1.15779210e-02,	/*   H(  7) = -1.15779210e-02 =  H( 41) */
	6.24762010e-03,	/*   H(  5) =  6.24762010e-03 =  H( 43) */
	-2.68876483e-03,	/*   H(  3) = -2.68876483e-03 =  H( 45) */
	8.75213591e-04,	/*   H(  1) =  8.75213591e-04 =  H( 47) */
#endif
};
/* 20 times interpolation using  8 original samples

   Mon Jan 22 13:19:52 2001
        Filter Length = 159
********* IMPULSE RESPONSE **********/
ftype is2cI[19][ 8]={
/* Start of coefficients for  1/20 interpolation */
{	3.39443982e-03,	/*   H( 19) =  3.39443982e-03 =  H(141) */
	-7.47941062e-03,	/*   H( 39) = -7.47941062e-03 =  H(121) */
	-4.72677639e-03,	/*   H( 59) = -4.72677639e-03 =  H(101) */
	9.53945220e-01,	/*   H( 79) =  9.53945220e-01 =  H( 81) */
	7.89186135e-02,	/*   H( 61) =  7.89186135e-02 =  H( 99) */
	-3.13919224e-02,	/*   H( 41) = -3.13919224e-02 =  H(119) */
	9.04073752e-03,	/*   H( 21) =  9.04073752e-03 =  H(139) */
	-2.12436169e-03,	/*   H(  1) = -2.12436169e-03 =  H(159) */
},/* Start of coefficients for  2/20 interpolation */
{	1.06672023e-03,	/*   H( 18) =  1.06672023e-03 =  H(142) */
	2.73180893e-03,	/*   H( 38) =  2.73180893e-03 =  H(122) */
	-3.90473157e-02,	/*   H( 58) = -3.90473157e-02 =  H(102) */
	9.41873729e-01,	/*   H( 78) =  9.41873729e-01 =  H( 82) */
	1.27669692e-01,	/*   H( 62) =  1.27669692e-01 =  H( 98) */
	-4.46800031e-02,	/*   H( 42) = -4.46800031e-02 =  H(118) */
	1.22673567e-02,	/*   H( 22) =  1.22673567e-02 =  H(138) */
	-1.50557060e-03,	/*   H(  2) = -1.50557060e-03 =  H(158) */
},/* Start of coefficients for  3/20 interpolation */
{	-8.97647638e-04,	/*   H( 17) = -8.97647638e-04 =  H(143) */
	1.15931705e-02,	/*   H( 37) =  1.15931705e-02 =  H(123) */
	-6.81654587e-02,	/*   H( 57) = -6.81654587e-02 =  H(103) */
	9.21985090e-01,	/*   H( 77) =  9.21985090e-01 =  H( 83) */
	1.80476859e-01,	/*   H( 63) =  1.80476859e-01 =  H( 97) */
	-5.84908836e-02,	/*   H( 43) = -5.84908836e-02 =  H(117) */
	1.56717729e-02,	/*   H( 23) =  1.56717729e-02 =  H(137) */
	-1.98004348e-03,	/*   H(  3) = -1.98004348e-03 =  H(157) */
},/* Start of coefficients for  4/20 interpolation */
{	-2.49319524e-03,	/*   H( 16) = -2.49319524e-03 =  H(144) */
	1.90289281e-02,	/*   H( 36) =  1.90289281e-02 =  H(124) */
	-9.20635462e-02,	/*   H( 56) = -9.20635462e-02 =  H(104) */
	8.94620001e-01,	/*   H( 76) =  8.94620001e-01 =  H( 84) */
	2.36788452e-01,	/*   H( 64) =  2.36788452e-01 =  H( 96) */
	-7.24865869e-02,	/*   H( 44) = -7.24865869e-02 =  H(116) */
	1.91630404e-02,	/*   H( 24) =  1.91630404e-02 =  H(136) */
	-2.50134687e-03,	/*   H(  4) = -2.50134687e-03 =  H(156) */
},/* Start of coefficients for  5/20 interpolation */
{	-3.72756622e-03,	/*   H( 15) = -3.72756622e-03 =  H(145) */
	2.50113048e-02,	/*   H( 35) =  2.50113048e-02 =  H(125) */
	-1.10823952e-01,	/*   H( 55) = -1.10823952e-01 =  H(105) */
	8.60244393e-01,	/*   H( 75) =  8.60244393e-01 =  H( 85) */
	2.95957297e-01,	/*   H( 65) =  2.95957297e-01 =  H( 95) */
	-8.62829685e-02,	/*   H( 45) = -8.62829685e-02 =  H(115) */
	2.26327982e-02,	/*   H( 25) =  2.26327982e-02 =  H(135) */
	-3.05396295e-03,	/*   H(  5) = -3.05396295e-03 =  H(155) */
},/* Start of coefficients for  6/20 interpolation */
{	-4.61932970e-03,	/*   H( 14) = -4.61932970e-03 =  H(146) */
	2.95562688e-02,	/*   H( 34) =  2.95562688e-02 =  H(126) */
	-1.24622084e-01,	/*   H( 54) = -1.24622084e-01 =  H(106) */
	8.19439113e-01,	/*   H( 74) =  8.19439113e-01 =  H( 86) */
	3.57250899e-01,	/*   H( 66) =  3.57250899e-01 =  H( 94) */
	-9.94540974e-02,	/*   H( 46) = -9.94540974e-02 =  H(114) */
	2.59566940e-02,	/*   H( 26) =  2.59566940e-02 =  H(134) */
	-3.61745874e-03,	/*   H(  6) = -3.61745874e-03 =  H(154) */
},/* Start of coefficients for  7/20 interpolation */
{	-5.19614015e-03,	/*   H( 13) = -5.19614015e-03 =  H(147) */
	3.27178054e-02,	/*   H( 33) =  3.27178054e-02 =  H(127) */
	-1.33715853e-01,	/*   H( 53) = -1.33715853e-01 =  H(107) */
	7.72887051e-01,	/*   H( 73) =  7.72887051e-01 =  H( 87) */
	4.19863969e-01,	/*   H( 67) =  4.19863969e-01 =  H( 93) */
	-1.11539721e-01,	/*   H( 47) = -1.11539721e-01 =  H(113) */
	2.89965551e-02,	/*   H( 27) =  2.89965551e-02 =  H(133) */
	-4.16626362e-03,	/*   H(  7) = -4.16626362e-03 =  H(153) */
},/* Start of coefficients for  8/20 interpolation */
{	-5.49218571e-03,	/*   H( 12) = -5.49218571e-03 =  H(148) */
	3.45824845e-02,	/*   H( 32) =  3.45824845e-02 =  H(128) */
	-1.38437182e-01,	/*   H( 52) = -1.38437182e-01 =  H(108) */
	7.21358061e-01,	/*   H( 72) =  7.21358061e-01 =  H( 88) */
	4.82931882e-01,	/*   H( 68) =  4.82931882e-01 =  H( 92) */
	-1.22055210e-01,	/*   H( 48) = -1.22055210e-01 =  H(112) */
	3.16035897e-02,	/*   H( 28) =  3.16035897e-02 =  H(132) */
	-4.66988748e-03,	/*   H(  8) = -4.66988748e-03 =  H(152) */
},/* Start of coefficients for  9/20 interpolation */
{	-5.54615725e-03,	/*   H( 11) = -5.54615725e-03 =  H(149) */
	3.52633893e-02,	/*   H( 31) =  3.52633893e-02 =  H(129) */
	-1.39176115e-01,	/*   H( 51) = -1.39176115e-01 =  H(109) */
	6.65692031e-01,	/*   H( 71) =  6.65692031e-01 =  H( 89) */
	5.45547903e-01,	/*   H( 69) =  5.45547903e-01 =  H( 91) */
	-1.30501211e-01,	/*   H( 49) = -1.30501211e-01 =  H(111) */
	3.36221270e-02,	/*   H( 29) =  3.36221270e-02 =  H(131) */
	-5.09371655e-03,	/*   H(  9) = -5.09371655e-03 =  H(151) */
},/* Start of coefficients for 10/20 interpolation */
{	-5.39937709e-03,	/*   H( 10) = -5.39937709e-03 =  H(150) */
	3.48940268e-02,	/*   H( 30) =  3.48940268e-02 =  H(130) */
	-1.36373043e-01,	/*   H( 50) = -1.36373043e-01 =  H(110) */
	6.06780469e-01,	/*   H( 70) =  6.06780469e-01 =  H( 90) */
	6.06780469e-01,	/*   H( 70) =  6.06780469e-01 =  H( 90) */
	-1.36373043e-01,	/*   H( 50) = -1.36373043e-01 =  H(110) */
	3.48940268e-02,	/*   H( 30) =  3.48940268e-02 =  H(130) */
	-5.39937709e-03,	/*   H( 10) = -5.39937709e-03 =  H(150) */
},/* Start of coefficients for 11/20 interpolation */
{	-5.09371655e-03,	/*   H(  9) = -5.09371655e-03 =  H(151) */
	3.36221270e-02,	/*   H( 29) =  3.36221270e-02 =  H(131) */
	-1.30501211e-01,	/*   H( 49) = -1.30501211e-01 =  H(111) */
	5.45547903e-01,	/*   H( 69) =  5.45547903e-01 =  H( 91) */
	6.65692031e-01,	/*   H( 71) =  6.65692031e-01 =  H( 89) */
	-1.39176115e-01,	/*   H( 51) = -1.39176115e-01 =  H(109) */
	3.52633893e-02,	/*   H( 31) =  3.52633893e-02 =  H(129) */
	-5.54615725e-03,	/*   H( 11) = -5.54615725e-03 =  H(149) */
},/* Start of coefficients for 12/20 interpolation */
{	-4.66988748e-03,	/*   H(  8) = -4.66988748e-03 =  H(152) */
	3.16035897e-02,	/*   H( 28) =  3.16035897e-02 =  H(132) */
	-1.22055210e-01,	/*   H( 48) = -1.22055210e-01 =  H(112) */
	4.82931882e-01,	/*   H( 68) =  4.82931882e-01 =  H( 92) */
	7.21358061e-01,	/*   H( 72) =  7.21358061e-01 =  H( 88) */
	-1.38437182e-01,	/*   H( 52) = -1.38437182e-01 =  H(108) */
	3.45824845e-02,	/*   H( 32) =  3.45824845e-02 =  H(128) */
	-5.49218571e-03,	/*   H( 12) = -5.49218571e-03 =  H(148) */
},/* Start of coefficients for 13/20 interpolation */
{	-4.16626362e-03,	/*   H(  7) = -4.16626362e-03 =  H(153) */
	2.89965551e-02,	/*   H( 27) =  2.89965551e-02 =  H(133) */
	-1.11539721e-01,	/*   H( 47) = -1.11539721e-01 =  H(113) */
	4.19863969e-01,	/*   H( 67) =  4.19863969e-01 =  H( 93) */
	7.72887051e-01,	/*   H( 73) =  7.72887051e-01 =  H( 87) */
	-1.33715853e-01,	/*   H( 53) = -1.33715853e-01 =  H(107) */
	3.27178054e-02,	/*   H( 33) =  3.27178054e-02 =  H(127) */
	-5.19614015e-03,	/*   H( 13) = -5.19614015e-03 =  H(147) */
},/* Start of coefficients for 14/20 interpolation */
{	-3.61745874e-03,	/*   H(  6) = -3.61745874e-03 =  H(154) */
	2.59566940e-02,	/*   H( 26) =  2.59566940e-02 =  H(134) */
	-9.94540974e-02,	/*   H( 46) = -9.94540974e-02 =  H(114) */
	3.57250899e-01,	/*   H( 66) =  3.57250899e-01 =  H( 94) */
	8.19439113e-01,	/*   H( 74) =  8.19439113e-01 =  H( 86) */
	-1.24622084e-01,	/*   H( 54) = -1.24622084e-01 =  H(106) */
	2.95562688e-02,	/*   H( 34) =  2.95562688e-02 =  H(126) */
	-4.61932970e-03,	/*   H( 14) = -4.61932970e-03 =  H(146) */
},/* Start of coefficients for 15/20 interpolation */
{	-3.05396295e-03,	/*   H(  5) = -3.05396295e-03 =  H(155) */
	2.26327982e-02,	/*   H( 25) =  2.26327982e-02 =  H(135) */
	-8.62829685e-02,	/*   H( 45) = -8.62829685e-02 =  H(115) */
	2.95957297e-01,	/*   H( 65) =  2.95957297e-01 =  H( 95) */
	8.60244393e-01,	/*   H( 75) =  8.60244393e-01 =  H( 85) */
	-1.10823952e-01,	/*   H( 55) = -1.10823952e-01 =  H(105) */
	2.50113048e-02,	/*   H( 35) =  2.50113048e-02 =  H(125) */
	-3.72756622e-03,	/*   H( 15) = -3.72756622e-03 =  H(145) */
},/* Start of coefficients for 16/20 interpolation */
{	-2.50134687e-03,	/*   H(  4) = -2.50134687e-03 =  H(156) */
	1.91630404e-02,	/*   H( 24) =  1.91630404e-02 =  H(136) */
	-7.24865869e-02,	/*   H( 44) = -7.24865869e-02 =  H(116) */
	2.36788452e-01,	/*   H( 64) =  2.36788452e-01 =  H( 96) */
	8.94620001e-01,	/*   H( 76) =  8.94620001e-01 =  H( 84) */
	-9.20635462e-02,	/*   H( 56) = -9.20635462e-02 =  H(104) */
	1.90289281e-02,	/*   H( 36) =  1.90289281e-02 =  H(124) */
	-2.49319524e-03,	/*   H( 16) = -2.49319524e-03 =  H(144) */
},/* Start of coefficients for 17/20 interpolation */
{	-1.98004348e-03,	/*   H(  3) = -1.98004348e-03 =  H(157) */
	1.56717729e-02,	/*   H( 23) =  1.56717729e-02 =  H(137) */
	-5.84908836e-02,	/*   H( 43) = -5.84908836e-02 =  H(117) */
	1.80476859e-01,	/*   H( 63) =  1.80476859e-01 =  H( 97) */
	9.21985090e-01,	/*   H( 77) =  9.21985090e-01 =  H( 83) */
	-6.81654587e-02,	/*   H( 57) = -6.81654587e-02 =  H(103) */
	1.15931705e-02,	/*   H( 37) =  1.15931705e-02 =  H(123) */
	-8.97647638e-04,	/*   H( 17) = -8.97647638e-04 =  H(143) */
},/* Start of coefficients for 18/20 interpolation */
{	-1.50557060e-03,	/*   H(  2) = -1.50557060e-03 =  H(158) */
	1.22673567e-02,	/*   H( 22) =  1.22673567e-02 =  H(138) */
	-4.46800031e-02,	/*   H( 42) = -4.46800031e-02 =  H(118) */
	1.27669692e-01,	/*   H( 62) =  1.27669692e-01 =  H( 98) */
	9.41873729e-01,	/*   H( 78) =  9.41873729e-01 =  H( 82) */
	-3.90473157e-02,	/*   H( 58) = -3.90473157e-02 =  H(102) */
	2.73180893e-03,	/*   H( 38) =  2.73180893e-03 =  H(122) */
	1.06672023e-03,	/*   H( 18) =  1.06672023e-03 =  H(142) */
},/* Start of coefficients for 19/20 interpolation */
{	-2.12436169e-03,	/*   H(  1) = -2.12436169e-03 =  H(159) */
	9.04073752e-03,	/*   H( 21) =  9.04073752e-03 =  H(139) */
	-3.13919224e-02,	/*   H( 41) = -3.13919224e-02 =  H(119) */
	7.89186135e-02,	/*   H( 61) =  7.89186135e-02 =  H( 99) */
	9.53945220e-01,	/*   H( 79) =  9.53945220e-01 =  H( 81) */
	-4.72677639e-03,	/*   H( 59) = -4.72677639e-03 =  H(101) */
	-7.47941062e-03,	/*   H( 39) = -7.47941062e-03 =  H(121) */
	3.39443982e-03,	/*   H( 19) =  3.39443982e-03 =  H(141) */
}};



/*
;Useage: I2xd16PI(signed int far *in, signed int far *out, unsigned int num_sp)
;	in	points to reference sample value (within an array of sample
;		values).  First output value = reference sample value
;	out	points to initial result storage location (within an array)
;	num_sp	is an integer euqal to half the number of output samples
;		generated
*/
void	I2xxPI(ftype *in, ftype *out, unsigned int num_sp)
{/** Note that in[-(lcstg1 -1)] ... in[num_sp+(lcstg1)] must be valid
	in order for valid results to be produced **/
 ftype		acum;
 signed int register	si;
 unsigned int	ic, oc;

 oc=0;
 for (ic=0; ic<num_sp; ic++)
 {out[oc]=in[ic]; oc++; acum=0.0;
  for (si=0; si<lcstg1; si++)
  {acum+= is1cI[lcstg1-si-1]* (in[ic-si] + in[ic+si+1]);}
  out[oc]=acum; oc++;
 }
}

ftype	I20xxPI(unsigned int frac, ftype *in)

{/**
; REV 13 May  1998	 Function to generate specified one of 19 intermediate
;			 points between specified sample point and the
;			 next sample point

;Useage: var=I20xxPI(unsigned int frac, far *in);
;	frac	is an integer [0, 19], indicating for which intermediate point
;		an interpolated value is desired (0 => no interpolation needed)
;	in	points to reference sample value (within an array of sample
;		values)

; NOTE that x20xcI is a doubly subscripted array containing coefficients
;		x20xcI[0 through xfstg2-2][0 through num_sptstg2-1]
  NOTE that in[-((num_sptstg2/2)-1)] and its num_sptstg2-1 successors
		must be valid in order for valid results to be generated
**/
 ftype		acum;
 signed int register	si, di;
 unsigned int	cf;

 if (frac == 0) return in[0];
 cf= frac-1;
 si= -(da_offstg2); acum=0.0;
 for (di=0; di<is2len; di++) {acum += in[si] * is2cI[cf][di]; si++;}
 return acum;
}

/* NOTE: In order to generate the 1st stage_2 interpolated point after the
	reference point the stage_1 result array (interpolated by a factor
	of 2) must contain {num_sptstg2/2 -1} entries prior to the
	reference point. (which is 2 for num_sptstg2=6) 

	Thus, for the num_sptstg2=6 case, the 1st entry in the stage_1 result
	array (TEMP) will be in[-1]; TEMP[0]=in[-1]; TEMP[2]=in[0] 

	In order to generate the stage_2 interpolated points after the last
	specified "in" value {num_sptstg2/2} entries are needed in the TEMP
	array following the entry for in[last].  TEMP[2*last]= in[last];
	TEMP[2*last+2]= in[last+1]; TEMP[2*last+4]= in[last+2].  Using the
	stage_1 (2X) interpolation routine to generate TEMP means that an 
	extra  entry will be generated for TEMP
*/
/**
;Useage: I10xd16PI(signed int far *in, signed int far *out, unsigned int num_sp)
;	in	points to reference sample value (within an array of sample
;		values).  First output value = reference sample value.
;	out	points to initial result storage location (within an array)
;	num_sp	is an integer indicating for how many sample points, 12
;		interpolated values are desired
;
;	2*(num_sp + num_sptstg2/2) WORDS are stored on stack from operation
;	of I2XD16PI.  These are used as input for DI10X8PI to generate the
;	10X (with respect to the original samples) values.
;	let TEMPx2 be the name of this array on the stack.
;	j=0
;	for (i=0; i<num_sp; i+=2)
;	{out[j]=TEMPx2[num_sptstg2/2 +i]; j+=1
;	 for (k=1; k<xfstg2; k++)
;	 {out[j]=DI05X8PI(k, &TEMPx2[num_sptstg2/2 +i]); j+=1}
;	 out[j]=TEMPx2[num_sptstg2/2 +i +1]; j+=1
;	 for (k=1; k<xfstg2; k++)
;	 {out[j]=DI05X8PI(k, &TEMPx2[num_sptstg2/2 +i +1]); j+=1}
;	}
**/
signed int	I40xxPI(ftype *in, ftype *out, unsigned int num_sp)
{/**	in	points to reference sample value (within an array of sample
;		values).  First output value = reference sample value.
;	out	points to initial result storage location (within an array)
;	num_sp	is an integer indicating for how many sample points, 40
;		interpolated values are desired
	returns -1 if space could not be allocated for 2_times array,
		0 otherwise
**/
 unsigned int	oc, si, frac, d1;
 ftype 	*x2; 	/* points to 2_times array */

 x2=malloc((size_t)(num_sp+pre_len+suf_len)*is1inf*(sizeof(ftype)));
 if (x2==NULL) return -1;
 I2xxPI(&in[-pre_len], x2, num_sp+pre_len+suf_len);
 oc=0;
 for (si=0; si<num_sp; si++)
 {d1=is1inf*(pre_len+si);
  out[oc]=x2[d1] ; oc++;
  for (frac=1; frac<is2inf; frac++) {out[oc]=I20xxPI(frac, &x2[d1]); oc++;}
  d1++;
  out[oc]=x2[d1]; oc++;
  for (frac=1; frac<is2inf; frac++) {out[oc]=I20xxPI(frac, &x2[d1]); oc++;}
 }
 free(x2);
 return 0;
}
