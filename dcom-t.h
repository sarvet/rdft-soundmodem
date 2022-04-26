/*
    definitions used by Wyman1x digital communications programs
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

/* Rev 02 Jun 2001      definitions to use interpolation filter on narrowband
	filtered data */
/* Rev 21 Jan 2001	filters for 230 Hz spaced sub-carriers at ~115 baud */
/* Rev 16 Dec 2000	narrow band filter for clock rate of 11,025 Hz introduced */
/* Rev 11 Oct 2000	wider filters, and clock rate of 11,025 Hz */
/* Rev 09 Oct 2000	frevt 23, Qt 12 introduced, for use with phase
	modulated version */

#include "dcom.h"


#if (frevt>=23)
#define	is1len		24	/*; number of sample points to use for each
				; stage 1 interpolated point*/
#define	is2len		8	/*; number of sample points to use for each
				; stage 2 interpolated point*/
#endif	/* frevt>=23 */


#define	lcstg1		(is1len/2)	/*; loop counter for stage 1 (X2)*/
#define	is1inf		2	/*; interpolation expansion factor, stage 1*/

#define	is2inf		20	/*; interpolation expansion factor, stage 2*/
#define	da_offstg2	((is2len/2)-1)	/*; in[da_offstg2] corresponds
				; to coefficient k20xcI[frac-1][0]*/

#define	ds1dcf		5
#define	ds2dcf		4
#define	ds3dcf		2


#define	lp1len		33
#define	lp1dly		((lp1len-1)/2)

#if (frevt>=22)

#if (Qt==12)
#define	ds1len		17
#define	ds2len		19
#define	ds3len		29
#endif	/* Qt==12 */

#endif	/* frevt>=22 */

#define	ds1dly		((ds1len-1)/2)
#define	ds2dly		((ds2len-1)/2)
#define	ds3dly		((ds3len-1)/2)

#define	DCF		(ds1dcf*ds2dcf*ds3dcf)
#define	INFx		(is1inf*is2inf)

#define	fap		920		/* number of final amplitude points */
#define	ipts		(fap/INFx)	/* # of points to generate
					interpolated values for */

#define pre_len	(is2len/(2*is1inf))
#define suf_len	(pre_len+1)


/* following is valid only when is1len/2 points prior to 1st desired
	amplitude point at stage 3 (of decimation) output are needed for
	subsequent interpolation (depends on is2len) */
#define	s3blen		(ipts + is1len + pre_len + suf_len)

#define	is1off		((is1len/2) + pre_len)	/* 98dec03  */
/* interpolation stage 1 offset + extra sample(s) needed for stage 1 interpolation */

#define	s2blen		((s3blen-1)*ds3dcf + ds3len)
#define	s1blen		((s2blen-1)*ds2dcf + ds2len)
#define	s0blen		((s1blen-1)*ds1dcf + ds1len)

#define	dtdly		((ds3dly*ds2dcf + ds2dly)*ds1dcf +ds1dly)
			/* total delay of 3 stages of decimation    */
#define	s0off		(dtdly + is1off*DCF)
			/* offset into initial array of reference point   */
#define	s1off		((s0off - ds1dly)/ds1dcf)	/* offset into
				stage 1 array of reference point */
#define	s2off		((s1off - ds2dly)/ds2dcf)	/* offset into
				stage 2 array of reference point */
#define	s3off		((s2off - ds3dly)/ds3dcf)	/* offset into
				stage 3 array of reference point */

#define	ds1nb2dcf		25
#define	ds2nb2dcf		4
#define	ds3nb2dcf		2
#define	DCFnb2			(ds1nb2dcf*ds2nb2dcf*ds3nb2dcf)
#define	ds1nb2len		59
#define	ds2nb2len		13
#define	ds3nb2len		23
#define	ds1nb2dly		((ds1nb2len-1)/2)
#define	ds2nb2dly		((ds2nb2len-1)/2)
#define	ds3nb2dly		((ds3nb2len-1)/2)

#define	fapnb2		880		/* number of final amplitude points */
#define	iptsnb2		(fapnb2/INFx)	/* # of points to generate
					interpolated values for */

/* following is valid only when is1len/2 points prior to 1st desired
	amplitude point at stage 3 (of decimation) output are needed for
	subsequent interpolation (depends on is2len) */
#define	s3nb2blen		(iptsnb2 + is1len + pre_len + suf_len)

#define	s2nb2blen		((s3nb2blen-1)*ds3nb2dcf + ds3nb2len)
#define	s1nb2blen		((s2nb2blen-1)*ds2nb2dcf + ds2nb2len)
#define	s0nb2blen		((s1nb2blen-1)*ds1nb2dcf + ds1nb2len)

#define	dtdlynb2	((ds3nb2dly*ds2nb2dcf + ds2nb2dly)*ds1nb2dcf +ds1nb2dly)
			/* total delay of 3 stages of decimation    */
#define	s0nb2off	(dtdlynb2 + is1off*DCFnb2)
			/* offset into initial array of reference point   */
#define	s1nb2off	((s0nb2off - ds1nb2dly)/ds1nb2dcf)	/* offset into
				stage 1 array of reference point */
#define	s2nb2off	((s1nb2off - ds2nb2dly)/ds2nb2dcf)	/* offset into
				stage 2 array of reference point */
#define	s3nb2off	((s2nb2off - ds3nb2dly)/ds3nb2dcf)	/* offset into
				stage 3 array of reference point */
