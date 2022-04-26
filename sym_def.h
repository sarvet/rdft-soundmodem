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
/* Rev 25 Nov  2001	version 1.1.0, for Wyman1x initial public release */

#define	num_sc_i	4	/* number of information sub-carriers */
#define	num_sc_r	4	/* number if redundancy sub-carriers */
#define	num_sc_t	(num_sc_i + num_sc_r)
				/* total number of sub_carriers */
#define sub_sym_type	unsigned short int
#define	ch_sym_type	unsigned int

#include "code1_def.h"

#define	soh	0x100	/* symbol # for start of header */
#define	sod	0x101	/* symbol # for start of data */
#define	dummy1	0x102	/* symbol # for dummy data (filler) */
#define	app1	0x103	/* symbol # for app_symbol1 */
#define	app2	0x104	/* symbol # for app_symbol2 */

#define	phased2	72.2	/* square of divisor determining unit phase step */

/* NOTE: the #defines below specify where the phase change patterns
	for the indicated cases must be stored */

#define	no_mod	N1	/* symbol # for no modulation */
#define	frame1	(N1+1)	/* symbol # for frame 1 symbol */
#define	frame2	(N1+2)	/* symbol # for frame 2 symbol */
#define	frame3	(N1+3)	/* symbol # for frame 3 symbol */
#define	frame4	(N1+4)	/* symbol # for frame 3 symbol */
#define	frame1r	(N1+5)	/* symbol # for frame 1 redundancy symbol */
#define	frame2r	(N1+6)	/* symbol # for frame 2 redundancy symbol */
#define	frame3r	(N1+7)	/* symbol # for frame 3 redundancy symbol */
#define	frame4r	(N1+8)	/* symbol # for frame 3 redundancy symbol */

#define	id1x05	(N1+9)	/* symbol # for mode_id1, ~10% redundancy symbols */
#define	id2x05	(N1+10)	/* symbol # for mode_id2, ~10% redundancy symbols  */
#define	id1x10	(N1+11)	/* symbol # for mode_id1, ~20% redundancy symbols */
#define	id2x10	(N1+12)	/* symbol # for mode_id2, ~20% redundancy symbols */
#define	id1x20	(N1+13)	/* symbol # for mode_id1, ~40% redundancy symbols */
#define	id2x20	(N1+14)	/* symbol # for mode_id2, ~40% redundancy symbols */
#define	id1x40	(N1+15)	/* symbol # for mode_id1, ~70% redundancy symbols */
#define	id2x40	(N1+16)	/* symbol # for mode_id2, ~70% redundancy symbols */

#define lead_i	(N1+17)	/* symbol # for leader information symbols */
#define lead_r	(N1+18)	/* symbol # for leader redundancy symbols */
