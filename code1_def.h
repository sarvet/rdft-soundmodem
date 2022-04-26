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

#define	N1	307	/* number of elements in GF */
#define	n1	(N1-1)	/* number of symbols per codeblock */
#define	a1	5	/* primative root used to generate code */
#define	leader_len	74	/* # of symbols in leader (0.60408 seconds) */
#define trailer_len	37	/* # of symbols in trailer (0.30204 seconds) */
#define	mode_len	3	/* # of pairs of symbols in mode specifier */
#define	frame_len	4	/* # of symbols in frame sequence */
#define	blks_per_frame	1	/* # of code blocks between frame sequences */
#define	num_modes	4	/* number of modes supported by this program */
#define num_ec_lvls	num_modes	/* number of error corrections levels */
#define	num_freqs	num_sc_t	/* number of sub-carrier frequencies */
#define	cbl1		(N1-1)	/* # of symbols per code block */
#define	iidx_mx1	(cbl1-1)	/* information index max */

#define	k05	32
#define	k10	(2*k05)
#define k20	(2*k10)
#define k40	(214)

#define	sym_type1	unsigned short int
#define sym_type2	unsigned int
#define	sym_typei	int
#define	sym_type_ch	unsigned int
sym_typei		ch_sym_buf_i[num_sc_i];
sym_typei		ch_sym_buf_r[num_sc_r];
#define	bytes_per_ch_sym	(sizeof(sym_type2))
#define	fnlmax		80	/* max file name length */
#define	flmax		65535	/* 2^16 -1 max file length */

#define	N2	9	/* number of elements in GF */
#define	t2	2	/* max number of correctable errors */
#define	a2	3	/* primative root used to generate code */
#define	fc2	3	/* characteristic of field */
#define	n2	(N2-1)	/* number of symbols per codeword */
