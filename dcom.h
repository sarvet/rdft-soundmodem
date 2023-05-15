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
/* rev 25 Nov 2001	version 1.1.0, for initial public release */

#define frevt 23
#define Qt    12
#define frevr 21
#define Qr    4

#define ftype float

#define pm_step 230.0
/* rev 26 Feb 2001	pmc0 changed from 590 to 575 */
/* rev 23 Jan 2001	pmc0 changed from 600 to 590 */
#define pmc0    575.0
#define pmc1    (pmc0 + pm_step)
#define pmc2    (pmc1 + pm_step)
#define pmc3    (pmc2 + pm_step)
#define pmc4    (pmc3 + pm_step)
#define pmc5    (pmc4 + pm_step)
#define pmc6    (pmc5 + pm_step)
#define pmc7    (pmc6 + pm_step)

/* below are clock period offsets to effect staggering of phase transitions */
#define soff0 0
#define soff1 11
#define soff2 22
#define soff3 34
#define soff4 45
#define soff5 56
#define soff6 67
#define soff7 79

/* rev 15 Nov 2001	for systems that don't define PI and PI2 (in math.h) */
#ifndef PI
  #define PI 3.14159265358979323846 /* pi */
#endif
#ifndef PI2
  #define PI2 (PI + PI)
#endif
