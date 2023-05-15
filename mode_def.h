/*
    initialized arrays used by Wyman1x digital communications programs
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
/* Rev 12 May  2001  version 1.0.0 for pm7b initial release */

unsigned int t1, k1, kx[num_modes] = { k05, k10, k20, k40 };
unsigned int iidx_mn1, iidx_mnx[num_modes] = { k05, k10, k20, k40 };
