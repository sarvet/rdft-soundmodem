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

extern	void	ds1t(ftype *in, ftype *out, unsigned int count);
extern	void	ds2t(ftype *in, ftype *out, unsigned int count);
extern	void	ds3t(ftype *in, ftype *out, unsigned int count);

extern	void	ds1nb2(ftype *in, ftype *out, unsigned int count);
extern	void	ds2nb2(ftype *in, ftype *out, unsigned int count);
extern	void	ds3nb2(ftype *in, ftype *out, unsigned int count);

signed int	I40xxPI(ftype *in, ftype *out, unsigned int num_sp);

extern	void	lpf01(ftype *in, ftype *out, unsigned int count);
