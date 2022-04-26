/*
    Reed-Solomon (8,4) code over GF(9), used by Wyman1x digital communications programs
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
/* Rev 25 Nov 2001	version 1.1.0, for Wyman1x initial public release */

#define	N2	9	/* number of elements in GF */
#define	t2	2	/* max number of correctable errors */
#define	a2	3	/* primative root used to generate code */
#define	fc2	3	/* characteristic of field */
#define	n2	(N2-1)	/* number of symbols per codeword */

/* arithemtic tables for Galois Field of 9 elements */
/* rev 14 Dec 2000  gf9mul fixed */

int	gf9add[9][9]={
{0, 1, 2, 3, 4, 5, 6, 7, 8},
{1, 2, 0, 4, 5, 3, 7, 8, 6},
{2, 0, 1, 5, 3, 4, 8, 6, 7},
{3, 4, 5, 6, 7, 8, 0, 1, 2},
{4, 5, 3, 7, 8, 6, 1, 2, 0},
{5, 3, 4, 8, 6, 7, 2, 0, 1},
{6, 7, 8, 0, 1, 2, 3, 4, 5},
{7, 8, 6, 1, 2, 0, 4, 5, 3},
{8, 6, 7, 2, 0, 1, 5, 3, 4}};

int	gf9mul[9][9]={
{0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 1, 2, 3, 4, 5, 6, 7, 8},
{0, 2, 1, 6, 8, 7, 3, 5, 4},
{0, 3, 6, 7, 1, 4, 5, 8, 2},
{0, 4, 8, 1, 5, 6, 2, 3, 7},
{0, 5, 7, 4, 6, 2, 8, 1, 3},
{0, 6, 3, 5, 2, 8, 7, 4, 1},
{0, 7, 5, 8, 3, 1, 4, 2, 6},
{0, 8, 4, 2, 7, 3, 1, 6, 5}};

int	gf9addi[9]={0, 2, 1, 6, 8, 7, 3, 5, 4};
int	gf9muli[9]={0, 1, 2, 4, 3, 7, 8, 5, 6};

	/* powers of a2 */
int	expa2[9]={1, 3, 7, 8, 2, 6, 5, 4, 1 };
	/* multiplicative inverses of field elements */
int	inve2[9]={0, 1, 2, 4, 3, 7, 8, 5, 6 };
	/* multiplicative inverses of element a2 */
int	inva2[9]={0, 4, 5, 6, 2, 8, 7, 3, 1 };

int	rgm2[4][4]={
	{7, 6, 2, 7 },
	{2, 2, 2, 1 },
	{7, 8, 1, 6 },
	{4, 5, 2, 5 }
};

int	encode2(int *extra, int *info)
{/* extra	points to array where redundancy symbols of code word will
			be stored,
    info	points to array where information symbols have been stored

    returns:	 0 if no errors were encountered
		-1 if an illegal symbol is found in the "info" array */
 register int	si, di;
 int		tmpi;

 for (si=0; si<(n2-2*t2); si++) if ((info[si]<0) || (info[si]>n2)) return -1;
 for (di=0; di<(2*t2); di++)
 {tmpi=0;
  for (si=0; si<(n2-2*t2); si++)
  {tmpi = gf9add[tmpi][gf9mul[info[si]][rgm2[di][si]]];}
  extra[di]= tmpi;
 }
 return 0;
}

int	decode2(int *codeword, int rho, int *erase)
{/* codeword	points to the codeword to be decoded.  This is  in the
		  form [r r r r i i i i], where the r's are the 
		  redundancy symbols, and the i's are the information symbols,
		  both with least significant symbols first.
    rho		is the number of erasures
    erase	points to an array of erasure indexes.  Each index is
		  in the range 0 through n-1 (i.e. 0 through 8)

   returns	a positive number indicating how many symbols were corrected
		0 if no symbols were in error
		-1 if any codeword symbols are illegal
		-2 if rho > 4 (too many erasures)
		-3 if any erasure indexes are out of bounds

		if return value >= 0, then *codeword contains the corrected
		  values.
*/
 register int	si, di;
 int		V[n2], S[n2], lam[n2], B[n2], T[n2], del[N2], C[n2], c[n2];
 int		sh[n2];
 int		r, L, tmpi, tmpj, tmpk, texp, ec;

 for (si=0; si<n2; si++) if ((codeword[si] < 0) || (codeword[si]> n2)) return -1;
 if (rho > (2*t2)) return -2;
 if (rho > 0)
 {for (si=0; si< rho; si++)
  if ((erase[si] < 0) || (erase[si] > (n2-1))) return -3;
 }
 L=0; r=0; B[0]=1; lam[0]=1;
 for (si=1; si<n2; si++) {B[si]=0; lam[si]=0;}
/* Transform received codeword */
 for (di=0; di<n2; di++)
 {tmpi=0;
  for (si=0; si<n2; si++)
  {tmpi = gf9add[tmpi][gf9mul[expa2[(si*di)%n2]][codeword[si]]];} 
  V[di]= tmpi;
  S[di]=V[di];
 }
/* note that while doing transforms S indexes range from 0 through n2-1, but
	when revising values of S the index rnages from 1 to n2.  This is
	resolved by using 0 through n2-1 while transforming, and using
	1 through n2 (mod n2) when revising values of S.
*/
 while(r<rho)
 {r++;
  {for (si=n2-1; si>0; si--) sh[si]= gf9mul[expa2[erase[r-1]]][lam[si-1]];
   sh[0]=0;
   for (si=0; si<n2; si++) {lam[si]= gf9add[lam[si]][gf9addi[sh[si]]]; B[si]= lam[si];}
   L++;
  } 
 } /* end of (r<rho) */
 
 while(r < 2*t2)
 {r++;
  /* calculate del[r] */
  tmpi=0;
  for (si=0; si<=(r<n2 ? r: n2-1); si++) tmpi= gf9add[tmpi][gf9mul[lam[si]][S[(r-si)%n2]]];
  del[r%n2]= tmpi;
  while(del[r%n2] < 0) del[r%n2]+=N2;
  if (del[r] == 0)
  {for (si=n2-1; si>0; si--) B[si]= B[si-1];
   B[0]=0;
  } /* end of del[r]==0 */
  else
  {for (si=1; si<n2; si++) sh[si]= gf9mul[del[r%n2]][B[si-1]];
   sh[0]=0;
   for (si=0; si<n2; si++) T[si]= gf9add[lam[si]][gf9addi[sh[si]]];
   if (2*L > r+rho-1)
   {for (si=0; si<n2; si++) lam[si]= T[si];
    for (si=n2-1; si>0; si--) B[si]= B[si-1];
    B[0]=0;
   } /* end of no change in L */
   else
   {L= r-L+rho;
    tmpi= inve2[del[r]];
    for (si=0; si<n2; si++) B[si]= gf9mul[tmpi][lam[si]];
    for (si=0; si<n2; si++) lam[si]= T[si];
   } /* end of increase in L */
  } /* end of del[r]!=0 */
 } /* end of (r <= 2*t2) */

 while(r<n2)
 {r++;
  /* calculate del[r] */
  tmpi=0;
  for (si=0; si<=(r<n2 ? r: n2-1); si++) tmpi= gf9add[tmpi][gf9mul[lam[si]][S[(r-si)%n2]]];
  del[r%n2]= tmpi;
  while(del[r%n2] < 0) del[r%n2]+=N2;
  S[r%n2]= gf9add[S[r%n2]][gf9addi[del[r%n2]]];
 } /* end of (r<=n) */

/* correct transformed codeword */
 for (si=0; si<n2; si++) C[si]= gf9add[V[si]][gf9addi[S[si]]];
/* inverse transform corrected codeword */
 tmpi= inve2[n2%fc2];	/* fc2 is characteristic of field */
			/* characteristic is # of elements in smallest sub-field */
 for (si=0; si<n2; si++)
 {tmpj=0;
  for (di=0; di<n2; di++)
  {texp= (si*di)%n2;	 
   if (texp == 0) tmpk=1; else tmpk= inva2[texp] ;
   tmpj = gf9add[tmpj][gf9mul[tmpk][C[di]]];
  } /* end of inner inverse transform loop  */
  c[si]= gf9mul[tmpi][tmpj];
  while(c[si] < 0) c[si] += N2;
 } /* end of inverse transform */
/* correct and count wrong symbols */
 ec=0;
 for (si=0; si<n2; si++)
 {if (c[si] != codeword[si])  {ec++; codeword[si]=c[si];}
 }
 return ec;
}
