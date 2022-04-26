/*
    demodulator-decoder program for Wyman1x digital communications method
    Copyright (C) 2001, 2002, 2003  Barry Sanderson

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
/* Rev 02 Jun 2003	Transfered changes to implement Frequency Change
	Compensation from development version of the program to this file.
	Those changes are summarized below:
			Increased search range for trailer to
	-250 Hz - +250 Hz.
			 Fixed calculation of frequency offset
	of trailer (exidxW[]).
			Changed to new model of linear drift for
	local oscillator.
			Limited lead_x[0].sidx to be >= 0.
			Fixed initialization condition for pv[lo].
			Sorted some erased symbols, based on distance
	from a legal value.
			Added saving of distances for symbols needing
	3 subcarrier values changed in order to match a legal value.
			Changed threshold for framing sequence cross-correlation
	good/bad discrimination.
			Added limits on standard deviation of
	demodulated values for acceptance of a good framing sequence.
			Added second stage of frequency offset drift
	compensation.
			Added compensation for frequency offset drift.
	Barry Sanderson 2003 Jun 02 */
/* Rev 29 May 2003	Removed: some commented out code; 'fnok'.
	Added: printf of copyright statement and warranty disclamer. 
	Changed logic on use of	default file name for decoded file.  Now the
	default file name is used if the first block fails to decode properly.
	Updated verstr[]. Made the initialization of 'ilc' independent of
	value of 'dbgl'.
	Barry Sanderson 2003 May 28 */
/* Rev 20 Mar 2003	Restored lsi_centerW back to 47 (from 46).
	Added MCWIL for maximum code word index limit, and set its
	value to 261 (previously 260 was used).  Added 2002, and 2003
	to copyright statement.  Updated verstr[].
	Barry Sanderson 2003 Mar. 21 */
/* Rev 31 Dec 2002	Fixed limit on index when searching for leader
	for the case of a single code block with 3 second 2 tone signal
	and chirp signal.  Changed calculations used to determine best
	framing sequence for the case of only two framing sequences
	(only one code block). Removed "Hamming window and" from message
	identifying program used.  Removed code for "Rectangular window
	and ..." message. Updated verstr[].
	Barry Sanderson 2002 Dec. 31 */
/* Rev 08 Dec 2002	Changed criterion for detection of end of file
	(for the file being transfered) for the case of unsuccessful
	decoding of the first block. Introduced ermsg2, to send error
	message to stderr without changing the value of ilc.  Introduced
	a version string in the print out of the program name.
	Barry Sanderson 2002 Dec. 08 */
/* Rev 25 Oct 2002	added saving of each correctly decoded block in a
	separate file, and the deletion of these files when all blocks are
	decoded correctly;  removed constraint on file name extension; 
	removed default mode as optional command line parameter; changed text
	for error message number 27; added text for error message 28;
	corrected condition for ERN=0 after nrmex:
	Barry Sanderson 2002 Oct. 25 */
/* Rev 08 Oct 2002	limited dmd_fidx to non-negative values
	Barry Sanderson 2002 Oct. 08 */
/* Rev 24 Sep 2002	changed value of gap_offi due to some subcarriers
	being modulated during leader and trailer; relaxed criteria for
	"good" subcarrier in lm_match3; relaxed criteria for success in
	lm_match3 routine; relaxed criteria for having detected leader and
	trailer; fixed array of "good" subcarrier values passed to "cbfitx"
	(previously this was wrong unless all subcariers were "good").
	Barry Sanderson 2002 Sep. 24 */
/* Rev 31 Aug 2002	added: detection of wrong number of bytes written
	to disk files; freeing of dmamp[][]; closing of irf2 for certain
	error exits; removing scratch files. Barry Sanderson 2002 Aug. 31 */
/* Rev 10 Dec 2001	changed time offset between demodulated, filtered,
	angle differences, and demodulated, average amplitudes; dbgl expanded 
	dblg = 35 introduced, for control of information dumped to stderr;
	changed	dumping of demodulated signals so that 40< dbgl <50 dumps
	demodulated angle differences, and dbgl >= 50 dumps demodulated
	amplitudes */
/* Rev 09 Dec 2001	added dump of results from adjusting symbol positions
	based on amplitude information to stderr */
/* Rev 05 Dec 2001	percentage of reserved erasures (era_rsv[]) doubled */
/* Rev 25 Nov 2001	cosmetic changes, functionallity not changed */
/* Rev 21 Nov 2001	inner decoding processed augmented to take
	care of some cases that actually have more than 2 errors, but
	are "close enough" with only 1 or 2 subcarriers erased */
/* Rev 20 Nov 2001	fixes to addition of 15 Nov 2001 */
/* Rev 15 Nov 2001	addition of compensation for rapid variations
	in propogation delay, using variations in amplitude, 
	format of dbgl > 5 changed */
/* Rev 18 Aug 2001	debuging of new inner decoding proceedure */
/* Rev 05 Aug 2001	new decoding proceedure for inner code */
/* Rev 16 Jun 2001	low pass filter used instead of averaging, on
	normalized angle differences */
/* Rev 12 Jun 2001	use of actual framing sequence positions, and best
	fit of these positions, for best subcarrier added between 1st
	and 2nd pass through demodulator */
/* Rev 05 Jun 2001	fine tuning bypassed, 2 demodulation passes forced */
/* Rev 02 Jun 2001	introduction of :
	  interpolation of narrowband filtered data (to get fine tuning slope
	  and interecpt) 
          2nd pass through demodulator
        Replacement of numerically unstable correlation coefficient calculation
	  in determining which fine tuing data to use
	use of only "good enough" fine tuning data */
/* Rev 01 Jun 2001	changed version 1.0.0 of pm7b for Wyman1x (pm8a) case */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <malloc.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <stddef.h>
#include <sys/stat.h>
#include <errno.h>

#include "dcom-t.h"
#include "exfun-t.h"
#include "sym_code.h"
#include "pm8a-frame.h"
#include "mode_def.h"
#include "pm8a-ref-spec2.h"
#include "gf307-inva.h"

#define	dbgl	30
/* dbgl value controls what information is printed to stdout, and stderr

 (dbgl > 0)  =>  final tuning, final cppblk, mode weights, and error counts
		 written
 (dbgl > 10)  =>  leader & trailer information and initial framing sequence
		 search information written
 (dbgl > 20)  =>  final framing sequence search information written
 (dbgl > 30)  =>  positions for all framing sequences written, and
		   adjustments (at framing sequence times) to compensate
		   for changes in frequency written
 (dbgl > 35) => symbol position adjustment information, based on
		 amplitude, written to stderr
 (dbgl > 40)  =>  gnuplot style files containing demodulation results
		 for each subcarrier written; 40 < dbgl < 50 dumps
		 filtered, normalized angle differences; dbgl >= 50
		 dumps average amplitudes
 (dbgl > 50)  =>  symbol numbers for both codes for each block written
*/
#define	in_type	signed short int
#define	dftsiz		4096
#define	dftsizW		(dftsiz/4)
#define	dmoff		(dftsiz/2)
#define	dmoffW		(dftsizW/2)
#define	wavhdrlenb	44L		/* .wav header length to skip */
#define	clk_freq_nom	11025.0
#define	dumd		1	/* dummy data value */
#define	delay_nom	90	/* clock periods per symbol */
#define	ave_nom		21	/* number of values to average together */

float		t_step; /* nominally 90.702947e-6 for 11025 Hz */

#define	NASN	0x200 	/* for 1 byte per symbol case */
#define	cppsym		delay_nom	/* clock periods per symbol, nominal */
#define	cppblk_nom	((n1+frame_len)*cppsym) /* clock periods per block, nominal */

#define num_modes_l     1       /* number of modes (different leaders) supported by this program */
#define num_files       10

/* for + and - 250 Hz around nominal frequencies, 1024 point DFT */
#define	tsi_start	23	/* 500 - 250 Hz */
#define	tsi_count	47	/* 506 Hz range */
/* for + and - 200 Hz around nominal frequencies */
#define	lsi_startW	27	/* 500 - 200 Hz */
#define lsi_countW	37	/* 400 Hz range */
#define	lsi_centerW	47	/* index of start of pattern, for no shift */


#define	l1_found	0x01
#define	t1_found	0x02

#define	frsq_srsym	60	/* # of symbol periods over which to search for framing sequence */
#define	frsq_srsam	(frsq_srsym*cppsym)
	/* # of sample periods over which to search for framing sequence */
#define	frosrng		42	/* one sided acceptance range for initial
	framing sequence */
#define	frosftr		34	/* one sided search range for fine tuning
	of framing sequence ranges */
#define	frosftr2	72	/* one sided search range for fine tuning
	of framing sequence ranges, version 2 */
#define	cppblk_ossr	26	/* one sided search range for fine tuning
	of cppblk */

#define	gap_offi	8	/* number of 4096 DFT points from center to
 				measurement point between sub-carriers */
#define	fdthld_nom	0.01	/* threshold for sub-carrier spacing ratio */

#define	xcthld		0.6	/* minimum cross-correlation value to accept */
#define nsc_min		5	/* minimum number of subcarriers required to be present */

#define	ossr2		6	/* one sided search range around calculated
				positions of local maximums, from cross-correlation */

#define	mode1dthld	13.000		/* distance threshold for mode 1, no erasures */
#define	mode2dthld	13.60147	/* distance threshold for mode 2, no erasures */


#define	slope_max	1.01	/* 10,000 parts per million */
#define slope_min	0.99
#define	off_hz_max	400.
#define	off_hz_min	-400.
#define	cppblk_max	28179	/* 10,000 parts per million */
#define	cppblk_min	27621
#define	cppblk_off_lim	278
#define	max_usr_sym	0xff

#define	mseavexf	1.0	/* factor to multiply average mean squared error
	by, to get threshold for good/bad decision on individual framing sequences */

#define	ossi_amp1	18	/* one sided search range for initial symbol
				position refinement search */
#define	thld_amp1	0.6	/* threshold value for cross-correlation value
				in initial symbol position refinement */
#define	ossi_amp2	18	/* one sided search range for 2nd pass symbol
				position refinement search */
#define	thld_amp2	0.5	/* threshold value for cross-correlation value
				in 2nd pass symbol position refinement */
#define	nfmin		2	/* number of framing sequences where a
				different equation is used in finding best */

#define	MCWIL		261	/* code word index limit to pass to
				   inner decoding routines */
#define	frxcthld	0.96	/* framing sequence cross-correlation threshold
				value, for "good"/"bad" determination */
#define	pm8a_frame_ave	1.993029	/* average value of pm8a_framing
					sequence samples (phase change units) */
#define pm8a_frame_sd_min	0.97	/* minimum standard deviation for
				"good" framing sequence */
#define pm8a_frame_sd_max	1.51	/* maximum standard deviation for
				"good" framing sequence */

/*** external functions associated with inner decoding process ***/
extern int	decode2ed1mode(ftype *data, ftype *dist, ftype dthld);
extern int	decode2ed2mode(ftype *data, ftype *dist, ftype dthld);
extern void	init_code2(void);

/* start of added subroutines, 29jul2001 */
extern int	decode2wne2(ftype *data, int *scc, unsigned short int *ob_cnt,
	unsigned short int *ob_mask);
extern int	decode2ed2w1e(ftype *data, ftype *dist, int mcwi, int *era, int *scc,
	unsigned short int ob_mask, unsigned short int *bec_mask);
extern int	decode2ed2w2e(ftype *data, ftype *dist, int mcwi, int *era, int *scc,
	unsigned short int ob_mask, unsigned short int *bec_mask);
extern int	decode2ed2w3e(ftype *data, ftype *dist, int mcwi, int *era, int *scc,
	unsigned short int ob_mask, unsigned short int *bec_mask);
extern int	decode2ed2w4e(ftype *data, ftype *dist, int mcwi, int *era, int *scc,
	unsigned short int ob_mask, unsigned short int *bec_mask);
/* end of added subroutines, 29jul2001 */

/*** external functions associated with DFT ***/
extern double	init_dft3(unsigned int dft_siz, float t_step, float *snqt);
extern void	dfft3f(float *inr, float *ini, float *outr, float *outi,
		float  *sint, unsigned int N, unsigned int inv,
		unsigned int real_in);


extern	int	errno;

int	soff[num_sc_t]={soff0, soff1, soff2, soff3,
			soff4, soff5, soff6, soff7};

unsigned short int	eraii[num_sc_t];

#if (s0blen > s0nb2blen)
in_type		samsi[s0blen];
#else
in_type		samsi[s0nb2blen];
#endif
in_type		datarI[dftsiz];
float 		inr[dftsiz], ini[dftsiz];
float		outr[dftsiz], outi[dftsiz];
float 		mag[3][dftsiz],  wind[dftsiz], windW[dftsizW];
float		snqt[(dftsiz/4)+1];
float		snqtW[(dftsizW/4)+1];

int		val_read, val_to_write;
unsigned long	sidx, hdr_lenb, sstepW, xctstart;
int		eoff, win_flg;
FILE		*inf, *irf, *irf2;
unsigned int	rsmode, bbs;

ch_sym_type	ch_sym_buf[cbl1];	/* output symbol buffer */
unsigned int	cb[2][cbl1];	/* code block buffers */

float	ct_slope, ct_intercept, f_slope, f_intercept;

/* following is array of phase modulation center frequencies */
double pm_nom_ctr[num_sc_t]={pmc0, pmc1, pmc2, pmc3, pmc4, pmc5, pmc6, pmc7};

#define	num_sc_l	12
/* following is array of nominal frequencies for leader and trailer */
double pm_nom_lt[num_sc_l]={pmc0, pmc1,
	1006.167, 1128.684, 1236.168, 1358.666,
	1401.326, 1523.835, 1631.327, 1753.832,
	pmc6, pmc7};
 

#if (fap > fapnb2)
ftype		Ssbc[fap], Ssbs[fap], Amp[fap], buf2[fap];
#else
ftype		Ssbc[fapnb2], Ssbs[fapnb2], Amp[fapnb2], buf2[fapnb2];
#endif

#if (s0blen > s0nb2blen)
ftype  		Sambuf[s0blen], Lo1c[s0blen], Lo1s[s0blen];
#else
ftype  		Sambuf[s0nb2blen], Lo1c[s0nb2blen], Lo1s[s0nb2blen];
#endif
#if (s1blen > s1nb2blen)
ftype  		Lo1cf[s1blen], Lo1sf[s1blen];
ftype		Lo2cc[s1blen], Lo2cs[s1blen], Lo2ss[s1blen], Lo2sc[s1blen];
#else
ftype  		Lo1cf[s1nb2blen], Lo1sf[s1nb2blen];
ftype		Lo2cc[s1nb2blen], Lo2cs[s1nb2blen], Lo2ss[s1nb2blen], Lo2sc[s1nb2blen];
#endif
#if (s2blen > s2nb2blen)
ftype		Lo2cs1[s2blen], Lo2ss1[s2blen];
#else
ftype		Lo2cs1[s2nb2blen], Lo2ss1[s2nb2blen];
#endif
#if (s3blen > s3nb2blen)
ftype		Lo2cs2[s3blen], Lo2ss2[s3blen], Ampb[s3blen];
#else
ftype		Lo2cs2[s3nb2blen], Lo2ss2[s3nb2blen], Ampb[s3nb2blen];
#endif

/* following are file names used for demodulated output values for each subcarrier */
unsigned char	sb0db[]="lo0dif23n.gnu";
unsigned char	sb1db[]="lo1dif23n.gnu";
unsigned char	sb2db[]="lo2dif23n.gnu";
unsigned char	sb3db[]="lo3dif23n.gnu";

unsigned char	sb4db[]="lo4dif23n.gnu";
unsigned char	sb5db[]="lo5dif23n.gnu";
unsigned char	sb6db[]="lo6dif23n.gnu";
unsigned char	sb7db[]="lo7dif23n.gnu";

unsigned char	sb8db[]="lo8dif23n.gnu";
unsigned char	sb9db[]="lo9dif23n.gnu";
unsigned char	*lodif7nfn[num_files]={sb0db, sb1db, sb2db, sb3db, sb4db,
		 sb5db, sb6db, sb7db, sb8db, sb9db};

ftype		*dif7n[num_sc_t], *dmamp[num_sc_t];
ftype		ec_mode_dist[num_ec_lvls];

/*** NOTE the two arrays below MUST be set to agree with the
	symbol values actually used for these modes
***/
sym_typei	ec_mode_id1xxsym[num_ec_lvls]={12,  6, 10,  3};
sym_typei	ec_mode_id2xxsym[num_ec_lvls]={ 6,  5,  9, 10};

int		ec_mode;
ftype		ref_a;

int		exidx[num_sc_l], exidxW[num_sc_l];

signed int	frsqbioff[num_sc_t];
float		frsqbmse[num_sc_t];

double		f_stepd, f_stepdW;

struct	lxcd2{
	unsigned long 	sidx;
	float 		r53;
	float 		off_hz;
	float 		slope;
	unsigned int	state;
	unsigned int	mode;
	signed int	offi;
	unsigned int	nlm;
	float		xcv; 
	unsigned int	lmbm;
	float		lmfreq[num_sc_l];
	float		lmqual[num_sc_l];
	unsigned int	lmidx[num_sc_l];
	unsigned int	lmgoodi[num_sc_l];
	signed int	lmcount; };

signed int	midx[3];
struct lxcd2	lead_x[3];


int do_dft(long idx, int win_flg)
{/* read dftsiz values from inf, starting at idx; apply proper
	window function; take DFT;
  Return:
	-1 if seek failed
	 1 if end of file was reached, while trying to read
	 0 otherwise
 */

 int	val_read, si, eoff;
 
 eoff=0;
 if (fseek(inf, hdr_lenb+ (idx>=0 ? idx*sizeof(in_type): 0L), SEEK_SET) == 0)
 {val_read=fread(&datarI[0], sizeof(in_type), dftsiz, inf);
  if (val_read != dftsiz)
  {eoff=1; for (si=val_read; si<dftsiz; si++) datarI[si]=0.0;}
  if (win_flg==1)
   for (si=0; si<dftsiz; si++) inr[si]= (float)(datarI[si]*wind[si]);
  else
   for (si=0; si<dftsiz; si++) inr[si]= (float)(datarI[si]);
  dfft3f(inr, ini, outr, outi, snqt, dftsiz-1, 0, 1);
  return(eoff);
 }	/* end of successful seek */
 else  return(-1);   /* seek failed */
}

int do_dftW(long idx, int win_flg)
{/* read dftsizW values from inf, starting at idx; apply proper
	window function; take DFT;
  Return:
	-1 if seek failed
	 1 if end of file was reached, while trying to read
	 0 otherwise
 */

 int	val_read, si, eoff;
 
 eoff=0;
 if (fseek(inf, hdr_lenb+ idx*sizeof(in_type), SEEK_SET) == 0)
 {val_read=fread(&datarI[0], sizeof(in_type), dftsizW, inf);
  if (val_read != dftsizW)
  {eoff=1; for (si=val_read; si<dftsizW; si++) datarI[si]=0.0;}
  if (win_flg==1)
   for (si=0; si<dftsizW; si++) inr[si]= (float)(datarI[si]*windW[si]);
  else
   for (si=0; si<dftsizW; si++) inr[si]= (float)(datarI[si]);
  dfft3f(inr, ini, outr, outi, snqtW, dftsizW-1, 0, 1);
  return(eoff);
 }	/* end of successful seek */
 else  return(-1);   /* seek failed */
}

/* function to calculate linear (magnitude) of complex sequence */

void	set_linmag(float *inc, float *ins, float *mag, float *arg, unsigned int count)
{double	dumy1, dummy2, dummy3;
 unsigned int i;

 for (i=0; i<count; i++)
 {dumy1=(double)(*(inc+i));  dummy2=(double)(*(ins+i));
  dummy3=sqrt(dumy1*dumy1 + dummy2*dummy2);
  dummy3= (dummy3>1e-10?dummy3:1e-3);
  *(mag+i)=(float)(dummy3);
 }
}


int	decode1(int *codeword, int rho, int *erase)
{/* codeword	points to the codeword to be decoded.  This is  in the
		  form [r r r r ... i i i i i i], where the r's are the 
		  redundancy symbols, and the i's are the information symbols,
		  both with least significant symbols first.
    rho		is the number of erasures
    erase	points to an array of erasure indexes.  Each index is
		  in the range 0 through n1-1

   returns	a positive number indicating how many symbols were corrected
		0 if no symbols were in error
		-1 if any codeword symbols are illegal
		-2 if rho > 2*t (too many erasures)
		-3 if any erasure indexes are out of bounds

		if return value >= 0, then *codeword contains the corrected
		  values.
*/
 register int	si, di;
 int		V[n1], S[n1], lam[n1], B[n1], T[n1], del[N1], C[n1], c[n1];
 int		sh[n1];
 int		r, L, tmpi, tmpj, tmpk, texp, ec;

 for (si=0; si<n1; si++) if ((codeword[si] < 0) || (codeword[si]> n1)) return -1;
 if (rho > (2*t1)) return -2;
 if (rho > 0)
 {for (si=0; si< rho; si++)
  if ((erase[si] < 0) || (erase[si] > (n1-1))) return -3;
 }
 L=0; r=0; B[0]=1; lam[0]=1;
 for (si=1; si<n1; si++) {B[si]=0; lam[si]=0;}
/* Transform received codeword */
 for (di=0; di<n1; di++)
 {tmpi=0;
  for (si=0; si<n1; si++) tmpi += expa1[(si*di)%n1] * codeword[si];
  V[di]= tmpi%N1;
  S[di]=V[di];
 }
/* note that while doing transforms S indexes range from 0 through n1-1, but
	when revising values of S the index ranges from 1 to n1.  This is
	resolved by using 0 through n1-1 while transforming, and using
	1 through n1 (mod n1) when revising values of S.
*/
 while(r<rho)
 {r++;
  {for (si=n1-1; si>0; si--) sh[si]= (expa1[erase[r-1]]*lam[si-1])%N1;
   sh[0]=0;
   for (si=0; si<n1; si++) {lam[si]= (lam[si] - sh[si])%N1; B[si]= lam[si];}
   L++;
  } 
 } /* end of (r<rho) */
 
 while(r < 2*t1)
 {r++;
  /* calculate del[r] */
  tmpi=0;
  for (si=0; si<=(r<n1 ? r: n1-1); si++) tmpi+= lam[si]*S[(r-si)%n1];
  del[r%n1]= tmpi%N1;
  while(del[r%n1] < 0) del[r%n1]+=N1;
  if (del[r] == 0)
  {for (si=n1-1; si>0; si--) B[si]= B[si-1];
   B[0]=0;
  } /* end of del[r]==0 */
  else
  {for (si=1; si<n1; si++) sh[si]= (del[r%n1]*B[si-1])%N1;
   sh[0]=0;
   for (si=0; si<n1; si++) T[si]= (lam[si]-sh[si])%N1;
   if (2*L > r+rho-1)
   {for (si=0; si<n1; si++) lam[si]= T[si];
    for (si=n1-1; si>0; si--) B[si]= B[si-1];
    B[0]=0;
   } /* end of no change in L */
   else
   {L= r-L+rho;
    tmpi= inve1[del[r]];
    for (si=0; si<n1; si++) B[si]= (tmpi*lam[si])%N1;
    for (si=0; si<n1; si++) lam[si]= T[si];
   } /* end of increase in L */
  } /* end of del[r]!=0 */
 } /* end of (r <= 2*t1) */

 while(r<n1)
 {r++;
  /* calculate del[r] */
  tmpi=0;
  for (si=0; si<=(r<n1 ? r: n1-1); si++) tmpi+= lam[si]*S[(r-si)%n1];
  del[r%n1]= tmpi%N1;
  while(del[r%n1] < 0) del[r%n1]+=N1;
  S[r%n1]= (S[r%n1]-del[r%n1])%N1;
 } /* end of (r<=n1) */

/* correct transformed codeword */
 for (si=0; si<n1; si++) C[si]= (V[si]-S[si])%N1;
/* inverse transform corrected codeword */
 tmpi= inve1[n1];
 for (si=0; si<n1; si++)
 {tmpj=0;
  for (di=0; di<n1; di++)
  {/* tmpk = a^(-si*di) */
   texp= (si*di)%n1;
   if (texp == 0) tmpk=1; else tmpk= inva1[texp] ;
   tmpj+= (tmpk*C[di])%N1;
  } /* end of inner inverse transform loop  */
  c[si]= (tmpi*tmpj)%N1;
  while(c[si] < 0) c[si] += N1;
 } /* end of inverse transform */
/* correct and count wrong symbols */
 ec=0;
 for (si=0; si<n1; si++)
 {if (c[si] != codeword[si])  {ec++; codeword[si]=c[si];}
 }
 return ec;
}

int	lm_match3(float *lmfreq, int *lmidx, int count, int *lmgoodi, double fdthld)
{/* check for spacing of lmfreq[] compared to ideal values
	return number of "good" cases found, or -1 if an error is detected */
 float	md, nd;
 int	si, di, max_cnt, init_flg, gcnt, lc[num_sc_l];

 init_flg=1; max_cnt=0;
 for (si=0; si<count; si++) lc[si]=0;
 for (di=0; di<count; di++)
 {for (si=0; si<count; si++)
  {if (si==di) continue;
   md= lmfreq[di] - lmfreq[si];
   nd= pm_nom_lt[lmidx[di]] - pm_nom_lt[lmidx[si]];
   if (fabs((double)(md/nd -1.0)) < fdthld)
   {lc[di]++;
    if (init_flg==1) {init_flg=0; max_cnt=1;}
    else
    {if (lc[di]>max_cnt) max_cnt= lc[di];}
   } /* end of match found */
  } /* end of for si */
 } /* end of for di */
 gcnt=0;
 for (si=0; si<count; si++)
 {if (lc[si]>= (count-3)) lmgoodi[gcnt++]= lmidx[si];
 }
 if(gcnt>= nsc_min) return(gcnt);
 else return(-1);	/* error case */
}

int	find_local_max3(float *mag, float *lmfreq, float *lmqual,
	 int *lmidx, unsigned int *bm)
{/* find local max values near indexes listed in exidx[],
    return number of local max values found */

 int	count, max_mask, lm_map;
 int	si, di, bi, ni, ai;
 float	acumn, acumd, blm, tmpa, dummyf;

 max_mask=0x0001; count=0; lm_map=0;
 for (si=0; si< num_sc_l; si++)
 {ni=exidx[si];
  for (di= -ossr2; di<= ossr2; di++)
  {if (di== -ossr2) {blm= mag[ni+di]; bi=di+ni;}
   if ((mag[ni+di-2] < mag[ni+di-1]) &&
       (mag[ni+di-1] < mag[ni+di])   &&
       (mag[ni+di]   > mag[ni+di+1]) &&
       (mag[ni+di+1] > mag[ni+di+2]))
   {
    if (mag[ni+di] > blm) {blm= mag[ni+di]; bi=di+ni;}
    lm_map |= max_mask;
   } /* end of local max found */
  } /* end of for di */
  if ((lm_map & max_mask) == max_mask)
  {lmidx[count]= si;
   /* calculate weighted average for this local max */
   acumn=0.; acumd=0.;
   for (ai=-1; ai<2; ai++)
   {tmpa= mag[bi+ai]; dummyf= (float)(bi+ai)*(float)(f_stepd);
    acumn += tmpa*dummyf;
    acumd += tmpa;
   }
   lmfreq[count]= acumn/acumd;
   /* calculate quality number for this local max */
   acumd=0.;
   for (ai=-1; ai<2; ai++)
   {acumd += mag[bi+ gap_offi+ai] + mag[bi- gap_offi+ai];
   }
   acumd= (acumd > 1e-3 ? acumd : 1e-3);
   lmqual[count]= mag[bi]*6.0/acumd;
   count++;
  }/* end of bit set for local max found */
  max_mask= max_mask << 1;
 } /* end of for si */
 bm[0]= lm_map;
 return(count);
}

int	find_local_max3W(float *mag, float *lmfreq, float *lmqual,
	 int *lmidx, unsigned int *bm)
{/* find local max values near indexes listed in exidxW[],
    return number of local max values found */

 int	count, max_mask, lm_map;
 int	si, di, bi, ni, ai;
 float	acumn, acumd, blm, tmpa, dummyf;

 max_mask=0x0001; count=0; lm_map=0;
 for (si=0; si< num_sc_l; si++)
 {ni=exidxW[si];
  for (di= -2; di<= 2; di++)
  {if (di== -2) {blm= mag[ni+di]; bi=di+ni;}
   if ((mag[ni+di-1] < mag[ni+di])   &&
       (mag[ni+di]   > mag[ni+di+1]))
   {
    if (mag[ni+di] > blm) {blm= mag[ni+di]; bi=di+ni;}
    lm_map |= max_mask;
   } /* end of local max found */
  } /* end of for di */
  if ((lm_map & max_mask) == max_mask)
  {lmidx[count]= si;
   /* calculate weighted average for this local max */
   acumn=0.; acumd=0.;
   for (ai=-1; ai<2; ai++)
   {tmpa= mag[bi+ai]; dummyf= (float)(bi+ai)*(float)(f_stepdW);
    acumn += tmpa*dummyf;
    acumd += tmpa;
   }
   lmfreq[count]= acumn/acumd;
   /* calculate quality number for this local max */
   acumd=0.;
   for (ai=-1; ai<2; ai++)
   {acumd += mag[bi+ gap_offi/4 +ai] + mag[bi- gap_offi/4 +ai];
   }
   acumd= (acumd > 1e-3 ? acumd : 1e-3);
   lmqual[count]= mag[bi]*6.0/acumd;
   count++;
  }/* end of bit set for local max found */
  max_mask= max_mask << 1;
 } /* end of for si */
 bm[0]= lm_map;
 return(count);
}


float	xcor4b(float *pattern, float *data, unsigned int patlen,
		unsigned int count, unsigned int *countidx)
/* Floating point correlation calculation
   pattern	points to 1st element of pattern array
   data		points to 1st element of data array
   patlen	is length of pattern
   count	is number of offsets into data array at which to calculate
  		correlation coefficient
   countidx	points to where "count" value for greatest correlation coefficient
  		value is to be stored

   returns 	correlation coefficient between pattern[] and data[countidx]
		(which is greatest value found) */
{
 float		sump, sump2, sumd, sumd2, sumpd, n, dp, dd;
 float		patd, tmpd, ccmax;
 unsigned int register	si, di, cm;

 n=(float)(patlen);
 sump=0.0; sump2=0.0;
 for (si=0; si<patlen; si++) {dp=pattern[si]; sump+=dp; sump2+= dp*dp;}
 patd= n*sump2 - (sump*sump);
 for (di=0; di<count; di++)
 {sumd=0.0; sumd2=0.0; sumpd=0.0;
  for (si=0; si<patlen; si++)
  {dp=pattern[si]; dd=data[si+di];
   sumd+=dd; sumd2+= dd*dd; sumpd+= dd*dp;
  } /* end of for si */
  tmpd=(float)((double)(n*sumpd - sump*sumd)/(sqrt((double)(patd*(n*sumd2 - (sumd*sumd))))));
  if (di==0) {cm=di; ccmax=tmpd;}
  else
  {if (tmpd > ccmax) {cm=di; ccmax=tmpd;}
  }
 } /* end of for di */
 *countidx= cm;
 return ccmax;
}

float	xcor4c(float *pattern, float *data, unsigned int patlen,
		unsigned int count, unsigned int *countidx,
		float *ave, float *sd)
/* Floating point correlation calculation
   pattern	points to 1st element of pattern array
   data		points to 1st element of data array
   patlen	is length of pattern
   count	is number of offsets into data array at which to calculate
  		correlation coefficient
   countidx	points to where "count" value for greatest correlation coefficient
  		value is to be stored
   ave		points to where average of patlen values, starting at
		data[countidx] is stored
   sd		points to where standard deviation of patlen values, starting at
		data[countidx] is stored

   returns 	correlation coefficient between pattern[] and data[countidx]
		(which is greatest value found) */
{
 float		sump, sump2, sumd, sumd2, sumpd, n, dp, dd;
 float		patd, tmpd, ccmax;
 unsigned int register	si, di, cm;

 n=(float)(patlen);
 sump=0.0; sump2=0.0;
 for (si=0; si<patlen; si++) {dp=pattern[si]; sump+=dp; sump2+= dp*dp;}
 patd= n*sump2 - (sump*sump);
 for (di=0; di<count; di++)
 {sumd=0.0; sumd2=0.0; sumpd=0.0;
  for (si=0; si<patlen; si++)
  {dp=pattern[si]; dd=data[si+di];
   sumd+=dd; sumd2+= dd*dd; sumpd+= dd*dp;
  } /* end of for si */
  tmpd=(float)((double)(n*sumpd - sump*sumd)/(sqrt((double)(patd*(n*sumd2 - (sumd*sumd))))));
  if (di==0) {cm=di; ccmax=tmpd;}
  else
  {if (tmpd > ccmax) {cm=di; ccmax=tmpd;}
  }
 } /* end of for di */
 sumd= 0.0;
 for (si=cm; si<(cm+patlen); si++) sumd += data[si];
 sumd= sumd/n;
 sumd2=0.0;
 for (si=cm; si<(cm+patlen); si++)
 {dd= data[si]- sumd;
  sumd2= sumd2 + dd*dd;
 }
 *ave= sumd;
 *sd= (float)(sqrt((double)(sumd2/(n-1.0))));
 *countidx= cm;
 return ccmax;
}

int	cbfitx(float *lmfreq, int *lmidx, int count, float *slope, float *intercept)
{/* calculate least squares best fit between contents of lmfreq[] and
	the count relevant entries in pm_nom_lt[]
    return 0 on success, -1 on error */

 int	di;
 ftype	sumxy, sumx, sumy, x2sum, cx;

/* least squares fit to coarse tuning errors */
 sumx=0.0; sumy=0.0; sumxy=0.0; x2sum=0.0;
 for (di=0; di<count; di++)
 {cx= pm_nom_lt[lmidx[di]];
  sumx += cx; x2sum += cx*cx;
  sumxy += cx*lmfreq[di];
  sumy += lmfreq[di];
 }
 *slope = ((ftype)(count)*sumxy - sumx*sumy)/(count*x2sum - sumx*sumx);
 *intercept= sumy/(ftype)(count) - *slope*sumx/(ftype)(count);
 return(0);
}

int	find_local_max_amp(float *data, int ossi, int *lm_off, float *lm_value)
/* Search no more than ossi locations away from the starting location
   for a local maximum in the cross -correlation between the data[] values
   and the am0_7ref[] pattern.
   Return 1 and set lm_off[0] and lm_value[] if local max found,
   otherwise return 0 */
{float	cval, pval;
 int	idx, istep, sii, cmx;

 pval= xcor4b(&am_ref[0], &data[0], am_ref_len, 1, &cmx);
 cval= xcor4b(&am_ref[0], &data[1], am_ref_len, 1, &cmx);
 if (pval==cval) return(0); /* fail on adjacent equal values */
 if (cval > pval)
 {istep=1; idx=1;} else {istep=-1; idx=0;}
 for (sii=0; sii<ossi; sii++)
 {pval=cval; idx += istep;
  cval= xcor4b(&am_ref[0], &data[idx], am_ref_len, 1, &cmx);
  if (cval==pval) return(0);
  if (cval < pval)
  {lm_off[0]= idx-istep; lm_value[0]=pval;
   return(1);
  }
 }/* end of for sii */
 return(0);
}

int	write_code_block(unsigned char *data, int dcnt, unsigned char *fns,
	 int fncnt, int cbi)
{/* data	points to buffer of data to be written
    dcnt	is number of unsigned characters to be written
    fns		points to buffer containing file name base string
    fncnt	is number of characters to preserve in file name
    cbi		is index of code block to be written

    returns	-1 if file can't be opened
		-2 if error in writing file happens
		 0 on success
 */
 FILE	*cbfp;
 int	fwr;

 sprintf(&fns[fncnt], "%03u", cbi);
 cbfp=fopen(&fns[0], "wb");
 if (cbfp==NULL) return -1;
 fwr= fwrite(data, sizeof(unsigned char), dcnt, cbfp);
 fclose(cbfp);
 if (fwr<0) return -2; else return 0;
}

 char	err00[]="Successful termination";
 char	ins[]="Usage:\n  Wyman1x-demod-decode datafile [sample_freq(Hz)]";
 char	err02[]="malloc could not allocate enough memory.";
 char	err03[]="Unknown error";
 char	err04[]="fseek failed";
 char	err05[]="WRONG number of bytes written to disk file";
 char	err06[]="stat failed";
 char	err07[]="Couldn't open file for input";
 char	err08[]="Couldn't open file for output";
 char	err09[]="";
 char	err10[]="The leader was not found in the first half of the input file";
 char	err11[]="The trailer was not found in the last half of the input file";
 char	err12[]="leader too corrupted for fine tuning operations";
 char	err13[]="error trying to read from scratch file";
 char	err14[]="logic error in step 3 of symbol position refinement";
 char	err15[]="logic error in step 4 of symbol position refinement";
 char	err16[]="";
 char	err17[]="";
 char	err18[]="";
 char	err19[]="Start of Header symbol not in proper place";
 char	err20[]="recovered file name length is too long";
 char	err21[]="recovered file name length is too short";
 char	err22[]="recovered file length is too long";
 char	err23[]="illegal symbol in file name";
 char	err24[]="Start of Data symbol not in proper place";
 char	err25[]="illegal symbol in body of file";
 char	err26[]="error while writing output file";
 char	err27[]="last codeblock(s) ignored:\n\tpossibly due to incorrectly recovered file length\n\tor extraneous trailer or header near end of file";
 char   err28[]="error trying to save decoded block in separate file";
 char   err29[]="";
 char	err30[]="";
 char	err31[]="Sample Frequency is out of range\n";
 char	err32[]="Digital filter could not allocate memory\n";
 char	err33[]="";
 char		*erexmsg[]={err00, ins, err02, err03, err04, err05, err06,
	err07, err08, err09, err10, err11, err12, err13, err14, err15,
	err16, err17, err18, err19, err20, err21, err22, err23, err24,
	err25, err26, err27, err28, err29, err30, err31, err32, err33,
			 };

unsigned char	ilc;

void	ermsg(int ern)
{fprintf(stderr,"%s\n", erexmsg[ern]);
 ilc='*';
}

void	ermsg2(int ern)
{fprintf(stderr,"%s\n", erexmsg[ern]);
}

int main(int argc, char *argv[], char *env[])

{
 int		ERN,  err, j, ki, eci, wcbr;
 int register	si, di;

 struct stat	statb, statir;
 int		statr, statirr;
 FILE		*otf;
 char		*ofn;
 unsigned char	lss;
 
 double		t_step1st, dummyd1;
 float		r53[3];

 int		mi,  lsi, lo, fi;
 float		cx, sumx, x2sum, sumxy, sumy, sumx2;
 double		pm_ft_ctr[num_sc_t];
 float		f_slope, f_intercept;
 double		clk_period, ca, astep, ang_ref, ang_off;
 double		bias[num_sc_l], pv[num_sc_l];
 double		astepd, astepdot, drift_rate_hps;

 int		end_of_file, ir;
 int		val_to_read, dum_len, val_to_write;
 size_t		val_read, val_written;
 unsigned long	flen, fbidx, numel, flen_dmd, numel_dmd;
 unsigned char	scratch[]="pm8-scratch";
 unsigned char	scratch2[]="pm8-scratch2";
 unsigned char	default_out_name[]="Wyman1x-decode-default.out";
 ftype		*fbufin, *fbufot;
 int		delay, ns, dmd_cnt;
 double		NS, acum0d;
 ftype		acum0, acum1, dummyf1;

 int		ioffW, dmd_lidx;
 float		cppblk, cppblk_fr, cppblk_ft;
 int		num_frsq, dmd_fidx, fridxoff[num_sc_t], bfridx;
 int		*fridxadj[num_sc_t], bk;
 int		fridxoffsr[num_sc_t], slope_valid[num_sc_t], scv;
 double		frsq_mmse[num_sc_t], bfrsq_mmse;

 double		dcidx_step[num_sc_t];
 int		cbidx, cbidxm, eracnt1, scc, scca, eracnt1s;
 int		eracnt1ss, eracnt1sss, eracnt1ssss;
 int		erasi1[2][cbl1], eccnt1[2];
 int		erasi1s[2][cbl1], hec, sec, erag;
 int		erasi1ss[2][cbl1], ssec;
 int		erasi1sss[2][cbl1], sssec, ssssec;
 ftype		distg, distg_last;
 int		sccx[num_sc_t];
 ftype		e1distgsss[2][cbl1];
 int		ssorts[5]={1,4,13,40,121};

 
 unsigned char	fotbuf[cbl1+cbl1];	/* buffer for output file */
 sub_sym_type	*cbmerged, cbmd[cbl1+cbl1], cbmr[cbl1+cbl1];
 int		tmpi, fbtow, fbw, fwr, fnlen, flenb;
 unsigned char	fnbuf[fnlmax+1], fnbbuf[fnlmax+4];	/* buffers for file name */
 ftype		c2phbuf[num_sc_t];

 int		bfridx0;


 int		fmatch, fmatchi[num_sc_l], ji, jfm;
 float		as_off;
 float		*frsqmse1;
 int		frsqmse1i;

 float		frsqmseave[num_sc_t];
 float		frsqmsesd[num_sc_t];

 unsigned char	*fridxbetter[num_sc_t];
 float		cppblk_slope[num_sc_l], cppblk_off[num_sc_l];
 int		era_rsv[num_ec_lvls]={8, 12, 20, 48};
 int		bad_blks;
 double		Ampd;
 double		NORM;
 int		demod_lcnt;

 unsigned short int	ob_cnt, ob_mask, bec_mask;
 int		di_last, scc_last;

 double		sym_pos[num_sc_t][n1+frame_len], ar_step0;
 unsigned char	ma1[num_sc_t][n1+frame_len];
 int		sas, spni, lm_off_amp, ar_state, aridx;
 int		arldx, arzdx, arudx, nzidxl, nzidxr, ii0;
 float		lm_value_amp; 
 int		blk0_ok, eotff;
 float		*froffphu[num_sc_t], adjslope, adjy, ey, ly;
 int		ex, lx, tostate, jref;
 float		frame_ave, frame_sd;
 unsigned long	sidx_dmo, sidxnz;
 double		cx2, cx1, cx0, xidx;

 unsigned char	verstr[]="2003Jun02";
 unsigned char	progname[]="Wyman1x-demod-decode";
 unsigned char	cww[]="Copyright (C) 2003 Barry Sanderson
There is ABSOLUTELY NO WARRANTY  for Wyman1x-demod-decode.
Wyman1x-demod-decode is covered by the the GNU General Public License.
The file \"COPYING.txt\" states the conditions under which you may legally:
copy, modify, or distribute Wyman1x-demod-decode.";


 printf("\n%s version %s\n%s\n\n", progname, verstr, cww);


 if ((argc < 2)) {fprintf(stderr,"%s\n", ins);  exit(-1);}
 if (argc>2) 
 {clk_period= 1.0/(atof(argv[2]));
  if ((clk_period > 500e-6) || (clk_period < 10e-6 ))
  {ERN=31 ; goto ex14;}
 }
 else
 clk_period=1.0/clk_freq_nom ; /* 11025.0 */


 t_step= (float)(clk_period);
 
 statr= stat(argv[1],  &statb);
 if (statr!=0) {ERN=6 ; goto ex14;}
 flen=(unsigned long)(statb.st_size);
 numel= (flen-hdr_lenb)/(sizeof(in_type));

/* finish initialization of distance threshold arrays */
 init_code2();

/* open data file */
 inf=fopen(argv[1], "rb");
 if (inf==NULL) {fprintf(stderr,"Could not open file %s\n", argv[1]); exit(-7);}
 ofn=argv[1];

 sidx=0L; lss=0x00; ; sstepW=750L; hdr_lenb=wavhdrlenb;
 ERN=3; 

 t_step1st=(double)(t_step);
 f_stepd=init_dft3(dftsiz, t_step, snqt);
 f_stepdW=init_dft3(dftsizW, t_step, snqtW);
 dummyd1=PI2/(double)((dftsiz-1));
 for (si=0; si< dftsiz; si++) wind[si]= (float)(0.54 - 0.46*cos(dummyd1*si));
 dummyd1=PI2/(double)((dftsizW-1));
 for (si=0; si< dftsizW; si++) windW[si]= (float)(0.54 - 0.46*cos(dummyd1*si));
 win_flg=1;

 fprintf(stdout,"For data file: %s\n", argv[1]);
 fprintf(stdout,"\tusing version %s of %s program\n", &verstr[0], argv[0]);

/* find leader, restrict search to first half of file */
 lsi=2;
 lead_x[0].state=0;
 lead_x[0].mode= -1 ; /* an illegal value */
 while((sidx < (numel > 100000 ? (numel/2) : (numel*9/10)) ) && (lss==0x00))
 {eoff= do_dftW(sidx, win_flg);
  if (eoff<0) {ERN=4; goto exf;}	/* seek failed */
  set_linmag(outr, outi, &mag[lsi][0], &mag[lsi][0], dftsizW);
  lead_x[lsi].xcv= xcor4b(&trail_spec_nom[0], &mag[lsi][lsi_startW], trail_spec_len_nom,
		lsi_countW, &lead_x[lsi].offi);
  mi=0;
  while ((mi<num_modes_l) && (lss==0x00))
  {if (lead_x[lsi].xcv > xcthld) 
   {as_off= (lead_x[lsi].offi + lsi_startW - lsi_centerW)*(float)(f_stepdW);
    for (fi=0; fi< num_sc_l; fi++)
     exidx[fi]= (int)(floor((double)((pm_nom_lt[fi] + as_off)/(float)(f_stepd)+ 0.5)));

    eoff= do_dft(sidx+dftsizW/2 - dftsiz/2, win_flg);
    if (eoff<0) {ERN=4; goto exf;}	/* seek failed */
    set_linmag(outr, outi, &mag[0][0], &mag[0][0], dftsiz);
    lead_x[0].nlm= find_local_max3(&mag[0][0], &lead_x[0].lmfreq[0], &lead_x[0].lmqual[0], &lead_x[0].lmidx[0], &lead_x[0].lmbm);
    if (lead_x[0].nlm >= nsc_min)
    {fmatch= lm_match3(&lead_x[0].lmfreq[0], &lead_x[0].lmidx[0], lead_x[0].nlm,
	 &fmatchi[0], 2.0*fdthld_nom);
     dummyf1= 0.0;
     for (ji=0; ji< lead_x[0].nlm; ji++) dummyf1 += lead_x[0].lmqual[ji];
     r53[lsi]= (float)(fmatch)*dummyf1/(float)(lead_x[0].nlm);
     if (fmatch >= nsc_min)
     {jfm=0;
      for (ji=0; ji<lead_x[0].nlm; ji++)
      {if(lead_x[0].lmidx[ji] == fmatchi[jfm])
        lead_x[0].lmfreq[jfm++]= lead_x[0].lmfreq[ji];
      }
      sidxnz= (sidx+dftsizW/2 < dftsiz/2 ? 0 : sidx+dftsizW/2 - dftsiz/2);
      switch(lead_x[0].state)
       {case 0:
	 err= cbfitx(&lead_x[0].lmfreq[0], &fmatchi[0], fmatch, &lead_x[0].slope, &lead_x[0].off_hz);
	 if (err <0) {ERN=10; goto exf;}
	 lead_x[0].state=1;
	 lead_x[0].sidx= sidxnz;
	 lead_x[0].r53= r53[lsi];
	 lead_x[0].mode= mi;
	 lead_x[0].lmcount=fmatch;
	 for (ji=0; ji<fmatch; ji++) lead_x[0].lmgoodi[ji]= fmatchi[ji];
	 break;
        case 1:
	 if (r53[lsi] > lead_x[0].r53)
	 {err= cbfitx(&lead_x[0].lmfreq[0], &fmatchi[0], fmatch, &lead_x[0].slope, &lead_x[0].off_hz);
	  if (err <0) {ERN=10; goto exf;}
	  lead_x[0].sidx= sidxnz;
	  lead_x[0].r53= r53[lsi];
	  lead_x[0].mode= mi;
          lead_x[0].lmcount=fmatch;
          for (ji=0; ji<fmatch; ji++) lead_x[0].lmgoodi[ji]= fmatchi[ji];
	 }
	 else
	 {lead_x[0].state=2; 
	 }
	 break;
        case 2:
	 if (r53[lsi] > lead_x[0].r53)
	 {err= cbfitx(&lead_x[0].lmfreq[0], &fmatchi[0], fmatch, &lead_x[0].slope, &lead_x[0].off_hz);
	  if (err <0) {ERN=10; goto exf;}
	  lead_x[0].sidx= sidxnz;
	  lead_x[0].r53= r53[lsi];
	  lead_x[0].mode= mi;
          lead_x[0].lmcount=fmatch;
          for (ji=0; ji<fmatch; ji++) lead_x[0].lmgoodi[ji]= fmatchi[ji];
	 }
	 else
	 {lead_x[0].state=3; 
	  lss |= l1_found ;
	 }
	 break;
       } /* end of switch lead_dat.state */
     }  /* end of fmatch large enough */
    }	/* end of at least nsc_min loacl maximums present */
   } /* end of cross-corelation value large enough */   
   else
   {switch(lead_x[0].state)
    {case 1:
     case 2:
     lead_x[0].state=3;
     lss |= l1_found ;
    } /* end of switch */
   } /* end of cross-correlation value below threshold */
   mi++;
  } /* end of while mi < num_modes_l */
  sidx+=sstepW;
 }	/* end of while (looking for leader) */

 if(lead_x[0].slope > slope_max) lead_x[0].slope = slope_max ;
 else
 {if (lead_x[0].slope < slope_min) lead_x[0].slope= slope_min;
 }  
 if(lead_x[0].off_hz > off_hz_max) lead_x[0].off_hz = off_hz_max ;
 else
 {if (lead_x[0].off_hz < off_hz_min) lead_x[0].off_hz = off_hz_min;
 }  

#if(dbgl>10)
/* for now, print results to stdout */
 fprintf(stdout,"\nSample index  ---------- leader --------------------\n");
 fprintf(stdout,"              r53 val    offset (Hz)     freq. scale   mode index\n");
 fprintf(stdout,"%8ld     %7.4f     %7.2f           %9.6f      %2d\n",
		 lead_x[0].sidx, lead_x[0].r53, lead_x[0].off_hz, lead_x[0].slope,
		 lead_x[0].mode);
#endif	/* dbgl>10 */



 if(lss== 0x00) {ERN=10; goto exf;}
/* find trailer, restrict search to last half of file */
 lead_x[1].state=0;
 lead_x[1].mode= -1 ; /* an illegal value */
 sidx= numel-dftsizW;
 ioffW=(int)((double)(lead_x[0].off_hz)/f_stepdW);
 while((sidx > numel/2) && (lss==0x01))
 {eoff= do_dftW(sidx, win_flg);
  if (eoff<0) {ERN=4; goto exf;}	/* seek failed */
  set_linmag(outr, outi, &mag[1][0], &mag[1][0], dftsizW);
  lsi=1;
  lead_x[lsi].xcv= xcor4b(&trail_spec_nom[0], &mag[lsi][tsi_start+ioffW], trail_spec_len_nom,
		tsi_count, &lead_x[lsi].offi);

  mi=0;
  while ((mi<num_modes_l) && (lss==0x01))
  {
   if (lead_x[lsi].xcv > xcthld) 
   { for (si=0; si< num_sc_l; si++)
      exidxW[si]= (int)(floor((double)((pm_nom_lt[si] + lead_x[0].off_hz)/(float)(f_stepdW)+ 0.5))) +
	tsi_start +lead_x[lsi].offi -lsi_centerW ;

    lead_x[lsi].nlm= find_local_max3W(&mag[lsi][0], &lead_x[lsi].lmfreq[0], &lead_x[lsi].lmqual[0], &lead_x[lsi].lmidx[0], &lead_x[lsi].lmbm);
    if (lead_x[lsi].nlm >= nsc_min)
    {fmatch= lm_match3(&lead_x[lsi].lmfreq[0], &lead_x[lsi].lmidx[0], lead_x[lsi].nlm,
	 &fmatchi[0], 8.0*fdthld_nom);
     dummyf1= 0.0;
     for (ji=0; ji< lead_x[lsi].nlm; ji++) dummyf1 += lead_x[lsi].lmqual[ji];
     r53[lsi]= (float)(fmatch)*dummyf1/(float)(lead_x[lsi].nlm);
     if (fmatch>= nsc_min)
     {jfm=0;
      for (ji=0; ji<lead_x[lsi].nlm; ji++)
      {if(lead_x[lsi].lmidx[ji] == fmatchi[jfm])
        lead_x[lsi].lmfreq[jfm++]= lead_x[lsi].lmfreq[ji];
      }
      switch(lead_x[lsi].state)
       {case 0:
         err= cbfitx(&lead_x[lsi].lmfreq[0], &fmatchi[0], fmatch, &lead_x[lsi].slope, &lead_x[lsi].off_hz);
	 if (err <0) {ERN=10; goto exf;}
	 lead_x[lsi].state=1;
	 lead_x[lsi].sidx= sidx;
	 lead_x[lsi].r53= r53[lsi];
	 lead_x[lsi].mode= mi;
         lead_x[lsi].lmcount=fmatch;
         for (ji=0; ji<fmatch; ji++) lead_x[lsi].lmgoodi[ji]= fmatchi[ji];
	 break;
        case 1:
	 if (r53[lsi] > lead_x[lsi].r53)
	 {
           err= cbfitx(&lead_x[lsi].lmfreq[0], &fmatchi[0], fmatch, &lead_x[lsi].slope, &lead_x[lsi].off_hz);
	  if (err <0) {ERN=10; goto exf;}
	  lead_x[lsi].sidx= sidx;
	  lead_x[lsi].r53= r53[lsi];
	  lead_x[lsi].mode= mi;
          lead_x[lsi].lmcount=fmatch;
          for (ji=0; ji<fmatch; ji++) lead_x[lsi].lmgoodi[ji]= fmatchi[ji];
	 }
	 else
	 {lead_x[lsi].state=2; 
	 }
	 break;
        case 2:
	 if (r53[lsi] > lead_x[lsi].r53)
	 {
           err= cbfitx(&lead_x[lsi].lmfreq[0], &fmatchi[0], fmatch, &lead_x[lsi].slope, &lead_x[lsi].off_hz);
	  if (err <0) {ERN=10; goto exf;}
	  lead_x[lsi].sidx= sidx;
	  lead_x[lsi].r53= r53[lsi];
	  lead_x[lsi].mode= mi;
          lead_x[lsi].lmcount=fmatch;
          for (ji=0; ji<fmatch; ji++) lead_x[lsi].lmgoodi[ji]= fmatchi[ji];
	 }
	 else
	 {lead_x[lsi].state=3; 
	  lss |= t1_found ;
	 }
	 break;
       } /* end of switch lead_dat.state */
     } /* end of fmatch great enough */
    } /* enough local maximums */
   } /* end of cross-correlation value great enough */   
   else
   {lsi=1;
    switch(lead_x[lsi].state)
    {case 1:
     case 2:
	lead_x[lsi].state=3;
	lss |= t1_found ;
    } /* end of switch */
   } /* end of cross-correlation value below threshold */
   mi++;
  } /* end of while mi < num_modes_l */
  sidx-=sstepW;
 }	/* end of while (looking for trailer) */
  
#if(dbgl>10)
/* for now, print results to stdout */
 fprintf(stdout,"\nSample index  ---------- trailer --------------------\n");
 fprintf(stdout,"              r53 val    offset (Hz)     freq. scale   mode index\n");
 fprintf(stdout,"%8ld     %7.4f     %7.2f           %9.6f      %2d\n",
		 lead_x[1].sidx, lead_x[1].r53, lead_x[1].off_hz, lead_x[1].slope,
		 lead_x[1].mode);
#endif	/* dbgl>10 */

 if(lss!= 0x03) {ERN=11; goto exf;}

 f_slope= lead_x[0].slope;
 f_intercept= lead_x[0].off_hz;
 demod_lcnt=0;
 delay= delay_nom;
 ns= (int)(ave_nom);
 NS=(double)(ave_nom);

 drift_rate_hps= (lead_x[1].off_hz - lead_x[0].off_hz)/((lead_x[1].sidx - lead_x[0].sidx)*clk_period);
 astepdot= drift_rate_hps*PI2*clk_period*clk_period;
 cx2= astepdot/2.0;

 cppblk= cppblk_nom/f_slope;	/* set clock periods per block */

/* found end of message */
 dmd_lidx= lead_x[1].sidx - (dmoff + lead_x[0].sidx)- (delay+ns) + dmoffW;

/* calculate expected number of framing sequences */
 num_frsq= (int)(ceil((double)(dmd_lidx/cppblk)));

/* allocate memory to store framing sequence mean squared errors */
  frsqmse1=malloc((size_t)(num_frsq*sizeof(float)));
  if(frsqmse1==NULL) {ERN=2; goto ex14;}

/* allocate memory to store positions for each framing sequence */
  for (j=0; j<num_sc_t; j++)
  {fridxadj[j]=malloc((size_t)(num_frsq*sizeof(int)));
   if (fridxadj[j]==NULL)
   {j--; while(j>=0) free(fridxadj[j--]);
     free(frsqmse1); ERN=2; goto ex14;
   }
  }

/* allocate memory to store "good/bad" indication for each framing sequence */
 for (j=0; j<num_sc_t; j++)
 {fridxbetter[j]=malloc((size_t)(num_frsq*sizeof(unsigned char)));
  if (fridxbetter[j]==NULL)
  {j--; while(j>=0) free(fridxbetter[j--]);
    free(frsqmse1); ERN=2; goto exf1;
  }
 }

/* allocate memory to store offset information for each framing sequence */
 for (j=0; j<num_sc_t; j++)
 {froffphu[j]=malloc((size_t)(num_frsq*sizeof(float)));
  if (froffphu[j]==NULL)
  {j--; while(j>=0) free(froffphu[j--]);
    free(frsqmse1); ERN=2; goto exf2;
  }
 }


DEMOD:
/* set local oscillator frequencies */
 for (lo=0; lo<num_sc_t; lo++) pm_ft_ctr[lo]= pm_nom_ctr[lo]*f_slope + f_intercept;
/* set origin for local oscillator angle calculation */
 sidx_dmo= dmoff + lead_x[0].sidx -s0off;

 cppblk= cppblk_nom/f_slope;	/* set clock periods per block */
#if(dbgl > 0)
 fprintf(stdout,"\tclock periods per block is: %10.3f\n", cppblk);
#endif


/* demodulate each subcarrier */
 for (lo=0; lo<num_sc_t; lo++)
 {irf=fopen(&scratch[0], "wb");
  if (irf==NULL)
  {fprintf(stderr,"Could not open file %s\n", &scratch[0]); ERN=8; goto exf3;}
  irf2=fopen(&scratch2[0], "wb");
  if (irf2==NULL)
  {fprintf(stderr,"Could not open file %s\n", &scratch2[0]);
   ERN=8; fclose(irf); remove(&scratch[0]); goto exf3;
  }
  sidx=dmoff + lead_x[0].sidx;
  end_of_file=0;  ca= (cx0=0.0);
  astep= (pm_ft_ctr[lo]-(drift_rate_hps*clk_period*s0off))*PI2*clk_period;
  astepd= astep;
  cx1= astep + cx2;
  bias[lo]=0.0;

  while (end_of_file == 0)
  {val_to_read=s0blen; dum_len=0;

   if (sidx < s0off)	/* this case is not correct 2003 Apr 09 */
   {dum_len=(int)(s0off-sidx);
    for (si=0; si < dum_len; si++) samsi[si]= dumd;
    val_to_read -= dum_len;
    if (fseek(inf, hdr_lenb+0L, SEEK_SET)!=0)
    {ERN=4; fclose(irf); fclose(irf2); goto exj;}
    ang_ref= 0.0;
    ang_off= astep*(double)(dum_len) + astepdot*(double)(dum_len*(dum_len-1))/2.0;
    ca= ang_ref-ang_off;
    ca= fmod(ca, PI2);
    if (ca < 0.) ca += PI2; else {if (ca > PI2) ca -= PI2;}
   }
   else
   {if (fseek(inf, hdr_lenb+(sidx-s0off)*sizeof(in_type), SEEK_SET)!=0)
    {ERN=4; fclose(irf); fclose(irf2); goto exj;}
   }
   val_to_write= fap;
   val_read=fread(&samsi[dum_len], sizeof(in_type), val_to_read, inf);
   if (val_read < val_to_read)
   {if (val_read < s0off) {end_of_file=1; val_to_write=0;}
    else
    {val_to_write= (val_read > (s0off+fap) ? fap: (val_read - s0off));}
    for (si=val_read; si< val_to_read; si++) samsi[si]= dumd;
   }
/* translate to baseband, for data in samsi */
/* multiply by sin and cos of local oscillator */
   for (si=0; si< s0blen; si++)
   {Lo1c[si]=(ftype)((double)(samsi[si])*cos(ca));
    Lo1s[si]=(ftype)((double)(samsi[si])*sin(ca));
    astepd += astepdot;
    ca+= astepd;    if (ca > PI2) ca-= PI2; else {if (ca < -PI2) ca+= PI2;}
   }

/* do decimation */
   ds1t(&Lo1c[ds1dly], Lo1cf, s1blen);
   ds2t(&Lo1cf[ds2dly], Lo2cs1, s2blen);
   ds3t(&Lo2cs1[ds3dly], Lo2cs2, s3blen);
   ds1t(&Lo1s[ds1dly], Lo1sf, s1blen);
   ds2t(&Lo1sf[ds2dly], Lo2ss1, s2blen);
   ds3t(&Lo2ss1[ds3dly], Lo2ss2, s3blen);

/* do interpolation */
   di=I40xxPI(&Lo2cs2[is1off], Ssbc, ipts);
   if (di <0) {ERN=32; fclose(irf); fclose(irf2); goto exj;}
   di=I40xxPI(&Lo2ss2[is1off], Ssbs, ipts);
   if (di <0) {ERN=32; fclose(irf); fclose(irf2); goto exj;}

/* calculate angle, and remove discontinuities */
   for (si=0; si< fap; si++)
   {Ampd= -atan2((double)(Ssbs[si]), (double)(Ssbc[si]));
    if ((sidx== (dmoff + lead_x[0].sidx))&& (si==0)) {pv[lo]=Ampd;}
    if ((Ampd+bias[lo] -pv[lo]) < -PI) bias[lo] += PI2;
    else {if ((Ampd+bias[lo]-pv[lo]) > PI) bias[lo] -= PI2;}
    pv[lo]= Ampd + bias[lo];
    Amp[si] = (ftype)(pv[lo]);
   }
   val_written= fwrite(&Amp[0], sizeof(ftype), val_to_write, irf);
   if (val_written != (size_t)(val_to_write))
   {ERN=5; fclose(irf); fclose(irf2); goto exj;}
/* calculate Amplitude */
   for (si=0; si< fap; si++)
   {buf2[si]=(ftype)(sqrt((double)(Ssbs[si]*Ssbs[si] + Ssbc[si]*Ssbc[si])));}
   val_written= fwrite(&buf2[0], sizeof(ftype), val_to_write, irf2);
   if (val_written != (size_t)(val_to_write))
   {ERN=5; fclose(irf); fclose(irf2); goto exj;}   
   sidx+= fap;

/* adjust angle of local oscillators */
   xidx= (double)(sidx -sidx_dmo -s0off);
   astepd= astep + astepdot*(xidx+1.0);
   ang_off= cx2*xidx*xidx +cx1*xidx +cx0;
   ca= fmod(ang_off, PI2);
   if (ca < 0.) ca += PI2; else {if (ca > PI2) ca-= PI2;}
/***99feb03	ca2 advances fap steps per loop with no overlap or gap
		from loop to loop (unlike ca, which has a lot of overlap)
99feb03***/

  } /* end of while end_of_file */
  fflush(irf); fclose(irf); fclose(irf2);


 NORM= PI2/sqrt(phased2);
/* calculate and store (in a separate buffer for each lo) the average difference
	data for the pm8-scratch file just written */
  if(lo==0)	/* if first pass */
  {statirr= stat(&scratch[0],  &statir);
   if (statirr!=0) {ERN=6 ; goto exj;}
   flen_dmd=(unsigned long)(statir.st_size);
   numel_dmd= flen_dmd/(sizeof(ftype));

   if(demod_lcnt==0)
   {
/* allocate buffers */
    for (j=0; j<num_sc_t; j++)
    {dif7n[j]=malloc((size_t)(numel_dmd*sizeof(ftype)));
     if (dif7n[j]==NULL)
     {j--; while(j>=0) free(dif7n[j--]);
      ERN=2; goto exj;
     }
    }
    for (j=0; j<num_sc_t; j++)
    {dmamp[j]=malloc((size_t)(numel_dmd*sizeof(ftype)));
     if (dmamp[j]==NULL)
     {j--; while(j>=0) free(dmamp[j--]);
      for (j=0; j<num_sc_t; j++) free(dif7n[j]);
      ERN=2; goto exj;
     }
    }
   }

   ERN=2;
   if ((fbufin=malloc((size_t)((numel_dmd)* sizeof(ftype)))) == NULL) goto exk;
   if ((fbufot=malloc((size_t)((numel_dmd)* sizeof(ftype) ))) == NULL) {free(fbufin); goto exk;}
  } /* end of lo==0 */

  ERN=7;
  irf=fopen(&scratch[0], "rb");
  if (irf==NULL)
  {fprintf(stderr,"Could not open file %s\n", scratch);
   free(fbufin); free(fbufot); goto exk;
  }
  irf2=fopen(&scratch2[0], "rb");
  if (irf2==NULL)
  {fprintf(stderr,"Could not open file %s\n", scratch2);
   fclose(irf); free(fbufin); free(fbufot); goto exk;
  }

  ir= fread(&fbufin[0], sizeof(ftype), numel_dmd, irf);
  if (ir!=numel_dmd)
  {ERN=13 ; fclose(irf); free(fbufin); free(fbufot); goto exk;}
  ir= fread(&fbufot[0], sizeof(ftype), numel_dmd, irf2);
  if (ir!=numel_dmd)
  {ERN=13 ; fclose(irf2); fclose(irf); free(fbufin); free(fbufot); goto exk;}

  dmd_cnt=0;
  for (fbidx=delay+ns; fbidx<numel_dmd; fbidx++)
  {acum1=0.0;
   for (si=1-ns; si<1; si++) acum1+= fbufot[fbidx+si];
   dmamp[lo][dmd_cnt++]= acum1/NS ;
  } 


  dmd_cnt=0;
  for (fbidx=delay-ns+lp1dly; fbidx<numel_dmd; fbidx++)
   fbufot[dmd_cnt++]= (fbufin[fbidx] - fbufin[fbidx-delay])/NORM;
  lpf01(&fbufot[lp1dly], &dif7n[lo][0], dmd_cnt-lp1len+1);

  fclose(irf); fclose(irf2);
 } /* end of for lo (demodulation) */

 free(fbufin); free(fbufot);
 dmd_cnt = dmd_cnt - lp1len +1;

#if(dbgl > 30)
 fprintf(stdout," %8d demodulated values calculated.  End of message at index %8d\n",
	dmd_cnt, dmd_lidx);
#endif


/* calculate starting point for search */
 dmd_fidx= dmd_lidx -(pm8frame_len -1) -(int)((float)((num_frsq-1))*cppblk) -(frsq_srsam-1);
 if (dmd_fidx < 0) dmd_fidx=0;
#if(dbgl > 10)
 fprintf(stdout,"Expecting to find %3d framing sequences\n", num_frsq);
 fprintf(stdout," Starting search at index %6d\n", dmd_fidx);
#endif


/* search for framing sequences */
 for (lo=0; lo<num_sc_t; lo++)
 {for (si=dmd_fidx; si<dmd_fidx+frsq_srsam; si++)
  {acum0d=1.0; sumx=0.0;
   for (j=0; j<num_frsq; j++)
   {acum1= 0.0; 
    di= si+ (int)((float)(j)*cppblk);
    for (ki=0; ki<pm8frame_len; ki++)
    {dummyf1=dif7n[lo][di+ki]- pm8_frame[ki];
     acum1 += dummyf1*dummyf1;
    } /* end of for ki */
    frsqmse1[j] = acum1/(float)(pm8frame_len);
    sumx += frsqmse1[j];
    if(j==0) {sumy= frsqmse1[j]; frsqmse1i=j;}
    else
    if (frsqmse1[j] > sumy) {sumy= frsqmse1[j]; frsqmse1i=j;}
   } /* end of for j */
   /* calculate product of mean and standard deviation,
	 leaving out largest value, if num_frsq > 2 */
   sumx = sumx - (num_frsq > nfmin ? frsqmse1[frsqmse1i] : 0);
   sumx = sumx/(float)(num_frsq - (num_frsq > nfmin ? 1 : 0));
   sumx2=0.0;
   for (j=0; j<num_frsq; j++)
   {if((j==frsqmse1i) && (num_frsq > nfmin)) continue;
    sumx2 = sumx2 + (frsqmse1[j] - sumx) * (frsqmse1[j] - sumx);
   } /* end of for j */
   acum0d= (double)(sumx)* sqrt((double)(sumx2/(double)(num_frsq- (num_frsq > nfmin ? 2 : 1) )));
   if (si== dmd_fidx)
   {frsq_mmse[lo]= acum0d; fridxoff[lo]= si-dmd_fidx;
    frsqmsesd[lo]= (float)(acum0d)/sumx;
    frsqmseave[lo]= sumx; 
   }
   else
   {if (acum0d < frsq_mmse[lo])
    {frsq_mmse[lo]= acum0d; fridxoff[lo]= si-dmd_fidx;
     frsqmsesd[lo]= (float)(acum0d)/sumx;
     frsqmseave[lo]= sumx; 
    }
   }
  } /* end of for si */
 } /* end of for lo searching for framing sequences */

#if(dbgl > 10)
 fprintf(stdout,"Results of initial search for framing sequences\n");
 fprintf(stdout,"  sub-carrier  offset   min product\n");
 for (lo=0; lo<num_sc_t; lo++)
  fprintf(stdout,"     %1d      %6d      %10.5g\n", lo, fridxoff[lo], frsq_mmse[lo]);
#endif


/* determine clock periods per block, from best framing sequence data*/
/* don't use data from subcarrier 0, due to non-linear phase filters in
	radios */
/* set index of best framing sequence */

/* use best (except for lowest freq) overall subcarrier for next step */
 bfridx= 1;
 bfrsq_mmse= frsq_mmse[1];
 for (lo=2; lo<num_sc_t; lo++)
 {if (frsq_mmse[lo]< bfrsq_mmse)
  {bfrsq_mmse= frsq_mmse[lo]; bfridx=lo;}
 }

/* use best framing sequence for this subcarrier as the reference */
 for (j=0; j<num_frsq; j++)
 {acum1= 0.0;
  di= fridxoff[bfridx]+dmd_fidx + (int)((float)(j)*cppblk);
  for (ki=0; ki<pm8frame_len; ki++)
  {dummyf1=dif7n[bfridx][di+ki]- pm8_frame[ki];
   acum1 += dummyf1*dummyf1;
  } /* end of for ki */
  acum1 = acum1/(float)(pm8frame_len);
  if(j==0) {sumy= acum1; bfridx0=j;}
  else
  if (acum1 < sumy) {sumy= acum1; bfridx0=j;}
 } /* end of for j */

 fridxadj[bfridx][0]= fridxoff[bfridx] + dmd_fidx;
 for (j=1; j<num_frsq; j++)
  fridxadj[bfridx][j]= fridxadj[bfridx][0]+ (int)(floor((double)((float)(j)*cppblk + 0.5)));

/* search for a better value for cppblk */
 lo= bfridx;
 for (si=-cppblk_ossr; si<=cppblk_ossr; si++)
 {acum0d=1.0; sumx=0.0;
  for (j=0; j<num_frsq; j++)
  {acum0= (ftype)(j-bfridx0)*(cppblk+(float)(si));
   bk= fridxadj[bfridx][bfridx0] + (int)(floor((double)(acum0+0.5)));
   acum1= 0.0;
   for (ki=0; ki<pm8frame_len; ki++)
   {dummyf1=dif7n[lo][bk+ki]- pm8_frame[ki];
    acum1 += dummyf1*dummyf1;
   } /* end of for ki */
   frsqmse1[j] = acum1/(float)(pm8frame_len);
   sumx += frsqmse1[j];
   if(j==0) {sumy= frsqmse1[j]; frsqmse1i=j;}
   else
   if (frsqmse1[j] > sumy) {sumy= frsqmse1[j]; frsqmse1i=j;}
  } /* end of for j */

/* calculate product of mean and standard deviation,
	 leaving out largest value, if num_frsq > 2 */
   sumx = sumx - (num_frsq > nfmin ? frsqmse1[frsqmse1i] : 0);
   sumx = sumx/(float)(num_frsq - (num_frsq > nfmin ? 1 : 0));
   sumx2=0.0;
   for (j=0; j<num_frsq; j++)
   {if((j==frsqmse1i) && (num_frsq > nfmin)) continue;
    sumx2 = sumx2 + (frsqmse1[j] - sumx) * (frsqmse1[j] - sumx);
   } /* end of for j */
   acum0d= (double)(sumx)* sqrt((double)(sumx2/(double)(num_frsq- (num_frsq > nfmin ? 2 : 1) )));
  if (si== -cppblk_ossr)
  {frsq_mmse[lo]= acum0d; fridxoffsr[lo]= si;}
  else
  {if (acum0d < frsq_mmse[lo]) {frsq_mmse[lo]= acum0d; fridxoffsr[lo]= si;}
  }
 } /* end of for si */
 if(fridxoffsr[lo]==0) cppblk_fr= cppblk;
 else cppblk_fr= cppblk+(float)(fridxoffsr[lo]);
 
#if(dbgl>20)
 fprintf(stdout,"Best min product found is: %10.5g, using an offset of %3d\n", frsq_mmse[lo], fridxoffsr[lo]);
 fprintf(stdout,"Now using %10.3f clock periods per block\n", cppblk_fr);
#endif /* dbgl>20 */

/* refine framing sequence positions, based on actual data, for best subcarrier */
 for (j=0; j<num_frsq; j++)
 {acum0= xcor4c(&pm8_frame[0], &dif7n[bfridx][fridxadj[bfridx][j]-frosftr2],
		 pm8frame_len, 2*frosftr2 +1, &ki, &frame_ave, &frame_sd);
  if((acum0 > frxcthld) &&
     (frame_sd >= pm8a_frame_sd_min) &&
     (frame_sd <= pm8a_frame_sd_max)
    )
  {bk= (int)(ki)-frosftr2;
   fridxadj[bfridx][j] += bk;
   fridxbetter[bfridx][j]= 0x01;
  }
  else {fridxbetter[bfridx][j]= 0x00;}
 } /* end of for j */

/* calculate least squares best fit straight line through the "better than
	 average" framing sequence positions */
 sumx=0.0; sumy=0.0; sumxy=0.0; x2sum=0.0; ki=0;
 for (j=0; j<num_frsq; j++)
 {if (fridxbetter[bfridx][j] == 0x01)
  {cx=(float)(j);
   sumx += cx; x2sum += cx*cx;
   sumxy += cx*(float)(fridxadj[bfridx][j]);
   sumy += (float)(fridxadj[bfridx][j]);
   ki++;
  }
 }
 if (ki > 2)
 {cppblk_slope[bfridx] = ((float)(ki)*sumxy - sumx*sumy)/((float)(ki)*x2sum - sumx*sumx);
  cppblk_off[bfridx]= sumy/(float)(ki) - cppblk_slope[bfridx]*sumx/(float)(ki);
 }
 else
 {cppblk_slope[bfridx]= cppblk_fr;
  cppblk_off[bfridx]= dmd_fidx + fridxoff[bfridx];
 }

/* if offset is different from 0, and this is the first pass,
   demodulate again */
 if (fridxoffsr[lo] != 0)
 {if (demod_lcnt<1)
  {f_slope= (float)(cppblk_nom)/cppblk_slope[bfridx]; /*(float)(cppblk_nom)/cppblk_fr;*/
#if(dbgl >20)
   fprintf(stdout,"Starting pass number %1d through the demodulator.\n", demod_lcnt+2);
#endif
   demod_lcnt++;
   goto DEMOD;
  }
 }


/* calculate expected position for all framing sequences for all subcarriers,
	based on best framing sequence for best subcarrier */
 for (j=0; j<num_frsq; j++)
 {acum0= (ftype)(j-bfridx0)*cppblk_fr;
  bk= fridxadj[bfridx][bfridx0] - soff[bfridx] + (int)(floor((double)(acum0+0.5)));
  for (lo=0; lo<num_sc_t; lo++) fridxadj[lo][j] = bk+ soff[lo];
 }

/* take care of lowest frequency subcarrier as a special case */
 lo=0;
 for (si=-36; si< 48; si++)
 {acum0d=1.0; sumx=0.0;
  for (j=0; j<num_frsq; j++)
  {acum1= 0.0;
   di= si+ fridxadj[lo][j];
   for (ki=0; ki<pm8frame_len; ki++)
   {dummyf1=dif7n[lo][di+ki]- pm8_frame[ki];
    acum1 += dummyf1*dummyf1;
   } /* end of for ki */
   frsqmse1[j] = acum1/(float)(pm8frame_len);
   sumx += frsqmse1[j];
   if(j==0) {sumy= frsqmse1[j]; frsqmse1i=j;}
   else
   if (frsqmse1[j] > sumy) {sumy= frsqmse1[j]; frsqmse1i=j;}
  } /* end of for j */

/* calculate product of mean and standard deviation,
	leaving out largest value, if num_frsq > 2 */
  sumx = sumx - (num_frsq > nfmin ? frsqmse1[frsqmse1i] : 0);
  sumx = sumx/(float)(num_frsq - (num_frsq > nfmin ? 1 : 0));
  sumx2=0.0;
  for (j=0; j<num_frsq; j++)
  {if((j==frsqmse1i) && (num_frsq > nfmin)) continue;
   sumx2 = sumx2 + (frsqmse1[j] - sumx) * (frsqmse1[j] - sumx);
  } /* end of for j */
  acum0d= (double)(sumx)* sqrt((double)(sumx2/(double)(num_frsq- (num_frsq > nfmin ? 2 : 1) )));
  if (si== -36)
  {frsq_mmse[lo]= acum0d; fridxoffsr[lo]= si;}
  else
  {if (acum0d < frsq_mmse[lo]) {frsq_mmse[lo]= acum0d; fridxoffsr[lo]= si;}
  }
 } /* end of for si */
 if(fridxoffsr[lo] != 0)
 {for (j=0; j<num_frsq; j++) fridxadj[lo][j] += fridxoffsr[lo];}


/* refine framing sequence positions, based on actual data */
 for (j=0; j<num_frsq; j++)
 {for (lo=0; lo<num_sc_t; lo++)
  {acum0= xcor4c(&pm8_frame[0], &dif7n[lo][fridxadj[lo][j]-frosftr2],
		pm8frame_len, 2*frosftr2 +1, &ki, &frame_ave, &frame_sd);
   if ((acum0 > frxcthld) &&
       (frame_sd >= pm8a_frame_sd_min) &&
       (frame_sd <= pm8a_frame_sd_max)
      )
   {bk= (int)(ki)-frosftr2;
    fridxadj[lo][j] += bk;
    fridxbetter[lo][j]= 0x01;
    acum1= 0.0;
    for (di=0; di< pm8frame_len; di++) acum1 += dif7n[lo][fridxadj[lo][j]+di];
    froffphu[lo][j]= acum1/(float)(pm8frame_len) - pm8a_frame_ave;
   }
   else
   {fridxbetter[lo][j]= 0x00;
    froffphu[lo][j]= 0.0;
   }
  } /* end of for lo */
 } /* end of for j */
 free(frsqmse1);

#if (dbgl>20)
 fprintf(stdout,"Good cases => 1,  Bad cases => 0\n");
 fprintf(stdout," framing --------- sub-carrier index----------\n");
 fprintf(stdout,"  index   0    1    2    3    4    5    6    7\n");
 for (j=0; j<num_frsq; j++)
 {fprintf(stdout,"   %3d  ", j);
  for (lo=0; lo<num_sc_t; lo++)
   fprintf(stdout,"  %1c  ", fridxbetter[lo][j]+ 0x30);
  fprintf(stdout,"\n");
 }
#endif

/* calculate least squares best fit straight line through the "better than
	 average" framing sequence positions */
 scv=0;	/* to count number of "valid subcarriers" */
 for (lo=0; lo<num_sc_t; lo++)
 {sumx=0.0; sumy=0.0; sumxy=0.0; x2sum=0.0; ki=0;
  for (j=0; j<num_frsq; j++)
  {if (fridxbetter[lo][j] == 0x01)
   {cx=(float)(j);
    sumx += cx; x2sum += cx*cx;
    sumxy += cx*(float)(fridxadj[lo][j]);
    sumy += (float)(fridxadj[lo][j]);
    ki++;
   }
  }
  if (ki > 2)
  {slope_valid[scv++]= lo;
   cppblk_slope[lo] = ((float)(ki)*sumxy - sumx*sumy)/((float)(ki)*x2sum - sumx*sumx);
   cppblk_off[lo]= sumy/(float)(ki) - cppblk_slope[lo]*sumx/(float)(ki);
  }
  else
  {cppblk_slope[lo]= 0.0;
   cppblk_off[lo]= dmd_fidx + fridxoff[lo];
  }
 } /* end of for lo */


#if(dbgl > 20)
 fprintf(stdout,"Results of best fit straight lines through better points\n");
 fprintf(stdout,"  sub-carrier    offset         slope\n");
 for (lo=0; lo<num_sc_t; lo++)
  fprintf(stdout,"     %1d   %12.7g      %12.7g\n", lo, cppblk_off[lo], cppblk_slope[lo]);
#endif

/* leave out highest and lowest, use average of remaining slopes for subsequent calculations */
 if (scv > 2)
 {sumy=(sumx= cppblk_slope[slope_valid[0]]);
	/* find sum and max */
  for (j=1; j<scv; j++)
  {sumx+= cppblk_slope[slope_valid[j]];
   if (cppblk_slope[slope_valid[j]] > sumy) sumy= cppblk_slope[slope_valid[j]];
  }
  sumx = sumx - sumy;
	/* find min */
  sumy= cppblk_slope[slope_valid[0]];
  for (j=1; j<scv; j++)
  {if (cppblk_slope[slope_valid[j]] < sumy) sumy= cppblk_slope[slope_valid[j]];
  }
  sumx = sumx - sumy;
  cppblk_ft= sumx/(float)(scv -2);
 }
 else
 {sumx= cppblk;
  if (scv > 0)
  {for (j=0; j< scv; j++) sumx += cppblk_slope[slope_valid[j]];
   cppblk_ft= sumx/(float)(scv+1);
  }
  else cppblk_ft= cppblk;
 }

/* limit cppblk_ft */
 if (cppblk_ft > cppblk_max) cppblk_ft= cppblk_max;
 else {if (cppblk_ft < cppblk_min) cppblk_ft= cppblk_min;}

#if(dbgl > 0)
 fprintf(stdout,"Now using %12.7g clock periods per block\n", cppblk_ft);
#endif 

/* limit cppblk_off values */
 for (lo=0; lo<num_sc_t; lo++)
 {if (lo==bfridx) continue;
  if (((cppblk_off[lo] - dmd_fidx) > (fridxoff[lo] + cppblk_off_lim)) ||
      ((cppblk_off[lo] - dmd_fidx) < (fridxoff[lo] - cppblk_off_lim)))
  {cppblk_off[lo]= cppblk_off[bfridx] - soff[bfridx] + soff[lo];
  }
 }

/* substitute for the average or worse cases, framing sequence positions
	calculated from the best fit straight lines */
 for (lo=0; lo<num_sc_t; lo++)
 {for (j=0; j<num_frsq; j++)
  {if (fridxbetter[lo][j]== 0x00)
   fridxadj[lo][j]= (int)(floor((double)(cppblk_off[lo] + cppblk_ft*(float)(j)) + 0.5));
  }
 } /* end of for lo */

#if (dbgl>30)
 fprintf(stdout,"Framing Sequence Positions\n");
 fprintf(stdout," framing   --------------------- sub-carrier index------------------\n");
 fprintf(stdout,"  index    0       1       2       3       4       5       6       7\n");
 for (j=0; j<num_frsq; j++)
 {fprintf(stdout,"   %3d  ", j);
  for (lo=0; lo<num_sc_t; lo++)
   fprintf(stdout,"%6d  ", fridxadj[lo][j]);
  fprintf(stdout,"\n");
 }
#endif

#if (dbgl>30)
 fprintf(stdout,"Framing Sequence Offsets\n");
 fprintf(stdout,"framing ------------------------- sub-carrier index ---------------------------\n");
 fprintf(stdout," index      0       1        2        3        4        5        6        7\n");
 for (j=0; j<num_frsq; j++)
 {fprintf(stdout,"  %3d  ", j);
  for (lo=0; lo<num_sc_t; lo++)
   fprintf(stdout," %7.3f ", froffphu[lo][j]);
  fprintf(stdout,"\n");
 }
#endif

#if(dbgl>40)
/* for testing purposes dump demodulated signals */
 for (lo=0; lo<num_sc_t; lo++)
 {otf=fopen(lodif7nfn[lo], "wt");
  if (otf==NULL)
  {fprintf(stderr,"Could not open file %s\n", lodif7nfn[lo]);
   goto exk;
  }
#if(dbgl <50)
  for (si=fridxadj[2][0]-10*cppsym; si<dmd_lidx; si++)  fprintf(otf,"%8.4f\n", dif7n[lo][si]);
#else
  for (si=fridxadj[2][0]-10*cppsym; si<dmd_lidx; si++)  fprintf(otf,"%8.4f\n", dmamp[lo][si]);
#endif /* dbgl <50 */
  fclose(otf);
 }
#endif /* dbgl>40 */

/* adjust for tuning offset changes */
 for (lo=0; lo<num_sc_t; lo++)
 {ex= fridxadj[lo][0];
  ey= froffphu[lo][0];
  jref=0;
  tostate= (fridxbetter[lo][0] == 0x01 ? 1 : 0);
  for (j=1; j<num_frsq; j++)
  {switch(tostate)
   {case 0:
	if (fridxbetter[lo][j] == 0x01)
	{tostate=1; goto doupdates;}
	break;
    case 1:
	if (fridxbetter[lo][j] == 0x01) goto doupdates;
	else tostate=2;
	break;
    case 2:
	if (fridxbetter[lo][j] == 0x01)
	{tostate=1; goto doupdates;}
	else
	{if (j== (num_frsq -1)) goto doupdates;}
	break;
   } /* end of switch(tostate) */
   continue;
doupdates:
   lx= fridxadj[lo][j];
   ly= froffphu[lo][j];
   adjslope= (ly-ey)/((float)(lx-ex));	/* don't divide by 0 */
   adjy= ey;
   for (si=0; si< (lx-ex); si++)
   {dif7n[lo][fridxadj[lo][jref]+ si] -= adjy;
    adjy += adjslope;
   } /* end of for si */
   ex= lx;
   ey= ly;
   jref=j;
  } /* end of for j */
 } /* end of for lo */


/* determine error correction mode (# of redundancy symbols per block) */
 for (si=0; si<num_ec_lvls; si++) ec_mode_dist[si]=0.0;

/* check three pair of symbols preceding 1st framing sequence */
 for (si= -6*cppsym; si<0; si+= 2*cppsym)
 {for (lo=0; lo<num_sc_t; lo++)
  {c2phbuf[lo]= dif7n[lo][si+fridxadj[lo][0]];
  }
  di= decode2ed1mode(&c2phbuf[0], &distg, mode1dthld);
  if (di >= 0) ec_mode_dist[di] += 1.0/(distg > 1e-3 ? distg : 1e-3);
  for (lo=0; lo<num_sc_t; lo++)
  {c2phbuf[lo]= dif7n[lo][si+cppsym+fridxadj[lo][0]];
  }
  di= decode2ed2mode(&c2phbuf[0], &distg, mode2dthld);
  if (di >= 0) ec_mode_dist[di] += 1.0/(distg > 1e-3 ? distg : 1e-3);
 } /* end of for si, checking initial mode id symbols */

/* check three pair of symbols following last framing sequence */
 for (si= frame_len*cppsym; si<=(frame_len+4)*cppsym; si+= 2*cppsym)
 {for (lo=0; lo<num_sc_t; lo++)
  {c2phbuf[lo]= dif7n[lo][si+fridxadj[lo][num_frsq-1]];
  }
  di= decode2ed1mode(&c2phbuf[0], &distg, mode1dthld);
  if (di >= 0) ec_mode_dist[di] += 1.0/(distg > 1e-3 ? distg : 1e-3);
  for (lo=0; lo<num_sc_t; lo++)
  {c2phbuf[lo]= dif7n[lo][si+cppsym+fridxadj[lo][num_frsq-1]];
  }
  di= decode2ed2mode(&c2phbuf[0], &distg, mode2dthld);
  if (di >= 0) ec_mode_dist[di] += 1.0/(distg > 1e-3 ? distg : 1e-3);
 } /* end of for si, checking final mode id symbols */


#if(dbgl > 0)
 fprintf(stdout,"Error Correction Mode Determination Weights\n");
 fprintf(stdout,"  Mode Index    Weights\n");
 for (si=0; si<num_ec_lvls; si++)
  fprintf(stdout,"     %2d       %8.4f\n", si, ec_mode_dist[si]);
#endif

 dummyf1=ec_mode_dist[0]; ec_mode=0;
 for (ki=1; ki<num_ec_lvls; ki++)
 {if (ec_mode_dist[ki] > dummyf1) {dummyf1= ec_mode_dist[ki]; ec_mode= ki;}
 }
 iidx_mn1=iidx_mnx[ec_mode];
 k1=kx[ec_mode];
 t1=k1/2;
 bad_blks=0; ERN=0;

/* for each block 
	decode inner code from demodulated phase changes to channel symbol numbers
	decode array of channel symbol numbers
	output erasure and error correction information
	if 1st block
		extract filename, and file length, open file, store
		 remaining data in file
	else
		store all "non-dummy" data in file
*/
#if (dbgl > 0)
/*  fprintf(stdout,"Blk # Code2_Chngs Erasures Corrections Made & Used \n");*/
fprintf(stdout,"Block inner  |      outer code      | outer   |       inner code\n");
fprintf(stdout,"  #   code   |       erasures       | code    |         changes\n");
fprintf(stdout,"      changes|1st   2nd   3rd   4th | changes |    0   1   2   3   4\n");
fprintf(stdout,"----  ------ |---   ---   ---   --- |   ---   |  --- --- --- --- ---\n");
#endif  /* dbgl > 0 */

 for (j=0; j< num_frsq-1; j++)
 {
/* refine symbol positions, based on amplitude variation information:

   1. Fill sym_pos[][] with nominal symbol positions, based on surrounding
      framing sequence positions.  Fill ma1[][] with 0x00, to indicate
      nominal calculated positions.
   2. For all symbol locations with 0x00 in ma1[][], if local maximum of
      cross-correlation between amplitude reference pattern and actual
      amplitude values is greater than thld_amp1, and is within ossi_amp1
      index values of the nominal position, then use the location of this
      local maximum as the position of that symbol, and set the
      corresponding location in ma1[][] to 0x01.
   3. Extend, form both ends, each run of 0x01 values in ma1[][].  A run 
      is extended by adjacent local maximum cross-correlation values being
      > thld_amp2 and within ossi_amp2 locations of the, f_slope adjusted,
      position from the previous symbol's position.  Set the corresponding
      locations in ma1[][] to 0x02, for those symbols whose position was
      changed in this step.
   4. Replace symbol positions corresponding to remaining runs of 0x00 in
      ma1[][] with linearly interpolated values, based on the positions 
      corresponding to the non-0x00 entries in ma1[][] adjacent to the
      runs of 0x00.  Leave any end symbol positions with 0x00 values in
      ma1[][] as they are.  Set locations in ma1[][] corresponding to
      symbol positions changed in this step to 0x03.
*/

/* 1. Fill sym_pos[][] with nominal symbol positions, based on surrounding
      framing sequence positions.  Fill ma1[][] with 0x00, to indicate
      nominal calculated positions. */

  for (lo=0; lo<num_sc_t; lo++)
  {sym_pos[lo][0]= (double)(fridxadj[lo][j]);
   dcidx_step[lo]= ((double)(fridxadj[lo][j+1]) - sym_pos[lo][0])/ (double)(n1 + frame_len);
   for (si=1; si< n1+frame_len; si++) 
   {sym_pos[lo][si]= sym_pos[lo][si-1] + dcidx_step[lo];}
   for (si=0; si< n1+frame_len; si++) ma1[lo][si]=0x00;
  }

#if(dbgl > 35)
fprintf(stderr,"#BLOCK: %2d  Step 1\n", j);
for (si=0; si<n1+frame_len; si++)
{for (lo=0; lo<num_sc_t; lo++) fprintf(stderr,"%7.0f ", sym_pos[lo][si]);
 for (lo=0; lo<num_sc_t; lo++) fprintf(stderr,"%1d ", (int)(ma1[lo][si]));
 fprintf(stderr,"\n");
}
#endif
/* 2. For all symbol locations with 0x00 in ma1[][], if local maximum of
      cross-correlation between amplitude reference pattern and actual
      amplitude values is greater thanthld_amp1, and is within ossi_amp1
      index values of the nominal position, then use the location of this
      local maximum as the position of that symbol, and set the
      corresponding location in ma1[][] to 0x01. */

  
  for (lo=0; lo<num_sc_t; lo++)
  {for (si=0; si< n1+frame_len; si++)
   {spni=(int)(floor(sym_pos[lo][si]+0.5));
    sas= find_local_max_amp(&dmamp[lo][spni], ossi_amp1, &lm_off_amp, &lm_value_amp);
    if ((sas==1) && (lm_value_amp > thld_amp1))
    {sym_pos[lo][si] += lm_off_amp;
     ma1[lo][si]= 0x01;
    }
   } /* end of for si */
  } /* end of for lo */

#if(dbgl > 35)
fprintf(stderr,"#BLOCK: %2d  Step 2\n", j);
for (si=0; si<n1+frame_len; si++)
{for (lo=0; lo<num_sc_t; lo++) fprintf(stderr,"%7.0f ", sym_pos[lo][si]);
 for (lo=0; lo<num_sc_t; lo++) fprintf(stderr,"%1d ", (int)(ma1[lo][si]));
 fprintf(stderr,"\n");
}
#endif
/* 3. Extend, form both ends, each run of 0x01 values in ma1[][].  A run 
      is extended by adjacent local maximum cross-correlation values being
      > thld_amp2 and within ossi_amp2 locations of the, f_slope adjusted,
      position from the previous symbol's position.  Set the corresponding
      locations in ma1[][] to 0x02, for those symbols whose position was
      changed in this step. */
  for (lo=0; lo<num_sc_t; lo++)
  {ar_state=0; aridx=0;
   for (si=0; si< n1+frame_len; si++)
   {switch(ar_state)
    {case 0:
      if (ma1[lo][aridx]==0x00) {arzdx= aridx; ar_state=1;}
      else {arldx= aridx; ar_state=2;}
      break;
     case 1:
      if (ma1[lo][++aridx]==0x01) ar_state=3;
      break;
     case 2:
      if (ma1[lo][++aridx]==0x00) ar_state=4;
      break;
     case 3:
      for (arudx= aridx-1; arudx >= arzdx; arudx--)
      {spni=(int)(floor(sym_pos[lo][arudx]+0.5));
       sas= find_local_max_amp(&dmamp[lo][spni], ossi_amp2, &lm_off_amp, &lm_value_amp);
       if ((sas==1) && (lm_value_amp > thld_amp2))
       {sym_pos[lo][arudx] += lm_off_amp;
        ma1[lo][arudx]= 0x02;
       }
       else break;
      } /* end of for arudx */
      arldx= aridx;
      ar_state=2;
      break;
     case 4:
      spni=(int)(floor(sym_pos[lo][aridx]+0.5));
       sas= find_local_max_amp(&dmamp[lo][spni], ossi_amp2, &lm_off_amp, &lm_value_amp);
       if ((sas==1) && (lm_value_amp > thld_amp2))
       {sym_pos[lo][aridx] += lm_off_amp;
        ma1[lo][aridx]= 0x02;
        if (aridx < n1+frame_len -1) ar_state=2; else ar_state=6;
       }
       else
       {if (aridx < (n1+frame_len -1)) {arzdx= aridx; ar_state=5;}
        else ar_state=6;
       }
       break;
      case 5:
       if (ma1[lo][++aridx]==0x00)
       {if (aridx < (n1+frame_len-1)) ar_state =1; else ar_state=6;}
       else ar_state=3;
       break;
      case 6:
       break;
      default:
       ermsg(ERN=14);
    } /* end of switch(ar_state) */
   } /* end of for si */
  } /* end of for lo */

#if(dbgl > 35)
fprintf(stderr,"#BLOCK: %2d  Step 3\n", j);
for (si=0; si<n1+frame_len; si++)
{for (lo=0; lo<num_sc_t; lo++) fprintf(stderr,"%7.0f ", sym_pos[lo][si]);
 for (lo=0; lo<num_sc_t; lo++) fprintf(stderr,"%1d ", (int)(ma1[lo][si]));
 fprintf(stderr,"\n");
}
#endif
/* 4. Replace symbol positions corresponding to remaining runs of 0x00 in
      ma1[][] with linearly interpolated values, based on the positions 
      corresponding to the non-0x00 entries in ma1[][] adjacent to the
      runs of 0x00.  Leave any end symbol positions with 0x00 values in
      ma1[][] as they are.  Set locations in ma1[][] corresponding to
      symbol positions changed in this step to 0x03. */
  for (lo=0; lo < num_sc_t; lo++)
  {ar_state=0; aridx=0; nzidxl=0; nzidxr=0;
   for (si=0; si<(n1+frame_len); si++)
   {switch(ar_state)
    {case 0:
      if (ma1[lo][aridx]==0x00) ar_state=2; else ar_state=1;
      break;
     case 1:
      if (ma1[lo][++aridx]==0x00)
      {nzidxl= aridx-1; ar_state=2;}
      else
      {if (si==(n1+frame_len-1)) ar_state=3;}
      break;
     case 2:
      if ((ma1[lo][++aridx]!=0x00) || (si==(n1+frame_len-1)))
      {nzidxr= aridx;
       dummyd1= sym_pos[lo][nzidxr] - sym_pos[lo][nzidxl];
       ar_step0= dummyd1/(double)(nzidxr-nzidxl);
       for (ii0=(nzidxl+1); ii0< nzidxr; ii0++)
       {sym_pos[lo][ii0]= sym_pos[lo][ii0-1] + ar_step0;
        ma1[lo][ii0]= 0x03;
       }
       if (si==(n1+frame_len-1)) ar_state=3; else ar_state=1;
      }
      break;
     case 3:
      break;
     default:
       ermsg(ERN=15);
    } /* end of switch(ar_state) */
   } /* end of for si */
  } /* end of for lo */

#if(dbgl > 35)
fprintf(stderr,"#BLOCK: %2d  Step 4\n", j);
for (si=0; si<n1+frame_len; si++)
{for (lo=0; lo<num_sc_t; lo++) fprintf(stderr,"%7.0f ", sym_pos[lo][si]);
 for (lo=0; lo<num_sc_t; lo++) fprintf(stderr,"%1d ", (int)(ma1[lo][si]));
 fprintf(stderr,"\n");
}
#endif
  eracnt1=0; scca=0; eracnt1s=0; eracnt1ss=0; eracnt1sss=0; eracnt1ssss=0;
  for (si=0; si<num_sc_t; si++) sccx[si]=0;
  for (si= 0; si< n1; si++)
  {for (lo=0; lo<num_sc_t; lo++)
   {c2phbuf[lo]= dif7n[lo][(int)(floor(sym_pos[lo][si+frame_len]+0.5))];
   } /* end of for lo */


   di= decode2wne2(&c2phbuf[0], &scc, &ob_cnt, &ob_mask);
#if (dbgl > 50)
   fprintf(stdout,"%3d: %3d_%1d  ", si, di, scc);
   for (lo=0; lo<num_sc_t; lo++) fprintf(stdout,"%5.2f  ", c2phbuf[lo]);
#endif
   if ((di<0) || (scc > 0))
   {if(ob_cnt > num_sc_r)	/* too many out of bounds subcarrier values */
    {erasi1[0][eracnt1++]= iidx_mx1-si;	/* hard erasure */
     if (di >= 0) cb[0][iidx_mx1-si]=(unsigned int)(di);
    }
    else
    {if (di < 0)
     {switch(ob_cnt)
      {case 0:
       case 1: di_last= decode2ed2w1e(&c2phbuf[0], &distg, 
		(si<(n1-iidx_mn1)?MCWIL:n1), &erag, &scc_last, ob_mask, &bec_mask);
               break;
       case 2: di_last= decode2ed2w2e(&c2phbuf[0], &distg, 
		(si<(n1-iidx_mn1)?MCWIL:n1), &erag, &scc_last, ob_mask, &bec_mask);
               break;
       case 3: di_last= decode2ed2w3e(&c2phbuf[0], &distg, 
		(si<(n1-iidx_mn1)?MCWIL:n1), &erag, &scc_last, ob_mask, &bec_mask);
               break;
       case 4: di_last= decode2ed2w4e(&c2phbuf[0], &distg, 
		(si<(n1-iidx_mn1)?MCWIL:n1), &erag, &scc_last, ob_mask, &bec_mask);
               if(erag==1)
               {erasi1[0][eracnt1++]= iidx_mx1-si; /* hard erasure */}
               else
               {erasi1s[0][eracnt1s++]= iidx_mx1-si; /* soft erasure */}
               if (di_last >= 0) cb[0][iidx_mx1-si]=(unsigned int)(di_last);
               break;
      } /* end of switch ob_cnt */
     } /* end of di < 0 */
     else
     {switch(scc)	/* # of valid cases depends on amount of redundancy in inner code */
      {case 1: di_last= decode2ed2w2e(&c2phbuf[0], &distg, 
	         (si<(n1-iidx_mn1)?MCWIL:n1), &erag, &scc_last, ob_mask, &bec_mask);
               break;
       case 2: di_last= decode2ed2w3e(&c2phbuf[0], &distg, 
	         (si<(n1-iidx_mn1)?MCWIL:n1), &erag, &scc_last, ob_mask, &bec_mask);
               break;
      } /* end of switch scc */
     } /* end of di >= 0 */
     while((scc_last < num_sc_r) && (di != di_last))
     {di= di_last; /* save previous value */
      switch(scc_last)
      {case 1: di_last= decode2ed2w2e(&c2phbuf[0], &distg, 
		(si<(n1-iidx_mn1)?MCWIL:n1), &erag, &scc_last, ob_mask, &bec_mask);
               break;
       case 2: di_last= decode2ed2w3e(&c2phbuf[0], &distg, 
		(si<(n1-iidx_mn1)?MCWIL:n1), &erag, &scc_last, ob_mask, &bec_mask);
               break;
       case 3: di_last= decode2ed2w4e(&c2phbuf[0], &distg, 
		(si<(n1-iidx_mn1)?MCWIL:n1), &erag, &scc_last, ob_mask, &bec_mask);
               if (di != di_last)
               {if (erag==1) erasi1[0][eracnt1++]= iidx_mx1 -si; /* hard erasure */
                else erasi1ss[0][eracnt1ss++]= iidx_mx1-si; /* soft erasure */
                if (di_last >= 0) cb[0][iidx_mx1-si]=(unsigned int)(di_last);
               }
               else /* count this as successful decode with 3 erasures */
               {e1distgsss[0][eracnt1sss]= distg;
                erasi1sss[0][eracnt1sss++]= iidx_mx1-si;
               }
      } /* end of switch scc_last */
     } /* end of while */
     if (di == di_last)
     {/* start of addition of 21 Nov 2001 */
      if (scc_last==2)
      {distg_last= distg;
       di_last= decode2ed2w3e(&c2phbuf[0], &distg, 
		(si<(n1-iidx_mn1)?MCWIL:n1), &erag, &scc_last, ob_mask, &bec_mask);
       if (di == di_last)  
       {cb[0][iidx_mx1-si]=(unsigned int)(di);
        scca += scc_last -2;
        sccx[scc_last-2]++;
       }
       else /* di != di_last */
       {di=di_last;
        di_last= decode2ed2w4e(&c2phbuf[0], &distg, 
		(si<(n1-iidx_mn1)?MCWIL:n1), &erag, &scc_last, ob_mask, &bec_mask);
        if (di != di_last)
        {if (erag==1) erasi1[0][eracnt1++]= iidx_mx1 -si; /* hard erasure */
         else erasi1ss[0][eracnt1ss++]= iidx_mx1-si; /* soft erasure */
         if (di_last >= 0) cb[0][iidx_mx1-si]=(unsigned int)(di_last);
        }
        else /* count this as successful decode with 3 erasures */
        {e1distgsss[0][eracnt1sss]= distg;
         erasi1sss[0][eracnt1sss++]= iidx_mx1-si;
         cb[0][iidx_mx1-si]=(unsigned int)(di);
         scca += scc_last -1;
         sccx[scc_last-1]++;
        }
       }
      }/* end of scc_last == 2 */
      else /* end of addition of 21 Nov 2001, except for braces enclosing next three lines */
      {cb[0][iidx_mx1-si]=(unsigned int)(di);
       scca += scc_last -1;
       sccx[scc_last-1]++;
      }
     } /* end of di == di_last  and scc_last != 2 */
#if (dbgl > 50)
   fprintf(stdout,"%3d_%1d_%1d d:%5.2f  ", di_last, scc_last-1, erag, distg);
#endif
    } /* end of not too many out of bounds subcarrier values */
   } /* end of di<0 || scc>0 */
   else	/* initial decode2 result >= 0 AND no subcarriers changed */
   {cb[0][iidx_mx1-si]= (unsigned int)(di);
    scca += scc;
    sccx[scc]++;
   }
#if (dbgl > 50)
   fprintf(stdout,"\n");
#endif	/* dbgl >50 */
  } /* end of for si */  

/* merge raw blocks before decoding */
  fbtow=0;
  cbidx=iidx_mx1;
  while(cbidx >= iidx_mn1)
  {cbmr[fbtow++]= cb[0][cbidx--];
  }
  hec= eracnt1; sec=eracnt1s;	/* save for diognostic output */
  ssec=eracnt1ss; sssec=eracnt1sss; ssssec=eracnt1ssss;

/* use hard erasures first, then soft erasures, until maximum number
	of allowed erasures are used up, or all hard and soft erasures
	are used up
	Soft erasures should really be sorted first, so that "most likely
	to be wrong" symbols are erased first */

  ir=0;
  while((eracnt1<k1) && (ir<eracnt1s))
  {erasi1[0][eracnt1++]= erasi1s[0][ir++];
  }
  ir=0;
  while((eracnt1<(k1- era_rsv[ec_mode])) && (ir<eracnt1ss))
  {erasi1[0][eracnt1++]= erasi1ss[0][ir++];
  }
  if ((eracnt1 + eracnt1sss) > (k1-era_rsv[ec_mode]))	/* sort */
  {if(eracnt1sss < 2) goto sss;
   for (ji=0; ji<5; ji++)
   {if(ssorts[ji] > eracnt1sss) {lsi= ssorts[ji-1]; goto sset;}
   }
   lsi=ssorts[4];
sset:
   do
   {lsi /=3;
    for (ji= lsi; ji< eracnt1sss; ji++)
    {cx= e1distgsss[0][ji];
     ki= erasi1sss[0][ji];
     ir= ji;
     while (e1distgsss[0][ir-lsi] < cx )
     {e1distgsss[0][ir] = e1distgsss[0][ir-lsi];
      erasi1sss[0][ir] = erasi1sss[0][ir-lsi];
      ir -= lsi;
      if (ir < lsi) break;
     } /* end of while */
     e1distgsss[0][ir] = cx;
     erasi1sss[0][ir] = ki;
    }
   } while(lsi > 0);
  }
sss:
  ir=0;
  while((eracnt1<(k1- era_rsv[ec_mode])) && (ir<eracnt1sss))
  {erasi1[0][eracnt1++]= erasi1sss[0][ir++];
  }
#if (dbgl > 50)
  for (si=0; si<=iidx_mx1; si++) fprintf(stdout,"%3d ", cb[0][si]);
  fprintf(stdout,"\nerasure indexes: ");
  for (si=0; si<eracnt1; si++) fprintf(stdout,"%3d ", erasi1[0][si]);
  fprintf(stdout,"\n");
#endif	/* dbgl > 50 */

/* decode code block [0] */
  eccnt1[0]= decode1(&cb[0][0], eracnt1, &erasi1[0][0]);

#if (dbgl > 50)
  fprintf(stdout,"outer code correction count: %3d\n", eccnt1[0]);
  for (si=0; si<=iidx_mx1; si++) fprintf(stdout,"%3d ", cb[0][si]);
  fprintf(stdout,"\n");
#endif	/* dbgl > 50 */

/* choose which merged blocks to use */
  eci= k1-eracnt1; cbmerged= &cbmd[0];
  if (eci<0) cbmerged= &cbmr[0];
  else
  {if ((eci/2 + eracnt1) < eccnt1[0]) cbmerged= &cbmr[0];
  }

  if(cbmerged == &cbmd[0])
  {
/* merge decoded blocks */
   fbtow=0;
   cbidx=iidx_mx1;
   while(cbidx >= iidx_mn1)
   {cbmd[fbtow++]= cb[0][cbidx--];
   }
  }

  if(cbmerged == &cbmd[0]) ilc= '+'; else ilc= '-';
#if (dbgl > 0)
  fprintf(stdout,"%3u    %4u   %3uh  %3us  %3ut  %3uf    %3d ",
	 j, scca, hec, sec, ssec, sssec,  eccnt1[0]);
  fprintf(stdout," || ");
  for (si=0; si<=num_sc_r; si++) fprintf(stdout," %3d", sccx[si]);
#endif  /* dbgl > 0 */

  if (j==0)	/* 1st block is special case */
/* extract header
	get file name, and file length
	open file, write data from first block
*/
  {cbidxm=0 ;
   if (cbmerged[cbidxm++] != soh)  {ermsg(ERN=19);}
   fnlen= cbmerged[cbidxm++] ; 
   if (fnlen > fnlmax) {ermsg(ERN=20); fnlen=fnlmax;}
   if (fnlen < 1)  {ermsg(ERN=21);}
   flenb= cbmerged[cbidxm++];
   flenb+= cbmerged[cbidxm++]<<8; 
   if (flenb > flmax)  {ermsg(ERN=22);}
   fbw=0; 

/* get file name */
   if (ilc == '+')	/* 1st block decoded properly */
   {for (si=0; si<fnlen; si++) 
    {tmpi= cbmerged[cbidxm++]; 
     if (tmpi== sod) {cbidxm--; goto eofn;}
     if (tmpi > max_usr_sym) {ermsg(ERN=23);}
     else fnbbuf[si]= (fnbuf[si]= (unsigned char)(tmpi));
    }
eofn:
    fnbbuf[si+3]= (fnbuf[si]= 0x00);  	/* for string termination */
    fnlen= si;		/* 25oct2002  in case of uncorrected fnlen */
   }
   else			/* 1st block failed to decode properly */
   {strcpy(&fnbuf[0], &default_out_name[0]);
    strcpy(&fnbbuf[0], &default_out_name[0]);
    fnlen= strlen(&default_out_name[0]);
    fnbbuf[fnlen+3]= (fnbuf[fnlen]= 0x00); /* for string termination */
	/* in case the sod symbol is correct */
    for (si=0; si<fnlmax; si++)
    {tmpi= cbmerged[cbidxm++];
     if (tmpi== sod) {cbidxm--; goto fsod;}
    }
fsod: ;
   }

/* for now, overwrite any existing file of the same name */
   otf=fopen(&fnbuf[0], "wb");
   if (otf==NULL) {ERN=8; goto exl;}

   if (cbmerged[cbidxm++] != sod) {ermsg(ERN=24);}
   fbtow=0;
   while((cbidxm < (n1-iidx_mn1)) &&(fbtow<flenb))
   {tmpi= cbmerged[cbidxm++];
    if (tmpi > max_usr_sym) {ermsg(ERN=25); fotbuf[fbtow++]= 0x00; }
    else fotbuf[fbtow++]= (unsigned char)(tmpi);
   } /* end of while cbidx >= iidx_mn1 */
   fwr= fwrite(fotbuf, sizeof(unsigned char), fbtow, otf);
   if (fwr!=fbtow) {ermsg2(ERN=26);}
   if (ilc != '+') {bad_blks++; blk0_ok=0;}
   else
   {wcbr= write_code_block(&fotbuf[0], fbtow, &fnbbuf[0], fnlen, j);
    if(wcbr<0) ermsg2(ERN=28);
    blk0_ok=1;
   }
   fbw+= fbtow;
   if (fbw >= flenb)
   {fclose(otf);
#if(dbgl > 0)
 fprintf(stdout,"   %c\n", ilc);
#endif
    goto nrmex;
   }
  } /* end of j==0 */
  else

/* process remaining blocks */
/* if block 0 was not decoded correctly, then flenb can not be trusted */
  {cbidxm=0;  fbtow=0;
   if (blk0_ok==1) eotff= (( (fbw+fbtow)<flenb) ? 0: 1);
   else eotff=0;
   while ((cbidxm < (n1-iidx_mn1)) && (eotff==0))
   {tmpi= cbmerged[cbidxm++];
    if (blk0_ok==1)
    {if (tmpi > max_usr_sym) {ermsg(ERN=25); fotbuf[fbtow++]=0x00;}
     else fotbuf[fbtow++]= (unsigned char)(tmpi);
     if ((fbw+fbtow)>= flenb) eotff=1;
    }
    else
    {if ((tmpi==dummy1) || ( (cbidxm==(n1-iidx_mn1+1)) && (j==(num_frsq-2)) ))
     {eotff=1;}
     else
     {if (tmpi > max_usr_sym) {ermsg(ERN=25); fotbuf[fbtow++]=0x00;}
      else fotbuf[fbtow++]= (unsigned char)(tmpi);
     }
    }
   } /* end of while cbidx >= iidx_mn1 */
   fwr= fwrite(fotbuf, sizeof(unsigned char), fbtow, otf);
   if (fwr!=fbtow) {ermsg2(ERN=26); }
   if (ilc != '+') bad_blks++;
   else
   {wcbr= write_code_block(&fotbuf[0], fbtow, &fnbbuf[0], fnlen, j);
    if(wcbr<0) ermsg2(ERN=28);
   }
   fbw+= fbtow;
   if ( (fbw >= flenb) || (eotff==1) )
   {fclose(otf);
#if(dbgl > 0)
 fprintf(stdout,"   %c\n", ilc);
#endif
    if (j!=(num_frsq-2))  ERN=27;
    goto nrmex;
   }
  }
#if(dbgl > 0)
 fprintf(stdout,"   %c\n", ilc);
#endif

 } /* end of for j (over each code block) */


nrmex:			/* normal exit */
 fprintf(stdout,"Decoded result is file: %s\n", &fnbuf[0]);
 
 fprintf(stdout,"%3d bad blocks detected, out of %3d total blocks.\n", bad_blks, num_frsq-1);
 if(j==(num_frsq-2)) ERN=0;
 if (bad_blks==0)
 {for (di=0; di<(num_frsq-1); di++)
  {sprintf(&fnbbuf[fnlen], "%03u", di);
   remove(&fnbbuf[0]);	/* remove separate files for each code block */
  }
 }

exl:
exk: for (j=0; j<num_sc_t; j++) {free(dif7n[j]); free(dmamp[j]);}
exj: remove(&scratch[0]); remove(&scratch2[0]);
exf3: for (j=0; j<num_sc_t; j++) free(froffphu[j]);
exf2: for (j=0; j<num_sc_t; j++) free(fridxbetter[j]);
exf1: for (j=0; j<num_sc_t; j++) free(fridxadj[j]);
exf:
ex14:
 fclose(inf);
 if (ERN !=0) fprintf(stderr,"%s\n", erexmsg[ERN]);
 if (ERN==0)
 {if (bad_blks !=0) exit(bad_blks) ; else exit(0);
 }
 else exit(-ERN);
}
