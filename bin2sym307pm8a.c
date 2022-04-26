/*
    bin2sym307pm8a produces encoded symbol number representing input file
    Copyright (C) 2001, 2003  Barry Sanderson

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
/* Rev 2003 May 29	If the compiler preprocessor symbol "PATH_DELIM_BACK"
  	is defined, then the '\' character will be used  as the path separator
	character; else, if the preprocessor symbol "PATH_DELIM_COLON" is
	defined, then the ':' character will be used as the path separator
	character; otherwise, the '/' character will be used as the path
	separator character.
			No changes to code, only this comment added.
	Barry Sanderson 2003 May 29 */
/* Rev 2003 May 14	Changed so that input file name can now be a
	path name, with just the basename being encoded in the output
	of this program; any leading directory names are removed from
	what is to be sent as the file name.
	Barry Sanderson 2003 May 14 */
/* Rev 25 Nov  2001	version 1.1.0, for initial public release */
/* This program requires a filename and a numerical value from the set
	{10, 20, 40, 70} as its first two arguments, in that order.  The
	output of this program is written to stdout, and is a sequence of
	symbol numbers, representing the encoded input file in the message
	format defined below. */
/* Message format:
	symbol	# of symbols		coments
	lead_mod     74		used to allow receiver to measure
				transmitter's sample clock rate
	id1xxx	      1		1st of two symbols used to specify mode
	id2xxx	      1		2nd of two symbols used to specify mode
	  .
	  .	      2		total of 3 [id1xxx id2xxx] pairs
	  .
	frame1	      1
	frame2	      1		4 symbol frame sequence delimits all code blocks
	frame3	      1
	frame4	      1
	soh	      1		start of header
	nnnnn	      1		length of file name, binary
	fffff	      2		length of file, binary LSB first
	xxxxx	   nnnnn	symbols for file name
	sod	      1		start of data
	sssss	   fffff	symbols for data in file
	  .
	  .			after every data block
	  .			the following 4 symbol sequence will be
				inserted:
				frame1 frame2 frame3 frame4
	ddddd	     yy		dummy symbols to fill last block, as
				required
	id1xxx	      1		1st of two symbols used to specify mode
	id2xxx	      1		2nd of two symbols used to specify mode
	  .
	  .	      2		total of 3 [id1xxx id2xxx] pairs
	  .
	lead_mod     37		trailer to mark end of file		
*/


#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>

#include "sym_code.h"
#include "rs8-4.h"
#include "mode_def.h"
#include "rs306-xxx-rgm.h"

unsigned int	rsmode, bbs;
unsigned int	kxx[num_modes]={10, 20, 40, 70};

sym_type1		frame_buf[frame_len]=
	{frame1, frame2, frame3, frame4};
sym_type1		framer_buf[frame_len]=
	{frame1r, frame2r, frame3r, frame4r};

sym_type1		mode_buf[mode_len+mode_len];
sym_type1		id1xx[num_modes]={id1x05, id1x10, id1x20, id1x40};
sym_type1		id2xx[num_modes]={id2x05, id2x10, id2x20, id2x40};

sym_typei	cb[cbl1];	/* code block buffer */
int		cbidx;
int		fnlenmod, ERN;
sym_type_ch	frame_buf_ch[frame_len];
sym_type_ch	mode_buf_ch[mode_len+mode_len];
sym_type_ch	leader_buf_ch[leader_len];
sym_type_ch	sym_buf_ot[cbl1];	/* output symbol buffer */

int	write_frame(void)
{int	fwr;
 fwr= fwrite(frame_buf_ch, bytes_per_ch_sym, frame_len, stdout);
 if (fwr!=frame_len) return -1;
 return 0;
}

 char	noerr[]="Successful termination";
 char   ins[]=
"Useage:\n  bin2sym307pm8a bin_input_file %_of_redundancy_symbols(10|20|40|70) [ > output_file ]\n";
 char	ukner[]="Unknown error\n";
 char	noif[]="Couldn't open input file\n";
 char	noof[]="call to stat() failed\n";
 char	em05[]="Invalid symbol (out of range)\n";
 char	em06[]="";
 char	em07[]="";
 char	em08[]="";
 char	em09[]="Input file is too long\n";
 char	em10[]="Memory allocation failed\n";
 char	em11[]="Error reading input file\n";
 char	em12[]="Error writing leader\n";
 char	em13[]="Error writing body of file\n";
 char	em14[]="";
 char	em15[]="";
 char	em16[]="";
 char	em17[]="an unsupported number of redundancy symbols was requested\n";
 char	em18[]="";
 char	em19[]="Error writing frame sequence\n";
 
 char		*erexmsg[]={noerr, ins, ukner, noif, noof, em05, em06, em07,
			 em08, em09, em10, em11, em12, em13, em14, em15, em16,
			 em17, em18, em19};

sym_type_ch	mk_ch_sym_num(sym_typei ra[], sym_typei ia[])
{/* calculate and return the integer corresponding to the
	concatenation of the base 9 values in ra[] and ia[]
*/
 int	si ;
 sym_type_ch	acum;

 acum=0;
 for (si=num_sc_i-1; si>=0; si--) acum = (acum+ ia[si])*N2;
 for (si=num_sc_r-1; si>0; si--) acum = (acum+ ra[si])*N2;
 return (acum+ ra[0]);
}

int	encode1(int *extra, int *info)
{/* extra	points to array where redundancy symbols of code word will
			be stored,
    info	points to array where information symbols have been stored

    returns:	 0 if no errors were encountered
		-1 if an illegal symbol is found in the "info" array */
 register int	si, di;
 int		tmpi;

 for (si=0; si<(n1-2*t1); si++) if ((info[si]<0) || (info[si]>n1)) return -1;
 for (di=0; di<(2*t1); di++)
 {tmpi=0;
  for (si=0; si<(n1-2*t1); si++) tmpi+= (info[si]*rgm1xx[rsmode][di][si])%N1;
  extra[di]= tmpi%N1;
 }
 return 0;
}

int do_block(void)
{/* encode1 contents of cb[],
    encode2 each symbol in cb[]
    write channel symbol numbers to sym_buf_ot,
    fwrite result to stdout 

    return 0 on success, -1 on failure */
 register int	j, p;
 int		i, er;

 
 i=encode1(&cb[0], &cb[iidx_mn1]);
 if(i!=0){ERN=5 ; return -1;}
 for (j=0; j< cbl1; j++)
 {for (p=0; p<num_sc_i; p++) ch_sym_buf_i[p]= sym_to_phi[cb[iidx_mx1-j]][p];
  er= encode2(&ch_sym_buf_r[0], &ch_sym_buf_i[0]);
  if(er !=0) {ERN=13; return -1 ;}
  sym_buf_ot[j]= mk_ch_sym_num(&ch_sym_buf_r[0], &ch_sym_buf_i[0]);
 } /* end of for j */
 i= fwrite(sym_buf_ot, bytes_per_ch_sym, cbl1, stdout);
 if(i!=cbl1) {ERN=13; return -1;}
 return 0;
}



int main(int argc, char *argv[], char *env[])
{

 register int	si, di;
 FILE		*inf;
 struct stat	statb;
 int		flen, bytes_used, bytes_left;
 int		statr, fnlen, ir, fwr, tl, sl, bl;
 int		bcnt;
 unsigned char	fnbuf[fnlmax+2], *fbuf;
 sym_type_ch	csn;

/* following is definition of character used as a path delimiter and
	a pointer to the basename, added 2003 May 14 */
#ifdef PATH_DELIM_BACK
 int	path_delim='\\';
#else
#ifdef PATH_DELIM_COLON
 int	path_delim=':';
#else
 int	path_delim='/';
#endif
#endif

 char		*basename;

 if (argc < 3) {ERN=1; goto exa;}
 statr= stat(argv[1],  &statb);
 if (statr!=0) {ERN=4 ; goto exa;}
 flen=(unsigned long)(statb.st_size);
 if (flen > flmax) {ERN=9; goto exa;}

/* start of modification #1 to get basename from path name, 2003 May 14 */
 basename= strrchr(argv[1], path_delim);
 if(basename==NULL)
  basename=argv[1]; 	/* no path delimeter found */
 else
  basename= &basename[1];	/* skip the last path delimeter */
 fnlen=strlen(basename);
/* end of modification #1 to get basename from path name, 2003 May 14 */

 if (fnlen > fnlmax)
 {fprintf(stderr,"WARNING: filename is being limited to %3u characters\n", fnlmax);
  fnlen=fnlmax;
 }

 si=atoi(argv[2]); di=0;
 while(di<num_modes)
 {if (si != kxx[di]) di++;
  else
  {rsmode=di; k1=kx[di]; t1= k1/2; iidx_mn1= iidx_mnx[rsmode];
  
   for (ir=0; ir< 2*mode_len; ir +=2) mode_buf[ir]=id1xx[di];
   for (ir=1; ir< 2*mode_len; ir +=2) mode_buf[ir]=id2xx[di];
   bbs= n1-k1;
   goto modeok;
  }
 }
 ERN=17; goto exa;

modeok:

/* initialize leader_buf_ch */
 for (si=0; si< num_sc_i; si++) ch_sym_buf_i[si]= sym_to_phi[lead_i][si];
 for (si=0; si< num_sc_r; si++) ch_sym_buf_r[si]= sym_to_phi[lead_r][si];
 csn= mk_ch_sym_num(&ch_sym_buf_r[0], &ch_sym_buf_i[0]);
 for (si=0; si<leader_len; si++) leader_buf_ch[si]=csn;

/* initialize frame_buf_ch */
 for (si=0; si< num_sc_i; si++) ch_sym_buf_i[si]= sym_to_phi[frame1][si];
 for (si=0; si< num_sc_r; si++) ch_sym_buf_r[si]= sym_to_phi[frame1r][si];
 frame_buf_ch[0]= mk_ch_sym_num(&ch_sym_buf_r[0], &ch_sym_buf_i[0]);
 for (si=0; si< num_sc_i; si++) ch_sym_buf_i[si]= sym_to_phi[frame2][si];
 for (si=0; si< num_sc_r; si++) ch_sym_buf_r[si]= sym_to_phi[frame2r][si];
 frame_buf_ch[1]= mk_ch_sym_num(&ch_sym_buf_r[0], &ch_sym_buf_i[0]);
 for (si=0; si< num_sc_i; si++) ch_sym_buf_i[si]= sym_to_phi[frame3][si];
 for (si=0; si< num_sc_r; si++) ch_sym_buf_r[si]= sym_to_phi[frame3r][si];
 frame_buf_ch[2]= mk_ch_sym_num(&ch_sym_buf_r[0], &ch_sym_buf_i[0]);
 for (si=0; si< num_sc_i; si++) ch_sym_buf_i[si]= sym_to_phi[frame4][si];
 for (si=0; si< num_sc_r; si++) ch_sym_buf_r[si]= sym_to_phi[frame4r][si];
 frame_buf_ch[3]= mk_ch_sym_num(&ch_sym_buf_r[0], &ch_sym_buf_i[0]);
 
/* initialize mode_buf_ch */
 for (si=0; si< num_sc_i; si++) ch_sym_buf_i[si]= sym_to_phi[mode_buf[0]][si];
 ir= encode2(&ch_sym_buf_r[0], &ch_sym_buf_i[0]);
 if(ir!=0){ERN=5 ; goto exa;}
 csn= mk_ch_sym_num(&ch_sym_buf_r[0], &ch_sym_buf_i[0]);
 for (si=0; si<mode_len+mode_len; si+=2) mode_buf_ch[si]=csn;
 for (si=0; si< num_sc_i; si++) ch_sym_buf_i[si]= sym_to_phi[mode_buf[1]][si];
 ir= encode2(&ch_sym_buf_r[0], &ch_sym_buf_i[0]);
 if(ir!=0){ERN=5 ; goto exa;}
 csn= mk_ch_sym_num(&ch_sym_buf_r[0], &ch_sym_buf_i[0]);
 for (si=1; si<mode_len+mode_len; si+=2) mode_buf_ch[si]=csn;

/* start of modification #2 to get basename from path name, 2003 May 14 */
 strncpy(fnbuf, basename, fnlen);
/* end of modification #2 to get basename from path name, 2003 May 14 */

 fnbuf[fnlen]=0x00;  fnbuf[fnlen+1]=0x00;
 fbuf=(unsigned char *)malloc((unsigned long)((flen+2)*sizeof(char)));
 if (fbuf==NULL) {ERN=10; goto exa;}
 inf=fopen(argv[1], "rb");
 if (inf==NULL) {ERN=3; goto exb;}
 ir= fread(fbuf, 1, flen, inf);
 if (ir!=flen) {ERN=11 ; goto exc;}
 if (ferror(inf) !=0) {ERN=11 ; goto exc;}

 fbuf[flen]=0x00; fbuf[flen+1]=0x00;

/* write leader */
 fwr= fwrite(leader_buf_ch, bytes_per_ch_sym, leader_len, stdout);
 if (fwr!=leader_len) {ERN=12; goto exc;}
 
/* put mode specifying sequence into output file */
 ir=mode_len+mode_len;
 fwr= fwrite(mode_buf_ch, bytes_per_ch_sym, ir, stdout);
 if (fwr!=ir) {ERN=12; goto exc;}

/* put fixed-length part of header into cb[] */
 fwr= write_frame();
 if (fwr !=0) {ERN=19; goto exc;}
 cbidx= iidx_mx1;
 cb[cbidx--]=(unsigned int)(soh);
 cb[cbidx--]=(unsigned int)(fnlen);
 cb[cbidx--]=(unsigned int)(flen & 0xff);
 cb[cbidx--]=(unsigned int)((flen>>8) & 0xff);

/* put file name into header */
 for (si=0; si<fnlen; si++) cb[cbidx--]=(unsigned int)(fnbuf[si]);
/* put start of data into code block */
 cb[cbidx--]=(unsigned int)(sod);

 bytes_left= flen; bytes_used=0; tl=0; bcnt=0;

/* do bulk of data from file */
 while(bytes_left > bbs)
 {sl= cbidx-iidx_mn1+1;
  bl= sl;
  for (si=0; si<bl; si++) cb[cbidx--]= (unsigned int)(fbuf[bytes_used+si]);
  if (do_block()!=0) goto exc;
  bytes_used +=(unsigned long)bl;
  bytes_left -=(unsigned long)bl;
  cbidx= iidx_mx1;
  bcnt++;
  if(bcnt==blks_per_frame)
  {bcnt=0;
   fwr= write_frame();
   if (fwr!=0) {ERN=19; goto exc;}
  }
 } /* end of while bytes_left > bbs */

/* do rest of data from file */
 while(bytes_left > 0)
 {sl= cbidx-iidx_mn1+1;
  bl=sl;
  if (bytes_left >= bl)
  {for (si=0; si<bl; si++) cb[cbidx--]= (unsigned int)(fbuf[bytes_used+si]);
   if (do_block()!=0) goto exc;
   bytes_used +=bl;
   bytes_left -=bl;
  } /* end of bytes_left >= bl */
  else
  {if (bytes_left > 0)
   {for (si=0; si<bytes_left; si++) cb[cbidx--]= (unsigned int)(fbuf[bytes_used+si]);
    bytes_used += bytes_left;
    bytes_left=0;
    sl= cbidx-iidx_mn1+1;
    for (si=0; si<sl; si++) cb[cbidx--]= (unsigned int)(dummy1);
    tl+= sl;
    if (do_block()!=0) goto exc;
   } /* end of bytes_left > 0 */
  } /* end of bytes_left < bl */
  cbidx= iidx_mx1;
  bcnt++;
  if(bcnt==blks_per_frame)
  {bcnt=0;
   fwr= write_frame();
   if (fwr!=0) {ERN=19; goto exc;}
  }
 } /* end of while bytes_left > 0 */


/* put mode specifying sequence into output file */
 ir=mode_len+mode_len;
 fwr= fwrite(mode_buf_ch, bytes_per_ch_sym, ir, stdout);
 if (fwr!=ir) {ERN=12; goto exc;}

/* do trailer symbols */
 fwr= fwrite(leader_buf_ch, bytes_per_ch_sym, trailer_len, stdout);
 if (fwr!=trailer_len) {ERN=12; goto exc;}



 ERN=0;

exc: fclose(inf);
exb: free(fbuf);
exa: if (ERN != 0) fprintf(stderr, erexmsg[ERN]);
 exit(ERN);
}
