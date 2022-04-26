/*
    post-modulation processing for digital communications Tx file
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
/* Rev 19 May 2003	added: printf of: program name and version string,
	copyright statement and warranty disclamer.
			updated copyright date. */
/* Rev 30 Aug  2002  added: detection of wrong number of bytes written
	to disk file */
/* Rev 23 Apr  2001  version 1.0.0 for initial release */

#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>

#define	in_type		float
#define	out_type	signed short int

#define	wav_hdr_len	44
#define	wav_len1_off	4
#define	wav_len2_off	40
#define	wav_len_diff	(wav_len2_off - wav_len1_off)

struct stat	statb;
int		statr, ir;
unsigned long	flen, fbidx, numel;
int		val_toread, val_read;
FILE		*inf, *otf;
in_type		*fbufin;
out_type	*fbufot;

unsigned char	wav_hdr_buf[wav_hdr_len]={
	0x52, 0x49, 0x46, 0x46, 0x00, 0x00, 0x00, 0x00,
	0x57, 0x41, 0x56, 0x45, 0x66, 0x6d, 0x74, 0x20,
	0x10, 0x00, 0x00, 0x00, 0x01, 0x00, 0x01, 0x00,
	0x11, 0x2b, 0x00, 0x00, 0x22, 0x56, 0x00, 0x00,
	0x02, 0x00, 0x10, 0x00, 0x64, 0x61, 0x74, 0x61,
	0x00, 0x00, 0x00, 0x00};

unsigned int	ERN;

 char	noerr[]="Successful termination";
 char   ins[]=
"Useage: flt-clip-scale-wav in-file out-file clipping_level scale_factor\n";
 char	ukner[]="Unknown error\n";
 char	em03[]="Couldn't open input file\n";
 char	em04[]="Couldn't open output file\n";
 char	em05[]="WRONG number of bytes written to disk file\n";
 char	em06[]="Error in reading body of input file\n";
 char	em07[]="\n";
 char	em08[]="Initial memory allocation failed\n";
 char	em09[]="\n";
 char	em10[]="\n";
 char	em11[]="\n";
 char	em12[]="\n";
 char	em13[]="stat failed\n";
 char		*erexmsg[]={noerr, ins, ukner, em03, em04, em05, em06, em07,
			    em08, em09, em10, em11, em12, em13};

int main(int argc, char *argv[], char *env[])
{

 double		tv, cl;
 in_type	scalef;
 unsigned int	wav_bc1, wav_bc2;
 size_t		val_written;

 unsigned char	verstr[]="2003May19";
 unsigned char	progname[]="flt-clip-scale-wav";
 unsigned char	cww[]="Copyright (C) 2003 Barry Sanderson
There is ABSOLUTELY NO WARRANTY  for flt-clip-scale-wav.
flt-clip-scale-wav is covered by the the GNU General Public License.
The file \"COPYING.txt\" states the conditions under which you may legally:
copy, modify, or distribute flt-clip-scale-wav.";


 printf("\n%s version %s\n\n%s\n", progname, verstr, cww);


 if (argc<5)
 {ERN=1; goto exa;}
 statr= stat(argv[1],  &statb);
 if (statr!=0) {ERN=13 ; goto exa;}
 flen=(unsigned long)(statb.st_size);
 numel= flen/(sizeof(in_type));
 cl=atof(argv[3]);
 scalef=(in_type)(atof(argv[4]));
 ERN=8;
 if ((fbufin=malloc((size_t)((numel)* sizeof(in_type)))) == NULL) goto exa;
 if ((fbufot=malloc((size_t)((numel)* sizeof(out_type) ))) == NULL) goto exb;
 ERN=3;
 inf=fopen(argv[1], "rb");
 if (inf==NULL)
 {fprintf(stderr,"Could not open file %s\n", argv[1]); goto exc;}
 ERN=4;
 otf=fopen(argv[2], "wb");
 if (otf==NULL)
 {fprintf(stderr,"Could not open file %s\n", argv[2]); goto exd;}
 ir= fread(&fbufin[0], sizeof(in_type), numel, inf);
 if (ir!=numel) {ERN=6 ; goto exe;}
 for (fbidx=0; fbidx<numel; fbidx++)
 {tv=fbufin[fbidx];
  if(tv>cl) tv=cl;
  else {if(tv<-cl) tv=-cl;}
  fbufot[fbidx]= (out_type)(scalef*tv);
 }

/* put proper lengths into wav_hdr_buf */
 wav_bc2= (unsigned int)((numel)* sizeof(out_type));
 wav_bc1= (unsigned int)((numel)* sizeof(out_type) + wav_len_diff);

/* write wav file header */
 ERN=5;
 val_written= fwrite(&wav_hdr_buf[0], 1, wav_len1_off, otf);
 if (val_written != (size_t)(wav_len1_off)) goto exe;
 val_written= fwrite(&wav_bc1, sizeof(unsigned int), 1, otf);
 if (val_written != (size_t)(1)) goto exe;
 val_written= fwrite(&wav_hdr_buf[wav_len1_off + sizeof(unsigned int)], 1, wav_len_diff-sizeof(unsigned int), otf);
 if (val_written != (size_t)(wav_len_diff-sizeof(unsigned int))) goto exe;
 val_written= fwrite(&wav_bc2, sizeof(unsigned int), 1, otf);
 if (val_written != (size_t)(1)) goto exe;

 val_written= fwrite(&fbufot[0], sizeof(out_type), numel, otf);
 if (val_written != (size_t)(numel)) goto exe;


 ERN=0;

exe: fclose(otf);
exd: fclose(inf);
exc: free(fbufot);
exb: free(fbufin);
exa: if (ERN != 0) fprintf(stderr, erexmsg[ERN]);
     exit(ERN);
}
