/*
    modulator program for Wyman1x digital communications programs
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
                        changed: error messages written to stdout so that they
        are now written to stderr.
                        updated copyright date. */
/* Rev 30 Aug 2002	added: detection of wrong number of bytes written
        to disk files; and deletion of temporary files */
/* Rev 25 Nov 2001	version 1.1.0, for initial public release */
/* Rev 15 Nov  2001  order of include files changed, for conditional
        defn of PI and PI2 */

#define dbgl 0
#define OS   linux

#include "dcom-t.h"
#include "exfun-t.h"
#include "sym_def.h"
#include <malloc.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#define ALLOC_MEM malloc

#define dum_sub 1400.0 /* dummy subcarrier frequency */
#define dumd    0.001  /* dummy data value */

#define num_freqs     num_sc_t /* number of sub-carrier frequencies */
#define num_phases    9        /* number of phases of each sub-carrier */
#define sym_num_limit 43046721 /* num_phases^num_freqs */

/* following is aray of phase modulation center frequencies */
double pmlo[num_sc_t] = { pmc0, pmc1, pmc2, pmc3, pmc4, pmc5, pmc6, pmc7 };

#define clk_per_sym  90 /* # of clock periods per symbol */
#define clk_per_edge 30 /* # of clock periods consumed by edge */
#define hammlim      1.0

int soff[num_sc_t] = { soff0, soff1, soff2, soff3, soff4, soff5, soff6, soff7 };

#define DS1  ds1t
#define DS2  ds2t
#define DS3  ds3t
#define I40x I40xxPI

double asteps[num_phases][clk_per_sym];
#define init_flag 1
double loai[num_sc_t] = { 0.0, 180.0, 0.0, 0.0, 90.0, 0.0, 180.0, 0.0 };

unsigned int sym_count, cfreq, nfreq, ERN;
signed int   ampi, dcf, dcc;

double astep, ca, astep2, ca2;
double clk_period; /* 11,025 KHz clock frequency */

void set_asteps(double start_val, double end_val, double *steps)
/* this subroutine sets up the steps[] array for the given starting
   and ending values*/
{
  double range, pv, nv;
  int register si;

  range = end_val - start_val;
  astep = hammlim * PI / ((double) (clk_per_edge - 1));
  pv    = start_val;
  if(range == 0.0)
  {
    for(si = 0; si < clk_per_sym; si++)
      steps[si] = 0.0;
  }
  else
  {
    if(range > 0.0)
    {
      ca = PI - PI * hammlim;
      for(si = 0; si < clk_per_edge; si++)
      {
        nv        = start_val + range * (0.5 - 0.5 * cos(ca));
        ca        += astep;
        steps[si] = nv - pv;
        pv        = nv;
      }
      steps[clk_per_edge] = end_val - pv;
      for(si = clk_per_edge + 1; si < clk_per_sym; si++)
        steps[si] = 0.0;
    } /* end of range > 0 */
    else
    {
      ca = PI;
      for(si = 0; si < clk_per_edge; si++)
      {
        nv        = end_val - range * (0.5 - 0.5 * cos(ca));
        ca        += astep;
        steps[si] = nv - pv;
        pv        = nv;
      }
      steps[clk_per_edge] = end_val - pv;
      for(si = clk_per_edge + 1; si < clk_per_sym; si++)
        steps[si] = 0.0;
    } /* end of range < 0 */
  }   /* end of range != 0.0 */
}

ch_sym_type base9exp[num_sc_t]
    = { 1, 9, 81, 729, 6561, 59049, 531441, 4782969 };

ch_sym_type sym_to_phi(ch_sym_type sym_num, int sc)
{ /* return modulation phase change index for given symbol and sub-carrier
   */

  ch_sym_type b9d[num_sc_t], acum;
  int         si;

  si      = num_sc_t - 1;
  acum    = sym_num;
  b9d[si] = sym_num / base9exp[si];
  si--;
  while(si >= sc)
  {
    acum    = acum - b9d[si + 1] * base9exp[si + 1];
    b9d[si] = acum / base9exp[si];
    si--;
  }
  return (b9d[++si]);
}

ch_sym_type cflags, csym, nsym;
FILE       *inf, *otf;
FILE       *lo_db[num_freqs];
FILE       *sb_db[num_freqs];

unsigned char lo0db[] = "lo0.flt";
unsigned char lo1db[] = "lo1.flt";
unsigned char lo2db[] = "lo2.flt";
unsigned char lo3db[] = "lo3.flt";

unsigned char lo4db[] = "lo4.flt";
unsigned char lo5db[] = "lo5.flt";
unsigned char lo6db[] = "lo6.flt";
unsigned char lo7db[] = "lo7.flt";

unsigned char  lo8db[] = "lo8.flt";
unsigned char  lo9db[] = "lo9.flt";
unsigned char *lodbfn[]
    = { lo0db, lo1db, lo2db, lo3db, lo4db, lo5db, lo6db, lo7db, lo8db, lo9db };

unsigned char sb0db[] = "pm0.flt";
unsigned char sb1db[] = "pm1.flt";
unsigned char sb2db[] = "pm2.flt";
unsigned char sb3db[] = "pm3.flt";

unsigned char sb4db[] = "pm4.flt";
unsigned char sb5db[] = "pm5.flt";
unsigned char sb6db[] = "pm6.flt";
unsigned char sb7db[] = "pm7.flt";

unsigned char  sb8db[] = "pm8.flt";
unsigned char  sb9db[] = "pm9.flt";
unsigned char *sbdbfn[]
    = { sb0db, sb1db, sb2db, sb3db, sb4db, sb5db, sb6db, sb7db, sb8db, sb9db };

ftype Sambuf[s0blen], Amp[fap];
ftype Ssbc[fap], Ssbs[fap];
ftype Lo1c[s0blen], Lo1s[s0blen];
ftype Lo1cf[s1blen], Lo1sf[s1blen];
ftype Lo2cc[s1blen], Lo2cs[s1blen], Lo2ss[s1blen], Lo2sc[s1blen];
ftype Lo2cs1[s2blen], Lo2ss1[s2blen];
ftype Lo2cs2[s3blen], Lo2ss2[s3blen], Ampb[s3blen];

long sidx;

#define dblen 256 /* double buffer length */

ftype  db0[dblen], db1[dblen], db2[dblen], db3[dblen];
ftype  db4[dblen], db5[dblen], db6[dblen], db7[dblen];
ftype  db8[dblen], db9[dblen];
ftype *dbb[] = { db0, db1, db2, db3, db4, db5, db6, db7, db8, db9 };

short signed int ampsib[dblen];

char noerr[]  = "Successful termination\n";
char ins[]    = "Invoke via: mod-pm8a input-file output-file[.flt]\n";
char noif[]   = "Couldn't open input file\n";
char noof[]   = "Couldn't open output file\n";
char em04[]   = "wrong number of bytes written to disk file\n";
char ukner[]  = "Unknown error\n";
char nenm[]   = "ALLOC_MEM could not allocate enough memory.\n";
char seeker[] = "fseek failed\n";
char em08[]   = "stat failed\n";
char em09[]   = "Error reading input file\n";
char em10[]   = "Digital filter failure\n";
char em11[]   = "Illegal Symbol Number in input file\n";

char *erexmsg[] = { noerr, ins,    noif, noof, em04, ukner,
                    nenm,  seeker, em08, em09, em10, em11 };

int main(int argc, char *argv[], char *env[])
{
  int register si, di;
  unsigned int lo;
  int          end_of_file;
  int          val_to_read, dum_len, val_to_write;
  size_t       val_read, val_written;

  struct stat   statb;
  int           statr, ir, ci;
  unsigned long flen, fbidx;
  ch_sym_type  *fbuf;
  double        ang_ref, ang_off, max_amp;
  unsigned char en;
  double        ma;
  ftype         tbuf[clk_per_sym], sumx, sumx2, mean, sd;
  long          scnt;
  double        phased;

  unsigned char verstr[]   = "2003May19";
  unsigned char progname[] = "mod-pm8a";
  unsigned char cww[]
      = "Copyright (C) 2003 Barry Sanderson There is ABSOLUTELY \
	 NO WARRANTY  for mod-pm8a. mod-pm8a is covered by the the GNU General\
	  Public License. The file \"COPYING.txt\" states the conditions under\
	   which you may legally: copy, modify, or distribute mod-pm8a.";

  printf("\n%s version %s\n\n%s\n\n", progname, verstr, cww);
  if(argc < 3)
  {
    ERN = 1;
    goto ex;
  }
  ERN        = 3;
  clk_period = 1.0 / 11025.0;

  dcf = 1;
  en  = 'h';

  /* initialize asteps[][]
   */
  phased = sqrt(phased2);
  for(si = -4; si <= 4; si++)
    set_asteps(0.0, ((double) (si) / phased) * PI2, &asteps[si + 4][0]);

/*for test purposes, print asteps[][] to stdout, in gnuplot format */
#if(dbgl == 1)
  for(si = 0; si < num_phases; si++)
  {
    for(di = 0; di < clk_per_sym; di++)
    {
      printf("%8.5f\n", asteps[si][di]);
    }
    printf("\n");
  }
  exit(0);
#endif

#if(dbgl == 2)
  for(si = 0; si < num_phases; si++)
  {
    ma = 0.0;
    for(di = 0; di < clk_per_sym; di++)
    {
      ma += asteps[si][di];
      printf("%8.5f\n", ma);
    }
    printf("\n");
  }
  exit(0);
#endif

  ca    = 0.0;
  astep = dum_sub * PI2 * clk_period; /* for dummy sub-carrier */
                                      /* convert from degrees to radians */
  for(lo = 0; lo < num_sc_t; lo++)
  {
    loai[lo] = loai[lo] * PI2 / 360.0;
  }

  statr = stat(argv[1], &statb);
  if(statr != 0)
  {
    ERN = 8;
    goto ex;
  }
  flen = (unsigned long) (statb.st_size);
  flen = flen / (sizeof(ch_sym_type));

  /* allocate space for data arrays */
  ERN = 5;
  if((fbuf = ALLOC_MEM((size_t) ((flen + 3) * sizeof(ch_sym_type)))) == NULL)
    goto exf;

  ERN = 1;
  inf = fopen(argv[1], "rb");
  if(inf == NULL)
  {
    fprintf(stderr, "Could not open file %s\n", argv[1]);
    goto exg;
  }
  ir = fread(&fbuf[1], sizeof(ch_sym_type), flen, inf);
  if(ir != flen)
  {
    ERN = 9;
    goto exh;
  }
  if(ferror(inf) != 0)
  {
    ERN = 9;
    goto exh;
  }
  fclose(inf);
  ERN = 2;

  for(lo = 0; lo < num_freqs; lo++)
  {
    ERN       = 3;
    fbidx     = 0;
    lo_db[lo] = fopen(lodbfn[lo], "wb");
    if(lo_db[lo] == NULL)
    {
      fprintf(stderr, "Could not open file %s\n", lodbfn[lo]);
      goto exh;
    }
    cflags = init_flag;
    ca     = 0.0;

    ma = 0.0;
    for(di = 0; di < flen; di++)
    {
      csym = fbuf[di + 1];
      if(csym >= sym_num_limit)
      {
        ERN = 11;
        fclose(lo_db[lo]);
        goto exh;
      }
      ci = sym_to_phi(csym, (int) (lo));
      for(si = 0; si < clk_per_sym; si++)
      {
        tbuf[si] = (ftype) (cos(ca + ma));
        ca       += astep;
        if(ca > PI2)
          ca -= PI2;
        ma += asteps[ci][si];
        if(ma > PI2)
          ma -= PI2;
        else
        {
          if(ma < -PI2)
            ma += PI2;
        }
      }
      val_written
          = fwrite(&tbuf[0], sizeof(ftype), (size_t) (clk_per_sym), lo_db[lo]);
      if(val_written != (size_t) (clk_per_sym))
      {
        ERN = 4;
        fclose(lo_db[lo]);
        remove(lodbfn[lo]);
        goto exh;
      }
    }

    fflush(lo_db[lo]);
    fclose(lo_db[lo]);
  } /* end of for lo */

  /* save and translate desired band of frequencies */
  astep = (dum_sub) *PI2 * clk_period; /* for dummy sub-carrier */
  for(lo = 0; lo < num_freqs; lo++)
  { /* open files for this lo */

    ERN       = 2;
    sb_db[lo] = fopen(sbdbfn[lo], "wb");
    if(sb_db[lo] == NULL)
    {
      fprintf(stderr, "Could not open file %s\n", sbdbfn[lo]);
      goto exh;
    }
    lo_db[lo] = fopen(lodbfn[lo], "rb");
    if(lo_db[lo] == NULL)
    {
      fprintf(stderr, "Could not open file %s\n", lodbfn[lo]);
      goto exi;
    }
    sidx        = 0;
    end_of_file = 0;
    ca          = 0.0;
    ca2         = loai[lo];
    astep2      = pmlo[lo] * PI2 * clk_period;

    while(end_of_file == 0)
    {
      val_to_read = s0blen;
      dum_len     = 0;
      if(sidx < s0off)
      {
        dum_len = (int) (s0off - sidx);
        for(si = 0; si < dum_len; si++)
          Sambuf[si] = dumd;
        val_to_read -= dum_len;
        if(fseek(lo_db[lo], 0L, SEEK_SET) != 0)
        {
          ERN = 7;
          goto exj;
        }
        ang_ref = 0.0;
        ang_off = astep * (double) (dum_len);
        ca      = ang_ref - ang_off;
        ca      = fmod(ca, PI2);
        if(ca < 0.)
          ca += PI2;
      }
      else
      {
        if(fseek(lo_db[lo], (sidx - s0off) * sizeof(ftype), SEEK_SET) != 0)
        {
          ERN = 7;
          goto exj;
        }
      }
      val_to_write = fap;
      val_read = fread(&Sambuf[dum_len], sizeof(ftype), val_to_read, lo_db[lo]);
      if(val_read < val_to_read)
      {
        if(val_read < s0off)
        {
          end_of_file  = 1;
          val_to_write = 0;
        }
        else
        {
          val_to_write = (val_read > (s0off + fap) ? fap : (val_read - s0off));
        }
        for(si = val_read; si < val_to_read; si++)
          Sambuf[si] = dumd;
      }
      /* generate SSB signal for data in Sambuf */
      /* multiply by sin and cos of dum_sub */
      for(si = 0; si < s0blen; si++)
      {
        Lo1c[si] = Sambuf[si] * (ftype) (cos(ca));
        Lo1s[si] = Sambuf[si] * (ftype) (sin(ca));
        ca       += astep;
        if(ca > PI2)
          ca -= PI2;
      }
      /* do decimation */
      DS1(&Lo1c[ds1dly], Lo1cf, s1blen);
      DS2(&Lo1cf[ds2dly], Lo2cs1, s2blen);
      DS3(&Lo2cs1[ds3dly], Lo2cs2, s3blen);
      DS1(&Lo1s[ds1dly], Lo1sf, s1blen);
      DS2(&Lo1sf[ds2dly], Lo2ss1, s2blen);
      DS3(&Lo2ss1[ds3dly], Lo2ss2, s3blen);
      /* do interpolation */
      di = I40x(&Lo2cs2[is1off], Ssbc, ipts);
      if(di < 0)
      {
        ERN = 10;
        goto exj;
      }
      di = I40x(&Lo2ss2[is1off], Ssbs, ipts);
      if(di < 0)
      {
        ERN = 10;
        goto exj;
      }
      /* combine sin and cos channels */
      for(si = 0; si < fap; si++)
      {
        Amp[si] = Ssbc[si] * (ftype) (cos(ca2)) + Ssbs[si] * (ftype) (sin(ca2));
        ca2     += astep2;
        if(ca2 > PI2)
          ca2 -= PI2;
      }
      /* write result for this lo */
      val_written = fwrite(&Amp[0], sizeof(ftype), val_to_write, sb_db[lo]);
      if(val_written != (size_t) (val_to_write))
      {
        ERN = 4;
        goto exj;
      }
      sidx    += fap;
      /* adjust angle of local oscillators */
      ang_ref = astep * (double) (sidx);
      /**   ang_off= astep*(double)(sidx >= s0off+s0off ? s0off: sidx-s0off);
       * **/
      ang_off = astep * (double) (sidx >= s0off + s0off ? s0off : s0off);
      ca      = ang_ref - ang_off;
      ca      = fmod(ca, PI2);
      if(ca < 0.)
        ca += PI2;
      /***99feb03	ca2 advances fap steps per loop with no overlap or gap
                      from loop to loop (unlike ca, which has a lot of overlap)
      99feb03***/
    } /* end of while end_of_file */
    fflush(sb_db[lo]);
    fclose(sb_db[lo]);
    fclose(lo_db[lo]);
    remove(lodbfn[lo]); /* this reduces the amount of disk space used */

  } /* end of for lo */
    /* combine all channels into final output */
  ERN = 2;
  otf = fopen(argv[2], "wb");
  if(otf == NULL)
  {
    fprintf(stderr, "Could not open file %s\n", argv[2]);
    goto exh;
  }
  for(lo = 0; lo < num_freqs; lo++)
  {
    sb_db[lo] = fopen(sbdbfn[lo], "rb");
    if(sb_db[lo] == NULL)
    {
      fprintf(stderr, "Could not open file %s\n", sbdbfn[lo]);
      goto exk;
    }
  }
  end_of_file = 0;
  max_amp     = 0.0;
  sumx        = 0;
  sumx2       = 0;
  scnt        = 0;

  /* stagger phase transitions */
  for(lo = 0; lo < num_freqs; lo++)
  {
    for(si = 0; si < dblen; si++)
      dbb[lo][si] = 0.0;
  }
  for(lo = 0; lo < num_freqs; lo++)
  {
    val_to_read = dblen - soff[lo];
    val_read = fread(&dbb[lo][soff[lo]], sizeof(ftype), val_to_read, sb_db[lo]);
  }
  for(si = 0; si < dblen; si++)
  {
    Amp[si] = 0.0;
    for(lo = 0; lo < num_freqs; lo++)
      Amp[si] += dbb[lo][si];
    if(fabs((double) (Amp[si])) > max_amp)
      max_amp = fabs((double) (Amp[si]));
    sumx  += Amp[si];
    sumx2 += Amp[si] * Amp[si];
    scnt++;
  }
  val_written = fwrite(&Amp[0], sizeof(ftype), dblen, otf);
  if(val_written != (size_t) (dblen))
  {
    ERN = 4;
    lo  = num_freqs;
    goto exk;
  }

  while(end_of_file < num_freqs)
  {
    val_to_read = dblen;
    for(lo = 0; lo < num_freqs; lo++)
    {
      val_read = fread(dbb[lo], sizeof(ftype), val_to_read, sb_db[lo]);
      if(val_read < val_to_read)
      {
        end_of_file++;
        si = val_read;
        while(si < dblen)
        {
          dbb[lo][si] = 0.0;
          si++;
        } /* fill with 0 */
        if(lo == num_freqs - 1)
          val_to_read = val_read;
      }
    }
    for(si = 0; si < val_to_read; si++)
    {
      Amp[si] = 0.0;
      for(lo = 0; lo < num_freqs; lo++)
        Amp[si] += dbb[lo][si];
      if(fabs((double) (Amp[si])) > max_amp)
        max_amp = fabs((double) (Amp[si]));
      sumx  += Amp[si];
      sumx2 += Amp[si] * Amp[si];
      scnt++;
    }
    val_written = fwrite(&Amp[0], sizeof(ftype), val_to_read, otf);
    if(val_written != (size_t) (val_to_read))
    {
      ERN = 4;
      lo  = num_freqs;
      goto exk;
    }
  }
  fflush(otf);
  fclose(otf);
  for(lo = 0; lo < num_freqs; lo++)
  {
    fclose(sb_db[lo]);
    remove(sbdbfn[lo]);
  }

  mean = sumx / (ftype) (scnt);
  sd   = (ftype) (sqrt(
      (double) ((sumx2 - (ftype) (scnt) *mean * mean) / ((ftype) (scnt - 1)))));

  printf("The standard deviation of the sample values is %7.4f\n", sd);

  ERN = 0;
  goto exh;

exk:
  fclose(otf);
  lo--;
  while(lo >= 0)
  {
    fclose(sb_db[lo]);
    lo--;
    remove(sbdbfn[lo]);
  }
  goto exh;
exj:
  fclose(lo_db[lo]);
  remove(lodbfn[lo]);
exi:
  fclose(sb_db[lo]);
  remove(sbdbfn[lo]);
exh:
exg:
  free(fbuf);
exf:
ex:
  if(ERN != 0)
    fprintf(stderr, "%s", erexmsg[ERN]);
  exit(ERN);
}
