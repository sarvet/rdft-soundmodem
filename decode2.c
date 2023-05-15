/*
    inner decoding routines used by Wyman1x digital communications programs
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
/* Rev 25 Nov  2001  version 1.1.0, for initial public release */
/* Rev 17 Nov  2001  fixed detection of out of bounds demodulated values */
/* Rev 29 July 2001  decode2ed2... subroutines introduced */
/* Rev 12 May  2001  version 1.0.0 for pm7b initial release */

#include "code2-vectors-bias.h"
#include "code2-we.h"
#include "code2.h"
#include "dcom.h"
#include "mode-id-vectors-bias.h"
#include "rs8-4.h"
#include "sym_def.h"
#include <math.h>

#define bias_ps                                                   \
  4                   /* bias due to translation of -4 through +4 \
range to 0 through  8 for 8 sub-carriers,  9 phases case */
#define max_dcphi_b 8 /* maximum legal decoded phase step + bias */
#define min_dcphi_b 0 /* minimum legal decoded phase step + bias */

int decode2lut(ch_sym_type key)
{ /* Binary search of code2_srt_sym_num[].  If key is found in table, entry
         from code2_srt_sym_num_idx[] at the corresponding location is
         returned, otherwise -1 is returned. */

  int         loidx, hiidx, tiidx;
  ch_sym_type tv;

  loidx = 0;
  hiidx = N1 - 1;
  while((hiidx - loidx) > 1)
  {
    tiidx = loidx + (hiidx - loidx) / 2;
    tv    = code2_srt_sym_num[tiidx];
    if(key > tv)
      loidx = tiidx; /* raise loidx */
    else
    {
      if(key < tv)
        hiidx = tiidx; /* lower hiidx */
      else
      {
        return (code2_srt_sym_num_idx[tiidx]);
      } /* match found */
    }
  } /* end of while index diff > 1 */
  if(key == code2_srt_sym_num[loidx])
    return (code2_srt_sym_num_idx[loidx]);
  else
  {
    if(key == code2_srt_sym_num[hiidx])
      return (code2_srt_sym_num_idx[hiidx]);
  }
  return -1; /* no match found */
}

sym_type_ch mk_ch_sym_num(sym_typei ra[], sym_typei ia[])
{ /* calculate and return the integer corresponding to the
         concatenation of the base 9 values in ra[] and ia[]
 */
  int         si;
  sym_type_ch acum;

  acum = 0;
  for(si = num_sc_i - 1; si >= 0; si--)
    acum = (acum + ia[si]) * N2;
  for(si = num_sc_r - 1; si > 0; si--)
    acum = (acum + ra[si]) * N2;
  return (acum + ra[0]);
}

int decode2ed1mode(ftype *data, /*ftype vectors[num_ec_lvls][num_sc_t],*/
                   ftype *dist, ftype dthld)
{ /* decode data[0...sum_sc_r-1, num_sc_r...num_sc_t-1], by finding
     vector closest to data values. Return its index, and set dist
     to distance from this vector, or return -1 if decoding fails */

  int   si, bi;
  ftype tc, bd, acum;
  int   erai[num_sc_t], ci, ec;

  for(si = 0; si < num_ec_lvls; si++)
  {
    acum = 0.0;
    for(ci = 0; ci < num_sc_t; ci++)
    {
      tc   = data[ci] - modeid1_vectors[si][ci];
      acum += tc * tc;
    }
    if(si > 0)
    {
      if(acum < bd)
      {
        bd = acum;
        bi = si;
      }
    }
    else
    {
      bd = acum;
      bi = si;
    }
  } /* end of for si */
  dist[0] = (ftype) (sqrt((double) (bd)));
  if(dist[0] < dthld)
    return (bi);
  /* otherwise  try cases with 1 erasure */
  for(ec = 0; ec < c8t1; ec++)
  {
    for(ci = 0; ci < num_sc_t; ci++)
      erai[ci] = 1;
    for(ci = 0; ci < 1; ci++)
      erai[code2_comb1e[ec][ci]] = 0;
    for(si = 0; si < num_ec_lvls; si++)
    {
      acum = 0.0;
      for(ci = 0; ci < num_sc_t; ci++)
      {
        if(erai[ci])
        {
          tc   = data[ci] - modeid1_vectors[si][ci];
          acum += tc * tc;
        }
      } /* end of for ci */
      if(si > 0)
      {
        if(acum < bd)
        {
          bd = acum;
          bi = si;
        }
      }
      else
      {
        bd = acum;
        bi = si;
      }
    } /* end of for si<num_ec_lvls */
    bd = (ftype) (sqrt((double) (bd)));
    if(bd < mode1_comb1edthld[ec]) /* successful decode with 1 erasures */
    {
      dist[0] = bd;
      return (bi);
    }
  } /* end of for ec<c8t1 */

  /* next try cases with 2 erasures */
  for(ec = 0; ec < c8t2; ec++)
  {
    for(ci = 0; ci < num_sc_t; ci++)
      erai[ci] = 1;
    for(ci = 0; ci < 2; ci++)
      erai[code2_comb2e[ec][ci]] = 0;
    for(si = 0; si < num_ec_lvls; si++)
    {
      acum = 0.0;
      for(ci = 0; ci < num_sc_t; ci++)
      {
        if(erai[ci])
        {
          tc   = data[ci] - modeid1_vectors[si][ci];
          acum += tc * tc;
        }
      } /* end of for ci */
      if(si > 0)
      {
        if(acum < bd)
        {
          bd = acum;
          bi = si;
        }
      }
      else
      {
        bd = acum;
        bi = si;
      }
    } /* end of for si<mcwi */
    bd = (ftype) (sqrt((double) (bd)));
    if(bd < mode1_comb2edthld[ec]) /* successful decode with 2 erasures */
    {
      dist[0] = bd;
      return (bi);
    }
  }            /* end of for ec<c8t2 */
  return (-1); /* indicate failure to decode */
}

int decode2ed2mode(ftype *data, /*ftype vectors[num_ec_lvls][num_sc_t],*/
                   ftype *dist, ftype dthld)
{ /* decode data[0...sum_sc_r-1, num_sc_r...num_sc_t-1], by finding
     vector closest to data values. Return its index, and set dist
     to distance from this vector, or return -1 if decoding fails */

  int   si, bi;
  ftype tc, bd, acum;
  int   erai[num_sc_t], ci, ec;

  for(si = 0; si < num_ec_lvls; si++)
  {
    acum = 0.0;
    for(ci = 0; ci < num_sc_t; ci++)
    {
      tc   = data[ci] - modeid2_vectors[si][ci];
      acum += tc * tc;
    }
    if(si > 0)
    {
      if(acum < bd)
      {
        bd = acum;
        bi = si;
      }
    }
    else
    {
      bd = acum;
      bi = si;
    }
  } /* end of for si */
  dist[0] = (ftype) (sqrt((double) (bd)));
  if(dist[0] < dthld)
    return (bi);
  /* otherwise  try cases with 1 erasure */
  for(ec = 0; ec < c8t1; ec++)
  {
    for(ci = 0; ci < num_sc_t; ci++)
      erai[ci] = 1;
    for(ci = 0; ci < 1; ci++)
      erai[code2_comb1e[ec][ci]] = 0;
    for(si = 0; si < num_ec_lvls; si++)
    {
      acum = 0.0;
      for(ci = 0; ci < num_sc_t; ci++)
      {
        if(erai[ci])
        {
          tc   = data[ci] - modeid2_vectors[si][ci];
          acum += tc * tc;
        }
      } /* end of for ci */
      if(si > 0)
      {
        if(acum < bd)
        {
          bd = acum;
          bi = si;
        }
      }
      else
      {
        bd = acum;
        bi = si;
      }
    } /* end of for si<num_ec_lvls */
    bd = (ftype) (sqrt((double) (bd)));
    if(bd < mode2_comb1edthld[ec]) /* successful decode with 1 erasures */
    {
      dist[0] = bd;
      return (bi);
    }
  } /* end of for ec<c8t1 */

  /* next try cases with 2 erasures */
  for(ec = 0; ec < c8t2; ec++)
  {
    for(ci = 0; ci < num_sc_t; ci++)
      erai[ci] = 1;
    for(ci = 0; ci < 2; ci++)
      erai[code2_comb2e[ec][ci]] = 0;
    for(si = 0; si < num_ec_lvls; si++)
    {
      acum = 0.0;
      for(ci = 0; ci < num_sc_t; ci++)
      {
        if(erai[ci])
        {
          tc   = data[ci] - modeid2_vectors[si][ci];
          acum += tc * tc;
        }
      } /* end of for ci */
      if(si > 0)
      {
        if(acum < bd)
        {
          bd = acum;
          bi = si;
        }
      }
      else
      {
        bd = acum;
        bi = si;
      }
    } /* end of for si<mcwi */
    bd = (ftype) (sqrt((double) (bd)));
    if(bd < mode2_comb2edthld[ec]) /* successful decode with 2 erasures */
    {
      dist[0] = bd;
      return (bi);
    }
  }            /* end of for ec<c8t2 */
  return (-1); /* indicate failure to decode */
}

int decode2wne2(ftype *data, int *scc, unsigned short int *ob_cnt,
                unsigned short int *ob_mask)
{ /* decode data[0...sum_sc_r-1, num_sc_r...num_sc_t-1],
         return -1 if too many erased components
                -2 if decoder fails
                -3 if decoder produces codeword in GF(9) that is not
                      in table of values used
         set scc to number of sub-carriers changed, when decoding does not fail
         set ob_cnt and ob_mask, based on out of bounds subcarrier values
 */

  int                erai[num_sc_t], eracnt, si, fi, ci;
  sym_typei          code2blk[num_sc_t], c2blk[num_sc_t];
  ch_sym_type        cwin;
  int                cwot, dcec;
  unsigned short int ob_sel;

  /* convert ftype input data to sym_typei */
  for(si = 0; si < num_sc_t; si++) /* changed upper limit on si 17nov2001 */
  {
    fi = (int) (floor((double) (data[si])));
    ci = (int) (ceil((double) (data[si])));
    code2blk[si]
        = bias_ps
        + (sym_typei) (((data[si] - (ftype) (fi)) < ((ftype) (ci) -data[si]))
                           ? fi
                           : ci);
  }

  /* check for any out of bounds values */
  eracnt     = 0;
  ob_cnt[0]  = 0;
  ob_mask[0] = 0x0;
  ob_sel     = 0x01;
  for(si = 0; si < num_sc_t; si++)
  {
    if(code2blk[si] > max_dcphi_b)
    {
      erai[eracnt++] = si;
      ob_cnt[0]++;
      ob_mask[0] = ob_mask[0] | ob_sel;
    }
    else
    {
      if(code2blk[si] < min_dcphi_b)
      {
        erai[eracnt++] = si;
        ob_cnt[0]++;
        ob_mask[0] = ob_mask[0] | ob_sel;
      }
    }
    ob_sel = ob_sel << 1;
  } /* end of for si */

  if(eracnt == 0)
  {
    cwin = mk_ch_sym_num(&code2blk[0], &code2blk[num_sc_r]);
    cwot = decode2lut(cwin);
    if(cwot >= 0)
    {
      *scc = 0;
      return (cwot);
    }
    else
    {
      for(si = 0; si < num_sc_t; si++)
        c2blk[si] = code2blk[si];
      dcec = decode2(&c2blk[0], eracnt, &erai[0]);
      if((dcec >= 0) && (dcec <= 2))
      { /* need to get decoded symbol number {0,...,N1-1} to return */
        cwin = mk_ch_sym_num(&c2blk[0], &c2blk[num_sc_r]);
        cwot = decode2lut(cwin);
        if(cwot >= 0)
        {
          *scc = dcec;
          return (cwot);
        }
        else
          return (-3);
      }
      else
        return (-2);
    }
  }
  else
    return (-1);
}

int decode2ed2w1e(ftype *data, ftype *dist, int mcwi, int *era, int *scc,
                  unsigned short int ob_mask, unsigned short int *bec_mask)
{ /* decode data[0...sum_sc_r-1, num_sc_r...num_sc_t-1], by finding
     codeword closest to data values, trying all combinations of 1 erased
     symbol, that are consistant with ob_mask. Return index of closest match
     found, set dist to distance from this codeword, set bec_mask[0] to mask for
     best match,
     and set era[0] to 1 if returned codeword is NOT close enough to contents
     of data[], otherwise, set era[0] to 0.
  mcwi is code word index limit
 */

  int   erai[num_sc_t], si, ci, bi, ec, bec, init_flg;
  ftype tc, bd, acum;

  init_flg = 1;
  /* cases with 1 erasure */
  for(ec = 0; ec < c8t1; ec++)
  {
    if(ob_mask != 0x0)
    {
      if((ob_mask & code2_comb1em[ec]) != ob_mask)
        continue;
    }
    for(ci = 0; ci < num_sc_t; ci++)
      erai[ci] = 1;
    for(ci = 0; ci < 1; ci++)
      erai[code2_comb1e[ec][ci]] = 0;
    for(si = 0; si < mcwi; si++)
    {
      acum = 0.0;
      for(ci = 0; ci < num_sc_t; ci++)
      {
        if(erai[ci])
        {
          tc   = data[ci] - code2_vectors[si][ci];
          acum += tc * tc;
        }
      } /* end of for ci */
      if(init_flg == 1)
      {
        init_flg = 0;
        bd       = acum;
        bi       = si;
        bec      = ec;
      }
      else
      {
        if(acum < bd)
        {
          bd  = acum;
          bi  = si;
          bec = ec;
        }
      }
    } /* end of for si<mcwi */
  }   /* end of for ec<c8t1 */
  bd          = (ftype) (sqrt((double) (bd)));
  scc[0]      = 1;
  dist[0]     = bd;
  bec_mask[0] = code2_comb1em[bec];
  if(bd < code2_comb1edthld[bec])
    era[0] = 0;
  else
    era[0] = 1;
  return (bi);
}

int decode2ed2w2e(ftype *data, ftype *dist, int mcwi, int *era, int *scc,
                  unsigned short int ob_mask, unsigned short int *bec_mask)
{ /* decode data[0...sum_sc_r-1, num_sc_r...num_sc_t-1], by finding
     codeword closest to data values, trying all combinations of 2 erased
     symbols, that are consistant with ob_mask. Return index of closest match
     found, set dist to distance from this codeword, set bec_mask[0] to mask for
     best match,
     and set era[0] to 1 if returned codeword is NOT close enough to contents
     of data[], otherwise, set era[0] to 0.
  mcwi is code word index limit
 */

  int   erai[num_sc_t], si, ci, bi, ec, bec, init_flg;
  ftype tc, bd, acum;

  init_flg = 1;
  /* cases with 2 erasures */
  for(ec = 0; ec < c8t2; ec++)
  {
    if(ob_mask != 0x0)
    {
      if((ob_mask & code2_comb2em[ec]) != ob_mask)
        continue;
    }
    for(ci = 0; ci < num_sc_t; ci++)
      erai[ci] = 1;
    for(ci = 0; ci < 2; ci++)
      erai[code2_comb2e[ec][ci]] = 0;
    for(si = 0; si < mcwi; si++)
    {
      acum = 0.0;
      for(ci = 0; ci < num_sc_t; ci++)
      {
        if(erai[ci])
        {
          tc   = data[ci] - code2_vectors[si][ci];
          acum += tc * tc;
        }
      } /* end of for ci */
      if(init_flg == 1)
      {
        init_flg = 0;
        bd       = acum;
        bi       = si;
        bec      = ec;
      }
      else
      {
        if(acum < bd)
        {
          bd  = acum;
          bi  = si;
          bec = ec;
        }
      }
    } /* end of for si<mcwi */
  }   /* end of for ec<c8t2 */
  bd          = (ftype) (sqrt((double) (bd)));
  scc[0]      = 2;
  dist[0]     = bd;
  bec_mask[0] = code2_comb2em[bec];
  if(bd < code2_comb2edthld[bec])
    era[0] = 0;
  else
    era[0] = 1;
  return (bi);
}

int decode2ed2w3e(ftype *data, ftype *dist, int mcwi, int *era, int *scc,
                  unsigned short int ob_mask, unsigned short int *bec_mask)
{ /* decode data[0...sum_sc_r-1, num_sc_r...num_sc_t-1], by finding
     codeword closest to data values, trying all combinations of 3 erased
     symbols, that are consistant with ob_mask. Return index of closest match
     found, set dist to distance from this codeword, set bec_mask[0] to mask for
     best match,
     and set era[0] to 1 if returned codeword is NOT close enough to contents
     of data[], otherwise, set era[0] to 0.
  mcwi is code word index limit
 */

  int   erai[num_sc_t], si, ci, bi, ec, bec, init_flg;
  ftype tc, bd, acum;

  init_flg = 1;
  /* cases with 3 erasures */
  for(ec = 0; ec < c8t3; ec++)
  {
    if(ob_mask != 0x0)
    {
      if((ob_mask & code2_comb3em[ec]) != ob_mask)
        continue;
    }
    for(ci = 0; ci < num_sc_t; ci++)
      erai[ci] = 1;
    for(ci = 0; ci < 3; ci++)
      erai[code2_comb3e[ec][ci]] = 0;
    for(si = 0; si < mcwi; si++)
    {
      acum = 0.0;
      for(ci = 0; ci < num_sc_t; ci++)
      {
        if(erai[ci])
        {
          tc   = data[ci] - code2_vectors[si][ci];
          acum += tc * tc;
        }
      } /* end of for ci */
      if(init_flg == 1)
      {
        init_flg = 0;
        bd       = acum;
        bi       = si;
        bec      = ec;
      }
      else
      {
        if(acum < bd)
        {
          bd  = acum;
          bi  = si;
          bec = ec;
        }
      }
    } /* end of for si<mcwi */
  }   /* end of for ec<c8t3 */
  bd          = (ftype) (sqrt((double) (bd)));
  scc[0]      = 3;
  dist[0]     = bd;
  bec_mask[0] = code2_comb3em[bec];
  if(bd < code2_comb3edthld[bec])
    era[0] = 0;
  else
    era[0] = 1;
  return (bi);
}

int decode2ed2w4e(ftype *data, ftype *dist, int mcwi, int *era, int *scc,
                  unsigned short int ob_mask, unsigned short int *bec_mask)
{ /* decode data[0...sum_sc_r-1, num_sc_r...num_sc_t-1], by finding
     codeword closest to data values, trying all combinations of 4 erased
     symbols, that are consistant with ob_mask. Return index of closest match
     found, set dist to distance from this codeword, set bec_mask[0] to mask for
     best match,
     and set era[0] to 1 if returned codeword is NOT close enough to contents
     of data[], otherwise, set era[0] to 0.
  mcwi is code word index limit
 */

  int   erai[num_sc_t], si, ci, bi, ec, bec, init_flg;
  ftype tc, bd, acum;

  init_flg = 1;
  /* cases with 4 erasures */
  for(ec = 0; ec < c8t4; ec++)
  {
    if(ob_mask != 0x0)
    {
      if((ob_mask & code2_comb4em[ec]) != ob_mask)
        continue;
    }
    for(ci = 0; ci < num_sc_t; ci++)
      erai[ci] = 1;
    for(ci = 0; ci < 4; ci++)
      erai[code2_comb4e[ec][ci]] = 0;
    for(si = 0; si < mcwi; si++)
    {
      acum = 0.0;
      for(ci = 0; ci < num_sc_t; ci++)
      {
        if(erai[ci])
        {
          tc   = data[ci] - code2_vectors[si][ci];
          acum += tc * tc;
        }
      } /* end of for ci */
      if(init_flg == 1)
      {
        init_flg = 0;
        bd       = acum;
        bi       = si;
        bec      = ec;
      }
      else
      {
        if(acum < bd)
        {
          bd  = acum;
          bi  = si;
          bec = ec;
        }
      }
    } /* end of for si<mcwi */
  }   /* end of for ec<c8t4 */
  bd          = (ftype) (sqrt((double) (bd)));
  scc[0]      = 4;
  dist[0]     = bd;
  bec_mask[0] = code2_comb4em[bec];
  if(bd < code2_comb4edthld[bec])
    era[0] = 0;
  else
    era[0] = 1;
  return (bi);
}

void init_code2(void)
{
  int si;
  /* finish initialization of distance threshold arrays */
  for(si = 0; si < c8t1; si++)
    code2_comb1edthld[si]
        = (ftype) (0.5 * sqrt((double) (code2_comb1edthld[si])));
  for(si = 0; si < c8t2; si++)
    code2_comb2edthld[si]
        = (ftype) (0.5 * sqrt((double) (code2_comb2edthld[si])));
  for(si = 0; si < c8t3; si++)
    code2_comb3edthld[si]
        = (ftype) (0.5 * sqrt((double) (code2_comb3edthld[si])));
  for(si = 0; si < c8t4; si++)
    code2_comb4edthld[si]
        = (ftype) (0.5 * sqrt((double) (code2_comb4edthld[si])));
};
