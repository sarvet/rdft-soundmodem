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

#define trail_spec_len_nom	168
ftype	trail_spec_nom[trail_spec_len_nom]={
    10034.6,    14893.6,    11735.5,     7700.2,     2704.0,
   167494.3,   927897.5,   794328.2,    81470.4,    21305.9,
    11181.5,     7735.7,     6471.4,    14808.1,     6295.1,
    14638.6,     9088.7,     9528.0,    12927.1,    11259.0,
     8600.0,    11053.5,     9109.6,     7038.8,     7353.6,
    16557.7,    31081.4,   643428.1,  1004615.8,   271019.2,
     2588.2,    10197.6,     9828.8,    15977.2,    10977.4,
     8035.3,    12302.7,    11091.7,     7753.5,    10471.3,
    12867.7,     7194.5,     7063.2,    14723.1,     7717.9,
   134276.5,   871966.9,   804451.7,   100230.5,    22856.0,
    14060.5,     7533.6,     9057.3,     6950.2,    20797.0,
    15470.3,     9026.1,   121618.6,   209411.2,    65013.0,
     5064.1,    11324.0,     4325.1,     8423.6,    11695.0,
     7970.8,    21355.0,   580096.2,   974989.6,   290736.8,
    12661.9,     8433.3,     6464.0,    20941.1,    12691.1,
    13106.9,    23469.3,     8346.4,    63023.1,   210862.8,
   126182.8,     1999.9,    66527.3,   210862.8,   120503.6,
     4221.8,    17885.5,    10303.9,    12545.8,    20797.0,
    16885.0,     5539.9,     7647.2,   100000.0,   800755.6,
   845278.8,   129270.7,     9057.3,    12459.5,     9862.8,
     5388.9,    12867.7,     5949.8,    26546.1,   177214.8,
   188799.1,    23577.6,    12119.9,    17498.5,     7816.3,
     8007.6,     8375.3,     6745.3,    16982.4,    13366.0,
   509917.6,   994260.1,   346337.9,     5610.5,    12274.4,
     7211.1,     3459.4,     9236.3,     7825.3,    10185.9,
     9817.5,     9931.2,     6760.8,     9268.3,    13152.2,
    12416.5,     8619.9,     4130.5,    89639.6,   806306.2,
   919390.5,   159037.7,    20820.9,    15940.4,     6998.4,
     7473.1,    14206.9,     5457.6,    15667.5,     6187.3,
    11681.5,    11155.8,     8780.1,     6067.4,     8423.6,
     8800.4,     7473.1,     3881.5,    11995.0,     5970.4,
   498310.5,  1045923.7,   402717.0,     3771.4,     6637.4,
     4639.8,    15275.7,     8550.7,     6538.8,     8560.5,
     8521.2,     4737.0,     9977.0};

#define	am_ref_len	181
ftype am_ref[am_ref_len]={
  4304.860287,  4298.982300,  4290.361988,  4279.179250,
  4265.187212,  4248.349662,  4228.365137,  4205.350375,
  4179.563088,  4150.645237,  4118.810588,  4083.983662,
  4046.193162,  4005.589300,  3961.876712,  3915.171987,
  3865.774650,  3813.519862,  3758.442888,  3700.968263,
  3641.160250,  3579.032412,  3514.727813,  3448.660462,
  3380.872575,  3311.608813,  3240.924737,  3169.212462,
  3097.070075,  3024.577525,  2952.280962,  2880.517338,
  2809.831237,  2740.790350,  2673.706837,  2609.210100,
  2548.111625,  2490.859675,  2438.108750,  2390.720862,
  2349.302350,  2314.469675,  2286.775013,  2266.855738,
  2254.990262,  2251.552238,  2256.567725,  2270.135175,
  2292.243812,  2322.624037,  2360.899013,  2406.618513,
  2459.268112,  2518.277175,  2582.973213,  2652.731112,
  2726.949125,  2804.924362,  2886.094038,  2969.811675,
  3055.528063,  3142.685462,  3230.825100,  3319.425362,
  3408.060263,  3496.241788,  3583.705625,  3670.010600,
  3754.856725,  3837.918550,  3918.857300,  3997.488800,
  4073.460800,  4146.570462,  4216.606400,  4283.368575,
  4346.655025,  4406.314513,  4462.200437,  4514.125075,
  4561.975475,  4605.577337,  4644.907087,  4679.773063,
  4710.180613,  4735.877000,  4757.041625,  4773.455687,
  4785.200925,  4792.157925,  4794.305663,  4791.721025,
  4784.300650,  4772.110850,  4755.164313,  4733.504587,
  4707.187200,  4676.288637,  4640.867900,  4601.025588,
  4556.833375,  4508.373400,  4455.827837,  4399.258725,
  4338.881812,  4274.699387,  4207.087913,  4136.101125,
  4062.054525,  3985.056725,  3905.378137,  3823.297050,
  3739.041713,  3652.922250,  3565.229325,  3476.305900,
  3386.529788,  3296.279800,  3205.941725,  3116.000637,
  3026.864525,  2939.100325,  2853.231300,  2769.792812,
  2689.378888,  2612.624087,  2540.219088,  2472.786437,
  2410.955475,  2355.415263,  2306.716412,  2265.386925,
  2231.976625,  2206.779062,  2190.154187,  2182.214888,
  2182.909913,  2192.111400,  2209.594500,  2234.801875,
  2267.224888,  2306.470425,  2351.840238,  2402.533537,
  2457.905500,  2517.395075,  2580.365800,  2646.157937,
  2713.891912,  2783.332850,  2853.752725,  2924.710538,
  2995.794375,  3066.493612,  3136.719900,  3206.094262,
  3274.136737,  3340.761275,  3405.811363,  3468.771475,
  3529.518350,  3588.074075,  3644.275312,  3697.802325,
  3748.690938,  3796.920075,  3842.556788,  3885.468087,
  3925.279662,  3962.338600,  3996.397338,  4027.512625,
  4055.712975,  4080.864550,  4103.330775,  4123.067000,
  4139.884750,  4154.076150,  4165.727275,  4174.612850,
  4180.862363};