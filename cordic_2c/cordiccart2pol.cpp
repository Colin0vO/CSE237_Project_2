#include "cordiccart2pol.h"
#include <ap_fixed.h>
#include <cmath>
#include <cstdio>

#ifndef NO_ITER
#define NO_ITER 16
#endif

#define W_TOTAL 16
#define I_ALL   3  
typedef ap_fixed<W_TOTAL, I_ALL> fx_all;   

#ifndef W_TBL
#define W_TBL 16
#endif
typedef ap_fixed<W_TBL, I_ALL> fx_tbl;    

#define USE_SHIFT_ONLY 1

static inline fx_tbl atan_2mi(int i){
  double a = std::atan(std::ldexp(1.0, -i));
  return fx_tbl(a);
}

void cordiccart2pol(data_t x_in, data_t y_in, data_t *r_out, data_t *theta_out)
{
  static bool once=false;
  if(!once){
    once=true;
    std::printf("[CFG-2c] W_TOTAL=16 I_ALL=3  W_TBL=%d  NO_ITER=%d  SHIFT=%d\n",
                W_TBL, NO_ITER, USE_SHIFT_ONLY);
  }

  fx_all x = fx_all(x_in);
  fx_all y = fx_all(y_in);
  fx_all z = 0;

  fx_tbl alpha[NO_ITER];
#pragma HLS array_partition variable=alpha complete
  for(int i=0;i<NO_ITER;++i){
#pragma HLS unroll
    alpha[i] = atan_2mi(i);
  }

  for(int i=0;i<NO_ITER;++i){
#pragma HLS pipeline II=1
#if USE_SHIFT_ONLY
    fx_all x_shift = (fx_all)( y >> i );
    fx_all y_shift = (fx_all)( x >> i );
#else
#endif
    if (y >= 0){
      x = (fx_all)(x + x_shift);
      y = (fx_all)(y - y_shift);
      z = (fx_all)(z + (fx_all)alpha[i]);
    }else{
      x = (fx_all)(x - x_shift);
      y = (fx_all)(y + y_shift);
      z = (fx_all)(z - (fx_all)alpha[i]);
    }
  }

  fx_tbl K = fx_tbl(1.0);
  for(int i=0;i<NO_ITER;++i){
#pragma HLS unroll
    double term = std::sqrt(1.0 + std::ldexp(1.0, -2*i));
    K = K / fx_tbl(term);
  }

  fx_all x_abs = (x < 0) ? (fx_all)(-x) : x;
  fx_all r = (fx_all)( fx_all(K) * x_abs );

  *r_out     = (data_t)r;
  *theta_out = (data_t)z;
}
