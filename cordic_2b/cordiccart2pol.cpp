#include "cordiccart2pol.h"
#include <ap_fixed.h>
#include <cmath>
#include <cstdio>


#ifndef NO_ITER
#define NO_ITER 16
#endif

#ifndef I_XY
#define I_XY 1      // x, y
#endif
#ifndef I_R
#define I_R  2      // r
#endif
#ifndef I_TH
#define I_TH 3      // theta
#endif

#ifndef W_TOTAL
#define W_TOTAL 32
#endif

#ifndef W_TBL
#define W_TBL 16
#endif

#define USE_SHIFT_ONLY 0   

typedef ap_fixed<W_TOTAL, 3>    fx_xy_core;                
typedef ap_fixed<W_TOTAL, I_XY> fx_xy_io;                   
typedef ap_fixed<W_TOTAL, I_R > fx_r_io;                   
typedef ap_fixed<W_TOTAL, I_TH> fx_th_io;                   
typedef ap_fixed<W_TBL,  I_TH>  fx_tbl;                     
typedef ap_fixed<W_TOTAL+2, 4>  fx_acc;                     
typedef ap_fixed<W_TBL, 1>      fx_pow2;                    

static inline fx_tbl atan_2mi(int i){
  double a = std::atan(std::ldexp(1.0, -i));
  return fx_tbl(a);
}

void cordiccart2pol(data_t x_in, data_t y_in, data_t *r_out, data_t *theta_out)
{
  static bool once=false;
  if(!once){
    once=true;
    std::printf("[CFG] W_TOTAL=%d W_TBL=%d NO_ITER=%d (2b multiply)\n",
                W_TOTAL, W_TBL, NO_ITER);
  }

  fx_xy_core x = fx_xy_io(x_in);
  fx_xy_core y = fx_xy_io(y_in);
  fx_th_io   z = 0;

  const fx_th_io HALF_PI = fx_th_io(1.5707963267948966);
  if (x < 0) {
    fx_xy_core xt, yt;
    if (y >= 0) { xt =  y; yt = -x; z += HALF_PI; }
    else        { xt = -y; yt =  x; z -= HALF_PI; }
    x = xt; y = yt;
  }

  fx_tbl alpha[NO_ITER];
#pragma HLS array_partition variable=alpha complete
  for (int i=0;i<NO_ITER;++i){
#pragma HLS unroll
    alpha[i] = atan_2mi(i);
  }

  fx_pow2 pow2_neg[NO_ITER];
#pragma HLS array_partition variable=pow2_neg complete
  for (int i=0;i<NO_ITER;++i){
#pragma HLS unroll
    pow2_neg[i] = fx_pow2(std::ldexp(1.0, -i));
  }

  for (int i=0;i<NO_ITER;++i){
#pragma HLS pipeline II=1
    fx_xy_core x_shift = (fx_xy_core)((fx_acc)y * (fx_acc)pow2_neg[i]);
    fx_xy_core y_shift = (fx_xy_core)((fx_acc)x * (fx_acc)pow2_neg[i]);

    fx_xy_core x_next, y_next;
    fx_th_io   z_next;

    if (y >= 0) {
      fx_acc tx = (fx_acc)x + (fx_acc)x_shift;
      fx_acc ty = (fx_acc)y - (fx_acc)y_shift;
      x_next = (fx_xy_core)tx;
      y_next = (fx_xy_core)ty;
      z_next = (fx_th_io)( z + (fx_th_io)alpha[i] );
    } else {
      fx_acc tx = (fx_acc)x - (fx_acc)x_shift;
      fx_acc ty = (fx_acc)y + (fx_acc)y_shift;
      x_next = (fx_xy_core)tx;
      y_next = (fx_xy_core)ty;
      z_next = (fx_th_io)( z - (fx_th_io)alpha[i] );
    }

    x = x_next; y = y_next; z = z_next;
  }

  fx_tbl K = fx_tbl(1.0);
  for (int i=0;i<NO_ITER;++i){
#pragma HLS unroll
    double term = std::sqrt(1.0 + std::ldexp(1.0, -2*i));
    K = K / fx_tbl(term);
  }

  fx_xy_core x_abs = x; if (x < 0) x_abs = -x;

  fx_r_io r = (fx_r_io)x_abs * (fx_r_io)K;
  *r_out     = (data_t)r;
  *theta_out = (data_t)z;
}
