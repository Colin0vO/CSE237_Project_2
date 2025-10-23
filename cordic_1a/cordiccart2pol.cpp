#include "cordiccart2pol.h"
#include <cmath>

data_t Kvalues[20] = {
  1.0f,
  0.5f,
  0.25f,
  0.125f,
  0.0625f,
  0.03125f,
  0.015625f,
  0.0078125f,
  0.00390625f,
  0.001953125f,
  0.0009765625f,
  0.00048828125f,
  0.000244140625f,
  0.0001220703125f,
  6.10351562500000000e-05f,
  3.05175781250000000e-05f,
  1.52587890625000000e-05f,
  7.62939453125000000e-06f,
  3.81469726562500000e-06f,
  1.90734863281250000e-06f
};

data_t angles[20] = {
  0.78539816339744828f,   // atan(2^-0)
  0.46364760900080609f,   // atan(2^-1)
  0.24497866312686414f,   // atan(2^-2)
  0.12435499454676144f,   // atan(2^-3)
  0.06241880999595735f,   // atan(2^-4)
  0.031239833430268277f,  // atan(2^-5)
  0.015623728620476831f,  // atan(2^-6)
  0.0078123410601011111f, // atan(2^-7)
  0.0039062301319669718f, // atan(2^-8)
  0.0019531225164788188f, // atan(2^-9)
  0.00097656218955931946f,// atan(2^-10)
  0.00048828121119489829f,// atan(2^-11)
  0.00024414062014936177f,// atan(2^-12)
  0.00012207031189367021f,// atan(2^-13)
  0.0000610351561742087726,// atan(2^-14)
  0.0000305175781155260957,// atan(2^-15)
  1.52587890613157615e-05f,// atan(2^-16)
  7.62939453110196998e-06f,// atan(2^-17)
  3.81469726560649614e-06f,// atan(2^-18)
  1.90734863281018696e-06f // atan(2^-19)
};

static const data_t K_inv = 0.6072529350088813f;
static const data_t HALF_PI = (data_t)1.5707963267948966f;
static const data_t PI      = (data_t)3.1415926535897932f;
static const data_t ZERO    = (data_t)0.0f;

void cordiccart2pol(data_t x, data_t y, data_t * r,  data_t * theta)
{
    data_t xi = x, yi = y, zi = ZERO;
    if (yi >= ZERO) {
        data_t x_old = xi;
        xi =  yi;
        yi = -x_old;
        zi -= HALF_PI;
    } else {
        data_t x_old = xi;
        xi = -yi;
        yi =  x_old;
        zi += HALF_PI;
    }

    for (int i = 0; i < NO_ITER; ++i) {
        data_t di = (yi >= ZERO) ? (data_t)1.0f : (data_t)-1.0f;
        data_t x_next = xi + di * yi * Kvalues[i];
        data_t y_next = yi - di * xi * Kvalues[i];
        data_t z_next = zi + di * angles[i];
        xi = x_next; yi = y_next; zi = z_next;
    }

    *r = K_inv * ((xi >= ZERO) ? xi : -xi);

    *theta = zi + ((y >= ZERO) ? PI : -PI);
}