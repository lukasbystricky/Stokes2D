#ifndef SQ_COMMON
#define SQ_COMMON

#include <string.h>
#include "mex.h"
#include "omp.h"
#include "Complex.h"

#define pi 3.1415926535897932385

extern const double W16[16];
extern const double W32[32];
extern const double IP1[128];
extern const double IP2[128];

void vandernewton(Complex *T, Complex *b, int N);
void vandernewtonT(double *T, double *b, int N);
void IPmultR(Complex* in,Complex* out);

#endif
