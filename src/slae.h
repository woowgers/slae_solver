#include <stdio.h>
#include <iso646.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdlib.h>
#include <limits.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#if defined(DEBUG)
#include <assert.h>
#endif
#include <float.h>
#include <stdarg.h>
#include <alloca.h>


#if !defined(ushort)
typedef unsigned short ushort;
#endif

#if !defined(uint)
typedef unsigned int uint;
#endif



#define max(a, b) ((a) > (b) ? (a) : (b))
#define max3(a, b, c) (max(a, b) > max(b, c) ? max(a, b) : max(b, c))
#define max4(a, b, c, d) (max(a, b) > max(c, d) ? max(a, b) : max(c, d))

#define DOUBLE_IS_ZERO(x) (fabsl(x) < DBL_MIN)
#define DOUBLE_EQUALS(x, y) (DOUBLE_IS_ZERO(x - y))

#define SLAE3_A(self, i, j) (self->A[3*i + j])
#define SLAE3_B(self, i) (self->b[i])
#define SLAE3_XI(self, i) (self->xi[i])

#define SLAE_A(self, i, j) (self->A[self->W * i + j])
#define SLAE_B(self, i) (self->b[i])
#define SLAE_XI(self, i) (self->xi[i])

#define MAT3X3_GET(self, i, j) (self[3*i + j])
#define MAT(self, dim, i, j) (self[dim*i + j])

#define likely(x) (__builtin_expect(x, 1))
#define unlikely(x) (__builtin_expect(x, 0))


enum SLAE_STATUS
{
  SLAE_ONE,
  SLAE_NO,
  SLAE_INF
};

struct SLAE
{
  double * A,
         * b;
  size_t * xi;
  ushort H, W;
  uint HW;
};


void vec_dot(double * dot, const double * vec, const double * mat, ushort dim);
void vec_subtract(double * self, const double * vec, ushort dim);
void vec_fprint(const double * vec, ushort dim, FILE * stream);

bool SLAE_create(struct SLAE * self, ushort w, ushort h);
void SLAE_init(struct SLAE * self, const double * A, const double * b);
bool SLAE_create_from_stream(struct SLAE * self, FILE * stream);
bool SLAE_create_from_stream_height_first(struct SLAE * self, FILE * stream);
bool SLAE_create_from_file(struct SLAE * self, const char * filename);
void SLAE_create_random(struct SLAE * self, size_t wdim, size_t hdim, uint modulo);
void SLAE_free(struct SLAE * self);
void SLAE_fprint(const struct SLAE * self, FILE * stream);
enum SLAE_STATUS SLAE_solve(struct SLAE * self, double * x);
void SLAE_remove_row(struct SLAE * self, size_t i);
