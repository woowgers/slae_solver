#include <inttypes.h>
#include <stdint.h>
#include <assert.h>
#include <stdio.h>
#include <stdbool.h>



struct rational
{
  int num, denom;
};


struct rational rational(int num, int denom);
struct rational rational_inverse(const struct rational * self);
void rational_multiply(struct rational * self, const struct rational * multiplier);
void rational_multiply_by_int(struct rational * self, int multiplier);
double rational_to_double(const struct rational * self);
void rational_normalize(struct rational * self);
void rational_add(struct rational * self, const struct rational * delta);
void rational_divide(struct rational * self, const struct rational * divisor);
void rational_fprint(const struct rational * self, FILE * stream);
bool rational_fread(struct rational * self, FILE * stream);
