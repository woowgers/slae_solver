#include "rational.h"


static int gcd(int a, int b)
{
  int temp;

  while (b != 0)
  {
    temp = a % b;
    a = b;
    b = temp;
  }

  return a;
}


struct rational rational(int num, int denom)
{
  int _gcd = gcd(num, denom);
  return (struct rational){ num/_gcd, denom/_gcd };
}

struct rational rational_inverse(const struct rational * self)
{
  struct rational inverse = (struct rational) { self->denom, self->num };
  rational_normalize(&inverse);
  return inverse;
}


void rational_multiply(struct rational * self, const struct rational * multiplier)
{
  self->num *= multiplier->num;
  self->denom *= multiplier->denom;
  rational_normalize(self);
}

void rational_multiply_by_int(struct rational * self, int multiplier)
{
  self->num *= multiplier;
  rational_normalize(self);
}

double rational_to_double(const struct rational * self)
{
  return (double)self->num / self->denom;
}

void rational_normalize(struct rational * self)
{
  int _gcd = gcd(self->num, self->denom);
  self->num /= _gcd;
  self->denom /= _gcd;
}

void rational_add(struct rational * self, const struct rational * delta)
{
  struct rational _delta = *delta;
  int _gcd = gcd(self->denom, delta->denom),
      m1 = _gcd / self->denom,
      m2 = _gcd / delta->denom;
  rational_multiply_by_int(self, m1);
  rational_multiply_by_int(&_delta, m2);
}

void rational_divide(struct rational * self, const struct rational * divisor)
{
  self->num *= divisor->denom;
  self->denom *= divisor->num;
}


void rational_fprint(const struct rational * self, FILE * stream)
{
  fprintf(stream, "%d/%d", self->num, self->denom);
}

bool rational_fread(struct rational * self, FILE * stream)
{
  if (fscanf(stream, "%d%d", &self->num, &self->denom) != 2)
    return false;
  rational_normalize(self);
  return true;
}
