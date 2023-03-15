#include <criterion/criterion.h>
#include <sys/time.h>

#include "slae.h"


void set_random_seed()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  srand(tv.tv_usec);
}


TestSuite(one_simple_solution);


Test(one_simple_solution, one_dimensional0)
{
  enum { w = 1, h = 1 };
  const double A[w*h] = {1},
               b[h] = {0};
  double x[w] = {};
  struct SLAE slae;
  enum SLAE_STATUS status;

  SLAE_create(&slae, w, h);
  SLAE_init(&slae, A, b);
  // SLAE_solve(&slae, x);
  SLAE_free(&slae);
  // cr_assert(x[0] = 0);
}

Test(one_simple_solution, two_dimensional0)
{
  enum { w = 1, h = 1 };
  const double A[w*h] = {1},
               b[h] = {0};
  double x[w] = {};
  struct SLAE slae;
  enum SLAE_STATUS status;

  SLAE_create(&slae, w, h);
  SLAE_init(&slae, A, b);
  // SLAE_solve(&slae, x);
  SLAE_free(&slae);
  // cr_assert(x[0] = 0);
}

Test(one_simple_solution, three_dimensional0)
{
  enum { w = 3, h = 3 };
  const double A[w*h] = {
    1, 0, 2,
    3, 2, 1,
    1, 2, 0,
  },
               b[h] = {1, 2, 1};
  double x[w] = {};
  struct SLAE slae;
  enum SLAE_STATUS status;

  SLAE_create(&slae, w, h);
  SLAE_init(&slae, A, b);
  SLAE_solve(&slae, x);
  SLAE_free(&slae);
  // cr_assert(x[0] = 0);
}
