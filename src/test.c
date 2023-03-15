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
  const double A[] = {1},
               b[] = {0};
  double x[] = {};
  const size_t w = 1, h = 1;
  struct SLAE slae;
  enum SLAE_STATUS status;

  SLAE_create(&slae, w, h);
  SLAE_init(&slae, A, b);
  // SLAE_solve(&slae, x);
  SLAE_free(&slae);
  // cr_assert(x[0] = 0);
}

Test(one_simple_solution, test1)
{
  const double A[] = {1},
               b[] = {0};
  double x[] = {};
  const size_t w = 1, h = 1;
  struct SLAE slae;
  enum SLAE_STATUS status;

  SLAE_create(&slae, w, h);
  SLAE_init(&slae, A, b);
  // SLAE_solve(&slae, x);
  SLAE_free(&slae);
  // cr_assert(x[0] == 0);
}
