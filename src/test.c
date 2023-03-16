#include <criterion/criterion.h>
#include <sys/time.h>

#include "slae.h"


void set_random_seed()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  srand(tv.tv_usec);
}


TestSuite(one_solution);


Test(one_solution, one_dim0)
{
  enum { w = 1, h = 1 };
  const double A[w*h] = {1},
               b[h] = {0};
  double x[w] = {};
  struct SLAE slae;
  enum SLAE_STATUS status;

  SLAE_create(&slae, w, h);
  SLAE_init(&slae, A, b);
  status = SLAE_solve(&slae, x);
  SLAE_free(&slae);

  cr_assert_eq(status, SLAE_ONE);
  cr_assert_float_eq(x[0], 0, DBL_MIN);
}

Test(one_solution, two_dim0)
{
  enum { w = 1, h = 1 };
  const double A[w*h] = {1},
               b[h] = {0};
  double x[w] = {};
  struct SLAE slae;
  enum SLAE_STATUS status;

  SLAE_create(&slae, w, h);
  SLAE_init(&slae, A, b);
  status = SLAE_solve(&slae, x);
  SLAE_free(&slae);

  cr_assert_eq(status, SLAE_ONE);
  cr_assert_float_eq(x[0], 0, DBL_MIN);
}

Test(one_solution, three_dim0)
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
  status = SLAE_solve(&slae, x);
  SLAE_free(&slae);

  cr_assert_eq(status, SLAE_ONE);
  for (size_t i = 0; i < w; i++)
    cr_assert_float_eq(x[i], (double)1/3, DBL_MIN);
}

Test(one_solution, threeXfour_dim0)
{
  enum { w = 3, h = 4 };
  const double A[w*h] = {
    1, 0, 0,
    0, 1, 0,
    0, 0, 1,
    0, 0, 1,
  },
               b[h] = {1, 1, 1, 1};
  double x[w] = {};
  struct SLAE slae;
  enum SLAE_STATUS status;

  SLAE_create(&slae, w, h);
  SLAE_init(&slae, A, b);
  status = SLAE_solve(&slae, x);
  SLAE_free(&slae);

  cr_assert_eq(status, SLAE_ONE);
  for (size_t i = 0; i < w; i++)
    cr_assert_float_eq(x[i], 1, DBL_MIN);
}


TestSuite(no);


Test(no, one_dim0)
{
  enum { w = 1, h = 1 };
  const double A[w*h] = {0},
               b[h] = {1};
  double x[w] = {};
  struct SLAE slae;
  enum SLAE_STATUS status;

  SLAE_create(&slae, w, h);
  SLAE_init(&slae, A, b);
  status = SLAE_solve(&slae, x);
  SLAE_free(&slae);

  cr_assert_eq(status, SLAE_NO);
}

Test(no, two_dim0)
{
  enum { w = 2, h = 2 };
  const double A[w*h] = {
    0, 0,
    0, 0
  },
               b[h] = {1, 1};
  double x[w] = {};
  struct SLAE slae;
  enum SLAE_STATUS status;

  SLAE_create(&slae, w, h);
  SLAE_init(&slae, A, b);
  status = SLAE_solve(&slae, x);
  SLAE_free(&slae);

  cr_assert_eq(status, SLAE_NO);
}

Test(no, three_dim0)
{
  enum { w = 3, h = 3 };
  const double A[w*h] = {
    0, 0, 0,
    0, 0, 0,
    0, 0, 0
  },
               b[h] = {1, 1, 1};
  double x[w] = {};
  struct SLAE slae;
  enum SLAE_STATUS status;

  SLAE_create(&slae, w, h);
  SLAE_init(&slae, A, b);
  status = SLAE_solve(&slae, x);
  SLAE_free(&slae);

  cr_assert_eq(status, SLAE_NO);
}

Test(no, threeXone_dim0)
{
  enum { w = 3, h = 1 };
  const double A[w*h] = {
    0, 0, 0,
  },
               b[h] = {1};
  double x[w] = {};
  struct SLAE slae;
  enum SLAE_STATUS status;

  SLAE_create(&slae, w, h);
  SLAE_init(&slae, A, b);
  status = SLAE_solve(&slae, x);
  SLAE_free(&slae);

  cr_assert_eq(status, SLAE_NO);
}

Test(no, threeXtwo_dim0)
{
  enum { w = 3, h = 2 };
  const double A[w*h] = {
    0, 0, 0,
    0, 0, 0
  },
               b[h] = {1, 1};
  double x[w] = {};
  struct SLAE slae;
  enum SLAE_STATUS status;

  SLAE_create(&slae, w, h);
  SLAE_init(&slae, A, b);
  status = SLAE_solve(&slae, x);
  SLAE_free(&slae);

  cr_assert_eq(status, SLAE_NO);
}


TestSuite(infinite);


Test(infinite, one_dim0)
{
  enum { w = 1, h = 1 };
  const double A[w*h] = {0},
               b[h] = {0};
  double x[w] = {};
  struct SLAE slae;
  enum SLAE_STATUS status;

  SLAE_create(&slae, w, h);
  SLAE_init(&slae, A, b);
  status = SLAE_solve(&slae, x);
  SLAE_free(&slae);

  cr_assert_eq(status, SLAE_INF);
}

Test(infinite, two_dim0)
{
  enum { w = 2, h = 2 };
  const double A[w*h] = {
    0, 0,
    0, 0
  },
               b[h] = {0, 0};
  double x[w] = {};
  struct SLAE slae;
  enum SLAE_STATUS status;

  SLAE_create(&slae, w, h);
  SLAE_init(&slae, A, b);
  status = SLAE_solve(&slae, x);
  SLAE_free(&slae);

  cr_assert_eq(status, SLAE_INF);
}

Test(infinite, three_dim0)
{
  enum { w = 3, h = 3 };
  const double A[w*h] = {
    0, 0, 0,
    0, 0, 0,
    0, 0, 0
  },
               b[h] = {0, 0, 0};
  double x[w] = {};
  struct SLAE slae;
  enum SLAE_STATUS status;

  SLAE_create(&slae, w, h);
  SLAE_init(&slae, A, b);
  status = SLAE_solve(&slae, x);
  SLAE_free(&slae);

  cr_assert_eq(status, SLAE_INF);
}

Test(infinite, threeXone_dim0)
{
  enum { w = 3, h = 1 };
  const double A[w*h] = {
    0, 0, 0,
  },
               b[h] = {0};
  double x[w] = {};
  struct SLAE slae;
  enum SLAE_STATUS status;

  SLAE_create(&slae, w, h);
  SLAE_init(&slae, A, b);
  status = SLAE_solve(&slae, x);
  SLAE_free(&slae);

  cr_assert_eq(status, SLAE_INF);
}

Test(infinite, threeXtwo_dim0)
{
  enum { w = 3, h = 2 };
  const double A[w*h] = {
    0, 0, 0,
    0, 0, 0
  },
               b[h] = {0, 0};
  double x[w] = {};
  struct SLAE slae;
  enum SLAE_STATUS status;

  SLAE_create(&slae, w, h);
  SLAE_init(&slae, A, b);
  status = SLAE_solve(&slae, x);
  SLAE_free(&slae);

  cr_assert_eq(status, SLAE_INF);
}
