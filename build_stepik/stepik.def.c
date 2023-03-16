#include "slae.h"


int main()
{
  size_t i;
  struct SLAE slae;
  double * x;

  SLAE_create_from_stream_height_first(&slae, stdin);
  x = alloca(slae.W * sizeof (double));

  switch (SLAE_solve(&slae, x))
  {
    case SLAE_INF: puts("INF"); break;
    case SLAE_NO: puts("NO"); break;
    case SLAE_ONE: puts("YES");
      for (i = 0; i < slae.W; i++)
        printf("%lf ", x[i]);
      putchar('\n');
    break;
  }

  SLAE_free(&slae);

  return 0;
}
