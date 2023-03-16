#include "slae.h"

#include <err.h>


void process_args(int argc, char * argv[], const char ** pfilename)
{
  if (argc != 1 and argc != 2)
    errx(EXIT_SUCCESS, "Usage: ./main [<filename>]");

  if (argc == 2)
    *pfilename = argv[1];
  else
    *pfilename = NULL;
}


int main(int argc, char * argv[])
{
  const char * filename;
  struct timeval tv;
  ushort wdim, hdim;
  struct SLAE slae;
  double * x, * diff;
  size_t i;

  process_args(argc, argv, &filename);

  if (filename != NULL)
  {
    if (not SLAE_create_from_file(&slae, filename))
      return 1;
  }
  else
  {
    wdim = rand() % 10;
    hdim = wdim + (rand() > RAND_MAX/2 ? 1 : -1) * rand()%2;
    SLAE_create_random(&slae, wdim, hdim, 100);
  }

  puts("Got SLAE:");
  SLAE_fprint(&slae, stdout);

  x = alloca(slae.W * sizeof (double));
  diff = alloca(slae.W * sizeof (double));

  switch (SLAE_solve(&slae, x))
  {
    case SLAE_NO:
      puts("\033[1;31mNo solution exists.\033[m");
      break;
    case SLAE_INF:
      puts("\033[1;31mInfinitely many solutions exist.\033[m");
      break;
    case SLAE_ONE:
      puts("\033[1;32mSolution:");
      vec_fprint(x, slae.W, stdout);
      printf("\033[m");
      fflush(stdout);
      vec_dot(diff, x, slae.A, slae.W);
      vec_subtract(diff, slae.b, slae.W);
      puts("Difference vector:");
      vec_fprint(diff, slae.W, stdout);
      break;
  }

  return 0;
}

