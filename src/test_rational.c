#include "rational.h"

#include <stdio.h>
#include <stdlib.h>
#include <err.h>
#include <iso646.h>


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
  struct rational a, b;
  double a_d, b_d, division;

  rational_fread(&a, stdin);
  rational_fread(&b, stdin);

  a_d = rational_to_double(&a);
  b_d = rational_to_double(&b);

  rational_fprint(&a, stdout);
  printf("~= %lf\n", a_d);
  rational_fprint(&b, stdout);
  printf("~= %lf\n", b_d);

  division = a_d / b_d;
  rational_divide(&a, &b);
  printf("a/b = ");
  rational_fprint(&a, stdout);
  printf(" = %lf\n", rational_to_double(&a));
  printf("double division: %lf\n", division);

  return 0;
}

