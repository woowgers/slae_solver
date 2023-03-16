#include "slae.h"


int main()
{
  const char * filename;
  struct SLAE slae;
  ushort wdim, hdim;
  uint modulo;

  srand(time(NULL));

  modulo = 100;
  wdim = 1 + rand()%10;
  hdim = 1 + rand()%10;
  SLAE_create_random(&slae, wdim, hdim, 100);

  puts("Got SLAE:");
  SLAE_fprint(&slae, stdout);
  SLAE_remove_row(&slae, 0);
  puts("After removing 0 row:");
  SLAE_fprint(&slae, stdout);

  SLAE_free(&slae);

  return 0;
}
