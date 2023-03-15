void print_slae_solution_and_diff_vector(const struct SLAE * slae, const long double * x)
{
  long double diff[3];

  puts("Solution:");
  vec_fprint(x, slae->D, stdout);

  vec_dot(diff, x, slae->A, slae->D);
  puts("A*x:");
  vec_fprint(diff, slae->D, stdout);

  vec_subtract(diff, slae->b, slae->D);
  puts("Difference vector:");
  vec_fprint(diff, slae->D, stdout);
}

void test_slae(int argc, char * argv[])
{
  const char * filename;
  struct SLAE slae;
  enum SLAE_STATUS status;
  long double x[3];
  ushort dim;

  process_args(argc, argv, &filename);
  if (filename != NULL)
    SLAE_create_from_file(&slae, "slae.txt");
  else
  {
    dim = 1 + rand() % 5;
    SLAE_create_random(&slae, dim, 10000);
    // puts("Enter the SLAE:");
    // SLAE_init_from_stream(&SLAE, stdin);
  }

  puts("Got SLAE:");
  SLAE_fprint(&slae, stdout);

  switch (SLAE_solve(&slae, x))
  {
    case SLAE_ONE:
      print_slae_solution_and_diff_vector(&slae, x);
      break;
    case SLAE_NO:
      puts("No solution exists.");
      break;
    case SLAE_INF:
      puts("Infinite number of solutions exists.");
      break;
  }
  SLAE_fprint(&slae, stdout);
}
