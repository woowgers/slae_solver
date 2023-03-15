#include "slae.h"


__attribute__((always_inline))
static inline void double_swap(double * a, double * b)
{
  double tmp;
  tmp = *a;
  *a = *b;
  *b = tmp;
}

__attribute__((always_inline))
static inline void size_swap(size_t * a, size_t * b)
{
  size_t tmp;
  tmp = *a;
  *a = *b;
  *b = tmp;
}

__attribute__((always_inline))
static inline short random_sign()
{
  return rand() > (RAND_MAX / 2) ? 1 : -1;
}


static void * _valloc(size_t n_bytes, const char * format_on_error, va_list va_args)
{
  void * ptr = malloc(n_bytes);
  if (ptr == NULL)
    vfprintf(stderr, format_on_error, va_args);
  return ptr;
}

static void * alloc(size_t n_bytes, const char * format_on_error, ...)
{
  void * ptr;
  va_list args;
  va_start(args, format_on_error);
  ptr = _valloc(n_bytes, format_on_error, args);
  va_end(args);
  return ptr;
}

static size_t ndigits(double value)
{
  size_t n_digits = 0, // 0.0
         value_int = value;

  if (value == -0.0f)
    return 2;

  if (value < 0)
  {
    n_digits += 1;
    value_int = -value_int;
  }

  while (value_int > 0)
  {
    n_digits++;
    value_int /= 10;
  }

  return n_digits;
}

static ushort vec_width(const double * vec, ushort dim)
{
  ushort width = ndigits(vec[0]),
         value;

  for (size_t i = 1; i < dim; i++)
  {
    value = ndigits(vec[i]);
    if (value > width)
      width = value;
  }

  return 3 + width;
}


void vec_dot(double * dot, const double * vec, const double * mat, ushort dim)
{
  size_t i, j;

  memset(dot, 0, dim * sizeof (double));
  for (i = 0; i < dim; i++)
    for (j = 0; j < dim; j++)
      dot[i] += vec[j] * MAT(mat, dim, i, j);
}

void vec_subtract(double * self, const double * vec, ushort dim)
{
  size_t i;

  for (i = 0; i < dim; i++)
    self[i] -= vec[i];
}

void vec_fprint(const double * vec, ushort dim, FILE * stream)
{
  size_t i;

  if (dim == 1)
  {
    fprintf(stream, "(%.2lf)\n", vec[0]);
    return;
  }

  fprintf(stream, "⎛%e⎞\n", vec[0]);
  for (i = 1; i + 1 < dim; i++)
    fprintf(stream, "⎜%e⎟\n", vec[i]);
  fprintf(stream, "⎝%e⎠\n", vec[dim-1]);
}

static ushort SLAE_col_width(const struct SLAE * self, size_t i_col)
{
  ushort width = ndigits(SLAE_A(self, 0, i_col)),
         value;

  for (size_t i = 1; i < self->H; i++)
  {
    value = ndigits(SLAE_A(self, i, i_col));
    if (value > width)
      width = value;
  }

  return 3 + width;
}

static ushort SLAE_row_width(const struct SLAE * self, size_t i_row)
{
  ushort width = ndigits(SLAE_A(self, i_row, 0)),
         value;

  for (size_t i = 1; i < self->W; i++)
  {
    value = ndigits(SLAE_A(self, i_row, i));
    if (value > width)
      width = value;
  }

  return 3 + width;
}

static ushort SLAE_b_width(const struct SLAE * self)
{
  ushort width = ndigits(self->b[0]),
         value;

  for (size_t i = 1; i < self->W; i++)
  {
    value = ndigits(self->b[i]);
    if (value > width)
      width = value;
  }

  return 3 + width;
}

static void SLAE_init_xi(struct SLAE * self)
{
  for (size_t i = 0; i < self->W; i++)
    self->xi[i] = i;
}

__attribute__((always_inline))
static inline void SLAE_fprint_first_row(const struct SLAE * self, FILE * stream, ushort * width)
{
  size_t j;

  fprintf(stream, "⎛%*.2lf", width[0], SLAE_A(self, 0, 0));
  for (j = 1; j < self->W; j++)
    fprintf(stream, " %*.2lf", width[j], SLAE_A(self, 0, j));
  fprintf(stream, "⎟%*.2lf⎞\n", width[self->W-1], SLAE_B(self, 0));

}

__attribute__((always_inline))
static inline void SLAE_fprint_last_row(const struct SLAE * self, FILE * stream, ushort * width)
{
  const size_t i = self->H - 1;
  size_t j;

  fprintf(stream, "⎝%*.2lf", width[0], SLAE_A(self, i, 0));
  for (j = 1; j < self->W; j++)
    fprintf(stream, " %*.2lf", width[j], SLAE_A(self, i, j));
  fprintf(stream, "⎟%*.2lf⎠\n", width[self->W - 1], SLAE_B(self, i));
}

__attribute__((always_inline))
static inline void SLAE_fprint_row(const struct SLAE * self, FILE * stream, size_t i_row, ushort * width)
{
  size_t j;

  fprintf(stream, "⎟%*.2lf", width[0], SLAE_A(self, i_row, 0));
  for (j = 1; j < self->W; j++)
    fprintf(stream, " %*.2lf", width[j], SLAE_A(self, i_row, j));
  fprintf(stream, "⎟%*.2lf⎟\n", width[self->W - 1], SLAE_B(self, i_row));
}

static void SLAE_fprint_x(const struct SLAE * self, FILE * stream)
{
  size_t i;

  fprintf(stream, "(x%zu; ", self->xi[0]);
  for (i = 1; i + 1 < self->W; i++)
    fprintf(stream, "x%zu; ", self->xi[i]);
  if (self->W != 1)
    fprintf(stream, "x%zu", self->xi[i]);
  fputs(")\n", stream);
}

static void SLAE_permute_rows(struct SLAE * self, size_t i, size_t j)
{
#if defined(DEBUG)
  assert(i < self->W and j < self->H);
#endif

  if (i == j)
    return;

  for (size_t k = 0; k < self->W; k++)
    double_swap(&SLAE_A(self, i, k), &SLAE_A(self, j, k));
  double_swap(&SLAE_B(self, i), &SLAE_B(self, j));
}

static bool SLAE_alloc(struct SLAE * self, ushort w, ushort h)
{
  self->W = w;
  self->H = h;
  self->HW = w * h;

  self->A = alloc(self->HW * sizeof (double), "Failed to allocate memory for SLAE %zux%zu matrix.\n", self->W, self->H);
  if (self->A == NULL)
    return false;

  self->b = alloc(self->H * sizeof (double), "Failed to allocate memory for SLAE %zu b-vector.\n", self->H);
  if (self->b == NULL)
  {
    free(self->A);
    return false;
  }

  self->xi = alloc(self->W * sizeof (size_t), "Failed to allocate memory for SLAE %zu xi-vector.\n", self->W);
  if (self->xi == NULL)
  {
    free(self->A);
    free(self->b);
    return false;
  }

  return true;
}

static void SLAE_remove_row(struct SLAE * self, size_t i)
{
#if defined(DEBUG)
  assert(i < self->H);
  assert(self->H != 0);
#endif

  SLAE_permute_rows(self, i, self->H - 1);
  self->H--;
}

static void SLAE_subtract_multiplied_row(
  struct SLAE * self,
  size_t i_minuend,
  size_t i_subtrahend,
  double coeff
)
{
  for (size_t i = 0; i < self->W; i++)
    SLAE_A(self, i_minuend, i) -= coeff * SLAE_A(self, i_subtrahend, i);
  SLAE_B(self, i_minuend) -= coeff * SLAE_B(self, i_subtrahend);
}

static bool SLAE_row_is_zeros(const struct SLAE * self, size_t i)
{
  static size_t j; 

  for (j = 0; j < self->W; j++)
    if (not DOUBLE_IS_ZERO(SLAE_A(self, i, j)))
      return false;
  return true;
}

static size_t SLAE_i_diag_minor_i_maxabs(struct SLAE * self, size_t i)
{
  double maxabs = fabsl(SLAE_A(self, i, i)),
         value;
  size_t i_maxabs = i;

  for (size_t j = i; j < self->H; j++)
  {
    value = fabsl(SLAE_A(self, j, i));
    if (value > maxabs)
    {
      maxabs = value;
      i_maxabs = j;
    }
  }

  return i_maxabs;
}

static enum SLAE_STATUS SLAE_gauss_forward_i(struct SLAE * self, size_t i)
{
  static size_t j, i_maxabs;
  static double k;
  bool row_is_zeros;

  i_maxabs = SLAE_i_diag_minor_i_maxabs(self, i);
  if (i_maxabs != i)
    SLAE_permute_rows(self, i, i_maxabs);

  if (SLAE_A(self, i, i) == 0)
  {
    row_is_zeros = SLAE_row_is_zeros(self, i);
    if (row_is_zeros and DOUBLE_IS_ZERO(SLAE_B(self, i)))
      SLAE_remove_row(self, i);
    else if (row_is_zeros)
      return SLAE_NO;
    return SLAE_INF;
  }

  for (j = i + 1; j < self->H; j++)
  {
    k = SLAE_A(self, j, i) / SLAE_A(self, i, i);
    SLAE_subtract_multiplied_row(self, j, i, k);
  }

  return SLAE_ONE;
}

static enum SLAE_STATUS SLAE_gauss_forward(struct SLAE * self)
{
  size_t i;
  enum SLAE_STATUS status;

  for (i = 0; i < self->H; i++)
    if ((status = SLAE_gauss_forward_i(self, i)) != SLAE_ONE)
      return status;

  return SLAE_ONE;
}

static enum SLAE_STATUS SLAE_gauss_backward_i(struct SLAE * self, size_t i)
{
  static double k;
  static size_t j;
  static bool row_is_zeros;
  
  for (j = i - 1; j + 1 >= 1; j--)
  {
    k = SLAE_A(self, j, i) / SLAE_A(self, i, i);
    SLAE_subtract_multiplied_row(self, j, i, k);
  }

  return SLAE_ONE;
}

static enum SLAE_STATUS SLAE_gauss_backward(struct SLAE * self)
{
  enum SLAE_STATUS status;
  size_t i;

  for (i = self->H-1; i + 1 >= 1; i--)
    if ((status = SLAE_gauss_backward_i(self, i)) != SLAE_ONE)
      return status;

  return SLAE_ONE;
}

static enum SLAE_STATUS SLAE_gauss_solve(struct SLAE * self, double * x)
{
  enum SLAE_STATUS status;
  size_t i;

  if ((status = SLAE_gauss_forward(self)) != SLAE_ONE)
    return status;

  if (self->H > self->W)
    return SLAE_NO;
  else if (self->W > self->H)
    return SLAE_INF;

  if ((status = SLAE_gauss_backward(self)) != SLAE_ONE)
    return status;

  for (i = 0; i < self->W; i++)
    x[self->xi[i]] = SLAE_B(self, i) / SLAE_A(self, i, i);

  return SLAE_ONE;
}


bool SLAE_create(struct SLAE * self, ushort w, ushort h)
{
  return SLAE_alloc(self, w, h);
}

void SLAE_init(struct SLAE * self, const double * A, const double * b)
{
  memcpy(self->A, A, self->HW * sizeof (double));
  memcpy(self->b, b, self->H * sizeof (double));
}

bool SLAE_create_from_stream(struct SLAE * self, FILE * stream)
{
  size_t i, j;

  if (fscanf(stream, "%hu %hu", &self->W, &self->H) != 2)
  {
    fprintf(stderr, "Failed to read SLAE size from stream\n");
    return false;
  }

  if (not SLAE_alloc(self, self->W, self->H))
    return false;
  SLAE_init_xi(self);

  for (i = 0; i < self->H; i++)
  {
    for (j = 0; j < self->W; j++)
      if (fscanf(stream, "%lf", &SLAE_A(self, i, j)) != 1)
      {
        fprintf(stderr, "Failed to read SLAE from stream.\n");
        return false;
      }
    if (fscanf(stream, "%lf", &SLAE_B(self, i)) != 1)
    {
      fprintf(stderr, "Failed to read SLAE from stream.\n");
      return false;
    }
  }

  return true;
}

bool SLAE_create_from_file(struct SLAE * self, const char * filename)
{
  FILE * fp;

  if ((fp = fopen(filename, "r")) == NULL)
  {
    fprintf(stderr, "Failed to read SLAE from file \"%s\"\n", filename);
    return false;
  }

  if (not SLAE_create_from_stream(self, fp))
  {
    fclose(fp);
    return false;
  }

  fclose(fp);
  return true;
}

void SLAE_create_random(struct SLAE * self, size_t wdim, size_t hdim, uint modulo)
{
  SLAE_alloc(self, wdim, hdim);
  SLAE_init_xi(self);

  for (size_t i = 0; i < self->HW; i++)
    self->A[i] = random_sign() * rand() % modulo;
  for (size_t i = 0; i < self->H; i++)
    self->b[i] = random_sign() * rand() % modulo;
}

void SLAE_free(struct SLAE * self)
{
  free(self->A);
  free(self->b);
  free(self->xi);
}

void SLAE_fprint(const struct SLAE * self, FILE * stream)
{
  size_t i, j;
  ushort * width;

  if (unlikely(self->H == 0))
  {
    fprintf(stream, "() ()\n");
    return;
  }
  else if (unlikely(self->H == 1))
  {
    fprintf(stream, "(%.2lf⎟%.2lf) (x0)\n", SLAE_A(self, 0, 0), SLAE_B(self, 0));
    return;
  }

  width = alloca((self->W + 1) * sizeof (ushort));

  for (i = 0; i < self->W; i++)
    width[i] = SLAE_col_width(self, i);
  width[i] = SLAE_b_width(self);

  SLAE_fprint_x(self, stream);
  SLAE_fprint_first_row(self, stream, width);
  for (i = 1; i + 1 < self->H; i++)
    SLAE_fprint_row(self, stream, i, width);
  SLAE_fprint_last_row(self, stream, width);
}

enum SLAE_STATUS SLAE_solve(struct SLAE * self, double * x)
{
  return SLAE_gauss_solve(self, x);
}

