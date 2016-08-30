# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "matrix_exponential.h"

/******************************************************************************/

double *expm11 ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    EXPM11 is essentially MATLAB's built-in matrix exponential algorithm.

  Discussion:

    The GCC compiler feels that the name "expm1" belongs to it, so I 
    give up and use the ridiculous alternative of "expm11".

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 December 2011

  Author:

    Cleve Moler, Charles Van Loan

  Reference:

    Cleve Moler, Charles VanLoan,
    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
    Twenty-Five Years Later,
    SIAM Review,
    Volume 45, Number 1, March 2003, pages 3-49.

  Parameters:

    Input, int N, the dimension of the matrix.

    Input, double A[N*N], the matrix.

    Output, double EXPM1[N*N], the estimate for exp(A).
*/
{
  double *a2;
  double a_norm;
  double c;
  double *d;
  double *e;
  int ee;
  int k;
  const double one = 1.0;
  int p;
  const int q = 6;
  int s;
  double t;
  double *x;

  a2 = r8mat_copy_new ( n, n, a );

  a_norm = r8mat_norm_li ( n, n, a2 );

  ee = ( int ) ( r8_log_2 ( a_norm ) ) + 1;
  
  s = i4_max ( 0, ee + 1 );

  t = 1.0 / pow ( 2.0, s );

  r8mat_scale ( n, n, t, a2 );

  x = r8mat_copy_new ( n, n, a2 );

  c = 0.5;

  e = r8mat_identity_new ( n );

  r8mat_add ( n, n, one, e, c, a2, e );

  d = r8mat_identity_new ( n );

  r8mat_add ( n, n, one, d, -c, a2, d );

  p = 1;

  for ( k = 2; k <= q; k++ )
  {
    c = c * ( double ) ( q - k + 1 ) / ( double ) ( k * ( 2 * q - k + 1 ) );

    r8mat_mm ( n, n, n, a2, x, x );

    r8mat_add ( n, n, c, x, one, e, e );

    if ( p )
    {
      r8mat_add ( n, n, c, x, one, d, d );
    }
    else
    {
      r8mat_add ( n, n, -c, x, one, d, d );
    }

    p = !p;
  }
/*
  E -> inverse(D) * E
*/
  r8mat_minvm ( n, n, d, e, e );
/*
  E -> E^(2*S)
*/
  for ( k = 1; k <= s; k++ )
  {
    r8mat_mm ( n, n, n, e, e, e );
  }

  free ( a2 );
  free ( d );
  free ( x );

  return e;
}
/******************************************************************************/

double *expm2 ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    EXPM2 uses the Taylor series for the matrix exponential.

  Discussion:

    Formally,

      exp ( A ) = I + A + 1/2 A^2 + 1/3! A^3 + ...

    This function sums the series until a tolerance is satisfied.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 December 2011

  Author:

    Cleve Moler, Charles Van Loan

  Reference:

    Cleve Moler, Charles VanLoan,
    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
    Twenty-Five Years Later,
    SIAM Review,
    Volume 45, Number 1, March 2003, pages 3-49.

  Parameters:

    Input, int N, the dimension of the matrix.

    Input, double A[N*N], the matrix.

    Output, double EXPM2[N*N], the estimate for exp(A).
*/
{
  double *e;
  double *f;
  int k;
  const double one = 1.0;
  double s;

  e = r8mat_zero_new ( n, n );

  f = r8mat_identity_new ( n );

  k = 1;

  while ( r8mat_significant ( n, n, e, f ) )
  {
    r8mat_add ( n, n, one, e, one, f, e );

    r8mat_mm ( n, n, n, a, f, f );

    s = 1.0 / ( double ) ( k );

    r8mat_scale ( n, n, s, f );

    k = k + 1;
  }

  free ( f );

  return e;
}
/******************************************************************************/

double *expm3 ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    EXPM3 approximates the matrix exponential using an eigenvalue approach.

  Discussion:

    exp(A) = V * D * V

    where V is the matrix of eigenvectors of A, and D is the diagonal matrix
    whose i-th diagonal entry is exp(lambda(i)), for lambda(i) an eigenvalue
    of A.

    This function is accurate for matrices which are symmetric, orthogonal,
    or normal.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 December 2011

  Author:

    Cleve Moler, Charles Van Loan

  Reference:

    Cleve Moler, Charles VanLoan,
    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
    Twenty-Five Years Later,
    SIAM Review,
    Volume 45, Number 1, March 2003, pages 3-49.

  Parameters:

    Input, int N, the dimension of the matrix.

    Input, double A[N*N], the matrix.

    Output, double EXPM3[N*N], the estimate for exp(A).
*/
{
  double *e = NULL;
/*
  [ V, D ] = eig ( A );
  E = V * diag ( exp ( diag ( D ) ) ) / V;
*/
  return e;
}
/******************************************************************************/

int i4_max ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MAX returns the maximum of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, are two integers to be compared.

    Output, int I4_MAX, the larger of I1 and I2.
*/
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
/******************************************************************************/

int i4_min ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MIN returns the smaller of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, two integers to be compared.

    Output, int I4_MIN, the smaller of I1 and I2.
*/
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
/******************************************************************************/

double *mexp_a ( int test, int n )

/******************************************************************************/
/*
  Purpose:

    MEXP_A returns the matrix for a given test.

  Discussion:

     1) Diagonal example
     2) Symmetric example
     3) Laub
     4) Moler and Van Loan
     5) Moler and Van Loan
     6) Moler and Van Loan
     7) Moler and Van Loan
     8) Wikipedia example
     9) NAG F01ECF
    10) Ward #1
    11) Ward #2
    12) Ward #3
    13) Ward #4

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 November 2011

  Author:

    John Burkardt

  Reference:

    Alan Laub,
    Review of "Linear System Theory" by Joao Hespanha,
    SIAM Review,
    Volume 52, Number 4, December 2010, page 779-781.

    Cleve Moler, Charles VanLoan,
    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
    Twenty-Five Years Later,
    SIAM Review,
    Volume 45, Number 1, March 2003, pages 3-49.

    Robert Ward,
    Numerical computation of the matrix exponential with accuracy estimate,
    SIAM Journal on Numerical Analysis,
    Volume 14, Number 4, September 1977, pages 600-610.

  Parameters:

    Input, int TEST, the index of the test case.

    Input, int N, the order of the matrix.

    Output, double MEXP_A[N*N], the matrix.
*/
{
  double *a;
  static double a01[2*2] = {
      1.0, 0.0, 
      0.0, 2.0 };
  static double a02[2*2] = {
      1.0, 3.0, 
      3.0, 2.0 };
  static double a03[2*2] = {
      0.0, -39.0, 
      1.0, -40.0 };
  static double a04[2*2] = {
      -49.0, -64.0, 
       24.0,  31.0 };
  static double a05[4*4] = {
      0.0, 0.0, 0.0, 0.0, 
      6.0, 0.0, 0.0, 0.0, 
      0.0, 6.0, 0.0, 0.0, 
      0.0, 0.0, 6.0, 0.0 };
  static double a06[2*2] = {
      1.0, 0.0, 
      1.0, 1.0 };
  static double a08[3*3] = {
      21.0,  -5.0,   4.0, 
      17.0,  -1.0,   4.0, 
       6.0,  -6.0,  16.0 };
  static double a09[4*4] = {
      1.0, 3.0, 3.0, 3.0, 
      2.0, 1.0, 2.0, 3.0, 
      2.0, 1.0, 1.0, 3.0, 
      2.0, 2.0, 2.0, 1.0 };
  static double a10[3*3] = {
      4.0, 1.0, 1.0, 
      2.0, 4.0, 1.0, 
      0.0, 1.0, 4.0 };
  static double a11[3*3] = {
      29.87942128909879, 
       0.7815750847907159, 
      -2.289519314033932, 
       0.7815750847907159, 
      25.72656945571064, 
       8.680737820540137, 
      -2.289519314033932, 
       8.680737820540137, 
      34.39400925519054 };
  static double a12[3*3] = {
      -131.0, -390.0, -387.0, 
        19.0,   56.0,   57.0, 
        18.0,   54.0,   52.0 };
  int i;
  int j;

  if ( test == 1 )
  {
    a = r8mat_copy_new ( n, n, a01 );
  }
  else if ( test == 2 )
  {
    a = r8mat_copy_new ( n, n, a02 );
  }
  else if ( test == 3 )
  {
    a = r8mat_copy_new ( n, n, a03 );
  }
  else if ( test == 4 )
  {
    a = r8mat_copy_new ( n, n, a04 );
  }
  else if ( test == 5 )
  {
    a = r8mat_copy_new ( n, n, a05 );
  }
  else if ( test == 6 )
  {
    a = r8mat_copy_new ( n, n, a06 );
  }
  else if ( test == 7 )
  {
    a = ( double * ) malloc ( 2 * 2 * sizeof ( double ) );
    a[0+0*2] = 1.0 + r8_epsilon ( );
    a[1+0*2] = 0.0;
    a[0+1*2] = 0.0;
    a[1+1*2] = 1.0 - r8_epsilon ( );
  }
  else if ( test == 8 )
  {
    a = r8mat_copy_new ( n, n, a08 );
  }
  else if ( test == 9 )
  {
    a = r8mat_copy_new ( n, n, a09 );
  }
  else if ( test == 10 )
  {
    a = r8mat_copy_new ( n, n, a10 );
  }
  else if ( test == 11 )
  {
    a = r8mat_copy_new ( n, n, a11 );
  }
  else if ( test == 12 )
  {
    a = r8mat_copy_new ( n, n, a12 );
  }
  else if ( test == 13 )
  {
    a = ( double * ) malloc ( n * n * sizeof ( double ) );

    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < n; i++ )
      {
        if ( j == i + 1 )
        {
          a[i+j*n] = 1.0;
        }
        else if ( i == n - 1 && j == 0 )
        {
          a[i+j*n] = 1.0E-10;
        }
        else
        {
          a[i+j*n] = 0.0;
        }
      }
    }
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "MEXP_A - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of TEST = %d\n", test );
    exit ( 1 );
  }
  return a;
}
/******************************************************************************/

double *mexp_expa ( int test, int n )

/******************************************************************************/
/*
  Purpose:

    MEXP_EXPA returns the "exact" exponential matrix for a given test.

  Discussion:

    In some cases, the "exact" value is given to six significant digits.

     1) Diagonal example
     2) Symmetric example
     3) Laub
     4) Moler and Van Loan
     5) Moler and Van Loan
     6) Moler and Van Loan
     7) Moler and Van Loan
     8) Wikipedia example
     9) NAG F01ECF
    10) Ward #1
    11) Ward #2
    12) Ward #3
    13) Ward #4

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 November 2011

  Author:

    John Burkardt

  Reference:

    Alan Laub,
    Review of "Linear System Theory" by Joao Hespanha,
    SIAM Review,
    Volume 52, Number 4, December 2010, page 779-781.

    Cleve Moler, Charles VanLoan,
    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
    Twenty-Five Years Later,
    SIAM Review,
    Volume 45, Number 1, March 2003, pages 3-49.

    Robert Ward,
    Numerical computation of the matrix exponential with accuracy estimate,
    SIAM Journal on Numerical Analysis,
    Volume 14, Number 4, September 1977, pages 600-610.

  Parameters:

    Input, int TEST, the index of the test case.

    Input, int N, the order of the matrix.

    Output, double MEXP_EXPA[N*N], the exponential of the test matrix.
*/
{
  double exp16;
  double exp4;
  double *expa;
  static double expa01[2*2] = {
      2.718281828459046, 0.0, 
      0.0,               7.389056098930650 };
  static double expa02[2*2] = {
      39.322809708033859,  46.166301438885753, 
      46.166301438885768,  54.711576854329110 };
  static double expa03[2*2] = {
      0.0,               1.154822E-17,   
      2.718281828459046, 2.718281828459046 };
  static double expa04[2*2] = {
      -0.735759, -1.471518, 
       0.551819,  1.103638 };
  static double expa05[4*4] = {
      1.0,  0.0, 0.0, 0.0, 
      6.0,  1.0, 0.0, 0.0, 
     18.0,  6.0, 1.0, 0.0, 
     36.0, 18.0, 6.0, 1.0 };
  static double expa06[2*2] = {
      2.718281828459046, 0.0, 
      2.718281828459046, 2.718281828459046 };
  static double expa07[2*2] = {
      2.718309, 0.0, 
      2.718282, 2.718255 };
  static double expa09[4*4] = {
      740.7038, 731.2510, 823.7630, 998.4355, 
      610.8500, 603.5524, 679.4257, 823.7630, 
      542.2743, 535.0884, 603.5524, 731.2510, 
      549.1753, 542.2743, 610.8500, 740.7038 };
  static double expa10[3*3] = {
      147.8666224463699, 
      127.7810855231823, 
      127.7810855231824, 
      183.7651386463682, 
      183.7651386463682, 
      163.6796017231806, 
      71.79703239999647, 
      91.88256932318415, 
     111.9681062463718 };
  static double expa11[3*3] = {
     5.496313853692378E+15, 
    -1.823188097200899E+16, 
    -3.047577080858001E+16, 
    -1.823188097200898E+16, 
     6.060522870222108E+16, 
     1.012918429302482E+17, 
    -3.047577080858001E+16, 
     1.012918429302482E+17, 
     1.692944112408493E+17 };
  static double expa12[3*3] = {
    -1.509644158793135, 
    -5.632570799891469, 
    -4.934938326088363, 
     0.3678794391096522, 
     1.471517758499875, 
     1.103638317328798, 
     0.1353352811751005, 
     0.4060058435250609, 
     0.5413411267617766 };
  int i;
  int j;

  if ( test == 1 )
  {
    expa = r8mat_copy_new ( n, n, expa01 );
  }
  else if ( test == 2 )
  {
    expa = r8mat_copy_new ( n, n, expa02 );
  }
  else if ( test == 3 )
  {
    expa = r8mat_copy_new ( n, n, expa03 );
  }
  else if ( test == 4 )
  {
    expa = r8mat_copy_new ( n, n, expa04 );
  }
  else if ( test == 5 )
  {
    expa = r8mat_copy_new ( n, n, expa05 );
  }
  else if ( test == 6 )
  {
    expa = r8mat_copy_new ( n, n, expa06 );
  }
  else if ( test == 7 )
  {
    expa = r8mat_copy_new ( n, n, expa07 );
  }
  else if ( test == 8 )
  {
    expa = ( double * ) malloc ( 3 * 3 * sizeof ( double ) );
    exp16 = exp ( 16.0 );
    exp4 = exp ( 4.0 );
    expa[0+0*3] = 0.25 * ( 13.0 * exp16 -       exp4 );
    expa[1+0*3] = 0.25 * ( -9.0 * exp16 +       exp4 );
    expa[2+0*3] = 0.25 * ( 16.0 * exp16 );
    expa[0+1*3] = 0.25 * ( 13.0 * exp16 - 5.0 * exp4 );
    expa[1+1*3] = 0.25 * ( -9.0 * exp16 + 5.0 * exp4 );
    expa[2+1*3] = 0.25 * ( 16.0 * exp16 );
    expa[0+2*3] = 0.25 * (  2.0 * exp16 - 2.0 * exp4 );
    expa[1+2*3] = 0.25 * ( -2.0 * exp16 + 2.0 * exp4 );
    expa[2+2*3] = 0.25 * (  4.0 * exp16 );
  }
  else if ( test == 9 )
  {
    expa = r8mat_copy_new ( n, n, expa09 );
  }
  else if ( test == 10 )
  {
    expa = r8mat_copy_new ( n, n, expa10 );
  }
  else if ( test == 11 )
  {
    expa = r8mat_copy_new ( n, n, expa11 );
  }
  else if ( test == 12 )
  {
    expa = r8mat_copy_new ( n, n, expa12 );
  }
  else if ( test == 13 )
  {
    expa = ( double * ) malloc ( n * n * sizeof ( double ) );

    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < n; i++ )
      {
        expa[i+j*n] = 0.0;
      }
    }
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "MEXP_EXPA - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of TEST = %d\n", test );
    exit ( 1 );
  }
  return expa;
}
/******************************************************************************/

int mexp_n ( int test )

/******************************************************************************/
/*
  Purpose:

    MEXP_N returns the matrix order for a given test.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, int TEST, the index of the test case.

    Output, int MEXP_N, the order of the matrix.
*/
{
  int n;

  if ( test == 1 )
  {
    n = 2;
  }
  else if ( test == 2 )
  {
    n = 2;
  }
  else if ( test == 3 )
  {
    n = 2;
  }
  else if ( test == 4 )
  {
    n = 2;
  }
  else if ( test == 5 )
  {
    n = 4;
  }
  else if ( test == 6 )
  {
    n = 2;
  }
  else if ( test == 7 )
  {
    n = 2;
  }
  else if ( test == 8 )
  {
    n = 3;
  }
  else if ( test == 9 )
  {
    n = 4;
  }
  else if ( test == 10 )
  {
    n = 3;
  }
  else if ( test == 11 )
  {
    n = 3;
  }
  else if ( test == 12 )
  {
    n = 3;
  }
  else if ( test == 13 )
  {
    n = 10;
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "MEXP_N - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of TEST = %d\n", test );
    exit ( 1 );
  }

  return n;
}
/******************************************************************************/

void mexp_story ( int test )

/******************************************************************************/
/*
  Purpose:

    MEXP_STORY prints explanatory text for each problem.

  Discussion:

     1) Diagonal example
     2) Symmetric example
     3) Laub
     4) Moler and Van Loan
     5) Moler and Van Loan
     6) Moler and Van Loan
     7) Moler and Van Loan
     8) Wikipedia example
     9) NAG F01ECF
    10) Ward #1
    11) Ward #2
    12) Ward #3
    13) Ward #4

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 November 2011

  Author:

    John Burkardt

  Reference:

    Alan Laub,
    Review of "Linear System Theory" by Joao Hespanha,
    SIAM Review,
    Volume 52, Number 4, December 2010, page 779-781.

    Cleve Moler, Charles VanLoan,
    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
    Twenty-Five Years Later,
    SIAM Review,
    Volume 45, Number 1, March 2003, pages 3-49.

    Robert Ward,
    Numerical computation of the matrix exponential with accuracy estimate,
    SIAM Journal on Numerical Analysis,
    Volume 14, Number 4, September 1977, pages 600-610.

  Parameters:

    Input, int TEST, the index of the test case.
*/
{
  if ( test == 1 )
  {
    printf ( "\n" );
    printf ( "  This matrix is diagonal.\n" );
    printf ( "  The calculation of the matrix exponential is simple.\n" );
  }
  else if ( test == 2 )
  {
    printf ( "\n" );
    printf ( "  This matrix is symmetric.\n" );
    printf ( "  The calculation of the matrix exponential is straightforward.\n" );
  }
  else if ( test == 3 )
  {
    printf ( "\n" );
    printf ( "  This example is due to Laub.\n" );
    printf ( "  This matrix is ill-suited for the Taylor series approach.\n" );
    printf ( "  As powers of A are computed, the entries blow up too quickly.\n" );
  }
  else if ( test == 4 )
  {
    printf ( "\n" );
    printf ( "  This example is due to Moler and Van Loan.\n" );
    printf ( "  The example will cause problems for the series summation approach,\n" );
    printf ( "  as well as for diagonal Pade approximations.\n" );
  }
  else if ( test == 5 )
  {
    printf ( "\n" );
    printf ( "  This example is due to Moler and Van Loan.\n" );
    printf ( "  This matrix is strictly upper triangular\n" );
    printf ( "  All powers of A are zero beyond some (low) limit.\n" );
    printf ( "  This example will cause problems for Pade approximations.\n" );
  }
  else if ( test == 6 )
  {
    printf ( "\n" );
    printf ( "  This example is due to Moler and Van Loan.\n" );
    printf ( "  This matrix does not have a complete set of eigenvectors.\n" );
    printf ( "  That means the eigenvector approach will fail.\n" );
  }
  else if ( test == 7 )
  {
    printf ( "\n" );
    printf ( "  This example is due to Moler and Van Loan.\n" );
    printf ( "  This matrix is very close to example 5.\n" );
    printf ( "  Mathematically, it has a complete set of eigenvectors.\n" );
    printf ( "  Numerically, however, the calculation will be suspect.\n" );
  }
  else if ( test == 8 )
  {
    printf ( "\n" );
    printf ( "  This matrix was an example in Wikipedia.\n" );
  }
  else if ( test == 9 )
  {
    printf ( "\n" );
    printf ( "  This matrix is due to the NAG Library.\n" );
    printf ( "  It is an example for function F01ECF.\n" );
  }
  else if ( test == 10 )
  {
    printf ( "\n" );
    printf ( "  This is Ward's example #1.\n" );
    printf ( "  It is defective and nonderogatory.\n" );
    printf ( "  The eigenvalues are 3, 3 and 6.\n" );
  }
  else if ( test == 11 )
  {
    printf ( "\n" );
    printf ( "  This is Ward's example #2.\n" );
    printf ( "  It is a symmetric matrix.\n" );
    printf ( "  The eigenvalues are 20, 30, 40.\n" );
  }
  else if ( test == 12 )
  {
    printf ( "\n" );
    printf ( "  This is Ward's example #3.\n" );
    printf ( "  Ward's algorithm has difficulty estimating the accuracy\n" );
    printf ( "  of its results.  The eigenvalues are -1, -2, -20.\n" );
  }
  else if ( test == 13 )
  {
    printf ( "\n" );
    printf ( "  This is Ward's example #4.\n" );
    printf ( "  This is a version of the Forsythe matrix.\n" );
    printf ( "  The eigenvector problem is badly conditioned.\n" );
    printf ( "  Ward's algorithm has difficulty estimating the accuracy\n" );
    printf ( "  of its results for this problem.\n" );
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "MEXP_STORY - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of TEST = %d\n", test );
    exit ( 1 );
  }
  return;
}
/******************************************************************************/

int mexp_test_num ( )

/******************************************************************************/
/*
  Purpose:

    MEXP_TEST_NUM returns the number of matrix exponential tests.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 November 2011

  Author:

    John Burkardt

  Parameters:

    Output, int MEXP_TEST_NUM, the number of tests.
*/
{
  int test_num;

  test_num = 13;

  return test_num;
}
/******************************************************************************/

double r8_abs ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_ABS returns the absolute value of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the quantity whose absolute value is desired.

    Output, double R8_ABS, the absolute value of X.
*/
{
  double value;

  if ( 0.0 <= x )
  {
    value = + x;
  }
  else
  {
    value = - x;
  }
  return value;
}
/******************************************************************************/

double r8_add ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_ADD returns the sum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the numbers to be added.

    Output, double R8_ADD, the sum of X and Y.
*/
{
  double value;

  value = x + y;

  return value;
}
/******************************************************************************/

double r8_epsilon ( void )

/******************************************************************************/
/*
  Purpose:

    R8_EPSILON returns the R8 round off unit.

  Discussion:

    R8_EPSILON is a number R which is a power of 2 with the property that,
    to the precision of the computer's arithmetic,
      1 < 1 + R
    but
      1 = ( 1 + R / 2 )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2010

  Author:

    John Burkardt

  Parameters:

    Output, double R8_EPSILON, the R8 round-off unit.
*/
{
  double one;
  double temp;
  double test;
  double value;

  one = ( double ) ( 1 );

  value = one;
  temp = value / 2.0;
  test = r8_add ( one, temp );

  while ( one < test )
  {
    value = temp;
    temp = value / 2.0;
    test = r8_add ( one, temp );
  }

  return value;
}
/******************************************************************************/

double r8_huge ( void )

/******************************************************************************/
/*
  Purpose:

    R8_HUGE returns a "huge" R8.

  Discussion:

    The value returned by this function is NOT required to be the
    maximum representable R8.  This value varies from machine to machine,
    from compiler to compiler, and may cause problems when being printed.
    We simply want a "very large" but non-infinite number.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 October 2007

  Author:

    John Burkardt

  Parameters:

    Output, double R8_HUGE, a "huge" R8 value.
*/
{
  double value;

  value = 1.0E+30;

  return value;
}
/******************************************************************************/

double r8_log_2 ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_LOG_2 returns the logarithm base 2 of an R8.

  Discussion:

    R8_LOG_2 ( X ) = Log ( |X| ) / Log ( 2.0 )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the number whose base 2 logarithm is desired.
    X should not be 0.

    Output, double R8_LOG_2, the logarithm base 2 of the absolute
    value of X.  It should be true that |X| = 2**R_LOG_2.
*/
{
  double value;

  if ( x == 0.0 )
  {
    value = - r8_huge ( );
  }
  else
  {
    value = log ( fabs ( x ) ) / log ( 2.0 );
  }

  return value;
}
/******************************************************************************/

double r8_max ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MAX returns the maximum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MAX, the maximum of X and Y.
*/
{
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}
/******************************************************************************/

void r8mat_add ( int m, int n, double alpha, double a[], double beta, 
  double b[], double c[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_ADD computes C = alpha * A + beta * B for R8MAT's.

  Discussion:

    An R8MAT is an array of R8 values.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 December 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double ALPHA, the multiplier for A.

    Input, double A[M*N], the first matrix.

    Input, double BETA, the multiplier for A.

    Input, double B[M*N], the second matrix.

    Output, double C[M*N], the sum of alpha*A+beta*B.
*/
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      c[i+j*m] = alpha * a[i+j*m] + beta * b[i+j*m];
    }
  }
  return;
}
/******************************************************************************/

void r8mat_copy ( int m, int n, double a1[], double a2[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_COPY copies one R8MAT to another.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8's, which
    may be stored as a vector in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A1[M*N], the matrix to be copied.

    Output, double A2[M*N], the copy of A1.
*/
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
  }

  return;
}
/******************************************************************************/

double *r8mat_copy_new ( int m, int n, double a1[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_COPY_NEW copies one R8MAT to a "new" R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8's, which
    may be stored as a vector in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A1[M*N], the matrix to be copied.

    Output, double R8MAT_COPY_NEW[M*N], the copy of A1.
*/
{
  double *a2;
  int i;
  int j;

  a2 = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
  }

  return a2;
}
/******************************************************************************/

double *r8mat_fss_new ( int n, double a[], int nb, double b[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_FSS_NEW factors and solves a system with multiple right hand sides.

  Discussion:

    This routine uses partial pivoting, but no pivot vector is required.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input/output, double A[N*N].
    On input, A is the coefficient matrix of the linear system.
    On output, A is in unit upper triangular form, and
    represents the U factor of an LU factorization of the
    original coefficient matrix.

    Input, int NB, the number of right hand sides.

    Input, double B[N*NB], the right hand sides of the linear systems.

    Output, double R8MAT_FSS_NEW[N*NB], the solutions of the linear systems.
*/
{
  int i;
  int ipiv;
  int j;
  int jcol;
  double piv;
  double t;
  double *x;

  x = ( double * ) malloc ( n * nb * sizeof ( double ) );

  for ( j = 0; j < nb; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+j*n] = b[i+j*n];
    }
  }
  for ( jcol = 1; jcol <= n; jcol++ )
  {
/*
  Find the maximum element in column I.
*/
    piv = r8_abs ( a[jcol-1+(jcol-1)*n] );
    ipiv = jcol;
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( piv < r8_abs ( a[i-1+(jcol-1)*n] ) )
      {
        piv = r8_abs ( a[i-1+(jcol-1)*n] );
        ipiv = i;
      }
    }

    if ( piv == 0.0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "R8MAT_FSS_NEW - Fatal error!\n" );
      fprintf ( stderr, "  Zero pivot on step %d\n", jcol );
      exit ( 1 );
    }
/*
  Switch rows JCOL and IPIV, and X.
*/
    if ( jcol != ipiv )
    {
      for ( j = 1; j <= n; j++ )
      {
        t                 = a[jcol-1+(j-1)*n];
        a[jcol-1+(j-1)*n] = a[ipiv-1+(j-1)*n];
        a[ipiv-1+(j-1)*n] = t;
      }
      for ( j = 0; j < nb; j++ )
      {
        t            = x[jcol-1+j*n];
        x[jcol-1+j*n] = x[ipiv-1+j*n];
        x[ipiv-1+j*n] = t;
      }
    }
/*
  Scale the pivot row.
*/
    t = a[jcol-1+(jcol-1)*n];
    a[jcol-1+(jcol-1)*n] = 1.0;
    for ( j = jcol+1; j <= n; j++ )
    {
      a[jcol-1+(j-1)*n] = a[jcol-1+(j-1)*n] / t;
    }
    for ( j = 0; j < nb; j++ )
    {
      x[jcol-1+j*n] = x[jcol-1+j*n] / t;
    }
/*
  Use the pivot row to eliminate lower entries in that column.
*/
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( a[i-1+(jcol-1)*n] != 0.0 )
      {
        t = - a[i-1+(jcol-1)*n];
        a[i-1+(jcol-1)*n] = 0.0;
        for ( j = jcol+1; j <= n; j++ )
        {
          a[i-1+(j-1)*n] = a[i-1+(j-1)*n] + t * a[jcol-1+(j-1)*n];
        }
        for ( j = 0; j < nb; j++ )
        {
          x[i-1+j*n] = x[i-1+j*n] + t * x[jcol-1+j*n];
        }
      }
    }
  }
/*
  Back solve.
*/
  for ( jcol = n; 2 <= jcol; jcol-- )
  {
    for ( i = 1; i < jcol; i++ )
    {
      for ( j = 0; j < nb; j++ )
      {
        x[i-1+j*n] = x[i-1+j*n] - a[i-1+(jcol-1)*n] * x[jcol-1+j*n];
      }
    }
  }

  return x;
}
/******************************************************************************/

double *r8mat_identity_new ( int n )

/******************************************************************************/
/*
  Purpose:

    R8MAT_IDENTITY_NEW sets the square matrix A to the identity.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8's, which
    may be stored as a vector in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 September 2005

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of A.

    Output, double A[N*N], the N by N identity matrix.
*/
{
  double *a;
  int i;
  int j;
  int k;

  a = ( double * ) malloc ( n * n * sizeof ( double ) );

  k = 0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a[k] = 1.0;
      }
      else
      {
        a[k] = 0.0;
      }
      k = k + 1;
    }
  }

  return a;
}
/******************************************************************************/

void r8mat_minvm ( int n1, int n2, double a[], double b[], double c[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MINVM computes inverse(A) * B for R8MAT's.

  Discussion:

    An R8MAT is an array of R8 values.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N1, N2, the order of the matrices.

    Input, double A[N1*N1], B[N1*N2], the matrices.

    Output, double C[N1*N2], the result, C = inverse(A) * B.
*/
{
  double *alu;
  double *d;

  alu = r8mat_copy_new ( n1, n1, a );

  d = r8mat_fss_new ( n1, alu, n2, b );
 
  r8mat_copy ( n1, n2, d, c );

  free ( alu );
  free ( d );

  return;
}
/******************************************************************************/

void r8mat_mm ( int n1, int n2, int n3, double a[], double b[], double c[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MM multiplies two matrices.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    For this routine, the result is returned as the function value.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 December 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N1, N2, N3, the order of the matrices.

    Input, double A[N1*N2], double B[N2*N3], the matrices to multiply.

    Output, double C[N1*N3], the product matrix C = A * B.
*/
{
  double *d;
  int i;
  int j;
  int k;

  d = ( double * ) malloc ( n1 * n3 * sizeof ( double ) );

  for ( i = 0; i < n1; i ++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      d[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        d[i+j*n1] = d[i+j*n1] + a[i+k*n1] * b[k+j*n2];
      }
    }
  }

 for ( i = 0; i < n1; i ++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = d[i+j*n1];
    }
  }

  free ( d );

  return;
}
/******************************************************************************/

double r8mat_norm_l1 ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_NORM_L1 returns the matrix L1 norm of an R8MAT.

  Discussion:

    An R8MAT is an array of R8 values.

    The matrix L1 norm is defined as:

      R8MAT_NORM_L1 = max ( 1 <= J <= N )
        sum ( 1 <= I <= M ) abs ( A(I,J) ).

    The matrix L1 norm is derived from the vector L1 norm, and
    satisifies:

      r8vec_norm_l1 ( A * x ) <= r8mat_norm_l1 ( A ) * r8vec_norm_l1 ( x ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 December 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A(M,N), the matrix whose L1 norm is desired.

    Output, double R8MAT_NORM_L1, the L1 norm of A.
*/
{
  double col_sum;
  int i;
  int j;
  double value;

  value = 0.0;

  for ( j = 0; j < n; j++ )
  {
    col_sum = 0.0;
    for ( i = 0; i < m; i++ )
    {
      col_sum = col_sum + r8_abs ( a[i+j*m] );
    }
    value = r8_max ( value, col_sum );
  }
  return value;
}
/******************************************************************************/

double r8mat_norm_li ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_NORM_LI returns the matrix L-oo norm of an R8MAT.

  Discussion:

    An R8MAT is an array of R8 values.

    The matrix L-oo norm is defined as:

      R8MAT_NORM_LI =  max ( 1 <= I <= M ) sum ( 1 <= J <= N ) abs ( A(I,J) ).

    The matrix L-oo norm is derived from the vector L-oo norm,
    and satisifies:

      r8vec_norm_li ( A * x ) <= r8mat_norm_li ( A ) * r8vec_norm_li ( x ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 December 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the matrix whose L-oo
    norm is desired.

    Output, double R8MAT_NORM_LI, the L-oo norm of A.
*/
{
  int i;
  int j;
  double row_sum;
  double value;

  value = 0.0;

  for ( i = 0; i < m; i++ )
  {
    row_sum = 0.0;
    for ( j = 0; j < n; j++ )
    {
      row_sum = row_sum + r8_abs ( a[i+j*m] );
    }
    value = r8_max ( value, row_sum );
  }
  return value;
}
/******************************************************************************/

void r8mat_print ( int m, int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PRINT prints an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8's, which
    may be stored as a vector in column-major order.

    Entry A(I,J) is stored as A[I+J*M]

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the M by N matrix.

    Input, char *TITLE, a title.
*/
{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PRINT_SOME prints some of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8's, which
    may be stored as a vector in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, double A[M*N], the matrix.

    Input, int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column to be printed.

    Input, char *TITLE, a title.
*/
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );

  if ( m <= 0 || n <= 0 )
  {
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  (None)\n" );
    return;
  }
/*
  Print the columns of the matrix, in strips of 5.
*/
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    fprintf ( stdout, "\n" );
/*
  For each column J in the current range...

  Write the header.
*/
    fprintf ( stdout, "  Col:  ");
    for ( j = j2lo; j <= j2hi; j++ )
    {
      fprintf ( stdout, "  %7d     ", j - 1 );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Row\n" );
    fprintf ( stdout, "\n" );
/*
  Determine the range of the rows in this strip.
*/
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
/*
  Print out (up to) 5 entries in row I, that lie in the current strip.
*/
      fprintf ( stdout, "%5d:", i - 1 );
      for ( j = j2lo; j <= j2hi; j++ )
      {
        fprintf ( stdout, "  %14f", a[i-1+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
}
/******************************************************************************/

void r8mat_scale ( int m, int n, double s, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_SCALE multiplies an R8MAT by a scalar.

  Discussion:

    An R8MAT is an array of R8 values.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 December 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double S, the scale factor.

    Input/output, double A[M*N], the matrix to be scaled.
*/
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = a[i+j*m] * s;
    }
  }
  return;
}
/******************************************************************************/

int r8mat_significant ( int m, int n, double r[], double s[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_SIGNIFICANT determines if an R8MAT is significant compared to another.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the dimension of the matrices.

    Input, double R[M*N], the vector to be compared against.

    Input, double S[M*N], the vector to be compared.

    Output, int R8MAT_SIGNIFICANT, is TRUE if S is significant
    compared to R.
*/
{
  int i;
  int j;
  double t;
  double tol;
  int value;

  value = 0;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      t = r[i+j*m] + s[i+j*m];
      tol = r8_epsilon ( ) * r8_abs ( r[i+j*m] );

      if ( tol < r8_abs ( r[i+j*m] - t ) )
      {
        value = 1;
        break;
      }
    }
  }
  return value;
}
/******************************************************************************/

double *r8mat_zero_new ( int m, int n )

/******************************************************************************/
/*
  Purpose:

    R8MAT_ZERO_NEW returns a new zeroed R8MAT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 September 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Output, double R8MAT_ZERO[M*N], the new zeroed matrix.
*/
{
  double *a;
  int i;
  int j;

  a = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = 0.0;
    }
  }
  return a;
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
