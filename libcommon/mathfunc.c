#include<mathfunc.h>
#include<time.h>
#include<stdlib.h>
#include<stdio.h>

#define M_PI           3.14159265358979323846  /* pi */
#define M_SQRT2    1.41421356237309504880   // sqrt(2)

int __maxidx(double* data, int size )
{
	double max = data[0];
	int v;
	int i;
	v = 0;

	for ( i = 1 ; i < size; i++ )
	{
		if ( max < data[i] )
		{
			v  = i;
			max = data[i];
		}
	}
	return v;
}

double __max(double* data, int size )
{
	double max = data[0];
	int i;

	for ( i = 1 ; i < size; i++ )
	{
		if ( max < data[i] )
		{
			max = data[i];
		}
	}
	return max;
}

double __min(double* data, int size )
{
	double min = data[0];
	int i;
	
	for ( i = 1 ; i < size; i++ )
	{
		if ( min > data[i] )
		{
			min = data[i];
		}
	}
	return min;
}




double ipow(double val, int expo)
{
    double result = 1.0;
    unsigned int e = expo;
    
    if ((int)e < 0) {
        e = -e;
        val = 1.0 / val;
    }
    
    while (1) {
        if (e & 1)          //Bitwise AND
            result *= val;
        if ((e >>= 1) == 0) //Bitwise right shift assignment
            break;
        val *= val;
    }
    
    return result;
}

/*
 * URL : http://www.rskey.org/gamma.htm
 *       http://www.american.edu/academic.depts/cas/econ/gaussres/pdf/pdf.htm
 *       http://www.physics.unlv.edu/~pang/cp_c.html
 */

double __Qs[] = { 1.000000000190015, 76.18009172947146, -86.50532032941677,
		    	  24.01409824083091, -1.231739572450155, 1.208650973866179e-3, 
				  -5.395239384953e-6};

double unifrnd()
{
	return rand()/(RAND_MAX+1.0);
}

/***********************************************************************
 * Gamma Function
 ***********************************************************************/

double gamma(double a)
{
	return exp(gammaln(a));
}

#define M_SQRT_2PI		2.5066282746310002416123552393401042  // sqrt(2pi)
#define M_LOG_SQRT_2PI 	0.9189385332046726695409688545623794  // log(sqrt(2pi))

/*
 * URL: http://www.hpcsoft.com/products/MathSoL/specialFunction/gammaIn.html
 * \int^{\infty}_{0} t^{a-1}e^{-t}\;dt
 */
double gammaln(double a)
{
	int n;
	double p = __Qs[0];
	double a_add_5p5 = a + 5.5;
	for ( n = 1; n <= 6; n++ ) p += __Qs[n]/(a+n);
	return (a+0.5)*log(a_add_5p5) - (a_add_5p5) + log(M_SQRT_2PI*p/a);
}

/*
 * For x < \alpha + 1
 *  gamma(x,\alpha) = e^{-x}x^{\alpha} \sum_{i=0}^{\infty} \frac{\Gamma(\alpha)}{\Gamma(\alpha+1+i)} x^{i}
 */

double incgammaln_lower(double x, double a)
{
	int i;
	double p = 1/a;
	double t = 1/a;
	for ( i = 1 ; i < 1000 ; i++ )
	{
		t *= x/(a+i);
		if ( t < __EPS__ ) break;
		p += t;
	}
	return i == 1000 ? gammaln(a) :  log(p) + a * log(x) - x ;
}

/* \alpha > 0
 * P(x,\alpha) = \frac{1}{\Gamma{\alpha}\int_0^{x} x^{\alpha-1}e^{-t}\;dt
 */
double gammaincln(double x, double a)
{
	return incgammaln_lower(x,a) - gammaln(a);	
}

double gammainc(double x, double a)
{
	return exp(gammaincln(x,a));
}

/***********************************************************************
 * Beta Function
 ***********************************************************************/
double beta(double alpha, double beta)
{
	return betaln(alpha,beta);	
}

double betaln(double alpha, double beta)
{
	return gammaln(alpha) + gammaln(beta) - gammaln(alpha+beta);
}

#define BETA_CF_MIN 	1e-100

double __beta_cf( double a, double b, double x )
{
	int m;
	double a_add_b = a + b;
	double a_add_1 = a + 1;
	double a_sub_1 = a - 1;
	double c = 1;
	double d = 1 - a_add_b*x/a_add_1;
	if ( fabs(d) < BETA_CF_MIN ) d = BETA_CF_MIN;
	d = 1/d;
	double h = d;
	double m2, aa, del;
	for ( m = 1 ; m <= 200 ; m++ )
	{
		m2 = 2 * m;

		aa = m*(b-m)*x/(a_sub_1+m2)/(a+m2);	

		d = 1 + aa*d;
		if ( fabs(d) < BETA_CF_MIN ) d = BETA_CF_MIN;
		c = 1 + aa/c;
		if ( fabs(c) < BETA_CF_MIN ) c = BETA_CF_MIN;
		d = 1/d;
		h *= d*c;

		aa = -(a+m)*(a_add_b+m)*x/(a+m2)/(a_add_1+m2);

		d = 1 + aa*d;
		if ( fabs(d) < BETA_CF_MIN ) d = BETA_CF_MIN;
		c = 1 + aa/c;
		if ( fabs(c) < BETA_CF_MIN ) c = BETA_CF_MIN;
		d = 1/d;
		del = d*c;
		h *= del;
		if ( fabs(del-1) < __EPS__ ) break;
	}
	return h;
}

double betainc(double x, double alpha, double beta)
{
	if ( x == 0 || x == 1 ) return 0;
	double bconst = alpha * log(x) + beta * log(1-x) - betaln(alpha,beta);
	return ( x < (alpha+1)/(alpha+beta+2.0) ? exp(bconst)* __beta_cf(alpha,beta,x)/alpha : 
		                                   1 - exp(bconst)*  __beta_cf(beta,alpha,1-x)/beta);
}

/***********************************************************************
 * Factorial, permutation, combination 
 ***********************************************************************/

double choose(double n, double r)
{
	return exp(chooseln(n,r));
}

double chooseln(double n, double r)
{
	if ( r == 0 || ( n == 0 && r == 0) ) return 0;
	else if ( n <= 0 || r <= 0 ) return log(0);
	return gammaln(n+1) - gammaln(r+1) - gammaln(n-r+1);
}

double factorial(double n)
{
	return exp(gammaln(n+1));
}

double factorialln(double n)
{
	return gammaln(n+1);
}

/***********************************************************************
 * Gamma Distribution
 ***********************************************************************/

double gampdf(double x, double alpha, double beta)
/* 
 * alpha : shape parameter
 * beta : scale parameter
 */
{
	if ( x < 0 ) return 0;	
	if ( alpha < 0 || beta < 0 ) return 0;
	return exp( -(alpha*log(beta) + gammaln(alpha)) + (alpha-1)*log(x) - x/beta);
}

/*
 * \frac{ * \Gamma(\alpha,\frac{x}{\beta} } { \Gamma(\alpha) } 
 */

double gamcdf(double x, double alpha, double beta )
{
	return gammainc(x/beta,alpha);
}

/***********************************************************************
 * Beta Distribution
 ***********************************************************************/

double betapdf(double x, double alpha, double beta)
{
	if ( x < 0 || x > 1 ) return 0;
	return exp(-betaln(alpha,beta) + (alpha-1)*log(x) + (beta-1)*log(1-x)); 
}

/***********************************************************************
 * Chi Square Distribution
 ***********************************************************************/

double chi2pdf(double x, int df)
{
	return gampdf(x,df/2.0,2);
}

double chi2cdf(double x, int df )
{
	return gamcdf(x,df/2.0,2);
}

/***********************************************************************
 * F Distribution
 ***********************************************************************/

double fpdf(double x, int df1, int df2 )
{
	double v1 = df1/2.0;
	double v2 = df2/2.0;
	double p = gammaln(v1+v2) + v1 * log(v1/v2) + (v1-1)*log(x) - gammaln(v1) - gammaln(v2)
	           - ( v1 + v2 ) * log( 1+ x * v1/v2);
	return exp(p);
}

double fcdf(double x, int df1, int df2 )
{
	return 1 - betainc(df2/(df2+df1*x),df2/2.0, df1/2.0);
}

/***********************************************************************
 * t Distribution
 ***********************************************************************/

double tpdf(double x, int df)
{
	double v = (df + 1)/2.0;
	return exp(gammaln(v) - gammaln(v-0.5) - log(df*M_PI)/2 - v*log(1+x*x/df));
}

double tcdf(double x, int df)
{
	double tc = betainc(df/(df+x*x), df/2.0, 0.5 )/2;
	return x > 0 ? 1 - tc : tc;
}

/***********************************************************************
 * Binomial Distribution
 ***********************************************************************/

double binopdf(int x, int n, double p)
{
	return exp(chooseln(n,x) + x*log(p) * (n-x)*log(1-p));
}

double binocdf(int x, int n, double p)
{
	return betainc(1-p,n-x,x+1);
}

int binornd(int n, double p)
{
	int i;
	int cnt = 0;
	for ( i = 0 ; i < n ; i++ )
	{
		cnt += unifrnd() <= p ? 1 : 0;
	}
	return cnt;
}

/***********************************************************************
 * Poisson Distribution
 ***********************************************************************/

double poisspdf(int x, double lambda)
{
	return exp( x*log(lambda) - gammaln(x+1) - lambda );
}

double poisscdf(int x, double lambda )
{
	int i = 0;
	double cdf = 0;
	double LL = 1;
	double LX = 1;
	double tmp;
	for ( i = 1 ; i <= x ; i++ )
	{
		LL *= lambda;
		LX *= i;
		tmp = LL/LX;
		cdf += tmp;
		if ( tmp < __EPS__ ) break;
	}
	return (1 + cdf )*exp(-lambda);
}

// http://www.mindspring.com/~hamill4/hamillnumerics/idl/random/knuthpoissonq.html

int poissrnd(double lambda)
{
	int n = 0;	
	double q = 1;
	double p = exp(-lambda);
	q = q * unifrnd();
	while( q >= p )
	{
		n++;
		q = q * unifrnd();
	}
	return n;
}

/***********************************************************************
 * Exponential Distribution
 ***********************************************************************/

double exppdf(double x, double mu )
{
	return exp(-x/mu)/mu;
}

double exprnd(double mu)
{
	return -log(unifrnd())*mu;
}

/***********************************************************************
 * Rayleight Distribution
 ***********************************************************************/

double raylpdf(double x, double b)
{
	return x/(b*b) * exp( -x*x/(b*b));
}

/***********************************************************************
 * Gaussian Distribution
 ***********************************************************************/

#define __SQRT_2_PI	2.506628274631000241612355

double normpdf(double x, double mu, double sigma)
{
	double r = (x-mu)*(x-mu)/(2*sigma*sigma);
	return exp(-r)/(sigma*__SQRT_2_PI);
}

double normrnd(double mu, double sigma)
{
/*
	int i = 0;
	double sum = 0;
	for ( i = 0 ; i < 4 ; i++ )
	{
		sum += unifrnd();
	}
	return  ( sum - 2 ) * sigma/0.57735026918963 + mu;
*/
	double x1,x2,w;
	do 
	{
		x1 = 2 * unifrnd() - 1;
		x2 = 2 * unifrnd() - 1;
		w = x1 * x1 + x2 * x2;
	}while( w >= 1.0 );
	w = sqrt( (-2* log(w))/w );
	return mu + x1*w * sigma;
}


/* 
Source http://home.online.no/~pjacklam/notes/invnorm/impl/lea/lea.c
*/

/*
 ** An implementation of adaptive, recursive Newton-Cotes integration.
 ** Based on the MATLAB implementation, but covered in a lot of books...
 **
 ** This only does integration over the standard normal PDF.  It's just
 ** here to check the error function approximations.
 **/
#define LEVMAX 10
double quad8_stdnormal_pdf(double a, double b, double Q )
{
 /* The magic Newton-Cotes weights */
 const int w[9] = {3956, 23552, -3712, 41984, -18160, 41984, -3712, 23552,
3956};
 const int dw = 14175;
 static int level = -1;
 static double tol = 1e-30;
 register double h, Q1 = 0.0, Q2 = 0.0;
 register int i;

 level++;
 h = (b-a)/16.0;
 for (i = 0; i < 9; i++) {
  Q1 += h*w[i]*normpdf(a+i*h,0,1)/dw;
  Q2 += h*w[i]*normpdf(a+(i+8)*h,0,1)/dw;
 };
 /* This is the adaptive recursive bit.  We only recurse if we can
 * improve... */
 if (fabs(Q1+Q2-Q) > tol*fabs(Q1+Q2) && level <= LEVMAX) {
  tol = tol/2;
  Q1 = quad8_stdnormal_pdf(a,(a+b)/2,Q1);
  Q2 = quad8_stdnormal_pdf((a+b)/2,b,Q2);
  tol = tol*2;
 }
 level--;
 return Q1 + Q2;
}

/*
 ** The standard normal CDF, for one random variable.
 **
 **   Author:  W. J. Cody
 **   URL:   http://www.netlib.org/specfun/erf
 **
 ** This is the erfc() routine only, adapted by the
 ** transform stdnormal_cdf(u)=(erfc(-u/sqrt(2))/2;
 **/
#define  M_1_SQRTPI		0.564189583547756279280349644978  // 1/sqrt(pi);
#define  M_SQRT2PI		2.506628274631000241612355239340  // sqrt(2*pi);

double stdnormal_cdf(double u)
{
 const double a[5] = {
  1.161110663653770e-002,3.951404679838207e-001,2.846603853776254e+001,
  1.887426188426510e+002,3.209377589138469e+003
 };
 const double b[5] = {
  1.767766952966369e-001,8.344316438579620e+000,1.725514762600375e+002,
  1.813893686502485e+003,8.044716608901563e+003
 };
 const double c[9] = {
  2.15311535474403846e-8,5.64188496988670089e-1,8.88314979438837594e00,
  6.61191906371416295e01,2.98635138197400131e02,8.81952221241769090e02,
  1.71204761263407058e03,2.05107837782607147e03,1.23033935479799725E03
 };
 const double d[9] = {
  1.00000000000000000e00,1.57449261107098347e01,1.17693950891312499e02,
  5.37181101862009858e02,1.62138957456669019e03,3.29079923573345963e03,
  4.36261909014324716e03,3.43936767414372164e03,1.23033935480374942e03
 };
 const double p[6] = {
  1.63153871373020978e-2,3.05326634961232344e-1,3.60344899949804439e-1,
  1.25781726111229246e-1,1.60837851487422766e-2,6.58749161529837803e-4
 };
 const double q[6] = {
  1.00000000000000000e00,2.56852019228982242e00,1.87295284992346047e00,
  5.27905102951428412e-1,6.05183413124413191e-2,2.33520497626869185e-3
 };
 register double y, z;

 /*
 if (_isnan(u))
  return _Nan._D;
 if (!_finite(u))
  return (u < 0 ? 0.0 : 1.0);
 */
 y = fabs(u);
    if (y <= 0.46875*M_SQRT2) {
  /* evaluate erf() for |u| <= sqrt(2)*0.46875 */
  z = y*y;
  y = u*((((a[0]*z+a[1])*z+a[2])*z+a[3])*z+a[4])
       /((((b[0]*z+b[1])*z+b[2])*z+b[3])*z+b[4]);
  return 0.5+y;
 }
 z = exp(-y*y/2)/2;
 if (y <= 4.0) {
  /* evaluate erfc() for sqrt(2)*0.46875 <= |u| <= sqrt(2)*4.0 */
  y = y/M_SQRT2;
  y =
((((((((c[0]*y+c[1])*y+c[2])*y+c[3])*y+c[4])*y+c[5])*y+c[6])*y+c[7])*y+c[8])


/((((((((d[0]*y+d[1])*y+d[2])*y+d[3])*y+d[4])*y+d[5])*y+d[6])*y+d[7])*y+d[8]);

  y = z*y;
    } else {
  /* evaluate erfc() for |u| > sqrt(2)*4.0 */
  z = z*M_SQRT2/y;
  y = 2/(y*y);
        y = y*(((((p[0]*y+p[1])*y+p[2])*y+p[3])*y+p[4])*y+p[5])
    /(((((q[0]*y+q[1])*y+q[2])*y+q[3])*y+q[4])*y+q[5]);
        y = z*(M_1_SQRTPI-y);
    }
 return (u < 0.0 ? y : 1-y);
};

/*
 *  * The inverse standard normal distribution.
 *   *
 *    *   Author:      Peter J. Acklam <pjacklam@online.no>
 *     *   URL:         http://home.online.no/~pjacklam
 *      *
 *       * This function is based on the MATLAB code from the address above,
 *        * translated to C, and adapted for our purposes.
 *         */
double stdnormal_inv(double p)
{
 const double a[6] = {
  -3.969683028665376e+01,  2.209460984245205e+02,
  -2.759285104469687e+02,  1.383577518672690e+02,
  -3.066479806614716e+01,  2.506628277459239e+00
 };
 const double b[5] = {
  -5.447609879822406e+01,  1.615858368580409e+02,
  -1.556989798598866e+02,  6.680131188771972e+01,
  -1.328068155288572e+01
 };
 const double c[6] = {
  -7.784894002430293e-03, -3.223964580411365e-01,
  -2.400758277161838e+00, -2.549732539343734e+00,
   4.374664141464968e+00,  2.938163982698783e+00
 };
 const double d[4] = {
   7.784695709041462e-03,  3.224671290700398e-01,
   2.445134137142996e+00,  3.754408661907416e+00
 };

 register double q, t, u;
/*
 if (_isnan(p) || p > 1.0 || p < 0.0)
  return _Nan._D;
 if (p == 0.0)
  return -_Inf._D;
 if (p == 1.0)
  return  _Inf._D;
*/
 q = MIN(p,1-p);
 if (q > 0.02425) {
  /* Rational approximation for central region. */
  u = q-0.5;
  t = u*u;
  u = u*(((((a[0]*t+a[1])*t+a[2])*t+a[3])*t+a[4])*t+a[5])
    /(((((b[0]*t+b[1])*t+b[2])*t+b[3])*t+b[4])*t+1);
 } else {
  /* Rational approximation for tail region. */
  t = sqrt(-2*log(q));
  u = (((((c[0]*t+c[1])*t+c[2])*t+c[3])*t+c[4])*t+c[5])
   /((((d[0]*t+d[1])*t+d[2])*t+d[3])*t+1);
 }
 /* The relative error of the approximation has absolute value less
 *     than 1.15e-9.  One iteration of Halley's rational method (third
 *         order) gives full machine precision... */
 t = stdnormal_cdf(u)-q;    /* error */
 t = t*M_SQRT2PI*exp(u*u/2);   /* f(u)/df(u) */
 u = u-t/(1+u*t/2);     /* Halley's method */

 return (p > 0.5 ? -u : u);
};


double normcdf(double x, double mu, double sigma)
{
	return stdnormal_cdf((x-mu)/sigma);
}

double norminv(double x, double mu, double sigma)
{
	return stdnormal_inv((x-mu)/sigma);
}

/***********************************************************************
 * ETC
 ***********************************************************************/

double mean(double* data, int size)
{
	int i;
	double sum = 0;
	for( i = 0 ; i < size; i++ ) sum += data[i];
	return sum/size;
}

double variance(double* data, int size)
{
	int i;
	double m = mean(data,size);
	double v = 0;
	for ( i = 0 ; i < size; i++ ) v += ( data[i] - m ) * (data[i] - m);
	return v/(size-1);
}

double pvalue(double v, double* conddist, int size )
{
	int from = 0;
	int to = size - 1;
	int mi;
	while(from < to)
	{
		mi = from + (to-from)/2;
		if ( conddist[mi] > v )
		{
			to = mi - 1;
		}
		else if ( conddist[mi] < v )
		{
			from = mi + 1;
		}
		else
		{
			for( from = mi - 1; from >= 0 && conddist[from] == v ; from-- );
			for( to = mi + 1; to < size && conddist[to] == v ; to++ );
			from++, to--;
			break;
		}
	}
	if ( from > to ) to = from;
	return (double)( from + (conddist[from] <= v ? 1 : 0) + (to-from)/2.0)/(double)size;
}

double summation(double* data, int size)
{
	int i;
	double s = 0;
	for ( i = 0 ; i < size ; i++ )
	{
		s += data[i];	
	}
	return s;
}

double* vector_fraction(double* data, int size, double denominator)
{
	int i;
	for ( i = 0 ; i < size ; i++ )
	{
		data[i] = data[i]/denominator;	
	}
	return data;
}

int comp_double (const void * elem1, const void * elem2) 
{
    double f = *((double*)elem1);
    double s = *((double*)elem2);
    return (f > s) - (f < s);
    return 0;
}
