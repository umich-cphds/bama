#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include <R_ext/BLAS.h>

#define GAMMA_SHAPE    2.0
#define GAMMA_SCALE_M0 1E-4
#define GAMMA_SCALE_M1 1.0
#define GAMMA_SCALE    1.0

#define BETA_M    0
#define R1        1
#define ALPHA_A   2
#define R3        3
#define BETA_A    4
#define PI_M      4
#define PI_A      6
#define SIGMA_M0  7
#define SIGMA_M1  8
#define SIGMA_MA0 9
#define SIGMA_MA1 10

struct dataframe {
    double *y;
    double *a;
    double *m;
    double *c1;
    double *c2;
    int n;
    int nm;
    int nc1;
    int nc2;
};

struct norm {
  double a;
  double *m;
  double *c1;
  double *c2;
};

struct sigma {
  double m0;
  double m1;
  double ma0;
  double ma1;
  double a;
  double e;
  double g;
};

struct residuals {
    double *y;
    double *m;
    double *mc;
};

struct coeffs {
    double beta_a;
    double *beta_m;
    double *beta_c1;
    double *alpha_a;
    double *alpha_c2;
};

/*
 * wapper functions around R functions that generate random numbers.
 * Get/PutRNGstate handle updating the seed for R's RNG
 */
int rbernoulli(double p)
{
    GetRNGstate();
    double random = runif(0.0, 1.0);
    PutRNGstate();
    return random <= p ? 1 : 0;
}

double rinvgamma(double alpha, double beta)
{
    GetRNGstate();
    double random = rgamma(alpha, beta);
    PutRNGstate();
    return 1.0 / random;
}

double rnormal(double mu, double var)
{
    GetRNGstate();
    double random = rnorm(mu, var);
    PutRNGstate();
    return random;
}

double runiform(double a, double b)
{
    GetRNGstate();
    double random = runif(a, b);
    PutRNGstate();
    return random;
}

void calculate_residuals(struct residuals r, struct coeffs *coeffs,
                             struct dataframe *df)
{
    int n = df->n;
    int q = df->nm;
    int w = df->nc2;
    /*
     * calculate residuals for the the outcome model. Start off by
     * residuals.y := y
     */
    memcpy(r.y, df->y, n * sizeof(double));
    /* calcluate residuals.y := residuals.y - beta_a * A  */
    int u = 1;
    double a = -1.0 * coeffs->beta_a;
    double b = 1.0;
    daxpy_(&n, &a, df->a, &u, r.y, &u);

    a = -1.0;
    /* residuals.y := - M * beta_m + residuals.y */
    dgemv_("N", &n, &q, &a, df->m, &n, coeffs->beta_m, &u, &b, r.y, &u);

    /* residuals.y := - C1 * beta_c1 + residuals.y */
    dgemv_("N", &n, &df->nc1, &a, df->c1, &n, coeffs->beta_c1, &u, &b, r.y, &u);

    /* calculate residuals for the the mediator model */
    memcpy(r.m, df->m, n * q * sizeof(double));

    /* residuals.m := - alpha_c * c + residuals.m */
    dgemm_("N", "N", &n, &q, &w, &a, df->c2, &n, coeffs->alpha_c2, &w, &b, r.m, &n);

    memcpy(r.mc, r.m, n * q * sizeof(double));

    /*  residuals.m := - A * alpha_a + residuals.m */
    dger_(&n, &q, &a, df->a, &u, coeffs->alpha_a, &u, r.m, &n);
}

void update_beta_m(double *beta_m, double *r1, struct residuals residuals,
                       double pi_m,  double *var_m0, double *var_m1,
                       struct dataframe *df, struct sigma *sigma,
                       double *norm2_m)
{
    int n = df->n;
    int q = df->nm;

    double *mu_beta_m0 = malloc(q * sizeof(double));
    double *mu_beta_m1 = malloc(q * sizeof(double));

    int u = 1;
    double a = 1.0;
    double b = 1.0;
    double z = 0.0;

    /* calculate mu_beta_m = M**T * (r.y + M * beta_m) */
    /* set tmp := residuals.y */

    /* residuals.y := residuals.y + M * beta_M */
    dgemv_("N", &n, &q, &a, df->m, &n, beta_m, &u, &b, residuals.y, &u);
    /* mu_m := M**T * tmp */
    dgemv_("T", &n, &q, &a, df->m, &n, residuals.y, &u, &z, mu_beta_m0, &u);

    memcpy(mu_beta_m1, mu_beta_m0, q * sizeof(double));
    /* generate new beta_m */
    for (int j = 0; j < q; ++j) {
        mu_beta_m0[j] /= (sigma->e / sigma->m0 + norm2_m[j]);
        mu_beta_m1[j] /= (sigma->e / sigma->m1 + norm2_m[j]);
        beta_m[j] = r1[j] * rnormal(mu_beta_m1[j], sqrt(var_m1[j])) +
                        (1.0 - r1[j]) * rnormal(mu_beta_m0[j], sqrt(var_m0[j]));
    }
    /* residuals.y := - M * beta_m + residuals.y */
    a = -1.0;
    dgemv_("N", &n, &q, &a, df->m, &n, beta_m, &u, &b, residuals.y, &u);

    /* update r1 */
    for (int j = 0; j < q; ++j) {
        double p = 0.5 * (mu_beta_m1[j] * mu_beta_m1[j] / var_m1[j] -
                        mu_beta_m0[j] * mu_beta_m0[j] / var_m0[j] +
                        log(sigma->m0 * var_m1[j] / (sigma->m1 * var_m0[j]))) +
                        log(pi_m / (1.0 - pi_m));
        if (p > 100.0)
            r1[j] = 1.0;
        else
            r1[j] = rbernoulli(exp(p) / (1.0 + exp(p)));
    }

    free(mu_beta_m0);
    free(mu_beta_m1);
}

void update_alpha_a(double *alpha_a, double *r3, struct residuals residuals,
                          double pi_a, struct dataframe *df, struct sigma
                          *sigma, double norm2_a)
{
    int q = df->nm;
    int n = df->n;

    double var_alpha_a0 = sigma->g / (sigma->g / sigma->ma0 + norm2_a);
    double var_alpha_a1 = sigma->g / (sigma->g / sigma->ma1 + norm2_a);

    double *mu_alpha_a0 = malloc(q * sizeof(double));
    double *mu_alpha_a1 = malloc(q * sizeof(double));

    int u = 1;
    double a = 1.0;
    double z = 0.0;

    /* Calculate mu_alpha_a := A**T * residuals.c2 */
    dgemv_("T", &n, &q, &a, residuals.mc, &n, df->a, &u, &z, mu_alpha_a0, &u);
    memcpy(mu_alpha_a1, mu_alpha_a0, q * sizeof(double));
    double *new_alpha_a = malloc(q * sizeof(double));
    for (int j = 0; j < q; ++j) {
        mu_alpha_a0[j] *= var_alpha_a0 / sigma->g;
        mu_alpha_a1[j] *= var_alpha_a1 / sigma->g;
        new_alpha_a[j] = r3[j] * rnormal(mu_alpha_a1[j], sqrt(var_alpha_a1)) +
            (1.0 - r3[j]) * rnormal(mu_alpha_a0[j], sqrt(var_alpha_a0));
    }
    for (int j = 0; j < df->nm; ++j)
        alpha_a[j] -= new_alpha_a[j];
    /* residuals.m := +  a * (old_alpha_a - new_alpha_a)**T + residuals.m */
    dger_(&n, &q, &a, df->a, &u, alpha_a, &u, residuals.m, &n);

    memcpy(alpha_a, new_alpha_a, q * sizeof(double));
    free(new_alpha_a);

    /* update r3 */
    for (int j = 0; j < q; ++j) {
        double p = 0.5 * (mu_alpha_a1[j] * mu_alpha_a1[j] / var_alpha_a1 -
            mu_alpha_a0[j] * mu_alpha_a0[j] / var_alpha_a0 +
            log(var_alpha_a1 * sigma->ma0 / (var_alpha_a0 * sigma->ma1))) +
            log(pi_a / (1.0 - pi_a));
        if (p > 100)
            r3[j]= 1.0;
        else
            r3[j] = rbernoulli(exp(p) / (1.0 + exp(p)));
    }

    free(mu_alpha_a0);
    free(mu_alpha_a1);
}

void update_alpha_c(double *alpha_c2, struct residuals r, struct dataframe
                        *df, double sigma_g, double *norm2_c)
{
    int n = df->n;
    int q = df->nm;
    int w = df->nc2;

    double a = 1.0;
    double b = 1.0;
    double z = 0.0;

    double *mu_alpha_c = malloc(q * w * sizeof(double));
    double *tmp        = malloc(n * q * sizeof(double));
    memcpy(tmp, r.m, n * q * sizeof(double));
    /* tmp := r.m + C2 * alpha_c */
    dgemm_("N", "N", &n, &q, &w, &a, df->c2, &n, alpha_c2, &w, &b, tmp, &n);
    /* mu_alpha_c = C2**T * tmp */
    dgemm_("T", "N", &w, &q, &n, &a, df->c2, &n, tmp, &n, &z, mu_alpha_c, &w);
    free(tmp);

    double *new_alpha_c = malloc(w * q * sizeof(double));
    for (int j = 0; j < q; ++j) {
        for (int k = 0; k < w; ++k) {
            new_alpha_c[w * j + k] = rnormal(mu_alpha_c[w * j + k] / norm2_c[k],
                                                sqrt(sigma_g / norm2_c[k]));
            alpha_c2[w * j + k] -= new_alpha_c[w * j + k];
        }
    }
    free(mu_alpha_c);

    /* r.m := C2 * (old_alpha_c - new_alpha_c) + r.m  */
    dgemm_("N", "N", &n, &q, &w, &a, df->c2, &n, alpha_c2, &w, &b, r.m, &n);

    /* r.mc := C2 * (old_alpha_c - new_alpha_c) + r.mc */
    dgemm_("N", "N", &n, &q, &w, &a, df->c2, &n, alpha_c2, &w, &b, r.mc, &n);

    memcpy(alpha_c2, new_alpha_c, w * q * sizeof(double));
    free(new_alpha_c);
}

void update_beta_c(double *beta_c1, struct residuals r, struct dataframe *df,
                       double sigma_e, double *norm2_c1)
{
    int n = df->n;
    int t = df->nc1;

    int u = 1;
    double a = 1.0;
    double b = 1.0;
    double z = 0.0;

    /* mu_beta_c := C1**T * (r.y + C1 * beta_c1) */
    double *mu_beta_c1 = malloc(t * sizeof(double));
    /* r.y := r.y + C1 * beta_c1 */
    dgemv_("N", &n, &t, &a, df->c1, &n, beta_c1, &u, &b, r.y, &u);
    /* mu_beta_c1 := C1**T * tmp */
    dgemv_("T", &n, &t, &a, df->c1, &n, r.y, &u, &z, mu_beta_c1, &u);

    for (int j = 0; j < t; ++j) {
          beta_c1[j] = rnormal(mu_beta_c1[j] / norm2_c1[j],
                               sqrt(sigma_e / norm2_c1[j]));
    }
    free(mu_beta_c1);
    a = -1.0;
    /* r.y  := - C1 * beta_c1 + r.y */
    dgemv_("N", &n, &t, &a, df->c1, &n, beta_c1, &u, &b, r.y, &u);
}

double update_beta_a(double beta_a, double *a, int n, double var_a, double *ry,
                       double sigma_e)
{
    Rprintf("var a: %f\n", var_a);
    double mu_a = 0.0;
    for (int i = 0; i < n; ++i)
        mu_a += a[i] * (ry[i] + beta_a * a[i]);
    double new_beta_a = rnormal(mu_a * var_a / sigma_e, sqrt(var_a));
    beta_a -= new_beta_a;
    for (int i = 0; i < n; ++i)
        ry[i] += beta_a * a[i];
    Rprintf("NEW BETA_A: %f\n", new_beta_a);
    return new_beta_a;
}

void update_sigma(struct sigma *sigma, struct coeffs *coeffs, double *r1,
                      double *r3, double q)
{
    double sum_r1 = 0.0;
    double norm2_beta_m_r1 = 0.0;
    double norm2_beta_m    = 0.0;
    double sum_r3 = 0.0;
    double norm2_alpha_a  = 0.0;
    double norm2_alpha_a_r3 = 0.0;
    for (int j = 0; j < q; ++j) {
        sum_r1 += r1[j];
        norm2_beta_m += coeffs->beta_m[j] * coeffs->beta_m[j];
        norm2_beta_m_r1 += coeffs->beta_m[j] * coeffs->beta_m[j] * r1[j];
        sum_r3 += r3[j];
        norm2_alpha_a_r3 = coeffs->alpha_a[j] * coeffs->alpha_a[j] * r3[j];
        norm2_alpha_a = coeffs->alpha_a[j] * coeffs->alpha_a[j];
    }
    sum_r1           /= 2.0;
    norm2_beta_m     /= 2.0;
    norm2_beta_m_r1  /= 2.0;
    sum_r3           /= 2.0;
    norm2_alpha_a_r3 /= 2.0;
    norm2_alpha_a    /= 2.0;
    q                /= 2.0;
    sigma->a = rinvgamma(0.5 + GAMMA_SHAPE, 1.0 / (0.5 * coeffs->beta_a *
                                              coeffs->beta_a + GAMMA_SCALE));

    sigma->m1 = rinvgamma(sum_r1 + GAMMA_SHAPE,
                    1.0 / (norm2_beta_m_r1 + GAMMA_SCALE_M1));

    sigma->ma1 = rinvgamma(sum_r3 + GAMMA_SHAPE, 1.0 / (norm2_alpha_a_r3 +
                                                 GAMMA_SCALE_M1));

    sigma->m0 = rinvgamma(q - sum_r1 + GAMMA_SHAPE, 1.0 / (norm2_beta_m -
                                                    norm2_beta_m_r1 +
                                                    GAMMA_SCALE_M0));
                                                    
    sigma->ma0 = rinvgamma(q - sum_r3 + GAMMA_SHAPE, 1.0 / (norm2_alpha_a -
                                                            norm2_alpha_a_r3 +
                                                            GAMMA_SCALE_M0));
}

double posterior_distribution(int q, double *r_, double m_pi)
{
    double sum = 0.0;
    for (int j = 0; j < q; ++j)
        sum += log(m_pi) * r_[j] + log(1 - m_pi) * (1 - r_[j]);
    return sum;
}

void iterate(struct dataframe *df, struct sigma *sigma, struct norm norm2,
                 struct residuals residuals, struct coeffs *coeffs,
                 double *r1, double *r3, double *pi_m, double *pi_a)
{
    int n = df->n;
    int q = df->nm;

    double *var_m0 = malloc(q * sizeof(double));
    double *var_m1 = malloc(q * sizeof(double));
    for (int j = 0; j < q; ++j) {
        var_m0[j] = 1.0 / (1.0 / sigma->m0 + norm2.m[j] / sigma->e);
        var_m1[j] = 1.0 / (1.0 / sigma->m1 + norm2.m[j] / sigma->e);
    }

    double le1 = 0.0;
    double lg1 = 0.0;
    for (int i = 0; i < n; ++i)
        le1 += residuals.y[i] * residuals.y[i];
    Rprintf("sum r.y^2 %f\n", le1);
    for (int i = 0; i < n * q; ++i)
        lg1 += residuals.m[i] * residuals.m[i];
    Rprintf("sum r.m^2 %f\n", lg1);
    le1 = 1.0 / (0.5 * le1 + GAMMA_SCALE);
    lg1 = 1.0 / (0.5 * lg1 + GAMMA_SCALE);

    sigma->e = rinvgamma(GAMMA_SHAPE + 0.5 * n, le1);
    sigma->g = rinvgamma(GAMMA_SHAPE + 0.5 * n * q, lg1);

    Rprintf("beta a %f\n", coeffs->beta_a);
    Rprintf("sigma a %f\n", sigma->a);
    Rprintf("sigma e %f\n", sigma->e);
    Rprintf("sigma g %f\n", sigma->g);
    Rprintf("sigma m0 %f\n", sigma->m0);
    Rprintf("sigma ma0 %f\n", sigma->ma0);
    Rprintf("sigma m1 %f\n", sigma->m1);
    Rprintf("sigma ma1 %f\n", sigma->ma1);
    update_beta_m(coeffs->beta_m, r1, residuals, *pi_m, var_m0, var_m1, df,
                       sigma, norm2.m);
    free(var_m0);
    free(var_m1);
    double var_a = sigma->e / (sigma->e / sigma->a + norm2.a);
    update_alpha_a(coeffs->alpha_a, r3, residuals, *pi_a, df, sigma, norm2.a);

    update_alpha_c(coeffs->alpha_c2, residuals, df, sigma->g, norm2.c2);

    update_beta_c(coeffs->beta_c1, residuals, df, sigma->e, norm2.c1);

    // coeffs->beta_a = update_beta_a(coeffs->beta_a, df->a, n, var_a, residuals.y, sigma->e);

    update_sigma(sigma, coeffs, r1, r3, q);

    double m_pi_m = fabs(*pi_m * exp(runiform(-0.01, 0.01)));
    double m_pi_a = fabs(*pi_a * exp(runiform(-0.01, 0.01)));
    if (m_pi_m  > 1.0)
        m_pi_m = 1.0 / m_pi_m;
    if (m_pi_m  < 1.0 / q)
        m_pi_m = 1.0 / (q * q * m_pi_m);
    if (m_pi_a  > 1.0)
        m_pi_a = 1.0 / m_pi_a;
    if (m_pi_a  < 1.0 / q)
        m_pi_a = 1.0 / (q * q * m_pi_a);

    double p = posterior_distribution(q, r3, m_pi_a) -
                   posterior_distribution(q, r3, *pi_a) +
                   posterior_distribution(q, r1, m_pi_m) -
                   posterior_distribution(q, r1, *pi_m);

    if (p > log(runiform(0.0, 1.0))) {
        *pi_a = m_pi_a;
        *pi_m = m_pi_m;
     }
}

/* generate initial values for the sigmas, which are assumed to have an inverse
 * gamma posterior distribution */
struct sigma initialize_sigma()
{
    struct sigma sigma;
    sigma.m0  = rinvgamma(GAMMA_SHAPE, 1.0 / GAMMA_SCALE_M0);
    sigma.m1  = rinvgamma(GAMMA_SHAPE, 1.0 / GAMMA_SCALE_M1);
    sigma.ma0 = rinvgamma(GAMMA_SHAPE, 1.0 / GAMMA_SCALE_M0);
    sigma.ma1 = rinvgamma(GAMMA_SHAPE, 1.0 / GAMMA_SCALE_M1);
    sigma.a   = rinvgamma(GAMMA_SHAPE, 1.0 / GAMMA_SCALE);
    sigma.e   = rinvgamma(GAMMA_SHAPE, 1.0 / GAMMA_SCALE);
    sigma.g   = rinvgamma(GAMMA_SHAPE, 1.0 / GAMMA_SCALE);
    return sigma;
}

/* calcluate the L2 norm for each column of data. */
struct norm calcluate_norms(struct dataframe *df)
{
    int n = df->n;
    int u = 1;

    struct norm norm2;
    norm2.a = ddot_(&n, df->a, &u, df->a, &u);
    norm2.m = malloc(df->nm * sizeof(double));
    for (int j = 0; j < df->nm; ++j)
        norm2.m[j] = ddot_(&n, df->m + n * j, &u, df->m + n * j, &u);

    norm2.c1 = malloc(df->nc1 * sizeof(double));
    for (int j = 0; j < df->nc1; ++j)
        norm2.c1[j] = ddot_(&n, df->c1 + n * j, &u, df->c1 + n * j, &u);

    norm2.c2 = malloc(df->nc2 * sizeof(double));
    for (int j = 0; j < df->nc2; ++j)
        norm2.c2[j] = ddot_(&n, df->c2 + n * j, &u, df->c2 + n * j, &u);

    return norm2;
}

void hdbm_mcmc(struct dataframe *df, struct coeffs *coeffs, double pi_m,
                   double pi_a, int burnin, int ns, SEXP output)
{
    double *r1 = calloc(df->nm, sizeof(double));
    double *r3 = calloc(df->nm, sizeof(double));
    struct norm norm2  = calcluate_norms(df);
    struct sigma sigma = initialize_sigma();
    coeffs->beta_a = runiform(0, sqrt(sigma.a));
    struct residuals r;
    r.y  = malloc(df->n * sizeof(double));
    r.m  = malloc(df->n * df->nm * sizeof(double));
    r.mc = malloc(df->n * df->nm * sizeof(double));
    calculate_residuals(r, coeffs, df);
    /* run mcmc for the number of specified burn in iterations */
    for (int i = 0; i < burnin; ++i)
        iterate(df, &sigma, norm2, r, coeffs, r1, r3, &pi_m, &pi_a);
    /* sample from the markov chain. */
    for (int i = 0; i < ns; ++i) {
        /* run 50 iterations between samples to reduce autocorrelation(?) */
        for (int j = 0; j < 50; ++j)
            iterate(df, &sigma, norm2, r, coeffs, r1, r3, &pi_m, &pi_a);

        /* fill in row i of output with the current values */
        for (int j = 0; j < df->nm; ++j) {
            REAL(VECTOR_ELT(output, BETA_M))[i + ns * j]  = coeffs->beta_m[j];
            REAL(VECTOR_ELT(output, R1))[i + ns * j]      = r1[j];
            REAL(VECTOR_ELT(output, ALPHA_A))[i + ns * j] = coeffs->alpha_a[j];
            REAL(VECTOR_ELT(output, R3))[i + ns * j]      = r3[j];
        }
        REAL(VECTOR_ELT(output, BETA_A))[i]    = coeffs->beta_a;
        REAL(VECTOR_ELT(output, PI_M))[i]      = pi_m;
        REAL(VECTOR_ELT(output, PI_A))[i]      = pi_a;
        REAL(VECTOR_ELT(output, SIGMA_M0))[i]  = sigma.m0;
        REAL(VECTOR_ELT(output, SIGMA_M1))[i]  = sigma.m1;
        REAL(VECTOR_ELT(output, SIGMA_MA0))[i] = sigma.ma0;
        REAL(VECTOR_ELT(output, SIGMA_MA1))[i] = sigma.ma1;
    }

    free(norm2.m);
    free(norm2.c1);
    free(norm2.c2);
    free(r.y);
    free(r.m);
    free(r.mc);
    free(r1);
    free(r3);
}

SEXP hdbm_rwrapper(SEXP y, SEXP a, SEXP m, SEXP c1, SEXP c2, SEXP alpha_a, SEXP
                       beta_m, SEXP pi_m, SEXP pi_a, SEXP burnin, SEXP nsamples)
{
    struct dataframe df;
    df.a   = REAL(a);
    df.y   = REAL(y);
    df.m   = REAL(m);
    df.c1  = REAL(c1);
    df.c2  = REAL(c2);
    df.n   = length(y);
    df.nm  = ncols(m);
    df.nc1 = ncols(c1);
    df.nc2 = ncols(c2);
    struct coeffs coeffs;

    coeffs.beta_m =  malloc(df.nm * sizeof(double));
    memcpy(coeffs.beta_m, REAL(beta_m), df.nm * sizeof(double));
    coeffs.alpha_a =  malloc(df.nm * sizeof(double));
    memcpy(coeffs.alpha_a, REAL(alpha_a), df.nm * sizeof(double));
    coeffs.beta_c1  = calloc(df.nm, sizeof(double));
    coeffs.alpha_c2 = calloc(df.nc2 * df.nm, sizeof(double));
    /* bogus value. set in hdbm_mcmc */
    coeffs.beta_a = 31426;

    int ns = asInteger(nsamples);
    /* allocate memory for the R output */
    SEXP output = PROTECT(allocVector(VECSXP, 11));
    SET_VECTOR_ELT(output, BETA_M,    allocMatrix(REALSXP, ns, df.nm));
    SET_VECTOR_ELT(output, R1,        allocMatrix(REALSXP, ns, df.nm));
    SET_VECTOR_ELT(output, ALPHA_A,   allocMatrix(REALSXP, ns, df.nm));
    SET_VECTOR_ELT(output, R3,        allocMatrix(REALSXP, ns, df.nm));
    SET_VECTOR_ELT(output, BETA_A,    allocVector(REALSXP, ns));
    SET_VECTOR_ELT(output, PI_M,      allocVector(REALSXP, ns));
    SET_VECTOR_ELT(output, PI_A,      allocVector(REALSXP, ns));
    SET_VECTOR_ELT(output, SIGMA_M0,  allocVector(REALSXP, ns));
    SET_VECTOR_ELT(output, SIGMA_M1,  allocVector(REALSXP, ns));
    SET_VECTOR_ELT(output, SIGMA_MA0, allocVector(REALSXP, ns));
    SET_VECTOR_ELT(output, SIGMA_MA1, allocVector(REALSXP, ns));

    hdbm_mcmc(&df, &coeffs, asReal(pi_m), asReal(pi_a), asInteger(burnin), ns,
                  output);

    free(coeffs.alpha_a);
    free(coeffs.alpha_c2);
    free(coeffs.beta_m);
    free(coeffs.beta_c1);
    UNPROTECT(1);
    return output;
}
