#include <R.h>
#include <RcppArmadillo.h>

// constants for generating random inverse gammas
#define K    2.0
#define L_M0 1E-4
#define L_M1 1.0
#define L    1.0

// returns 1 with probability p
// [[Rcpp::export]]
int rand_bernoulli(double p)
{
    return R::runif(0.0, 1.0) <= p ? 1 : 0;
}

// rand_invgamma with shape - scale parameterization is generated by
// 1 / gamma( shape, scale^-1) according to wikipedia
// [[Rcpp::export]]
double rand_invgamma(double shape, double scale)
{
    return 1.0 / Rcpp::rgamma(1, shape, 1.0 / scale)[0];
}

class hdbm_mcmc {
    public:

    double sigma_e;
    double sigma_g;
    double sigma_m0;
    double sigma_m1;
    double sigma_ma0;
    double sigma_ma1;
    double sigma_a;

    arma::vec rY;
    arma::mat rM;
    arma::mat rMC;

    arma::vec r1;
    arma::vec r3;

    double beta_a;
    arma::vec beta_m;
    arma::vec beta_c1;
    arma::vec alpha_a;
    arma::mat alpha_c2;

    double norm2_a;
    arma::vec norm2_c1;
    arma::vec norm2_c2;
    arma::vec norm2_m;

    double pi_m;
    double pi_a;

    void update_beta_m(arma::mat &M, arma::vec &var_m0, arma::vec &var_m1)
    {
        arma::vec mu_beta_m = diagmat(M.t() * M) * beta_m;
        for (arma::uword j = 0; j < M.n_cols; ++j) {
            double t = dot(M.col(j), rY);
            auto mbm0 = (mu_beta_m(j) + t) / (sigma_e / sigma_m0 + norm2_m[j]);
            auto mbm1 = (mu_beta_m(j) + t) / (sigma_e / sigma_m1 + norm2_m[j]);
            auto new_beta = r1[j] * R::rnorm(mbm1, sqrt(var_m1(j)))
                            + (1.0 - r1[j]) * R::rnorm(mbm0, sqrt(var_m0(j)));

            rY += M.col(j) * (beta_m[j] - new_beta);
            beta_m[j] = new_beta;

            // update r1
            double p = 0.5 * (mbm1 * mbm1 / var_m1(j) - mbm0 * mbm0 / var_m0(j)
                       + log(var_m1(j) / sigma_m1) - log(var_m0(j) / sigma_m0)
                       ) + log(pi_m / (1.0 - pi_m)
            );
            if (p > 30)
                r1[j] = 1;
             else
                r1[j] = rand_bernoulli(exp(p) / (1.0 + exp(p)));
        }
    }

    void update_alpha_a(arma::vec &A)
    {
        double var_aa0 = sigma_g / (sigma_g / sigma_ma0 + norm2_a);
        double var_aa1 = sigma_g / (sigma_g / sigma_ma1 + norm2_a);

        arma::vec maa0 = rMC.t() * A;
        arma::vec maa1 = maa0;

        arma::vec new_alpha_a = alpha_a;
        for (arma::uword j = 0; j < alpha_a.n_elem; ++j) {
            maa0[j] *= var_aa0 / sigma_g;
            maa1[j] *= var_aa1 / sigma_g;
            new_alpha_a[j] = r3[j] * R::rnorm(maa1[j], sqrt(var_aa1))
                             + (1.0 - r3[j]) * R::rnorm(maa0[j], sqrt(var_aa0));
        }
        rM += A * (alpha_a - new_alpha_a).t();
        alpha_a = new_alpha_a;

        // update r3
        for (arma::uword j = 0; j < r3.n_elem; ++j) {
            double p = 0.5 * (maa1[j] * maa1[j] / var_aa1
                - maa0[j] * maa0[j] / var_aa0 + log(var_aa1 / sigma_ma1)
                + log(var_aa0 / sigma_ma0)
            ) + log(pi_a / (1.0 - pi_a));
            if (p > 30)
                r3[j] = 1;
            else
                r3[j] = rand_bernoulli(exp(p) / (1.0 + exp(p)));
        }
    }

    void update_alpha_c(arma::mat &C2)
    {
        rM  += C2 * alpha_c2;
        rMC += C2 * alpha_c2;
        arma::mat mu_alpha_c = C2.t() * rM;
        for (arma::uword j = 0; j < rM.n_cols; ++j) {
            for (arma::uword k = 0; k < C2.n_cols; ++k) {
                alpha_c2(k, j) = R::rnorm(mu_alpha_c(k, j) / norm2_c2[k],
                                        sqrt(sigma_g / norm2_c2[k]));
            }
        }
        rM  -= C2 * alpha_c2;
        rMC -= C2 * alpha_c2;
    }

    void update_beta_c(arma::mat &C1)
    {
        arma::vec mu_beta_c = diagmat(C1.t() * C1) * beta_c1;
        for (arma::uword j = 0; j < C1.n_cols; ++j) {
           double t   = dot(C1.col(j), rY);
           double mbc = (mu_beta_c(j) + t) / norm2_c1(j);
           double new_beta_c1 = R::rnorm(mbc, sqrt(sigma_e / norm2_c1(j)));
           rY += C1.col(j) * (beta_c1(j) - new_beta_c1);
           beta_c1(j) = new_beta_c1;
        }
    }

    void update_beta_a(arma::vec &A)
    {
        double var_a = sigma_e / (sigma_e / sigma_a + norm2_a);
        double mu_a = arma::dot(rY + A * beta_a, A) * var_a / sigma_e;

        double new_beta_a = R::rnorm(mu_a, sqrt(var_a));
        rY += A * (beta_a - new_beta_a);
        beta_a = new_beta_a;
    }

    hdbm_mcmc(arma::vec Y, arma::vec A, arma::mat M, arma::mat C1, arma::mat C2,
                  arma::vec beta_m, arma::vec alpha_a, double pi_m, double pi_a)
    {
        sigma_m0  = rand_invgamma(K, L_M0);
        sigma_m1  = rand_invgamma(K, L_M1);
        sigma_ma0 = rand_invgamma(K, L_M0);
        sigma_ma1 = rand_invgamma(K, L_M1);
        sigma_a   = rand_invgamma(K, L);
        sigma_e   = rand_invgamma(K, L);
        sigma_g   = rand_invgamma(K, L);

        this->beta_m = beta_m;
        this->alpha_a = alpha_a;

        // beta_a must be set before initializing the residuals!
        beta_a   = R::runif(0, sqrt(sigma_a));
        beta_c1  = arma::vec(C1.n_cols, arma::fill::zeros);
        alpha_c2 = arma::mat(C2.n_cols, M.n_cols, arma::fill::zeros);

        // Initialize residuals
        rY  = Y - A * beta_a - M * beta_m;
        rM  = M - A * alpha_a.t();
        rMC = rM;

        r1 = arma::vec(beta_m.n_elem, arma::fill::zeros);
        r3 = arma::vec(alpha_a.n_elem, arma::fill::zeros);

        norm2_a = dot(A, A);
        norm2_m = arma::vec(M.n_cols);
        for (arma::uword j = 0; j < M.n_cols; ++j)
            norm2_m[j] = dot(M.col(j), M.col(j));

        norm2_c1 = arma::vec(C1.n_cols);
        for (arma::uword j = 0; j < C1.n_cols; ++j)
            norm2_c1[j] = dot(C1.col(j), C1.col(j));

        norm2_c2 = arma::vec(C2.n_cols);
        for (arma::uword j = 0; j < C2.n_cols; ++j)
            norm2_c2[j] = dot(C2.col(j), C2.col(j));

        this->pi_m = pi_m;
        this->pi_a = pi_a;
    }

    void iterate(arma::vec &A, arma::mat &M, arma::mat &C1, arma::mat &C2)
    {
        int n = rM.n_rows;
        int q = rM.n_cols;

        arma::vec var_m1 = arma::vec(q);
        arma::vec var_m0 = arma::vec(q);
        for (int j = 0; j < q; ++j) {
            var_m0(j) = sigma_e / (sigma_e / sigma_m0 + norm2_m(j));
            var_m1(j) = sigma_e / (sigma_e / sigma_m1 + norm2_m(j));
        }

        double le1 = dot(rY, rY);
        // square every element of rM, and sum
        double lg1 = arma::accu(rM % rM);

        sigma_e = rand_invgamma(K + 0.5 * n, 0.5 * le1 + L);
        sigma_g = rand_invgamma(K + 0.5 * n * q, 0.5 * lg1 + L);

        //Rcpp::Rcout << "rY^2 " << le1 << "\n";
        //Rcpp::Rcout << "rM^2 " << lg1 << "\n";
        //Rcpp::Rcout << "beta_a "    << beta_a    << "\n";
        //Rcpp::Rcout << "pi_m "      << pi_m      << "\n";
        //Rcpp::Rcout << "pi_a "      << pi_a      << "\n";
        //Rcpp::Rcout << "sigma_a "   << sigma_a   << "\n";
        //Rcpp::Rcout << "sigma_e "   << sigma_e   << "\n";
        //Rcpp::Rcout << "sigma_g "   << sigma_g   << "\n";
        //Rcpp::Rcout << "sigma_m0 "  << sigma_m0  << "\n";
        //Rcpp::Rcout << "sigma_ma0 " << sigma_ma0 << "\n";
        //Rcpp::Rcout << "sigma_m1 "  << sigma_m1  << "\n";
        //Rcpp::Rcout << "sigma_ma1 " << sigma_ma1 << "\n";

        update_beta_m(M, var_m0, var_m1);
        update_alpha_a(A);
        update_alpha_c(C2);
        update_beta_c(C1);
        update_beta_a(A);

        // "%" means element by element multiplication
        double c1 = 0.5 * sum(r1);
        double c2 = 0.5 * sum(r1 % beta_m % beta_m);
        double c3 = 0.5 * sum(r3);
        double c4 = 0.5 * sum(r3 % alpha_a % alpha_a);

        double c5 = 0.5 * sum(1 - r1);
        // % means element by element multiplication
        double c6 = 0.5 * sum((1 - r1) % beta_m % beta_m);
        double c7 = 0.5 * sum(1 - r3);
        double c8 = 0.5 * sum((1 - r3) % (alpha_a % alpha_a));

        sigma_a   = rand_invgamma(0.5 + K, 0.5 * beta_a * beta_a + L);
        sigma_m1  = rand_invgamma(c1 + K, c2 + L_M1);
        sigma_ma1 = rand_invgamma(c3 + K, c4 + L_M1);
        sigma_m0  = rand_invgamma(c5 + K, c6 + L_M0);
        sigma_ma0 = rand_invgamma(c7 + K, c8 + L_M0);

        double m_pi_m = std::abs(pi_m * exp(R::runif(-.01, 0.01)));
        double m_pi_a = std::abs(pi_a * exp(R::runif(-.01, 0.01)));

        if (m_pi_m  > 1.0)
            m_pi_m = 1.0 / m_pi_m;
        if (m_pi_m  < 1.0 / q)
            m_pi_m = 1.0 / (q * q * m_pi_m);
        if (m_pi_a  > 1.0)
            m_pi_a = 1.0 / m_pi_a;
        if (m_pi_a  < 1.0 / q)
            m_pi_a = 1.0 / (q * q * m_pi_a);

        //inline function
        auto post_dist = [](arma::vec r, double pi)
        {
            return sum((log(pi) * r) + (log(1.0 - pi) * (1.0 - r)));
        };

        double p = post_dist(r3, m_pi_a) - post_dist(r3, pi_a)
                   + post_dist(r1, m_pi_m) - post_dist(r1, pi_m);

        if (p > log(R::runif(0.0, 1.0))) {
            pi_a = m_pi_a;
            pi_m = m_pi_m;
         }
    }
};

// [[Rcpp::depends(RcppArmadillo)]]
//' @export
//' @useDynLib hdbm
//' @importFrom Rcpp evalCpp
// [[Rcpp::export]]
Rcpp::List run_hdbm_mcmc(arma::vec Y, arma::vec A, arma::mat M, arma::mat C1,
                          arma::mat C2, arma::vec beta_m_init, arma::vec
                          alpha_a_init, double pi_m_init, double pi_a_init,
                          int burnin, int nsamples)
{
    auto mcmc = hdbm_mcmc(Y, A, M, C1, C2, beta_m_init, alpha_a_init, pi_m_init,
                              pi_a_init);

    // run mcmc for the number of specified burn in iterations
    for (int i = 0; i < burnin; ++i)
         mcmc.iterate(A, M, C1, C2);

    auto beta_m    = Rcpp::NumericMatrix(nsamples, beta_m_init.n_elem);
    auto r1        = Rcpp::NumericMatrix(nsamples, beta_m_init.n_elem);
    auto alpha_a   = Rcpp::NumericMatrix(nsamples, alpha_a_init.n_elem);
    auto r3        = Rcpp::NumericMatrix(nsamples, alpha_a_init.n_elem);
    auto beta_a    = Rcpp::NumericVector(nsamples);
    auto pi_m      = Rcpp::NumericVector(nsamples);
    auto pi_a      = Rcpp::NumericVector(nsamples);
    auto sigma_m0  = Rcpp::NumericVector(nsamples);
    auto sigma_m1  = Rcpp::NumericVector(nsamples);
    auto sigma_ma0 = Rcpp::NumericVector(nsamples);
    auto sigma_ma1 = Rcpp::NumericVector(nsamples);


    for (int i = 0; i < nsamples; ++i) {
        Rcpp::NumericMatrix::Row bm_row = beta_m(i, Rcpp::_);
        bm_row = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(mcmc.beta_m));

        Rcpp::NumericMatrix::Row r1_row = r1(i, Rcpp::_);
        r1_row = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(mcmc.r1));

        Rcpp::NumericMatrix::Row aa_row = alpha_a(i, Rcpp::_);
        aa_row = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(mcmc.alpha_a));

        Rcpp::NumericMatrix::Row r3_row = r3(i, Rcpp::_);
        r3_row = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(mcmc.r3));

        beta_a[i]    = mcmc.beta_a;
        pi_a[i]      = mcmc.pi_a;
        pi_m[i]      = mcmc.pi_m;
        sigma_m0[i]  = mcmc.sigma_m0;
        sigma_m1[i]  = mcmc.sigma_m1;
        sigma_ma0[i] = mcmc.sigma_ma0;
        sigma_ma1[i] = mcmc.sigma_ma1;
        int j = 0;
        while (j++ < 50)
            mcmc.iterate(A, M, C1, C2);
    }
    // need to write the sampler
    return Rcpp::List::create(Rcpp::Named("beta.m") = beta_m,
        Rcpp::Named("r1") = r1, Rcpp::Named("alpha.a") = alpha_a,
        Rcpp::Named("r3") = r3, Rcpp::Named("beta.a") = beta_a,
        Rcpp::Named("pi.m") = pi_m, Rcpp::Named("pi.a") = pi_a,
        Rcpp::Named("sigma.m0") = sigma_m0, Rcpp::Named("sigma.m1") = sigma_m1,
        Rcpp::Named("sigma.ma0") = sigma_ma0,
        Rcpp::Named("sigma.ma1") = sigma_ma1
    );
}
