#include <R.h>
#include <RcppArmadillo.h>

// constants for generating random inverse gammas
#define K    2.0
#define L_M0 1E-4
#define L_M1 1.0
#define L    1.0

// returns 1 with probability p
int rand_bernoulli(double p)
{
    return R::runif(0.0, 1.0) <= p ? 1 : 0;
}

class hdbm_mcmc {
    private:

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
    arma::vec beta_c1;
    arma::mat alpha_c2;

    double norm2_a;
    arma::vec norm2_c1;
    arma::vec norm2_c2;
    arma::vec norm2_m;

    void update_beta_m(arma::vec &beta_m, arma::mat &M, arma::vec &var_m0,
                           arma::vec &var_m1, double pi_m)
    {
        arma::vec mu_beta_m = M.t() * M * beta_m;
        for (arma::uword j = 0; j < M.n_cols; ++j) {
            double t = dot(M.col(j), rY);
            auto mbm0 = (mu_beta_m[j] + t) / (sigma_e / sigma_m0 + norm2_m[j]);
            auto mbm1 = (mu_beta_m[j] + t) / (sigma_e / sigma_m1 + norm2_m[j]);
            auto new_beta = r1[j] * R::rnorm(mbm1, sqrt(var_m1(j)))
                            + (1.0 - r1[j]) * R::rnorm(mbm0, sqrt(var_m0(j)));

            rY += M.col(j) * (beta_m[j] - new_beta);
            beta_m[j] = new_beta;
            /* update r1 */

            double p = exp(0.5 * (mbm1 * mbm1 / var_m1(j)
                       - mbm0 * mbm0 / var_m0(j) + log(var_m1(j))
                       - log(var_m0(j)) + log(sigma_m0) - log(sigma_m1)
                       ) + log(pi_m / (1.0 - pi_m))
            );
            r1[j] = rand_bernoulli(p / (1.0 + p));
        }
    }

    void update_alpha_a(arma::vec &alpha_a, arma::vec &A, double pi_a)
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
            double p = exp(0.5 * (maa1[j] * maa1[j] / var_aa1
                - maa0[j] * maa0[j] / var_aa0 + log(var_aa1) - log(var_aa0)
                + log(sigma_ma0) - log(sigma_ma1)
                ) + log(pi_a / (1.0 - pi_a))
            );
            r3[j] = rand_bernoulli(p / (1.0 + p));
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
        // This is what Yanyi wrote but I think it is incorrect
        for (arma::uword j = 0; j < C1.n_cols; ++j) {
            double mu_cj = 0.0;
            double old = beta_c1[j];
            for (arma::uword i = 0; i < C1.n_rows; ++i) {
                mu_cj += C1(i, j) * ( rY[i] + C1(i, j) * beta_c1[j]);
                mu_cj /= norm2_c1[j];
                beta_c1[j] = R::rnorm(mu_cj, sqrt(sigma_e /norm2_c1[j]));
            }
            rY += (old - beta_c1(j)) * C1.col(j);
        }
        // What I wrote -- Doesn't fix the bug(s)
        //arma::vec mu_beta_c = C1.t() * C1 * beta_c1;
        //for (arma::uword j = 0; j < C1.n_cols; ++j) {
        //    double t   = dot(C1.col(j), rY);
        //    double mbc = (mu_beta_c(j) + t) / norm2_c1(j);
        //    double new_beta_c1 = R::rnorm(mbc, sqrt(sigma_e / norm2_c1(j)));
        //    rY += C1.col(j) * (beta_c1(j) - new_beta_c1);
        //    beta_c1(j) = new_beta_c1;
        //}
    }

    void update_beta_a(arma::vec &A)
    {
        double var_a = sigma_e / (sigma_e / sigma_a + norm2_a);
        double mu_a = arma::dot(rY + A * beta_a, A) * var_a / sigma_e;

        double new_beta_a = R::rnorm(mu_a, sqrt(var_a));
        rY += A * (beta_a - new_beta_a);
        beta_a = new_beta_a;
    }

    public:

    hdbm_mcmc(arma::vec Y, arma::vec A, arma::mat M, arma::mat C1, arma::mat C2,
                  arma::vec beta_m, arma::vec alpha_a)
    {
        sigma_m0  = 1.0 / R::rgamma(K, 1.0 / L_M0);
        sigma_m1  = 1.0 / R::rgamma(K, 1.0 / L_M1);
        sigma_ma0 = 1.0 / R::rgamma(K, 1.0 / L_M0);
        sigma_ma1 = 1.0 / R::rgamma(K, 1.0 / L_M1);
        sigma_a   = 1.0 / R::rgamma(K, 1.0 / L);
        sigma_e   = 1.0 / R::rgamma(K, 1.0 / L);
        sigma_g   = 1.0 / R::rgamma(K, 1.0 / L);

        // Initialize residuals
        rY  = Y - A * beta_a - M * beta_m;
        rM  = M - A * alpha_a.t();
        rMC = rM;

        r1 = arma::vec(beta_m.n_elem, arma::fill::zeros);
        r3 = arma::vec(alpha_a.n_elem, arma::fill::zeros);

        beta_a = R::runif(0, sqrt(sigma_a));
        beta_c1 = arma::vec(C1.n_cols, arma::fill::zeros);
        alpha_c2 = arma::mat(C2.n_cols, M.n_cols, arma::fill::zeros);

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

    }

    void iterate(arma::vec &A, arma::mat &M, arma::mat &C1, arma::mat &C2,
                     arma::vec &beta_m, arma::vec &alpha_a, double &pi_m,
                     double &pi_a)
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

        Rcpp::Rcout << "rY^2 " << le1 << "\n";
        Rcpp::Rcout << "rM^2 " << lg1 << "\n";

        sigma_e = 1.0 / R::rgamma(K + 0.5 * n, 1.0 / (0.5 * le1 + L));
        sigma_g = 1.0 / R::rgamma(K + 0.5 * n * q, 1.0 / (0.5 * lg1 + L));

        Rcpp::Rcout << "beta_a "    << beta_a    << "\n";
        Rcpp::Rcout << "sigma_a "   << sigma_a   << "\n";
        Rcpp::Rcout << "sigma_e "   << sigma_e   << "\n";
        Rcpp::Rcout << "sigma_g "   << sigma_g   << "\n";
        Rcpp::Rcout << "sigma_m0 "  << sigma_m0  << "\n";
        Rcpp::Rcout << "sigma_ma0 " << sigma_ma0 << "\n";
        Rcpp::Rcout << "sigma_m1 "  << sigma_m1  << "\n";
        Rcpp::Rcout << "sigma_ma1 " << sigma_ma1 << "\n";

        update_beta_m(beta_m, M, var_m0, var_m1, pi_m);
        update_alpha_a(alpha_a, A, pi_a);
        update_alpha_c(C2);
        update_beta_c(C1);
        update_beta_a(A);

        double c1 = 0.5 * sum(r1);
        double c2 = 0.5 * sum(r1 % beta_m % beta_m);
        double c3 = 0.5 * sum(r3);
        double c4 = 0.5 * sum(r3 % alpha_a % alpha_a);

        double c5 = 0.5 * sum(1 - r1);
        // % means element by element multiplication
        double c6 = 0.5 * sum((1 - r1) % beta_m % beta_m);
        double c7 = 0.5 * sum(1 - r3);
        double c8 = 0.5 * sum((1 - r3) % (alpha_a % alpha_a));

        sigma_a   = 1.0 / R::rgamma(0.5 + K, 1.0 / (0.5 * beta_a * beta_a + L));
        sigma_m1  = 1.0 / R::rgamma(c1 + K, 1.0 / (c2 + L_M1));
        sigma_ma1 = 1.0 / R::rgamma(c3 + K, 1.0 / (c4 + L_M1));
        sigma_m0  = 1.0 / R::rgamma(c5 + K, 1.0 / (c6 + L_M0));
        sigma_ma0 = 1.0 / R::rgamma(c7 + K, 1.0 / (c8 + L_M0));

        double m_pi_m = abs(pi_m * exp(R::runif(-0.01, 0.01)));
        double m_pi_a = abs(pi_a * exp(R::runif(-0.01, 0.01)));
        if (m_pi_m  > 1.0)
            m_pi_m = 1.0 / m_pi_m;
        if (m_pi_m  < 1.0 / q)
            m_pi_m = 1.0 / (q * q * m_pi_m);
        if (m_pi_a  > 1.0)
            m_pi_a = 1.0 / m_pi_a;
        if (m_pi_a  < 1.0 / q)
            m_pi_a = 1.0 / (q * q * m_pi_a);

        // inline function
        auto post_dist = [](arma::vec r, double pi)
        {
            return sum((log(pi) * r) % (log(1.0 - pi) * (1.0 - r)));
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
arma::vec run_hdbm_mcmc(arma::vec Y, arma::vec A, arma::mat M, arma::mat C1,
                          arma::mat C2, arma::vec beta_m, arma::vec alpha_a,
                          double pi_m, double pi_a, int burnin, int nsamples)
{
    auto mcmc = hdbm_mcmc(Y, A, M, C1, C2, beta_m, alpha_a);
    // run mcmc for the number of specified burn in iterations
    for (int i = 0; i < burnin; ++i)
         mcmc.iterate(A, M, C1, C2, beta_m, alpha_a, pi_m, pi_a);

    return A;
}
