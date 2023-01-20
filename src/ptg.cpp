
#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist, BH)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// [[Rcpp::export]]
void ptg(arma::vec Y,
         arma::vec A,
         arma::mat M,
         arma::mat C1,
         arma::mat C2,
         int burnin,
         int ndraws,
         int seed, 
         double lambda0,
         double lambda1,
         double lambda2,
         double ha,
         double la,
         double h1,
         double l1,
         double h2,
         double l2,
         double km,
         double lm,
         double kma,
         double lma,
         arma::mat& samples_betam,
         arma::mat& samples_alphaa,
         arma::mat& betam_member,
         arma::mat& alphaa_member,
         arma::mat& samples_betaa,
         arma::mat& samples_alphac,
         arma::mat& samples_betac,
         arma::mat& samples_sigmasqa,
         arma::mat& samples_sigmasqe,
         arma::mat& samples_sigmasqg
) {
  
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(std::floor(std::fabs(seed)));
  
  int n = M.n_rows;
  int p = M.n_cols; // 0 < j < p
  int q1 = C1.n_cols;
  int q2 = C2.n_cols;
  
  arma::vec betam_til(p, arma::fill::zeros);
  arma::vec alphaa_til(p, arma::fill::zeros);
  
  arma::vec betam(p, arma::fill::zeros);
  arma::vec alphaa(p, arma::fill::zeros);
  
  arma::vec betac(q1, arma::fill::zeros);
  arma::mat alphac(q2, p, arma::fill::zeros);
  
  double tausq_b = 1.0;
  double tausq_a = 1.0;
  
  double B1 = 1, B2 = 1, B3 = 1;
  double A1 = 1, A2 = 1, A3 = 1;
  
  double sigmasqe = 1.0;
  double sigmasqa = 1.0;
  double sigmasqg = 1.0;
  double sigmasq1 = 1.0;
  
  arma::vec sseq(3, arma::fill::zeros);
  sseq(0) = 1.0;
  sseq(1) = 2.0;
  sseq(2) = 3.0;
  
  double betaa = 0.0;
  
  arma::vec prop(3, arma::fill::zeros);
  arma::vec propA(3, arma::fill::zeros);
  int Bk = 0;
  int Ak = 0;
  
  double inf = std::numeric_limits<double>::infinity();
  
  int draw = 0;
  for(int it = 0; it < ndraws; it++) {
    // sample (beta_m)j
    double u_betam_tillj = 0.0;
    for(int nj = 0; nj < p; nj++) {
      if(alphaa_til(nj) == 0.0) {
        u_betam_tillj = lambda1;
      } else {
        arma::vec lambmin(2, arma::fill::zeros);
        lambmin(0) = lambda1;
        lambmin(1) = lambda0/std::abs(alphaa_til(nj));
        u_betam_tillj = arma::min(lambmin);
      }
      
      double mu_mj = 0.0;
      double sigsq_mj = 0.0;
      for(int ui = 0; ui < n; ui++) {
        mu_mj += as_scalar(M(ui, nj) * (Y(ui) - A(ui) * betaa - M.row(ui) * betam -
          C1.row(ui) * betac + M(ui, nj) * betam(nj)));
      }
      mu_mj /= ((sigmasqe / tausq_b) + accu(M.col(nj) % M.col(nj)));
      sigsq_mj = 1 / ((1 / tausq_b) + (accu(M.col(nj) % M.col(nj)) / sigmasqe));
      
      B1 = 1 - 2.0 * R::pnorm(-1.0 * (u_betam_tillj / std::sqrt(tausq_b)), 0.0, 1.0, true, false);
      B2 = std::log(1 - R::pnorm(u_betam_tillj / std::sqrt(sigsq_mj), mu_mj, 1.0, true, false));
      B3 = std::log(R::pnorm((-u_betam_tillj / std::sqrt(sigsq_mj)), mu_mj, 1.0, true, false));
      
      B1 = std::log(B1);
      B2 = ((mu_mj * mu_mj) / (2.0 * (sigsq_mj))) + 
        0.5 * log(std::sqrt(sigsq_mj)/std::sqrt(tausq_b)) + B2;
      B3 = ((mu_mj * mu_mj) / (2.0 * (sigsq_mj))) + 
        0.5 * log(std::sqrt(sigsq_mj)/std::sqrt(tausq_b)) + B3;
      
      double Bmax = std::max(B1, std::max(B2, B3));
      
      B1 = std::exp(B1 - Bmax);
      B2 = std::exp(B2 - Bmax);
      B3 = std::exp(B3 - Bmax);
      
      if(B1 > 10000) { B1 = 10000; }
      if(B2 > 10000) { B2 = 10000; }
      if(B3 > 10000) { B3 = 10000; }
      
      double sumB = 0;
      sumB = B1 + B2 + B3;
      prop(0) = B1 / sumB;
      prop(1) = B2 / sumB;
      prop(2) = B3 / sumB;
      arma::vec Bk_samp(1, arma::fill::zeros);
      
      Bk_samp = RcppArmadillo::sample(sseq, 1, true, prop);
      Bk = Bk_samp(0);
      
      // sample betam_til from truncated normal distribution
      //abs(betam_til(nj)) < u_betam_tillj
      //betam_til(nj) >= u_betam_tillj
      //
      if (Bk == 1) {
        betam_til(nj) = r_truncnorm(0, std::sqrt(tausq_b), -u_betam_tillj, u_betam_tillj);
      } else if (Bk == 2) {
        betam_til(nj) = r_truncnorm(mu_mj, std::sqrt(sigsq_mj), u_betam_tillj, inf);
        if(std::isinf(betam_til(nj))) {
          betam_til(nj) = u_betam_tillj;
        }
      } else if (Bk == 3) {
        betam_til(nj) = r_truncnorm(mu_mj, std::sqrt(sigsq_mj), -inf, -u_betam_tillj);
        if(std::isinf(betam_til(nj))) {
          betam_til(nj) = -u_betam_tillj;
        }
      }
      
      if(std::abs(betam_til(nj)) < u_betam_tillj) {
        betam(nj) = 0.0;
      } else if(betam_til(nj) >= u_betam_tillj) {
        betam(nj) = betam_til(nj);
      } else if(betam_til(nj) <= -u_betam_tillj){
        betam(nj) = betam_til(nj);
      }
      if(draw >= burnin) {
        samples_betam(draw - burnin, nj) = betam(nj);
        if(betam(nj) != 0.0) {
          betam_member(draw - burnin, nj) = 1;
        }
      }
    }
    // sample (alpha_a)j
    double u_alphaa_tillj = 0.0;
    for(int nj = 0; nj < p; nj++) {
      if(betam_til(nj) == 0.0) {
        u_alphaa_tillj = lambda2;
      } else {
        arma::vec lambmin(2, arma::fill::zeros);
        lambmin(0) = lambda2;
        lambmin(1) = lambda0/std::abs(betam_til(nj));
        u_alphaa_tillj = arma::min(lambmin);
      }
      
      double mu_aj = 0.0;
      double sigsq_aj = 0.0;
      for(int ui = 0; ui < n; ui++) {
        arma::mat CibyAlphac(1, p, arma::fill::zeros);
        CibyAlphac = C2.row(ui) * alphac;
        mu_aj += as_scalar(A(ui) * (M(ui, nj) - CibyAlphac.col(nj)));
      }
      mu_aj /= ((sigmasqg / tausq_a) + accu(A % A));
      sigsq_aj = 1 / ((1 / tausq_a) + (accu(A % A) / sigmasqg));
      
      A1 = 1 - 2.0 * R::pnorm(-1.0 * (u_alphaa_tillj / std::sqrt(tausq_a)), 0.0, 1.0, true, false);
      A2 = std::log(1.0 - R::pnorm(u_alphaa_tillj / std::sqrt(sigsq_aj), mu_aj, 1.0, true, false));
      A3 = std::log(R::pnorm((-u_alphaa_tillj / std::sqrt(sigsq_aj)), mu_aj, 1.0, true, false));
      
      A1 = std::log(A1);
      A2 = ((mu_aj * mu_aj) / (2.0 * (sigsq_aj))) + 
        0.5 * log(std::sqrt(sigsq_aj)/std::sqrt(tausq_a)) + A2;
      A3 = ((mu_aj * mu_aj) / (2.0 * (sigsq_aj))) + 
        0.5 * log(std::sqrt(sigsq_aj)/std::sqrt(tausq_a)) + A3;
      
      double Amax = std::max(A1, std::max(A2, A3));
      
      A1 = std::exp(A1 - Amax);
      A2 = std::exp(A2 - Amax);
      A3 = std::exp(A3 - Amax);
      
      if(A1 > 10000) { A1 = 10000; }
      if(A2 > 10000) { A2 = 10000; }
      if(A3 > 10000) { A3 = 10000; }
      
      double sumA = 0;
      sumA = A1 + A2 + A3;
      propA(0) = A1 / sumA;
      propA(1) = A2 / sumA;
      propA(2) = A3 / sumA;
      arma::vec Ak_samp(1, arma::fill::zeros);
      
      Ak_samp = RcppArmadillo::sample(sseq, 1, true, propA);
      Ak = Ak_samp(0);
      
      // sample alphaa_til from truncated normal distribution
      // abs(alphaa_til(nj)) < u_alphaa_tillj
      // alphaa_til(nj) >= u_alphaa_tillj
      if (Ak == 1) {
        alphaa_til(nj) = r_truncnorm(0, std::sqrt(tausq_a), -u_alphaa_tillj, u_alphaa_tillj);
      } else if (Ak == 2) {
        alphaa_til(nj) = r_truncnorm(mu_aj, std::sqrt(sigsq_aj), u_alphaa_tillj, inf);
        if(std::isinf(alphaa_til(nj))) {
          alphaa_til(nj) = u_alphaa_tillj;
        }
      } else if (Ak == 3) {
        alphaa_til(nj) = r_truncnorm(mu_aj, std::sqrt(sigsq_aj), -inf, -u_alphaa_tillj);
        if(std::isinf(alphaa_til(nj))) {
          alphaa_til(nj) = -u_alphaa_tillj;
        }
      }
      
      if (std::abs(alphaa_til(nj)) < u_alphaa_tillj) {
        alphaa(nj) = 0.0;
      } else if (alphaa_til(nj) >= u_alphaa_tillj) {
        alphaa(nj) = alphaa_til(nj);
      } else if (alphaa_til(nj) <= -u_alphaa_tillj){
        alphaa(nj) = alphaa_til(nj);
      }
      if(draw >= burnin) {
        samples_alphaa(draw - burnin, nj) = alphaa(nj);
        if(alphaa(nj) != 0.0) {
          alphaa_member(draw - burnin, nj) = 1;
        }
      }
    }
    
    // sample betaa
    double mu_betaa = 0.0;
    double sigsq_betaa = 0.0;
    for(int ui = 0; ui < n; ui++) {
      mu_betaa += as_scalar(A(ui) * (Y(ui) - M.row(ui) * betam - C1.row(ui) * betac));
    }
    mu_betaa /= ((sigmasq1 / sigmasqa) + accu(A % A));
    sigsq_betaa = 1 / ((1 / sigmasqa) + (accu(A % A) / sigmasq1));
    betaa = R::rnorm(mu_betaa, std::sqrt(sigsq_betaa));
    if(draw >= burnin) {
      samples_betaa(draw - burnin, 0) = betaa;
    }
    // sample sigsqa
    sigmasqa = 1.0 / Rcpp::rgamma(1, (0.5 + ha), 1.0 / (((betaa * betaa) / 2) + la))[0];
    if(draw >= burnin) {
      samples_sigmasqa(draw - burnin, 0) = sigmasqa;
    }
    
    // sample sigsqe
    double l_sigsqe = 0.0;
    for(int ui = 0; ui < n; ui++) {
      double lsigres = as_scalar((Y(ui) - M.row(ui) * betam - A(ui) * betaa - C1.row(ui) * betac));
      l_sigsqe += lsigres * lsigres;
    }
    sigmasqe = 1.0 / Rcpp::rgamma(1, (n / 2) + h1, 1.0 / ((l_sigsqe / 2) + l1))[0];
    if(draw >= burnin) {
      samples_sigmasqe(draw - burnin, 0) = sigmasqe;
    }
    
    // sample sigsqg
    double l_sigsqg = 0.0;
    for(int ui = 0; ui < n; ui++) {
      arma::rowvec lsigsqg(p, arma::fill::zeros);
      lsigsqg = (M.row(ui) - (A(ui) * alphaa).t() - C2.row(ui) * alphac);
      l_sigsqg += accu(lsigsqg % lsigsqg);
    }
    sigmasqg = 1.0 / Rcpp::rgamma(1, ((p * n) / 2) + h2, 1.0 / ((l_sigsqg / 2) + l2))[0];
    if(draw >= burnin) {
      samples_sigmasqg(draw - burnin, 0) = sigmasqg;
    }
    
    // sample tausq_b
    double sum_betamjsq = 0.0;
    for(int jp = 0; jp < p; jp++) {
      sum_betamjsq += betam_til(jp) * betam_til(jp);
    }
    tausq_b = 1.0 / Rcpp::rgamma(1, (p / 2) + km, 1.0 / ((sum_betamjsq / 2) + lm))[0];
    
    // sample tausq_a
    double sum_alphaajsq = 0.0;
    for(int jp = 0; jp < p; jp++) {
      sum_alphaajsq += alphaa_til(jp) * alphaa_til(jp);
    }
    tausq_a = 1.0 / Rcpp::rgamma(1, (p / 2) + kma, 1.0 / ((sum_alphaajsq / 2) + lma))[0];
    
    // sample beta_cw
    for(int w = 0; w < q1; w++) {
      double mu_betac_w = 0.0;
      for(int ni = 0; ni < n; ni++) {
        mu_betac_w += as_scalar(C1(ni, w) * (Y(ni) - A(ni) * betaa - M.row(ni) *
          betam - C1.row(ni) * betac + C1(ni, w) * betac(w)));
      }
      betac(w) = R::rnorm(mu_betac_w / accu(C1.col(w) % C1.col(w)),
            std::sqrt(sigmasqe / accu(C1.col(w) % C1.col(w))));
      if(draw >= burnin) {
        samples_betac(draw - burnin, w) = betac(w);
      }
    }
    
    // sample (alpha_cw)j
    for(int w = 0; w < q2; w++) {
      for(int nj = 0; nj < p; nj++) {
        double mu_alphac_wj = 0.0;
        for(int ni = 0; ni < n; ni++) {
          mu_alphac_wj += as_scalar(C2(ni, w) * (M(ni, nj) - A(ni) * alphaa(nj) -
            C2.row(ni) * alphac.col(nj) + C2(ni, w) * alphac(w, nj)));
        }
        mu_alphac_wj /= accu(C2.col(w) % C2.col(w));
        alphac(w, nj) = R::rnorm(mu_alphac_wj, 
               std::sqrt(sigmasqg / accu(C2.col(w) % C2.col(w))));
        if(draw >= burnin) {
          samples_alphac(draw - burnin, w * p + nj) = alphac(w, nj);
        }
      }
    }
    draw += 1;
  }
}


// [[Rcpp::export]]
arma::vec rdirichletcpp(arma::vec alpha_m) {
  int distribution_size = alpha_m.n_elem;
  // each element will be a draw from a Dirichlet
  arma::vec distribution(distribution_size, arma::fill::zeros);
  
  double sum_term = 0;
  // loop through the distribution and draw Gamma variables
  for (int j = 0; j < distribution_size; ++j) {
    double cur = R::rgamma(alpha_m[j], 1.0);
    distribution(j) = cur;
    sum_term += cur;
  }
  // now normalize
  for (int j = 0; j < distribution_size; ++j) {
    distribution(j) = distribution(j)/sum_term;
  }
  return(distribution);
}

// [[Rcpp::export]]
double rand_igamma(double shape, double scale){
  return 1.0 / Rcpp::rgamma(1, shape, 1.0 / scale)[0];
}
// 
// // [[Rcpp::export]]
// void gmm(arma::vec Y,
//          arma::vec A,
//          arma::mat M,
//          arma::mat C1,
//          arma::mat C2,
//          int burnin,
//          int ndraws,
//          int seed,
//          double phi0,
//          double phi1,
//          int a0,
//          int a1,
//          int a2,
//          int a3,
//          double h1,
//          double l1,
//          double h2,
//          double l2,
//          double ha,
//          double la,
//          arma::mat& samples_betam,
//          arma::mat& samples_alphaa,
//          arma::mat& betam_member,
//          arma::mat& alphaa_member,
//          arma::mat& samples_betaa,
//          arma::mat& samples_alphac,
//          arma::mat& samples_betac,
//          arma::mat& samples_sigmasqa,
//          arma::mat& samples_sigmasqe,
//          arma::mat& samples_sigmasqg
// ) {
//   
//   
//   
//   
//   Rcpp::Environment base_env("package:base");
//   Rcpp::Function set_seed_r = base_env["set.seed"];
//   set_seed_r(std::floor(std::fabs(seed)));
//   
//   int n = M.n_rows;
//   int p = M.n_cols; // 0 < j < p
//   int q1 = C1.n_cols;
//   int q2 = C2.n_cols;
//   
//   arma::vec pi(4, arma::fill::zeros);
//   pi(0) = 0.25; pi(1) = 0.25;
//   pi(2) = 0.25; pi(3) = 0.25;
//   
//   int k = 3;
//   int v = 4;
//   
//   arma::vec diagPhinot(2, arma::fill::ones);
//   arma::mat Phinot(2, 2, arma::fill::zeros);
//   Phinot(0, 0) = phi0;
//   Phinot(1, 1) = phi1;
//   
//   arma::mat Wjmat(2, 2, arma::fill::eye);
//   arma::mat lwjmat(2, 1, arma::fill::ones); // need this?
//   
//   arma::vec betam(p, arma::fill::zeros);
//   arma::vec alphaa(p, arma::fill::zeros);
//   
//   arma::vec betac(q1, arma::fill::zeros);
//   arma::mat alphac(q2, p, arma::fill::zeros);
//   
//   arma::vec gamk(p, arma::fill::zeros); // need this?
//   
//   arma::vec logpgam(4, arma::fill::ones);
//   
//   arma::mat Vk1(2, 2, arma::fill::zeros);
//   Vk1(0, 0) = 1;
//   Vk1(0, 1) = 0.5;
//   Vk1(1, 0) = 0.5;
//   Vk1(1, 1) = 0.9;
//   
//   arma::mat Vk2(2, 2, arma::fill::zeros);
//   Vk2(0, 0) = 1;
//   
//   arma::mat Vk3(2, 2, arma::fill::zeros);
//   Vk3(1, 1) = 1;
//   
//   double sigmasqe = 1.0;
//   double sigmasqa = 1.0;
//   double sigmasqg = 1.0;
//   
//   arma::vec sseq(4, arma::fill::zeros);
//   sseq(0) = 0;
//   sseq(1) = 1;
//   sseq(2) = 2;
//   sseq(3) = 3; // start at 1?
//   
//   double betaa = 0.0;
//   
//   arma::mat bmaa(2, p, arma::fill::zeros);
//   
//   int draw = 0;
//   for(int it = 0; it < ndraws; it++) {
//     
//     for(int j = 0; j < p; j++) {
//       // Wj
//       arma::vec mj = M.col(j);
//       
//       Wjmat(0, 0) = accu(mj % mj) / sigmasqe;
//      
//       Wjmat(1, 1) = accu(A % A) / sigmasqg;
//       
//       //lwj
//      
//       lwjmat(0) = accu( (Y - A * betaa - M * betam + mj * betam(j)) % mj ) / sigmasqe;
//       lwjmat(1) = accu(mj % A) / sigmasqg;
//      
//       
//       // logpgam
//       arma::vec lggam(1, arma::fill::zeros);
//       lggam = ( - 0.5 * std::log(det(Wjmat + Vk1.i())) - 0.5 * std::log(det(Vk1)) + 
//         0.5 * lwjmat.t() * ((Wjmat + Vk1.i()).i()) * lwjmat + std::log(pi(0) + 0.01));
//       logpgam(0) = lggam(0);
//       lggam = ( - 0.5 * std::log(abs(Wjmat(0, 0) + (1 / Vk2(0, 0)))) - 
//         0.5 * log(Vk2(0, 0)) + 0.5 * lwjmat(0) * (1 / (Wjmat(0, 0) + (1 / Vk2(0, 0)))) * lwjmat(0) + std::log(pi(1) + 0.01));
//       logpgam(1) = lggam(0);
//       lggam = ( - 0.5 * std::log(abs(Wjmat(1, 1) + (1 / Vk3(1, 1)))) - 
//         0.5 * log(Vk3(1, 1)) + 0.5 * lwjmat(1) * (1 / (Wjmat(1, 1) + (1 / Vk3(1, 1)))) * lwjmat(1) + std::log(pi(2) + 0.01));
//       logpgam(2) = lggam(0);
//       logpgam(3) = std::log(pi(3) + 0.01);
//       
//       // sample k
//       double max = -100.0;
//       for(int i = 0; i < 4; i++) {
//         if(logpgam(i) > 10000) logpgam(i) = 10000;
//         if(logpgam(i) < -10000) logpgam(i) = -10000;
//         if(logpgam(i) > max) max = logpgam(i);
//       }
//       arma::vec ep = exp(logpgam - max);
//       arma::vec eprop = ep / accu(ep);
//       
//       k = as_scalar( RcppArmadillo::sample_main(sseq, 1, true, eprop) );
//       
//       // sample betam and alphaa
//       arma::mat WjVkInvInv(2, 2, arma::fill::zeros);
//       
//       switch(k) {
//       case 0: WjVkInvInv = (Wjmat + Vk1.i()).i();
//         bmaa.col(j) = rmvnorm(1, WjVkInvInv * lwjmat, WjVkInvInv).t();
//         betam(j) = bmaa(0,j);
//         alphaa(j) = bmaa(1,j);
//         break;
//       case 1: WjVkInvInv(0, 0) = 1 / (Wjmat(0, 0) + (1 / Vk2(0, 0)));
//         bmaa(0,j) = R::rnorm(WjVkInvInv(0, 0) * lwjmat(0), sqrt(WjVkInvInv(0, 0)));
//         bmaa(1,j) = 0.0;
//         betam(j) = bmaa(0,j);
//         break;
//       case 2: WjVkInvInv(1, 1) = 1 / (Wjmat(1, 1) + (1 / Vk3(1, 1)));
//         bmaa(0,j) = 0.0;
//         bmaa(1,j) = R::rnorm(WjVkInvInv(1, 1) * lwjmat(1), sqrt(WjVkInvInv(1, 1)));
//         alphaa(j) = bmaa(1,j);
//         break;
//       case 3: bmaa(0,j) = 0.0;
//         bmaa(1,j) = 0.0;
//         betam(j) = 0.0;
//         alphaa(j) = 0.0;
//         break;
//       }
//       
//       if(draw >= burnin) {
//         samples_betam(draw - burnin, j) = betam(j);
//         samples_alphaa(draw - burnin, j) = alphaa(j);
//         if(betam(j) != 0.0) {
//           betam_member(draw - burnin, j) = 1;
//         }
//         if(alphaa(j) != 0.0) {
//           alphaa_member(draw - burnin, j) = 1;
//         }
//       }
//     }
//     
//     // sample pi
//     arma::vec a(4, arma::fill::zeros);
//     a(0) = a0; a(1) = a1; a(2) = a2; a(3) = a3;
//     for(int i = 0; i < p; i++) {
//       int k = gamk(i);
//       switch(k) {
//       case 0: a(0) += 1;
//         break;
//       case 1: a(1) += 1;
//         break;
//       case 2: a(2) += 1;
//         break;
//       case 3: a(3) += 1;
//         break;
//       }
//     }
//     
//     pi = rdirichletcpp(a);
//     
//     // sample the vk's
//     double sgam = accu(gamk);
//    
//     arma::mat gam_sum(2, 2, arma::fill::zeros);
//     for(int t = 0; t < p; t++) {
//       gam_sum += gamk(t) * (bmaa(0,p) * bmaa(0,p) + bmaa(1,p) * bmaa(1,p));
//     }
//     Vk1 = riwish(v + sgam, Phinot + gam_sum);
//     
//     double gam_sum1 = 0;
//     for(int t = 0; t < p; t++) {
//       gam_sum1 += gamk(t) * (bmaa(0,p) * bmaa(0,p));
//     }
//     Vk2(0, 0) = rand_igamma((v / 2) + sgam, (Phinot(0, 0) / 2) + gam_sum1);
//     
//     double gam_sum2 = 0;
//     for(int t = 0; t < p; t++) {
//       gam_sum2 += gamk(t) * (bmaa(1,p) * bmaa(1,p));
//     }
//     Vk3(1, 1) = rand_igamma((v / 2) + sgam, (Phinot(1, 1) / 2) + gam_sum2);
//     
//     // sample betaa
//     arma::vec mvar(2, arma::fill::zeros);
//     mvar(0) = accu(A % (Y - M * betam - C1 * betac)) / ((sigmasqe / sigmasqa) + accu(A % A));
//     mvar(1) = std::sqrt(1 / ((1 / sigmasqa) + (accu(A % A) / sigmasqe)));
//     betaa = R::rnorm(mvar(0), mvar(1));
//     
//     if(draw >= burnin) {
//       samples_betaa(draw - burnin, 0) = betaa;
//     }
//     
//     // sigmasqa
//     sigmasqa = rand_igamma(0.5 + ha, 0.5 * betaa * betaa + la);
//     
//     if(draw >= burnin) {
//       samples_sigmasqa(draw - burnin, 0) = sigmasqa;
//     }
//     
//     // sample sigmasqe
//     double sum_r = accu(Y - M * betam - A * betaa - C1 * betac);
//     sigmasqe = rand_igamma(0.5 * n + h1, 0.5 * sum_r * sum_r + l1);
//     
//     if(draw >= burnin) {
//       samples_sigmasqe(draw - burnin, 0) = sigmasqe;
//     }
//     
//     // sample sigmasqg
//     double resg = 0;
//     arma::vec resv(1, arma::fill::zeros);
//     for(int i = 0; i < n; i++) {
//       resv = (M.row(i) - A.row(i) * alphaa - C2.row(i) * alphac) *
//         (M.row(i) - A.row(i) * alphaa - C2.row(i) * alphac).t();
//       resg += resv(0);
//     }
//     
//     sigmasqg = rand_igamma((p * n / 2) + h2, (resg / 2) + l2);
//     
//     if(draw >= burnin) {
//       samples_sigmasqg(draw - burnin, 0) = sigmasqg;
//     }
//     
//     // sample betac
//     double mu = 0;
//     double sd = 0;
//     for(int w = 0; w < q1; w++) {
//       arma::vec Cw = C1.col(w);
//       double betacw = betac(w);
//       mu = accu(Cw % (Y - A * betaa - M * betam - C1 * betac + Cw * betac(w))) / accu(Cw % Cw);
//       sd = sqrt(sigmasqe / accu(Cw % Cw));
//       betac(w) = R::rnorm(mu, sd);
//       
//       if(draw >= burnin) {
//         samples_betac(draw - burnin, w) = betac(w);
//       }
//     }
//     
//     // sample alphac
//     mu = 0;
//     sd = 0;
//     for(int w = 0; w < q2; w++) {
//       arma::vec Cw = C2.col(w);
//       for(int j = 0; j < p; j++) {
//         mu = accu(Cw % (M.col(j) - A * alphaa(j) - C2 * alphac + Cw * alphac.col(j))) / accu(Cw % Cw);
//         sd = sqrt(sigmasqg / accu(Cw % Cw));
//         alphac(w, j) = R::rnorm(mu, sd);
//         
//         if(draw >= burnin) {
//           samples_alphac(draw - burnin, w * p + j) = alphac(w, j);
//         }
//       }
//     }
//     
//     draw += 1;
//   }
//   
// }

