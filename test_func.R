n_obs <- 1000
X <- rnorm(n_obs)
ps <- pnorm(0.5 * X)
TA <- rbinom(n_obs, 1, ps)

gpbalancer::bal_stats(X, TA==1, 'continuous')
gpbalancer::bal_table(as.data.frame(X), col_ind = 1, treat_ind = TA==1)
gpbalancer::ks_avg_test(X, TA)
outro_par <- gpbalancer::gpbal_fixed(TA, sqexp_par(as.matrix(X), theta = c(1,0.25)), ep_vers = 'p')
plot(ps, outro_par$ps)
abline(0,1)

outro_seq <- gpbalancer::gpbal_fixed(TA, sqexp_par(as.matrix(X), theta = c(1,0.25)), ep_vers = 'seq')

all.equal(outro_par, outro_seq)

outro_opt <- gpbalancer::gpbal(as.matrix(X), TA, sqexp_par, init_theta = c(1), verbose=T)
plot(ps, outro_opt$ps)
abline(0,1)
