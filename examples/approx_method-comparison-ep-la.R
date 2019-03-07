wtd_group_mean <- function(X, A, wts) {sum(A * X * w) / sum(A * wts)}

nobs <- 100
X <- rnorm(nobs)
ps <- pnorm( 0.5 * X )
A <- rbinom(nobs, 1, ps)
Y1 <- X^2 + 10 + rnorm(nobs)
Y0 <- X + rnorm(nobs)
Yobs <- ifelse(A==1, Y1, Y0)


t1 <- system.time(outro_ep <- gpbal(X, A, cov_function = sqexp_poly, init_theta = c(1,0.5)))
print(outro_ep$thetas)
t1 <- system.time(outro_ep <- gpbal(X, A, cov_function = sqexp_poly))
print(outro_ep$thetas)
t2 <- system.time(outro_la <- gpbal(X, A, cov_function = sqexp_poly, init_theta = c(1,0.5), approx_method = 'laplace'))
par(mfrow=c(1,2))
plot(ps,outro_ep$ps,
     main = paste('Min Bal:', round(outro_ep$minval,3),
                  ', Runtime:', round(t1[3], 2)))
abline(0,1)
plot(ps,outro_la$ps,
     main = paste('Min Bal:', round(outro_la$minval,3),
                  ', Runtime:', round(t2[3], 2)))
abline(0,1)
