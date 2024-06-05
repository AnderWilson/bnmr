
# simulate some data
set.seed(1)
dat <- simmonodata(n=1000, scenario="linear")

# fit the model
# in practice the number of iterations should be much larger
# you should also consider checking different order polynomials (M)
fit_bnmr <- bnmr(y=dat$y,
                 x=dat$x,
                 niter=1000, 
                 nburn=500, 
                 nthin=2, 
                 M=20, 
                 priors = list(mu0=.5, phi0=.25))

# summarize results
summary(fit_bnmr)

# predict new values on a grid
fit_bnmr_predict <- predict(fit_bnmr, newdata=dat$grideval$x )

# predict the residual
fit_bnmr_predict_deriv <- predict(fit_bnmr, newdata=dat$grideval$x, deriv=1 )

# plot estimated function
plot(fit_bnmr)
