# # Test: do unmeasured logistic risk vars => attenuation of the measured var coefs?
# n = 1e5
# coef.risk = rnorm(n)
# coef.unmeas = rnorm(n)

# mix = c(0, 2^seq(-4, 4, 0.5))

# truerisk = lapply(mix, function(m) coef.risk + m*coef.unmeas)
# obs = lapply(truerisk, function(tr) 1/(1+exp(-tr)) > runif(n))
# coefs = sapply(obs, function(obsi) coef(glm(obsi ~ coef.risk, family = binomial()))[2])

# plot(coefs ~ mix, log = "x")


sim_depletion = function(prevalence, sd)
{
    # Simulate mean shift due to depletion. Consider risk R
    # that translates to disease state D.  R ~ N(mu, sd), 
    # E[D] = prevalence.  Estimates and returns the unselected
    # expected risk E[R], and the expected risk in the disease-
    # free population, E_notD[R].
    # 
    # I can't find an analytical solution to this, so go numerical.
    # Suspect there is no such solution known, it's a logit-normal
    # distrib.

    pr_r = function(r, mu) dnorm(r, mean=mu, sd=sd)
    pr_aff_r = function(r) 1/(1+exp(-r))
    pr_r_notaff = function(r, mu) (1-pr_aff_r(r))*pr_r(r, mu) / (1-prevalence)      # Bayes rule

    # Solve for mu to give E[D] = prevalence.
    start_mu = log(prevalence) - log(1-prevalence)
    mu = uniroot(function(mu_test) integrate(function(r) pr_r(r, mu_test)*pr_aff_r(r), mu_test-sd*10, mu_test+sd*10)$value - prevalence, c(-start_mu-sd*10, start_mu+sd*10))$root

    # Then estimate E_notD[R]
    E_notD_R = integrate(function(r) r*pr_r_notaff(r, mu), mu-sd*10, mu+sd*10)$value

    plot.x = seq(mu-sd*10, mu+sd*10, length.out = 1000)

    list(expectations = list(E_R = mu, E_notD_R = E_notD_R), plots = list(d_R = cbind(x = plot.x, y = pr_r(plot.x, mu)), d_notD_R = cbind(x = plot.x, y = pr_r_notaff(plot.x, mu))))
}


sim_depletion.check = function(prevalence, sd, B = 1e5)
{
    draw_x = function(mu) rnorm(B, mean=mu, sd=sd)
    draw_disease_x = function(x) 1/(1+exp(-x)) > runif(length(x))
    frac_diseased = function(mu) mean(draw_disease_x(draw_x(mu)))

    start_mu = log(prevalence) - log(1-prevalence)
    test_mu = seq(start_mu - sd*5, start_mu + sd*5, length.out = 100)
    test_prev = sapply(test_mu, function(mu) frac_diseased(mu))
    plot(test_prev ~ test_mu, type = "o")
    abline(h = prevalence)
    mu = test_mu[which.min(abs(test_prev - prevalence))]
    abline(v = mu)

    x = draw_x(mu)
    diseased = draw_disease_x(x)

    list(expectations = list(E_R = mean(x), E_notD_R = mean(x[!diseased])), plots = list(d_R = density(x), d_notD_R = density(x[!diseased])))
}


# After filtering, 
# ProstateCancer:Hoffmann2:10.1158/2159-8290.CD-15-0315      11:68994667:A:G      0.19062036      0.52
# ProstateCancer:Hoffmann2:10.1158/2159-8290.CD-15-0315      12:53273904:G:A      0.152711877     0.15
# ProstateCancer:Hoffmann2:10.1158/2159-8290.CD-15-0315      17:36074979:G:A     -0.135395102     0.19
# ProstateCancer:Hoffmann2:10.1158/2159-8290.CD-15-0315      17:36101156:T:C      0.173953307     0.60
# ProstateCancer:Hoffmann2:10.1158/2159-8290.CD-15-0315      17:69108753:G:T     -0.122217633     0.52
# ProstateCancer:Hoffmann2:10.1158/2159-8290.CD-15-0315      2:173311553:A:G     -0.223143551     0.06
# ProstateCancer:Hoffmann2:10.1158/2159-8290.CD-15-0315      22:43500212:G:T     -0.131028262     0.49
# ProstateCancer:Hoffmann2:10.1158/2159-8290.CD-15-0315      4:106061534:C:A     -0.113328685     0.41
# ProstateCancer:Hoffmann2:10.1158/2159-8290.CD-15-0315      6:117210052:T:C     -0.122217633     0.29
# ProstateCancer:Hoffmann2:10.1158/2159-8290.CD-15-0315      6:160833664:C:T      0.09531018      0.28
# ProstateCancer:Hoffmann2:10.1158/2159-8290.CD-15-0315      6:41536427:C:T       0.09531018      0.27
# ProstateCancer:Hoffmann2:10.1158/2159-8290.CD-15-0315      7:27976563:G:A      -0.157003749     0.23
# ProstateCancer:Hoffmann2:10.1158/2159-8290.CD-15-0315      8:128011937:T:C     -0.157003749     0.29
# ProstateCancer:Hoffmann2:10.1158/2159-8290.CD-15-0315      8:128093297:C:T      0.246860078     0.21
# ProstateCancer:Hoffmann2:10.1158/2159-8290.CD-15-0315      8:128095156:A:G     -0.09531018      0.29
# ProstateCancer:Hoffmann2:10.1158/2159-8290.CD-15-0315      8:128106880:A:C      0.476234179     0.035
# ProstateCancer:Hoffmann2:10.1158/2159-8290.CD-15-0315      8:128413305:G:T     -0.223143551     0.49
# ProstateCancer:Hoffmann2:10.1158/2159-8290.CD-15-0315      19:51364623:A:G      0.182321557     0.85

coefs = c(0.19062036 , 0.152711877, -0.135395102, 0.173953307, -0.122217633, -0.223143551, -0.131028262, -0.113328685, -0.122217633, 0.09531018 , 0.09531018 , -0.157003749, -0.157003749, 0.246860078, -0.09531018 , 0.476234179, -0.223143551, 0.182321557)
vafs = c(0.52, 0.15, 0.19, 0.60, 0.52, 0.06, 0.49, 0.41, 0.29, 0.28, 0.27, 0.23, 0.29, 0.21, 0.29, 0.035, 0.49, 0.85)
sim.mean = mean(coefs*2*vafs)
sim.sd = sqrt(mean(coefs^2*2*vafs*(1-vafs)))
# sim.sd = 10
sim.prevalence = 0.1        # Prevalence by age 70 of ~ 10%

simdep = sim_depletion(sim.prevalence, sim.sd)
simdep.check = sim_depletion.check(sim.prevalence, sim.sd)

simdep$expectations
# $E_R
# [1] -2.20082
# 
# $E_notD_R
# [1] -2.20172


simdep$expectations$E_R - simdep$expectations$E_notD_R
# [1] 0.0008997179


par(mfrow = c(2, 1))
plot(simdep$plots$d_R, type = "l", col = "red", xlim = c(-100, 100), main = "Quadrature")
lines(simdep$plots$d_notD_R, col = "blue")
plot(simdep.check$plots$d_R, type = "l", col = "red", xlim = c(-100, 100), main = "Monte Carlo")
lines(simdep.check$plots$d_notD_R, col = "blue")
par(mfrow = c(1,1))




# Simulate the full Hoffmann score.  Observe OR of 6 in 10th decile vs 1 in 1st decile
# for full 105-locus score (Fig 4 in Hoffmann).
# Given popn rate of 10%, what does this translate to in terms of normal distribution params of x?
prev_for_mu_sd = function(mu, sd)
{
    integrate(function(x) dnorm(x, mean=mu, sd=sd)*1/(1+exp(-x)), qnorm(0.0001, mean=mu, sd=sd), qnorm(0.9999, mean=mu, sd=sd))$value
}

or_for_mu_sd = function(mu, sd)
{
    pr_diseased_lt10 = integrate(function(x) dnorm(x, mean=mu, sd=sd)*1/(1+exp(-x)), qnorm(0.0001, mean=mu, sd=sd), qnorm(0.1, mean=mu, sd=sd))$value / 0.1
    pr_diseased_gt90 = integrate(function(x) dnorm(x, mean=mu, sd=sd)*1/(1+exp(-x)), qnorm(0.9, mean=mu, sd=sd), qnorm(0.9999, mean=mu, sd=sd))$value / 0.1
    odds_diseased_lt10 = pr_diseased_lt10 / (1-pr_diseased_lt10)
    odds_diseased_gt90 = pr_diseased_gt90 / (1-pr_diseased_gt90)
    odds_diseased_gt90 / odds_diseased_lt10
}

musd = optim(c(-1, 0.1), function(musd) abs(prev_for_mu_sd(musd[1], musd[2]) - 0.1) + abs(or_for_mu_sd(musd[1], musd[2]) - 6))$par

musd
prev_for_mu_sd(musd[1], musd[2])
or_for_mu_sd(musd[1], musd[2])

simdep = sim_depletion(0.1, musd[2])
simdep$expectations
simdep$expectations$E_R - simdep$expectations$E_notD_R
plot(simdep$plots$d_R, type = "l", col = "red", xlim = c(-4, 0), main = "Hoffmann risk score distribution", xlab = "Score", ylab = "Density")
lines(simdep$plots$d_notD_R, col = "blue")
abline(v = simdep$expectations$E_R, col = "red", lty = "dotted")
abline(v = simdep$expectations$E_notD_R, col = "blue", lty = "dotted")
legend("topright", inset = 0.05, col = c("red", "blue"), legend = c("All men", "PCa-free men"), lty = "solid")
