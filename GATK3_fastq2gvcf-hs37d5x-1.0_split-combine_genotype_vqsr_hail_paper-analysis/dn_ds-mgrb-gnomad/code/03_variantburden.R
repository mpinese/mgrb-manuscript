library(boot)
library(ggplot2)
library(HDInterval)
library(modeest)
library(plyr)


VariantBurden = function(symbol, class, ac1, an1, aa1, ac2, an2, aa2, meta = NULL)
{
    # All parameters except meta are equal-length vectors.
    # If meta is not null, nrow(meta) must equal the length of all vectors.
    # Each entry in the vectors / row in meta corresponds to a single variant allele.
    # symbol: Character vector of gene symbols
    # class: Character vector of variant classes
    # ac1: Integer vector of variant allele count in cohort 1
    # an1: Integer vector of variant locus total allele count in cohort 1
    # aa1: Integer vector of variant allele homozygote count in cohort 1
    # ac2: Integer vector of variant allele count in cohort 2
    # an2: Integer vector of variant locus total allele count in cohort 2
    # aa2: Integer vector of variant allele homozygote count in cohort 2
    # meta: data frame of metadata to keep aligned with variants
    result = list()

    stopifnot(length(symbol) == length(class))
    stopifnot(length(symbol) == length(ac1))
    stopifnot(length(symbol) == length(an1))
    stopifnot(length(symbol) == length(aa1))
    stopifnot(length(symbol) == length(ac2))
    stopifnot(length(symbol) == length(an2))
    stopifnot(length(symbol) == length(aa2))
    stopifnot(is.null(meta) || length(symbol) == nrow(meta))

    result$meta = meta
    result$data = data.frame(
        symbol = symbol,
        class = class,
        ac1 = ac1,
        an1 = an1,
        aa1 = aa1,
        ac2 = ac2,
        an2 = an2,
        aa2 = aa2)
    result$permuted = FALSE

    class(result) = "VariantBurden"
    result
}

is.VariantBurden = function(x) inherits(x, "VariantBurden")
length.VariantBurden = function(x) nrow(x$data)
`[.VariantBurden` = function(x, i) {
    newvb = VariantBurden(
        x$data$symbol[i], x$data$class[i], 
        x$data$ac1[i], x$data$an1[i], x$data$aa1[i],
        x$data$ac2[i], x$data$an2[i], x$data$aa2[i],
        x$meta[i,,drop=FALSE])
    newvb$permuted = x$permuted
    newvb }
print.VariantBurden = function(x, ...) cat(format(x, ...), "\n")
format.VariantBurden = function(x, ...) {
    str = "VariantBurden object"
    if (x$permuted == TRUE) { str = paste(str, "(PERMUTED)") }
    str = paste(str, sprintf("\n  %d variants, %d classes, %d genes", nrow(x$data), length(unique(x$data$class)), length(unique(x$data$symbol))))
    str
}


permute.VariantBurden = function(vb)
{
    # Permute each row (ie variant) between cohorts 1 and 2.
    # Assume diploidy.  In cases where a variant is missing
    # from a cohort (ie AC == NA), impute AC by sampling from
    # the non-missing ACs of that cohort.

    an1 = vb$data$an1
    ac1 = vb$data$ac1
    aa1 = vb$data$aa1
    an2 = vb$data$an2
    ac2 = vb$data$ac2
    aa2 = vb$data$aa2

    # Each input is a numeric vector.  All vectors are of equal length.
    # an1: vector of AN for cohort 1
    # ac1: vector of AC for cohort 1
    # aa1: vector of genotype AA (alt/alt) for cohort 1
    # etc for cohort 2

    # Return the permuted values as list(an1p, ac1p, aa1p, an2p, ac2p, aa2p)

    # Impute AN for absent variants
    # Reference non-missing ANs from which to sample:
    an1_nonmissing = an1[!is.na(an1)]
    an2_nonmissing = an2[!is.na(an2)]
    # Impute:
    ac1[is.na(ac1)] = 0
    an1[is.na(an1)] = sample(an1_nonmissing, sum(is.na(an1)), replace = TRUE)
    aa1[is.na(aa1)] = 0
    ac2[is.na(ac2)] = 0
    an2[is.na(an2)] = sample(an2_nonmissing, sum(is.na(an2)), replace = TRUE)
    aa2[is.na(aa2)] = 0

    # Permute
    # Compute sample numbers per cohort (assume diploid)
    ns1 = an1/2
    ns2 = an2/2

    # Compute genotype counts per cohort
    ra1 = ac1 - 2*aa1
    rr1 = ns1 - aa1 - ra1
    ra2 = ac2 - 2*aa2
    rr2 = ns2 - aa2 - ra2

    # Combine genotypes across cohorts to give total genotype counts
    rr = rr1 + rr2
    ra = ra1 + ra2
    aa = aa1 + aa2
    ns = ns1 + ns2
    stopifnot(ns == rr + ra + aa)

    # Redistribute non-reference genotypes to the two cohorts
    # Two step selection process:
    #   1. Choose RA heterozygotes vs RR/AA homozygotes.
    #      White balls in urn:     RA hets (ra)
    #      Black balls in urn:     RR/RA homs (rr+aa)
    #      Number of balls drawn:  NS (ns)
    #      Returned value:         number of white balls (RA permuted)
    #   2. Choose AA homozygotes from RR/AA set
    #      (ie contingent on not being RA from step 1)
    #      White balls in urn:     AA (aa)
    #      Black balls in urn:     RR (rr)
    #      Number of balls drawn:  NS-RAp = ns-rap
    #      Returned value:         number of white balls (AA permuted)
    ra1p = rhyper(length(ns1), ra, rr+aa, ns1)
    aa1p = rhyper(length(ns1), aa, rr, ns1-ra1p)
    ra2p = ra - ra1p
    aa2p = aa - aa1p

    # Allocate reference genotypes to make up the correct totals
    rr1p = ns1 - ra1p - aa1p
    rr2p = ns2 - ra2p - aa2p

    # Convert back to AC and AN
    an1p = 2*(rr1p + ra1p + aa1p)
    an2p = 2*(rr2p + ra2p + aa2p)
    ac1p = ra1p + 2*aa1p
    ac2p = ra2p + 2*aa2p

    # Consistency checks
    stopifnot(an1p == an1)
    stopifnot(an2p == an2)
    stopifnot(ac1p + ac2p == ac1 + ac2)
    stopifnot(aa1p + aa2p == aa1 + aa2)
    stopifnot(all(rr1p >= 0))
    stopifnot(all(ra1p >= 0))
    stopifnot(all(aa1p >= 0))
    stopifnot(all(rr2p >= 0))
    stopifnot(all(ra2p >= 0))
    stopifnot(all(aa2p >= 0))

    vbp = VariantBurden(vb$data$symbol, vb$data$class, ac1p, an1p, aa1p, ac2p, an2p, aa2p, vb$meta)
    vbp$permuted = TRUE

    vbp
}


statistic_ns_tc.VariantBurden = function(vb, syn_regex, nonsyn_regex, control_symbols, test_symbols)
{
    # syn_regex        Regular expression matching vb$class values to be treated as synonymous
    # nonsyn_regex     Regular expression matching vb$class values to be treated as nonsynonymous
    # control_symbols  Vector of symbols to form the control set
    # test_symbols     Vector of symbols to form the test set
    # Variants for which class matches neither syn_regex or nonsyn_regex, or neither control or test
    # symbols, are ignored.

    nonsyn = grepl(nonsyn_regex, vb$data$class)
    syn = grepl(syn_regex, vb$data$class)
    stopifnot(!(nonsyn & syn))
    stopifnot(any(nonsyn) == TRUE)
    stopifnot(any(syn) == TRUE)

    test = vb$data$symbol %in% test_symbols
    control = vb$data$symbol %in% control_symbols
    stopifnot(!(test & control))
    stopifnot(any(test) == TRUE)
    stopifnot(any(control) == TRUE)

    keep = (nonsyn | syn) & (test | control)
    rm(syn, control)
    af1 = vb$data$ac1[keep] / vb$data$an1[keep]
    af2 = vb$data$ac2[keep] / vb$data$an2[keep]
    nonsyn = nonsyn[keep]
    test = test[keep]

    af1[is.na(af1)] = 0
    af2[is.na(af2)] = 0

    # Sum of NS AFs in test genes in cohort 1
    nt1 = sum(af1[nonsyn & test])
    # Sum of S AFs in test genes in cohort 1
    st1 = sum(af1[!nonsyn & test])
    # Sum of NS AFs in control genes in cohort 1
    nc1 = sum(af1[nonsyn & !test])
    # Sum of S AFs in control genes in cohort 1
    sc1 = sum(af1[!nonsyn & !test])
    stopifnot(abs(nt1 + st1 + nc1 + sc1 - sum(af1)) < 1e-4)

    # Sum of NS AFs in test genes in cohort 2
    nt2 = sum(af2[nonsyn & test])
    # Sum of S AFs in test genes in cohort 2
    st2 = sum(af2[!nonsyn & test])
    # Sum of NS AFs in control genes in cohort 2
    nc2 = sum(af2[nonsyn & !test])
    # Sum of S AFs in control genes in cohort 2
    sc2 = sum(af2[!nonsyn & !test])
    stopifnot(abs(nt2 + st2 + nc2 + sc2 - sum(af2)) < 1e-4)

    # Cohort 1 dN/dS in test genes
    rt1 = nt1 / st1
    # Cohort 1 dN/dS in control genes
    rc1 = nc1 / sc1
    # Cohort 2 dN/dS in test genes
    rt2 = nt2 / st2
    # Cohort 2 dN/dS in control genes
    rc2 = nc2 / sc2

    # Cohort 1 dN/dS in test genes relative to dN/dS in control genes
    r1 = rt1 / rc1
    # Cohort 2 dN/dS in test genes relative to dN/dS in control genes
    r2 = rt2 / rc2

    # Final statistic:dN/dS in test genes relative to control genes, Cohort 1 relative to Cohort 2
    r1 / r2
}


statistic_ns_t.VariantBurden = function(vb, syn_regex, nonsyn_regex, test_symbols)
{
    # syn_regex        Regular expression matching vb$class values to be treated as synonymous
    # nonsyn_regex     Regular expression matching vb$class values to be treated as nonsynonymous
    # test_symbols     Vector of symbols to form the test set
    # Variants for which class matches neither syn_regex or nonsyn_regex, or neither control or test
    # symbols, are ignored.

    nonsyn = grepl(nonsyn_regex, vb$data$class)
    syn = grepl(syn_regex, vb$data$class)
    stopifnot(!(nonsyn & syn))
    stopifnot(any(nonsyn) == TRUE)
    stopifnot(any(syn) == TRUE)

    test = vb$data$symbol %in% test_symbols
    stopifnot(any(test) == TRUE)

    keep = (nonsyn | syn) & test
    rm(syn)
    af1 = vb$data$ac1[keep] / vb$data$an1[keep]
    af2 = vb$data$ac2[keep] / vb$data$an2[keep]
    nonsyn = nonsyn[keep]
    test = test[keep]

    af1[is.na(af1)] = 0
    af2[is.na(af2)] = 0

    # Sum of NS AFs in test genes in cohort 1
    nt1 = sum(af1[nonsyn & test])
    # Sum of S AFs in test genes in cohort 1
    st1 = sum(af1[!nonsyn & test])
    stopifnot(abs(nt1 + st1 - sum(af1)) < 1e-4)

    # Sum of NS AFs in test genes in cohort 2
    nt2 = sum(af2[nonsyn & test])
    # Sum of S AFs in test genes in cohort 2
    st2 = sum(af2[!nonsyn & test])
    stopifnot(abs(nt2 + st2 - sum(af2)) < 1e-4)

    # Cohort 1 dN/dS in test genes
    rt1 = nt1 / st1
    # Cohort 2 dN/dS in test genes
    rt2 = nt2 / st2

    # Final statistic:dN/dS in test genes, Cohort 1 relative to Cohort 2
    rt1 / rt2
}


test_ns_tc.VariantBurden = function(vb, test_symbols, control_symbols, syn_regex = "^Synonymous$", nonsyn_regex = "^Nonsynonymous$", B = 1000, seed = NULL, ...)
{
    if (!is.null(seed))
    {
        oldseed = .Random.seed
        set.seed(seed)
    }

    vb2 = vb[vb$data$symbol %in% c(control_symbols, test_symbols)]
    ts = statistic_ns_tc.VariantBurden(vb2, syn_regex, nonsyn_regex, control_symbols, test_symbols)
    tsp = laply(1:B, function(i) statistic_ns_tc.VariantBurden(permute.VariantBurden(vb2), syn_regex, nonsyn_regex, control_symbols, test_symbols), ...)

    next_left = sum(tsp <= ts)
    next_right = sum(tsp >= ts)
    next_min = min(next_left, next_right)
    pval = min(1, 2*(next_min+1)/(B+1))

    if (!is.null(seed))
        .Random.seed = oldseed

    list(statistic = ts, null = tsp, p = pval)
}


test_ns_t.VariantBurden = function(vb, test_symbols, syn_regex = "^Synonymous$", nonsyn_regex = "^Nonsynonymous$", B = 1000, seed = NULL, ...)
{
    if (!is.null(seed))
    {
        oldseed = .Random.seed
        set.seed(seed)
    }

    vb2 = vb[vb$data$symbol %in% test_symbols]
    ts = statistic_ns_t.VariantBurden(vb2, syn_regex, nonsyn_regex, test_symbols)
    tsp = laply(1:B, function(i) statistic_ns_t.VariantBurden(permute.VariantBurden(vb2), syn_regex, nonsyn_regex, test_symbols), ...)

    next_left = sum(tsp <= ts)
    next_right = sum(tsp >= ts)
    next_min = min(next_left, next_right)
    pval = min(1, 2*(next_min+1)/(B+1))

    if (!is.null(seed))
        .Random.seed = oldseed

    list(statistic = ts, null = tsp, p = pval)
}


boot.ns.rangen.single_cohort = function(an, ac, aa)
{
    # Bootstrap resample genotypes and variants in a single cohort.
    # Assume diploid (ie NS = AN/2).
    # All parameters are equal-length vectors:
    #   an  Number of alleles sequenced
    #   ac  Count of alternate alleles
    #   aa  Count of homozygous alternate individuals (genotype AA)
    ns = an/2
    sel = !is.na(an)
    nv = sum(sel)
    anv = an[sel]
    acv = ac[sel]
    aav = aa[sel]
    nsv = ns[sel]
    rav = acv - 2*aav
    rrv = nsv - rav - aav

    # Multinomial sample
    rrp = rep(NA, length(ns))
    rap = rep(NA, length(ns))
    aap = rep(NA, length(ns))
    # Use a sequential binomial construction for speed: the 
    # obvious multinomial (code commented out below) is takes
    # ~ 10 X as long due to lack of C vectorization and GC
    # overhead.
    rbinom2 = function(n, size, prob)
    {
        result = rep(0, n)
        nz = size > 0
        result[nz] = rbinom(sum(nz), size[nz], prob[nz])
        result
    }
    aap[sel] = rbinom(nv, nsv, aav/nsv)
    rap[sel] = rbinom2(nv, nsv-aap[sel], rav/(rav+rrv))
    rrp[sel] = nsv - aap[sel] - rap[sel]
    # samp = mapply(function(rrvi, ravi, aavi) rmultinom(1, rrvi+ravi+aavi, c(rrvi, ravi, aavi)), rrv, rav, aav)
    # rrp[sel] = samp[1,]
    # rap[sel] = samp[2,]
    # aap[sel] = samp[3,]
    nsp = rap+aap+rrp

    # Regenerate counts
    anp = 2*nsp
    acp = 2*aap + rap

    # Consistency checks
    stopifnot(nsp[sel] == ns[sel])
    stopifnot(rrp[sel] >= 0)
    stopifnot(anp[sel] == an[sel])

    # Resample variants
    # varsel = sample.int(length(anp), replace = TRUE)
    # anp = anp[varsel]
    # acp = acp[varsel]
    # aap = aap[varsel]

    list(an = anp, ac = acp, aa = aap)
}


boot.ns.rangen = function(vb, mle)
{
    # Resample each variant in vb, within cohorts.
    # By necessity variants must be split, so this is not true
    # sample-level resampling, but rather a marginal approximation.

    perm1 = boot.ns.rangen.single_cohort(vb$data$an1, vb$data$ac1, vb$data$aa1)
    perm2 = boot.ns.rangen.single_cohort(vb$data$an2, vb$data$ac2, vb$data$aa2)

    vbp = VariantBurden(vb$data$symbol, vb$data$class, perm1$ac, perm1$an, perm1$aa, perm2$ac, perm2$an, perm2$aa, vb$meta)
    vbp$permuted = TRUE

    vbp
}


boot.ns.statistic = function(vb, test_symbols, control_symbols, syn_regex = "^Synonymous$", nonsyn_regex = "^Nonsynonymous$")
{
    nonsyn = grepl(nonsyn_regex, vb$data$class)
    syn = grepl(syn_regex, vb$data$class)
    stopifnot(!(nonsyn & syn))
    stopifnot(any(nonsyn) == TRUE)
    stopifnot(any(syn) == TRUE)

    test = vb$data$symbol %in% test_symbols
    control = vb$data$symbol %in% control_symbols
    stopifnot(!(test & control))
    stopifnot(any(test) == TRUE)
    stopifnot(any(control) == TRUE)

    keep = (nonsyn | syn) & (test | control)
    rm(syn, control)
    af1 = vb$data$ac1[keep] / vb$data$an1[keep]
    af2 = vb$data$ac2[keep] / vb$data$an2[keep]
    nonsyn = nonsyn[keep]
    test = test[keep]

    af1[is.na(af1)] = 0
    af2[is.na(af2)] = 0

    # Sum of NS AFs in test genes in cohort 1
    nt1 = sum(af1[nonsyn & test])
    # Sum of S AFs in test genes in cohort 1
    st1 = sum(af1[!nonsyn & test])
    # Sum of NS AFs in control genes in cohort 1
    nc1 = sum(af1[nonsyn & !test])
    # Sum of S AFs in control genes in cohort 1
    sc1 = sum(af1[!nonsyn & !test])
    stopifnot(abs(nt1 + st1 + nc1 + sc1 - sum(af1)) < 1e-4)

    # Sum of NS AFs in test genes in cohort 2
    nt2 = sum(af2[nonsyn & test])
    # Sum of S AFs in test genes in cohort 2
    st2 = sum(af2[!nonsyn & test])
    # Sum of NS AFs in control genes in cohort 2
    nc2 = sum(af2[nonsyn & !test])
    # Sum of S AFs in control genes in cohort 2
    sc2 = sum(af2[!nonsyn & !test])
    stopifnot(abs(nt2 + st2 + nc2 + sc2 - sum(af2)) < 1e-4)

    # Cohort 1 dN/dS in test genes
    rt1 = nt1 / st1
    # Cohort 1 dN/dS in control genes
    rc1 = nc1 / sc1
    # Cohort 2 dN/dS in test genes
    rt2 = nt2 / st2
    # Cohort 2 dN/dS in control genes
    rc2 = nc2 / sc2

    # Cohort 1 dN_test/dN_control
    rn1 = nt1 / nc1
    # Cohort 2 dN_test/dN_control
    rn2 = nt2 / nc2

    c(
        C1TN = nt1,
        C1TS = st1,
        C1CN = nc1,
        C1CS = sc1,

        C2TN = nt2,
        C2TS = st2,
        C2CN = nc2,
        C2CS = sc2,

        `C1TN/C1TS` = rt1,
        `C1CN/C1CS` = rc1,
        `C2TN/C2TS` = rt2,
        `C2CN/C2CS` = rc2,

        `C1TN/C1CN` = rn1,
        `C2TN/C2CN` = rn2,

        `(C1TN/C1TS) / (C2TN/C2TS)` = rt1 / rt2,
        `((C1TN/C1TS)/(C1CN/C1CS)) / ((C2TN/C2TS)/(C2CN/C2CS))` = (rt1 / rc1) / (rt2 / rc2),

        `(C1TN/C1TS) - (C2TN/C2TS)` = rt1 - rt2
    )
}


boot.VariantBurden = function(vb, test_symbols, control_symbols, syn_regex = "^Synonymous$", nonsyn_regex = "^Nonsynonymous$", seed = NULL, ...)
{
    stopifnot(require(boot))

    if (!is.null(seed))
    {
        oldseed = .Random.seed
        set.seed(seed)
    }

    result = boot(data = vb[vb$data$symbol %in% c(control_symbols, test_symbols)], statistic = boot.ns.statistic, sim = "parametric", 
        ran.gen = boot.ns.rangen, syn_regex = syn_regex, nonsyn_regex = nonsyn_regex, control_symbols = control_symbols, test_symbols = test_symbols, ...)

    if (!is.null(seed))
        .Random.seed = oldseed

    result
}


plotboot.VariantBurden = function(boot, vars)
{
    stopifnot(require(ggplot2))
    stopifnot(require(HDInterval))
    stopifnot(require(modeest))
    sel = match(vars, names(boot$t0))
    ts = boot$t[,sel,drop=FALSE]
    colnames(ts) = names(boot$t0)[sel]
    n = ncol(ts)
    hdis = apply(ts, 2, hdi)
    modes = apply(ts, 2, function(x) mlv(x, method = "parzen")$M)
    df = data.frame(variable = colnames(ts), mode = modes, lcl = hdis[1,], ucl = hdis[2,])
    ggplot(df, aes(x = variable, y = mode, ymin = lcl, ymax = ucl)) + geom_bar(stat = "identity") + geom_errorbar(width = 0.1)
}


simulate.VariantBurden = function(nvars = 500, nsamp1 = 3000, nsamp2 = 5000, rseq = 0.98, rfixed = 0.05, rvar = 0.0001)
{
    nsamp = nsamp1+nsamp2
    isfixed = runif(nvars) < rfixed
    rvar = rep(rvar, nvars)
    rvar[isfixed] = 1-rvar[isfixed]

    rbinom2 = function(n, size, prob)
    {
        result = rep(0, n)
        nz = size > 0
        result[nz] = rbinom(sum(nz), size[nz], prob[nz])
        result
    }

    ns1 = rbinom(nvars, nsamp1, rseq)
    aa1 = rbinom(nvars, ns1, rvar^2)
    rr1 = rbinom2(nvars, ns1-aa1, (1-rvar)^2 / (1-rvar^2))
    ra1 = ns1 - aa1 - rr1
    an1 = ns1*2
    ac1 = 2*aa1 + ra1

    ns2 = rbinom(nvars, nsamp2, rseq)
    aa2 = rbinom(nvars, ns2, rvar^2)
    rr2 = rbinom2(nvars, ns2-aa2, (1-rvar)^2 / (1-rvar^2))
    ra2 = ns2 - aa2 - rr2
    an2 = ns2*2
    ac2 = 2*aa2 + ra2

    an1[ac1 == 0] = NA
    aa1[ac1 == 0] = NA
    ac1[ac1 == 0] = NA
    an2[ac2 == 0] = NA
    aa2[ac2 == 0] = NA
    ac2[ac2 == 0] = NA

    symbol = sample(c("Test", "Control"), nvars, replace = TRUE)
    class = sample(c("Synonymous", "Nonsynonymous"), nvars, replace = TRUE)

    result = VariantBurden(symbol, class, ac1, an1, aa1, ac2, an2, aa2)

    result[!is.na(an1) | !is.na(an2)]
}


#set.seed(314159)
#test.sim.statistic_ns_tc = replicate(10000, statistic_ns_tc.VariantBurden(simulate.VariantBurden(), test_symbols = "Test", control_symbols = "Control", syn_regex = "^Synonymous$", nonsyn_regex = "^Nonsynonymous$"))
#set.seed(314159)
#test.sim.perm.statistic_ns_tc = replicate(10000, statistic_ns_tc.VariantBurden(permute.VariantBurden(simulate.VariantBurden()), test_symbols = "Test", control_symbols = "Control", syn_regex = "^Synonymous$", nonsyn_regex = "^Nonsynonymous$"))
#plot(test.sim.perm.statistic_ns_tc ~ test.sim.statistic_ns_tc)
#ks.test(test.sim.perm.statistic_ns_tc, test.sim.statistic_ns_tc)

