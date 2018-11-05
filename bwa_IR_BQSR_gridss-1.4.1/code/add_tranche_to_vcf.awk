## From the GRIDSS docs, access 3 Mar 2017:
##
## Quality score
##
## GRIDSS calculates quality scores according to the model outlined in [paper].
## As GRIDSS does not yet perform multiple test correction or score recalibration, 
## QUAL scores are vastly overestimated for all variants. As a rule of thumb, 
## variants with QUAL >= 1000 and have assembles from both sides of the breakpoint
## (AS > 0 & RAS > 0) are considered of high quality, variant with QUAL >= 500 but
## can only be assembled from one breakend (AS > 0 | RAS > 0) are considered of 
## intermediate quality, and variants with low QUAL score or lack any supporting 
## assemblies are considered the be of low quality.

BEGIN {
    FS="\t"
    OFS="\t"
}

/^#/ {
    print
    next
}

($7 != ".") {
    next
}

($7 == ".") {
    split($8, info_pairs, ";")
    delete info

    for (info_pair in info_pairs) {
        split(info_pairs[info_pair], info_pair_parts, "=")
        info[info_pair_parts[1]] = info_pair_parts[2]
    }

    if ($6 >= 1000 && info["AS"] > 0 && info["RAS"] > 0) {
        tranche="HIGH"
    } else if ($6 >= 500 && (info["AS"] > 0 || info["RAS"] > 0)) {
        tranche="INTERMEDIATE"
    } else {
        tranche="LOW"
    }

    $8=$8 ";TRANCHE=" tranche

    print
}
