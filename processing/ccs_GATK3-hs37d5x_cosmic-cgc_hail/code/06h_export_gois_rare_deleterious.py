#!./bin/pyhail.sh
import hail

hc = hail.HailContext(log = 'log/06h_export_gois.log', tmp_dir = 'tmp/hail')

vds = hc.read('../ccs.cosmic_cgc.minrep.locusannot.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.vds')


goi_file = open('tmp/gois.txt', 'rt')
next(goi_file)
target_sym_list = ['"{}"'.format(line.rstrip().split('\t')[2]) for line in goi_file]
print(', '.join(target_sym_list))

vds = (vds
    .filter_variants_table(hail.KeyTable.import_bed('tmp/gois.bed.gz'), keep=True)
    .filter_variants_expr('''
        let
            target_syms = [{target_sym_list}]
        in
            va.alleles.map(aa => aa.annotation.vep.cross_references.symbol).exists(vep_symbol => target_syms.exists(target_symbol => target_symbol == vep_symbol))'''.format(target_sym_list=', '.join(target_sym_list)))
    .split_multi()
    .variant_qc()
    .filter_variants_expr('''
        (
            va.alleles[va.aIndex-1].annotation.vep.consequences.impact_class == "HIGH" ||
            va.alleles[va.aIndex-1].annotation.vep.predictions.ClinVar.exists(clinsig => clinsig == "pathogenic" || clinsig == "likely_pathogenic") ||
            va.alleles[va.aIndex-1].annotation.clinvar.ClinSigSimple == 1 ||
            va.alleles[va.aIndex-1].annotation.predictions.Eigen.EigenPhred >= 20 || 
            va.alleles[va.aIndex-1].annotation.predictions.CATO.score >= 0.5
        ) && (
            orElse(toFloat(va.alleles[va.aIndex-1].annotation.vep.population_frequencies.MAX_AF), toFloat(0)) < 0.05 &&
            orElse(toFloat(va.alleles[va.aIndex-1].annotation.freqs.hrc.AF), toFloat(0)) < 0.05 &&
            orElse(toFloat(va.alleles[va.aIndex-1].annotation.freqs.gnomad.AF_AFR), toFloat(0)) < 0.05 &&
            orElse(toFloat(va.alleles[va.aIndex-1].annotation.freqs.gnomad.AF_AMR), toFloat(0)) < 0.05 &&
            orElse(toFloat(va.alleles[va.aIndex-1].annotation.freqs.gnomad.AF_ASJ), toFloat(0)) < 0.05 &&
            orElse(toFloat(va.alleles[va.aIndex-1].annotation.freqs.gnomad.AF_EAS), toFloat(0)) < 0.05 &&
            orElse(toFloat(va.alleles[va.aIndex-1].annotation.freqs.gnomad.AF_FIN), toFloat(0)) < 0.05 &&
            orElse(toFloat(va.alleles[va.aIndex-1].annotation.freqs.gnomad.AF_NFE), toFloat(0)) < 0.05 &&
            orElse(toFloat(va.alleles[va.aIndex-1].annotation.freqs.gnomad.AF_OTH), toFloat(0)) < 0.05 &&
            va.qc.AF < 0.1
        )''')
)


print(vds.summarize())

(vds
    .export_genotypes('../ccs.cosmic_cgc.minrep.locusannot.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.goi.genotypes.rare.deleterious.txt', expr='''
        SampleID = s,
        Variant = v,
        Variant.filters = va.variant.filters.mkString(","),
        Variant.locus_QC_tier = va.locus.tier,
        WARNING_FAKE_REF = g.fakeRef,
        Genotype = "" + g.gtj() + "/" + g.gtk(),
        Depth = g.dp,
        AlleleFreq1 = g.ad[g.gtj()]/g.dp,
        AlleleFreq2 = g.ad[g.gtk()]/g.dp,
        Symbol = va.alleles[va.aIndex-1].annotation.vep.cross_references.symbol,
        HGVSp = va.alleles[va.aIndex-1].annotation.vep.HGVS.HGVSp,
        VEP.feature = va.alleles[va.aIndex-1].annotation.vep.consequences.feature.id,
        VEP.feature_type = va.alleles[va.aIndex-1].annotation.vep.consequences.feature.type,
        VEP.biotype = va.alleles[va.aIndex-1].annotation.vep.consequences.feature.biotype,
        VEP.canonical = va.alleles[va.aIndex-1].annotation.vep.consequences.feature.canonical,
        VEP.consequences = va.alleles[va.aIndex-1].annotation.vep.consequences.consequences.mkString(","),
        VEP.class = va.alleles[va.aIndex-1].annotation.vep.consequences.impact_class,
        ClinVar.ClinicalSignificance.VEP = va.alleles[va.aIndex-1].annotation.vep.predictions.ClinVar.mkString(","),
        ClinVar.ClinicalSignificance = va.alleles[va.aIndex-1].annotation.clinvar.ClinicalSignificance.mkString(","),
        ClinVar.rsid = va.alleles[va.aIndex-1].annotation.clinvar.rsid,
        ClinVar.ReviewStatus = va.alleles[va.aIndex-1].annotation.clinvar.ReviewStatus,
        AF.GnomAD.global = va.alleles[va.aIndex-1].annotation.freqs.gnomad.AF,
        AF.GnomAD.AFR = va.alleles[va.aIndex-1].annotation.freqs.gnomad.AF_AFR,
        AF.GnomAD.AMR = va.alleles[va.aIndex-1].annotation.freqs.gnomad.AF_AMR,
        AF.GnomAD.ASJ = va.alleles[va.aIndex-1].annotation.freqs.gnomad.AF_ASJ,
        AF.GnomAD.EAS = va.alleles[va.aIndex-1].annotation.freqs.gnomad.AF_EAS,
        AF.GnomAD.FIN = va.alleles[va.aIndex-1].annotation.freqs.gnomad.AF_FIN,
        AF.GnomAD.NFE = va.alleles[va.aIndex-1].annotation.freqs.gnomad.AF_NFE,
        AF.GnomAD.OTH = va.alleles[va.aIndex-1].annotation.freqs.gnomad.AF_OTH,
        AF.VEP.max = va.alleles[va.aIndex-1].annotation.vep.population_frequencies.MAX_AF,
        AF.HRC = va.alleles[va.aIndex-1].annotation.freqs.hrc.AF,
        AF.MGRB.total = va.qc.AF,
        Predictions.CATO.score = va.alleles[va.aIndex-1].annotation.predictions.CATO.score,
        Predictions.CATO.motif = va.alleles[va.aIndex-1].annotation.predictions.CATO.motif,
        Predictions.CATO.celltypes = va.alleles[va.aIndex-1].annotation.predictions.CATO.celltypes.mkString(","),
        Predictions.Eigen.EigenRaw = va.alleles[va.aIndex-1].annotation.predictions.Eigen.EigenRaw,
        Predictions.Eigen.EigenPhred = va.alleles[va.aIndex-1].annotation.predictions.Eigen.EigenPhred,
        Predictions.PolyPhen.VEP = va.alleles[va.aIndex-1].annotation.vep.predictions.PolyPhen,
        Predictions.SIFT.VEP = va.alleles[va.aIndex-1].annotation.vep.predictions.SIFT''')
)
