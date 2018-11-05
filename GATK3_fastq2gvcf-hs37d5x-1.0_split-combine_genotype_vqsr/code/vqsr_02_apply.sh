#!/bin/bash
set -euo pipefail

GATK_PATH="./bin/GenomeAnalysisTK.jar"
GATK_RAM="240G"
GATK_VQSR_SNP="-resource:hapmap,known=false,training=true,truth=true,prior=15.0 ../../../resources/hs37d5x/gatk_bundle/hapmap_3.3.b37.vcf.gz -resource:omni,known=false,training=true,truth=true,prior=12.0 ../../../resources/hs37d5x/gatk_bundle/1000G_omni2.5.b37.vcf.gz -resource:1000G,known=false,training=true,truth=false,prior=10.0 ../../../resources/hs37d5x/gatk_bundle/1000G_phase1.snps.high_confidence.b37.vcf.gz -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ../../../resources/hs37d5x/gatk_bundle/dbsnp_138.b37.vcf.gz -an QD  -an FS -an MQRankSum -an ReadPosRankSum -mode SNP"
GATK_VQSR_INDEL="-resource:mills,known=true,training=true,truth=true,prior=12.0 ../../../resources/hs37d5x/gatk_bundle/Mills_and_1000G_gold_standard.indels.b37.vcf.gz -an FS -an MQRankSum -an ReadPosRankSum -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0"
GATK_REFERENCE="../../../resources/hs37d5x/reference_genome/hs37d5x.fa"
GATK_THREADS=14
GATK_TMP="./tmp"

mkdir -p "${GATK_TMP}"

INPUT_VCF_PATTERN="../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype/MGRB.phase2.tier12.match.*.vcf.gz"
OUTPUT_DIRECTORY="../"

start_timestamp=$(date +%Y%m%d%H%M)

# Background monitoring of disk and RAM usage
echo time temp $(free -w | head -n 1) | sed -E 's/ +/\t/g' > VQSR.${start_timestamp}.mem.log
while [ 1 ]; do
    sleep 30
    echo -e "$(date -Iseconds)\t$(du ${GATK_TMP} | cut -f 1)\t$(free -w | grep Mem | sed -E 's/Mem: +//; s/ +/\t/g')" >> VQSR.${start_timestamp}.mem.log
done &
monitor_psid=$!

function cleanup {
    kill $monitor_psid
    rm -f .current_chrom
    rm -f "${GATK_TMP}/org.broadinstitute.gatk.*"
}
trap cleanup EXIT

echo "#ApplyRecalibration SNP $(date +%Y%m%d%H%M)" >> VQSR.${start_timestamp}.mem.log
parallel java -Xmx2G -Djava.io.tmpdir=${GATK_TMP} -jar ${GATK_PATH} -T ApplyRecalibration -R ${GATK_REFERENCE} -input {} -mode SNP --ts_filter_level 99.9 -recalFile vqsr.snp.output.recal -tranchesFile vqsr.snp.output.tranches -o ../{/}.vqsr_snp.vcf.gz -XL NC_007605 -XL hs37d5 ::: ${INPUT_VCF_PATTERN} | tee VQSR.${start_timestamp}.snv.apply.log

echo "#ApplyRecalibration Indel $(date +%Y%m%d%H%M)" >> VQSR.${start_timestamp}.mem.log
parallel java -Xmx2G -Djava.io.tmpdir=${GATK_TMP} -jar ${GATK_PATH} -T ApplyRecalibration -R ${GATK_REFERENCE} -input {} -mode INDEL --ts_filter_level 99.9 -recalFile vqsr.indel.output.recal -tranchesFile vqsr.indel.output.tranches -o ../{/}.vqsr_indel.vcf.gz -XL NC_007605 -XL hs37d5 ::: ../*.vqsr_snp.vcf.gz | tee VQSR.${start_timestamp}.indel.apply.log

echo "#Done $(date +%Y%m%d%H%M)" >> VQSR.${start_timestamp}.mem.log

