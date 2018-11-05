1. Place gVCFs from fastq2gvcf AWS pipeline in GATK3_fastq2gvcf-hs37d5x-1.0/gvcfs/
2. Queue QC split-combine jobs:
     GATK3_fastq2gvcf-hs37d5x-1.0_qc_split-combine/code/queue_jobs.byshard.bygroup.sh
3. Queue QC genotype jobs:
     GATK3_fastq2gvcf-hs37d5x-1.0_qc_split-combine_genotype/code/queue_jobs.sh
4. Merge QC genotypes:
     GATK3_fastq2gvcf-hs37d5x-1.0_qc_split-combine_genotype/code/merge.sh
5. Calculate sample QC metrics, identifying outliers:
     GATK3_fastq2gvcf-hs37d5x-1.0_qc_split-combine_genotype_sampleqc/code/sampleqc.py
6. Queue final split-combine jobs:
     GATK3_fastq2gvcf-hs37d5x-1.0_split-combine/code/queue_jobs.sh
7. Queue final genotype jobs:
     GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype/code/queue_jobs.sh
   Some final genotype jobs will not complete at NCI (too slow).  For these, 
   copy the relevant g.vcf.gz shards to a local node:
     cd GATK3_fastq2gvcf-hs37d5x-1.0_split-combine
     mkdir -p tier1.match tier1.mismatch tier2.match tier2.mismatch tier3.match tier3.mismatch tier4.all
     parallel -j8 rsync --append -avP mxp569@r-dm.nci.org.au:/g/data3/wq2/results/phase2/hs37d5x/GATK3_fastq2gvcf-hs37d5x-1.0_split-combine/*/*.tier[12].*.shard{}.* . ::: <shard_list>
     mv *.tier1.match.* tier1.match/
     mv *.tier1.mismatch.* tier1.mismatch/
     mv *.tier2.match.* tier2.match/
     mv *.tier2.match.* tier2.mismatch/
     mv *.tier3.match.* tier3.match/
     mv *.tier3.match.* tier3.mismatch/
     mv *.tier4.all.* tier4.all/
     cd ..
   then genotype locally:
     GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype/code/run_locally.sh <shard_list>
   where <shard_list> is a space-separated list of shard IDs to run, for example:
     run_locally.sh 0027 0023 0004 0025 0044 0041 0009 0019 0048 0038 0037
   Ideally order <shard_list> by decreasing g.vcf.gz file size, using:
     cd GATK3_fastq2gvcf-hs37d5x-1.0_split-combine/tier1.match
     for shard in $(ls -1 *.g.vcf.gz | sed -E 's/.*(shard[0-9]+).*/\1/' | sort | uniq); do echo -e ${shard}\\t$(du -ch *.${shard}.* | grep total | cut -f1); done | sort -k2,2rh

