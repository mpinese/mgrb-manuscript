# Code supporting the MGRB manuscript

This repo contains partial code supporting the paper _"The Medical Genome Reference Bank contains whole genome and phenotype data of 2570 healthy elderly"_ by Pinese et al (https://doi.org/10.1038/s41467-019-14079-0).

The code is heavily customised to its execution platform, and uses restricted data files as input. As such, it won't run out of the box, and this repo serves more as a reference for the analytical steps taken for this paper, as well as a space to store additional result and intermediate files which might be of use to the genomics community.

## Structure

|Directory                    |Contents|
|-----------------------------|--------|
|processing/                  |low-level pipeline & processing code|
|manuscript/                  |analyses specific to the manuscript. This could do with some cleaning up; currently there is work included which was not included in the final form of the paper.|
|extra\_results/              |additional intermediate files which may be of use / interest.|
|extra\_results/genome\_tiers/|contains three BED files defining genome confidence tiers on GRCh37. These tiers were useful in the MGRB work to separate loci which were well-genotyped on the Illumina WGS platform used vs those with more dubious data quality. These intermediate files are released here in the hope they will be useful to other researchers working with Illumina WGS data (especially TruSeq Nano).|

