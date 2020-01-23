# Genome confidence tiers from the MGRB manuscript

Confidence tiers for WGS genotypes based on Illumina TruSeq Nano sequencing. We defined three tiers:

| Tier   | Interpretation |
|--------|----------------|
| 1      | Very high quality regions. Loci in these regions have excellent quality metrics in MGRB samples, and are also present in the Genome in a Bottle high-confidence regions |
| 2      | High quality regions. Loci in these regions have excellent quality metrics in MGRB samples, however are *not* present in the Genome in a Bottle high-confidence regions |
| 3      | Low quality regions. Loci in these regions are of suspect quality |

In our analyses we have found that most misbehaving loci (eg which have an allele frequency imbalance between batches) fall within Tier 3 regions.

BEDs are on GRCh37 coordinates.

For more information on the definition of the tiers, see https://doi.org/10.1038/s41467-019-14079-0
