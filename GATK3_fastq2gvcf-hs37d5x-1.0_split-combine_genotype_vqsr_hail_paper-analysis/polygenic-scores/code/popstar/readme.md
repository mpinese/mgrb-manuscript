# POPSTAR: Poly/Oligogenic Phenotype Score Testing And Resampling

*POPSTAR* is a suite of tools to compute and manipulate 'polygenic' scores, such as are produced by genotype-phenotype association studies.  *POPSTAR's* core functionality is similar to PLINK's --score command, but in addition *POPSTAR* includes modules to estimate the significance of scores by resampling, and convenience functions to automatically import scores from online GWAS databases.

## Installation

## Usage

*POPSTAR* operates on an allele dosage file (**dosages**) and polygenic trait specification file (**PTS**) to generate polygenic score outputs.

```bash
popstar [options] --dosages=samples.dosages --models=models.pts
```

The polygenic trait specification (PTS) file defines the polygenic scores to be calculated by *POPSTAR*.  The `ebi2pts.R` script in `utils/` is a simple script to automatically generate a PTS file from EBI GWAS database tables, available at https://www.ebi.ac.uk/gwas/docs/file-downloads.

## Formats

### Allele dosage file format (.dosages)

### Polygenic trait specification file format (.pts)

### Output format
