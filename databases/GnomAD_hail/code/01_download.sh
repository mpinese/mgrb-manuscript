#!/bin/bash
set -euo pipefail

#gsutil -m cp -r gs://gnomad-public/release-170228/gnomad.exomes.r2.0.1.sites.vds ../
#gsutil -m cp -r gs://gnomad-public/release-170228/gnomad.genomes.r2.0.1.sites.vds ../

gsutil -m cp -r gs://gnomad-public/release/2.0.1/vds/exomes/gnomad.exomes.r2.0.1.sites.X.vds ../
gsutil -m cp -r gs://gnomad-public/release/2.0.1/vds/exomes/gnomad.exomes.r2.0.1.sites.Y.vds ../
gsutil -m cp -r gs://gnomad-public/release/2.0.1/vds/exomes/gnomad.exomes.r2.0.1.sites.autosomes.vds ../

gsutil -m cp -r gs://gnomad-public/release/2.0.1/vds/genomes/gnomad.genomes.r2.0.1.sites.autosomes.vds ../
gsutil -m cp -r gs://gnomad-public/release/2.0.1/vds/genomes/gnomad.genomes.r2.0.1.sites.X.vds ../
