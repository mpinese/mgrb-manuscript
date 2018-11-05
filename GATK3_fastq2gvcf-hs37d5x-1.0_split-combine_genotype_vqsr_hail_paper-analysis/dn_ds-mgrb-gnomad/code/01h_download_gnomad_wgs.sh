#!/bin/bash
set -euo pipefail

gsutil -m cp -r gs://gnomad-public/release/2.0.1/vds/genomes/gnomad.genomes.r2.0.1.sites.autosomes.vds ../
