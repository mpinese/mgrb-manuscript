#!/bin/bash
set -euo pipefail

./bin/pyhail.sh ./05h_generate_vep_variant_annots.py

./bin/pyhail.sh ./05h_add_db_vep_variant_annots.py
