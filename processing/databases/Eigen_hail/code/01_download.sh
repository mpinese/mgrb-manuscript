#!/bin/bash
set -euo pipefail

mkdir -p ../source_data
cd ../source_data
wget -c https://xioniti01.u.hpc.mssm.edu/v1.1/Eigen_hg19_coding_annot_04092016.tab.bgz
cd -
