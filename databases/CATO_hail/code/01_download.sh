#!/bin/bash
set -euo pipefail

mkdir -p ../source_data
cd ../source_data
wget -c http://www.mauranolab.org/CATO/dbSNP142.CATO.V1.1.starch
cd -
