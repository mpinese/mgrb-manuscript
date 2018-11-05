#!/bin/bash
set -euo pipefail

seed=314159
kmin=2
kmax=12
B=28
N=200
inrds_single_hd_ld=../04_combined_burden.single.hdld.rds
inrds_single_ld=../04_combined_burden.single.ld.rds
inrds_agegrouped=../04_combined_burden.agegrouped.rds
outrds_stem=../04_fact_runs/04_combined_burden.fact

mkdir -p ../04_fact_runs

queue=hugemem
maxqueued=200
cores=28
mem=950G

if [ ! -e "${inrds_single_hd_ld}" ] || [ ! -e "${inrds_single_ld}" ] || [ ! -e "${inrds_agegrouped}" ]; then
  Rscript 04h_combined_factorise_prepdata.R "${inrds_single_hd_ld}" "${inrds_single_ld}" "${inrds_agegrouped}"
fi

exit

num_queued=$((maxqueued+1))

seq 1 $N | while read i; do
  if [ ${num_queued} -gt ${maxqueued} ]; then
    set +e
    num_queued=$(qstat -u `whoami` | grep "${queue:0:8}" | awk '($10 == "Q"){total+=1} END{print total+0}')
    set -e
  fi
  if [ ${num_queued} -gt ${maxqueued} ]; then
    echo -e "\e[31mTerminating: queue full\e[39m"
    break
  fi

  inrds="${inrds_single_ld}"
  outprefix="${outrds_stem}.singleld.${kmin}-${kmax}-${B}.${seed}-${i}-${N}"
  queuefile="${outprefix}.queued"
  lockfile="${outprefix}.lock"
  donefile="${outprefix}.done"

  if [ -e "${queuefile}" ]; then
    echo "${i} already queued"
  elif [ -e "${lockfile}" ]; then
    echo "${i} already running"
  elif [ -e "${donefile}" ]; then
    echo "${i} already done"
  else
    qsub -z \
      -v INRDS="${inrds}",OUTPREFIX="${outprefix}",SEED=${seed},KMIN=$kmin,KMAX=$kmax,B=$B,CORES=$cores,QUEUEFILE="${queuefile}",LOCKFILE="${lockfile}",DONEFILE="${donefile}" \
      -l ncpus="$cores" -N "nf${i}" \
      factorise.pbs && \
    touch "${queuefile}" && \
    echo "${i} queued" && \
    num_queued=$((num_queued+1)) && \
    sleep 1
  fi

  inrds="${inrds_agegrouped}"
  outprefix="${outrds_stem}.agegrouped.${kmin}-${kmax}-${B}.${seed}-${i}-${N}"
  queuefile="${outprefix}.queued"
  lockfile="${outprefix}.lock"
  donefile="${outprefix}.done"

  if [ -e "${queuefile}" ]; then
    echo "${i} already queued"
  elif [ -e "${lockfile}" ]; then
    echo "${i} already running"
  elif [ -e "${donefile}" ]; then
    echo "${i} already done"
  else
    qsub -z \
      -v INRDS="${inrds}",OUTPREFIX="${outprefix}",SEED=${seed},KMIN=$kmin,KMAX=$kmax,B=$B,CORES=$cores,QUEUEFILE="${queuefile}",LOCKFILE="${lockfile}",DONEFILE="${donefile}" \
      -l ncpus="$cores" -N "nf${i}" \
      factorise.pbs && \
    touch "${queuefile}" && \
    echo "${i} queued" && \
    num_queued=$((num_queued+1)) && \
    sleep 1
  fi

  seed=$((seed+1))
done
