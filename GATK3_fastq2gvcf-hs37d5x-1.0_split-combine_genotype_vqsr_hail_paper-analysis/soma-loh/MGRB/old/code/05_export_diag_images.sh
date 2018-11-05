#!/bin/bash

function export_images
{
  infile="$1"
  sampleid=$(basename "${infile}")
  sampleid="${sampleid%.pdf}"
  echo "${sampleid}"
  gs -q -dNOPAUSE -dBATCH -dSAFER -dFirstPage=1 -dLastPage=1 -sDEVICE=pnggray -sOutputFile="../${sampleid}.1.raw_vaf.png" -r150 "${infile}"
  gs -q -dNOPAUSE -dBATCH -dSAFER -dFirstPage=2 -dLastPage=2 -sDEVICE=pnggray -sOutputFile="../${sampleid}.2.filt_vaf.png" -r150 "${infile}"
  gs -q -dNOPAUSE -dBATCH -dSAFER -dFirstPage=3 -dLastPage=3 -sDEVICE=pnggray -sOutputFile="../${sampleid}.3.raw_full_trace.png" -r150 "${infile}"
  gs -q -dNOPAUSE -dBATCH -dSAFER -dFirstPage=4 -dLastPage=4 -sDEVICE=png256 -dTextAlphaBits=4 -sOutputFile="../${sampleid}.4.smoothed_trace.png" -r150 "${infile}"
}

export -f export_images

parallel export_images {} ::: ../*.aneu.diag.pdf
