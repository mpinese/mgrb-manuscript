#!/bin/bash

pdf2jpg()
{
  echo "$1"
  convert -density 300 "$1" -resize 50% -quality 85 "$2"
}

export -f pdf2jpg

./bin/parallel --gnu pdf2jpg {} {.}.jpg ::: ../*.pdf
