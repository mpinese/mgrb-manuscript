BEGIN {
  FS="\t"
  OFS="\t"
}

(NR == 1 || ( \
  $8 >= 2 && \
  $17 >= 0.9 && \
  ($7 == 0 || $18 >= 0.9) && \
  $15 >= 50 && \
  $16 >= 50 && \
  $13 / $8 >= 30 && \
  ($7 == 0 || $12 / $7 >= 30) && \
  $19 <= 20 && \
  $20 <= 20 && \
  $23 <= 20 && \
  $26 <= 20 && \
  $29 <= 20 && \
  $30 <= 20));

