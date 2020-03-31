#! /usr/bin/awk -f

$1!~/^#/{
  split($9,a,"_");
  split(a[2],b,";");
  split(b[2],c,"=");
  print $1"_"$7"_"$4"_"$5"_""orf-"b[1]"\t""\""c[2]"\""
  }
