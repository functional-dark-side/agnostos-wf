#! /usr/bin/awk -f

BEGIN {FS=OFS="\t";}
    {
      print NR,$1,(NF-1)
    }
