#!/bin/bash

set -e
set -x
THREADS=$2

# We need the following fields from the HHblits output:
# pfam_accession, pfam_length, cluster_id, cl_length, probability, pfam_aln_start, pfam_aln_end, cl_aln_start, cl_aln_end
# 1: cl_name |2: pfam_length|3: cluster_id|4: cl_length|5: probability|6: pfam_aln_start|7: pfam_aln_end|8: cl_aln_start|9: cl_aln_end
awk '{print $1,$12,$2,$11,$3,$9,$10,$7,$8}' "${1}" | sed 's/ /\t/g' |
  sort --parallel "${THREADS}" -k 3,3 -k 8n -k 9n |
  perl -e 'while(<>){chomp;@a=split(/\t/,$_);if(($a[-1]-$a[-2])>80){print $_,"\t",($a[-3]-$a[-4])/$a[1],"\n" if $a[4]>=90;}else{print $_,"\t",($a[-3]-$a[-4])/$a[1],"\n" if $a[4]>=90;}}' |
  awk '$NF > 0.4' |
  awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' |
  perl -e 'while(<>){chomp;@a=split;next if $a[-1]==$a[-2];push(@{$b{$a[2]}},$_);}foreach(sort keys %b){@a=@{$b{$_}};for($i=0;$i<$#a;$i++){@b=split(/\t/,$a[$i]);@c=split(/\t/,$a[$i+1]);$len1=$b[-1]-$b[-2];$len2=$c[-1]-$c[-2];$len3=$b[-1]-$c[-2];if($len3>0 and ($len3/$len1>0.5 or $len3/$len2>0.5)){if($b[4]>$c[4]){splice(@a,$i+1,1);}else{splice(@a,$i,1);}$i=$i-1;}}foreach(@a){print $_."\n";}}'
