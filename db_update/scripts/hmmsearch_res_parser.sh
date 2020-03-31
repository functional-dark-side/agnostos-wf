# Yanbin Yin
# 07/21/2015
# hmmscan output parser
# Usage: sh hmmsearch-parser.sh hmmsearch-output-file e-value cov_fract

# 1. take hmmer3 --domtblout output and extract necessary columns
# 2. sort on the protein position columns
# 3. remove overlapping hmm matches and keep the one with the lower e-values
# 4. calculate the covered fraction of hmm &
#    apply the E-value cutoff and the covered faction cutoff

if [ "$(uname -s)" == "Darwin" ]; then
   mysed=gsed
else
   mysed=sed
fi

awk '{print $4,$6,$1,$3,$13,$16,$17,$18,$19}' ${1} | $mysed 's/ /\t/g' | \
  $mysed 's/ /\t/g' | \
  sort --parallel 3 -k 3,3 -k 8n -k 9n | \
  perl -e 'while(<>){chomp;@a=split;next if $a[-1]==$a[-2];push(@{$b{$a[2]}},$_);}foreach(sort keys %b){@a=@{$b{$_}};for($i=0;$i<$#a;$i++){@b=split(/\t/,$a[$i]);@c=split(/\t/,$a[$i+1]);$len1=$b[-1]-$b[-2];$len2=$c[-1]-$c[-2];$len3=$b[-1]-$c[-2];if($len3>0 and ($len3/$len1>0.5 or $len3/$len2>0.5)){if($b[4]<$c[4]){splice(@a,$i+1,1);}else{splice(@a,$i,1);}$i=$i-1;}}foreach(@a){print $_."\n";}}' | \
  E=${2} perl -e 'while(<>){chomp;@a=split(/\t/,$_);if(($a[-1]-$a[-2])>80){print $_,"\t",($a[-3]-$a[-4])/$a[1],"\n" if $a[4]<$ENV{E};}else{print $_,"\t",($a[-3]-$a[-4])/$a[1],"\n" if $a[4]<$ENV{E};}}' | awk -v C=${3} '$NF>C' | sort --parallel 4 -k 3 -k 8,9g
