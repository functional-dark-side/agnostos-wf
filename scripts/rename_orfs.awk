
#Execution:
#awk -f ~/opt/scripts/rename_orfs.awk orfs_info.gff
#Header format
#c_scaffold|contig_strand_start_end_orf-num

{
    if($0 ~ /^>/){
      split($9,a,"_");
      split(a[2],b,";");
      sub("_[^_]*$","",$1);
    if ($7 == "1"){
        STRAND = "+"
    }else{
        STRAND = "-"
    }
    print $1"_"STRAND"_"$3"_"$5"_orf-"b[1]
    }else{
    O=gsub("\\*","",$0)
    print $O
    }
}
