
#Execution:
#awk -f ~/opt/scripts/rename_tara_orfs.awk orfs_info.gff
#Header format
#c_scaffold|contig_strand_start_end_orf-num

{
    if($0 ~ /^>/){
    split($1,a,"_")
    if ($7 == "1"){
        STRAND = "+"
    }else{
        STRAND = "-"
    }
    print a[1]"_"a[2]"_"STRAND"_"$3"_"$5"_orf-"a[3]
    }else{
    O=gsub("\\*","",$0)
    print $O
    }
}
