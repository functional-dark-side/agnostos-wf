
#Execution:
#awk -f ~/opt/scripts/rename_tara_orfs.awk TARA_004_DCM_0.22-1.6_scaffold.nt.fasta
#Header format
#>{TARASAMPLE+SCAFFOLDNUM}_{STRAND}_{START}_{END}_orf-{GENENUM}
#TARA_004_DCM_0.22-1.6_scaffold2_1_2_+_1334_1912_orf_2

{
    if($0 ~ /^>/){
    split($1,a,"_")
    if ($7 == "1"){
        STRAND = "+"
    }else{
        STRAND = "-"
    }
    print a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5]"_"a[6]"_"STRAND"_"$3"_"$5"_orf-"a[7]
    }else{
    O=gsub("\\*","",$0)
    print $O
    }
}
