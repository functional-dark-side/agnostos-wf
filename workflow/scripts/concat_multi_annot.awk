#! /usr/bin/awk -f

BEGIN { getline; id=$1; l1=$1;l2=$2;l3=$3;l4=$4;}

{ if ($1 != id)
  {
    print l1,l2,l3,l4;
    l1=$1;l2=$2;l3=$3;l4=$4;
  }
  else
  {
    l2=l2"|"$2;
    l3 =l3"|"$3;
    l4=l4"|"$4
  } id=$1;
}

END { print l1,l2,l3,l4; }
