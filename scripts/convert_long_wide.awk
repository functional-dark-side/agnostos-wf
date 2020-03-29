#! /usr/bin/awk -f

BEGIN { getline; id=$1; line=$0 }

{ if ($1 != id)
  {
    print line;
    line = $0;
  }
  else
  {
    line = line"\t"$2;
  }
    id=$1;
}

END { print line; }
