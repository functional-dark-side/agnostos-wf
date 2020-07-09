#!/bin/awk -f

{
 if ($1 in array) {
	A=array[$1]
	L=(log($11)/log(10));
	 if(L <= A) {
		print $0;
	}
}else{
	L=P*(log($11)/log(10));
	array[$1]=L;
	print $0;
     }
}
