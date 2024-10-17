#! /usr/bin/gawk -f
# this script removes space between SNP genotypes
# skipping "skip" fields in the entrance
# and formatting as UGA
# call as oneline2gs3.awk skip=2 exemple
BEGIN{skip=1
      missing=0}
{
    nsnp=(NF-skip)/2
    for (i =1; i<=skip; i++){
      printf( "%10s",$i)
    }
    printf( "%1s"," ")
    pos=0
    for (i =(skip+1); i<=NF; i=i+2){
      out=$i+$(i+1)-2
      if ($i==0 || $(i+1)==0) {out=5; missing++}
      pos +=1
      if(pos%4 != 0){
      	printf( "%1s",out)
      }	
    }
    printf("\n")
 }
END{}
# END{printf( "%10s\n",missing ) > "/dev/stderr"}


