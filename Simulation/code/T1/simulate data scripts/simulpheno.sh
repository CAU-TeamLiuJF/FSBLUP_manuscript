# extract qtl alleles in a readable format
./oneline2uga.awk r_datasim/p0_qtl_001.txt  >kk
# run r script to simulate traits m, n and y
R CMD BATCH --no-restore --no-save ./simulN.r
# check h2 of y with resulting output - should be close to broad sense h2
# /save/alegarra/progs/oneline2uga.awk p0_mrk_001.txt |awk 'NR>1' > genotypes
# awk '{print 1,NR,$1}' y.txt > pheno
# awk '{print NR,0,0}' y.txt > ped
# awk '{print NR,NR}' y.txt > genotypes_XrefID
# check variances are the same
#echo "var 100 first y"
#head -100 y.txt | stat2.awk
#echo "var all y"
#stat2.awk y.txt
echo "id ebv" > my_bv.txt
awk '{print NR,$1}' y.txt >> my_bv.txt
