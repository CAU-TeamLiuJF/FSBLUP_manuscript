nQTL=20000
# we know that we have 20 chromosomes with 1000 QTLs each + 1 id
#nQTLperM=20
#nQTLar=1000
#sigma2m=2
#h2m=0.75
#nM=500
#h2=0.20
#h2ar=0.04
#varalpha=0.1
#c2m=(h2-h2ar)/h2m

# Values from Weishaar et al paper, FCR
h2=0.16
h2ar=0.1236
h2m=0.28
# more realistic values, also use 5000 putative QTLs (1 every 4 markers) and 500 QTL per omic  and for the "ar" trait
nloci=20000
inter=4
nQTL=nloci/inter
nQTLperM=500
nQTLar=500
sigma2m=2
nM=1200
sigma2alpha=0.01
c2m=(h2-h2ar)/h2m   # 0.13

cat("c2m = ",c2m)
if(c2m<0 | c2m>1){
	cat("impossible values\n")
}

# from here we obtain
sigma2am=h2m*nM*sigma2m*sigma2alpha
sigma2em=(1-h2m)*nM*sigma2m*sigma2alpha
cat("sigma2am= ",sigma2am,"\n")
cat("sigma2em= ",sigma2em,"\n")
sigma2p=sigma2am+sigma2em
cat("sigma2p=sigma2em+sigma2am= ",sigma2p,"\n")
sigma2y=h2m*sigma2p/(h2-h2ar) 
cat("sigma2y= ",sigma2y,"\n")
sigma2ar=sigma2y*h2ar
cat("sigma2ar= ",sigma2ar,"\n")
sigma2epsilon=sigma2y*(1-c2m-h2ar)
cat("sigma2epsilon= ",sigma2epsilon,"\n")



# number of generation need for rescaling
ngen=c(rep(1,1100))
for (i in 2:11){ ngen=c(ngen,rep(i,2000)) }

a=scan("kk") # file with genotypes
geno=matrix(a,ncol=nloci+1,byrow=TRUE)
# reduce to matrix of QTL keeping every 4 locus
pos=which(1:nloci %% inter ==0)
# add 1st col to keep id
pos=c(1,pos+1)
geno=geno[,pos]

# we compute scaling factors using all markers for simplicity
# remmber there is the id in geno
freq=apply(geno[1:1100,-1],2,mean)/2
sc = sqrt(nQTLperM*2*sum(freq*(1-freq))/nQTL)
nanim=nrow(geno)

# now we find out which is the new batch of animals, so we don't re-compute anything for the old ones
# this is to avoid that samples are run differently for the same animals across the new generations
lastgen=ngen[nanim] # the last generation is the generation of the last animal
geno=geno[ngen[1:nanim]==lastgen,]
nanim=nrow(geno)
cat("lastgen ",lastgen, "with ",nanim,"animals \n")
set.seed(1234*lastgen) # so that same residuals are not sampled over and over


id=geno[,1]
geno=geno[,-1] # remove id
rm(a)
qtl=matrix(NA,ncol=nM,nrow=nQTLperM)
u  =matrix(NA,ncol=nM,nrow=nQTLperM)
if(lastgen==1){
	alpha=rnorm(nM,0,sqrt(sigma2alpha))
	cat("sigma2alpha sampled= ",var(alpha),"\n")
	write(alpha,file="alpha.txt",ncolumns=1)
	for (i in 1:nM){
		qtl[,i]=sample(1:nQTL,nQTLperM) # perhaps with overlap across Ms but we assume overlap across Ms to be small
		u[,i]=rnorm(nQTLperM)
	}
	write(qtl,file="qtl.txt",ncolumns=nM)
	write(u,file="u.txt",ncolumns=nM)
	#draw ar
	qtlar=sample(1:nQTL,nQTLar)
	uar=rnorm(nQTLar)
	write(qtl,file="qtl.txt",ncolumns=nM)
	write(qtlar,file='qtlar.txt',ncolumns=1)
	write(u,file="u.txt",ncolumns=nM)
	write(uar,file='uar.txt',ncolumns=1)
}else{
	alpha=scan("alpha.txt")
	qtl=as.matrix(read.table("qtl.txt"))
	qtlar=as.matrix(read.table("qtlar.txt"))
	u=as.matrix(read.table("u.txt"))
	uar=as.matrix(read.table("uar.txt"))
}


# 1st layer
G=matrix(NA,ncol=nM,nrow=nanim) # here I keep only the genetic part, what Ole calls G
E=matrix(NA,ncol=nM,nrow=nanim) 
for (i in 1:nM){
	W=geno[,qtl[,i]]
        G[,i]=W%*%u[,i]
	# rescale to h2m so that on output var(M)=sigma2m
	# scale constant is based on the 1st 1100 animals (generation 1) so scaling is constant
	# across all calls
	G[,i]=sqrt(sigma2m)*sqrt(h2m)*G[,i]/sc
	varG=var(G[,i])
	varE=sigma2m*(1-h2m)
	E[,i]=rnorm(nanim,0,sd=sqrt(varE))
	varE=var(E[,i])
	cat("nM=",i,"varG= ",varG,"varE= ",varE,"varM= ",varG+varE,"h2m ",varG/(varG+varE),"\n")
}
write(t(G+E),file="M.txt",ncolumns=nM,append=(lastgen!=1))

tbv_m = as.vector(G%*%alpha)
e_m   = as.vector(E%*%alpha)
p = tbv_m + e_m
cat("var(tbv_m) sampled ",var(tbv_m),"\n")
cat("var(e_m) sampled ",var(e_m),"\n")
cat("var(p) sampled ",var(p),"\n")



cat("computing required numbers\n")
s2p=var(p)
cat("s2p= ",s2p,"\n")
s2y=h2m*s2p/(h2-h2ar)
cat("s2y= ",s2y,"\n")
var_ar=h2ar*s2y
cat("var_ar= ",var_ar,"\n")
c2m=(h2-h2ar)/h2m
cat("c2m= ",c2m,"\n")
var_epsilon=(1-c2m-h2ar)*s2y
cat("var_epsilon= ",var_epsilon,"\n")
W=geno[,qtlar]
ar=as.vector(W%*%uar)
ar=sqrt(var_ar)*ar/sd(ar)

epsilon=rnorm(nanim,0,sqrt(var_epsilon))

y=tbv_m+e_m+ar+epsilon
tbv=tbv_m+ar
res=e_m+epsilon

cat("after simulation: \n")
cat("vary= ",var(y),"\n")
cat("c2m= ",var(tbv_m+e_m)/var(y),"\n")
cat("h2= ",var(tbv_m+ar)/var(y),"\n")
cat("h2ar= ",var(ar)/var(y),"\n")

# these are total genetic values , part of them are TBV because we don't multiply here
write(t(cbind(tbv_m,ar,tbv)),file="tbv.txt",ncolumns=3,append=(lastgen!=1))
write(t(cbind(tbv_m,e_m,p,ar,epsilon)),file="components.txt",ncolumns=5,append=(lastgen!=1))
write(y,file="y.txt",ncolumns=1,append=(lastgen!=1))



