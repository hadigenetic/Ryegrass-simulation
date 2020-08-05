# ---------------------------------
# Functions------Start-------------
# ---------------------------------

# Function for Recombination	
rec<-function(sire,linkage_map,chrom_length){

# Hatman kamel tozihat benevis
# This function is to recombintae a genome of an individual.
# input: Sire [This is a one row matrix or numeric where genotypes are as 11, ... 
#             where same allele is from sire and the other from dam]
#        linkage_map [This is output of functions in xbreed, three colomn data frame] 
#        chrom_length [length of chromosome in cM]

# Output: Recombinated genome as a Matrix with two rows.


li_map<-linkage_map
nchr<-length(unique(li_map[,2]))
cd<-as.vector(table(li_map[,2]))
a<-cumsum(cd)*2
a
p1<-c(1,a+1)
p1<-p1[1:(length(p1)-1)]
p2<-a

cd2<-cumsum(cd)
d_1<-c(1,cd2+1)
d_1<-d_1[1:(length(d_1)-1)]
d_2<-cd2

recombinated<-matrix(ncol=max(a)/2,nrow=2)

for (chr in 1:nchr){
	int_map<-subset(li_map,li_map[,2]==chr)
	int_map<-int_map[,3]
	s_chr_i<-sire[p1[chr]:p2[chr]]
	seq1<-seq(1,length(s_chr_i),2)
	seq2<-seq(2,length(s_chr_i),2)
	Str_sire1<-s_chr_i[seq1]
	Str_sire2<-s_chr_i[seq2]

	#recombination in sires
	s1<-Str_sire1
	s2<-Str_sire2
	nrec<-rpois(1,chrom_length/100)
	#-------------------
	if (nrec==0){
	nrec<-1
	}
	#-------------------
	pos_rec<-sample(int_map,nrec) 
	pos_rec<-sort(pos_rec)
	qq<-match(pos_rec,int_map)
	pos_rec<-qq
	s1new<-as.numeric()
	s2new<-as.numeric()
	if (nrec==0){
	s1new<-s1
	s2new<-s2
	} else if (nrec>0){
	s1new<-s1
	s2new<-s2
	for (y in 1:nrec){
	s1new[pos_rec[y]:length(s1)]<-s2[pos_rec[y]:length(s1)]
	s2new[pos_rec[y]:length(s1)]<-s1[pos_rec[y]:length(s1)]
	s1<-s1new
	s2<-s2new
	}
	}
	sire1<-as.numeric(s1)
	sire2<-as.numeric(s2)
	

	recombinated[1,d_1[chr]:d_2[chr]]<-sire1
	recombinated[2,d_1[chr]:d_2[chr]]<-sire2
}
	return(recombinated)
	
}



# Function for F-parent mating in Syn0
X_parent_mating<-function(size_polycross,individuals,indi_id,noff,link_map,chrom_length){


offspring<-matrix(nrow=noff*2,ncol=(dim(individuals)[2]/2))
seqof1<-seq(1,noff*2,2)
seqof2<-seq(2,noff*2,2)

	all_mats<-combn(1:size_polycross,2)
	dim(all_mats)
	
	whole_offsprings<-c()
	for (comb in 1:dim(all_mats)[2]){
		for (T in 1:noff){
		ali1<-rec(individuals[all_mats[1,comb],],link_map,chrom_length)
		ali2<-rec(individuals[all_mats[2,comb],],link_map,chrom_length)
		ali<-rbind(ali1,ali2)
			
			res1<-sample(c(1,2),1)
			res2<-sample(c(3,4),1)

			offspring[seqof1[T],]<-ali[res1,]
			offspring[seqof2[T],]<-ali[res2,]
		}
		whole_offsprings<-rbind(whole_offsprings,offspring)
	}
	
	
	# Make pedigree
	ped<-c()
		for (comb in 1:dim(all_mats)[2]){
		ul_1<-rep(indi_id[all_mats[1,comb]],noff)
		ul_2<-rep(indi_id[all_mats[2,comb]],noff)
		ul<-cbind(ul_1,ul_2)
		ped<-rbind(ped,ul)
		}
	ped<-as.matrix(ped)

	
	whole_offsprings<-as.matrix(whole_offsprings)
	out<-list()
	out[[1]]<-whole_offsprings
	out[[2]]<-ped
	
return(out)
}

	 
	bin_snp<-function(mat){
	s1<-seq(1,ncol(mat),2)
	s2<-seq(2,ncol(mat),2)
	a1<-mat[,s1]+mat[,s2]
	a1[a1==3]=1
	a1[a1==4]=0
	snp_code<-a1
	return(snp_code)
	 }
	 
	
	# Function for Freq Calc
	 
	Calc_freq<-function(Loci_Mat){
	freq1<-as.numeric()
	s1<-seq(1,length(Loci_Mat[1,]),2)
	s2<-seq(2,length(Loci_Mat[1,]),2)
	a1<-Loci_Mat[,s1]+Loci_Mat[,s2]
	a1[a1==3]=1
	a1[a1==4]=0
	Loci_Mat2<-a1

	dunkan<-function(vec){
		sum(vec)/(length(vec)*2)
		}

freq1<-apply(Loci_Mat2,2,dunkan)
return(freq1)
	}
	
	
	# Function for TBV

	calc_TBV<-function(mat,add_eff,dom_eff){
	freq1<-Calc_freq(mat)
	freq2<-1-freq1
	mat<-bin_snp(mat)
	add_eff_1<-add_eff
	add_eff_2<-add_eff_1*(-1)
	#Breeding Value
	xprogeny<-mat
	xprogeny<-as.matrix(xprogeny)
	q1<-xprogeny
	for (i in 1:length(xprogeny[,1])){
	ti<-xprogeny[i,]
	two<-which(ti==2)
	one<-which(ti==1)
	zero<-which(ti==0)
	q1[i,two]<-((freq1[two])*add_eff_1[two])+
			 ((freq2[two])*dom_eff[two])
	q1[i,one]<-((1/2*freq1[one])*add_eff_1[one])+
			 ((1/2*freq2[one])*dom_eff[one])+
			 ((1/2*freq1[one])*dom_eff[one])+
			 ((1/2*freq2[one])*add_eff_2[one])
	q1[i,zero]<-((freq2[zero])*add_eff_2[zero])+
			 ((freq1[zero])*dom_eff[zero])	
	}
	xprogeny<-q1
	tbv<- rowSums(xprogeny)
	return(tbv)
	}
	
	
	calc_TBV2<-function(mat,add_eff){

	mat<-bin_snp(mat)
	add_eff_1<-matrix(0,ncol=1,nrow=length(add_eff))
	add_eff_1[,1]<-add_eff

	#Breeding Value
	xprogeny<-mat
	xprogeny<-as.matrix(xprogeny)
    tbv<-xprogeny%*%add_eff_1
	return(tbv)
	}
	
	
	# Function for TGV
	calc_TGV<-function(mat,add_eff,dom_eff){
	add_eff_1<-add_eff
	add_eff_2<-add_eff_1*-1
	mat<-bin_snp(mat)
	xprogeny<-mat
	xprogeny<-as.matrix(xprogeny)
	q1<-xprogeny
	for (i in 1:length(xprogeny[,1])){
	ti<-xprogeny[i,]
	two<-which(ti==2)
	one<-which(ti==1)
	zero<-which(ti==0)
	q1[i,two]<-add_eff_1[two]
	q1[i,one]<-dom_eff[one]
	q1[i,zero]<-add_eff_2[zero]
	}
	xprogeny<-q1
	tgv<- rowSums(xprogeny)
	return(tgv)
	}
	
	
	
	# Function for GEBV
	calc_GEBV<-function(mat_own,ghat1,ghat2,dhat){
	
	g1<-ghat1
	g2<-ghat2
	d<-dhat
	
	freq1<-as.numeric()
	for (i in 1:length(mat_own[1,])){
	freq1[i]<-sum(mat_own[,i])/(length(mat_own[,i])*2)
	}
	freq2<-1-freq1
	#Breeding Value
	xprogeny<-mat_own
	xprogeny<-as.matrix(xprogeny)
	q1<-xprogeny
	for (i in 1:length(xprogeny[,1])){
	ti<-xprogeny[i,]
	two<-which(ti==2)
	one<-which(ti==1)
	zero<-which(ti==0)
	q1[i,two]<-((freq1[two])*g1[two])+
			 ((freq2[two])*d[two])
	q1[i,one]<-((1/2*freq1[one])*g1[one])+
			 ((1/2*freq2[one])*d[one])+
			 ((1/2*freq1[one])*d[one])+
			 ((1/2*freq2[one])*g2[one])
	q1[i,zero]<-((freq2[zero])*g2[zero])+
			 ((freq1[zero])*d[zero])	
	}
	xprogeny<-q1
	gebvs<- rowSums(xprogeny)
		return(gebvs)
	}



	# Function for GEBV
	calc_GEBV2<-function(mat_own,ghat1){
	
	mat_own<-bin_snp(mat_own)
	g1<-matrix(0,ncol=1,nrow=length(ghat1))
	g1[,1]<-ghat1

	#Breeding Value
	xprogeny<-mat_own
	xprogeny<-as.matrix(xprogeny)
    gebvs<-xprogeny%*%g1
	return(gebvs)
	}

	# Function for scaling add effects to get desired h2
	
scale_trait_effect<-function(h2,phen_var,nqtl,sampled_effects,qtlMatrix){


addTra<-TRUE

# Create Trait 1
cat('Simulating trait ...',fill=TRUE)
var_add<-h2*phen_var
if(addTra==FALSE){
var_dom	<-d2*phen_var_tr1
} else {
var_dom<-0
d2<-0
}
add_eff_1<-sampled_effects
add_eff_2<-add_eff_1*-1

dom_degree<-rnorm(nqtl,mean=0.5)
if(addTra==FALSE){
dom_eff<-dom_degree*abs(add_eff_1_tr1)
} else {
dom_eff<-rep(0,length(dom_degree))
}


#TBV based on allele frequency
# freq1<-Calc_freq(qtlLoci)
qtlLoci<-qtlMatrix
qtlLoci<-as.matrix(qtlLoci)

freq1<-Calc_freq(qtlLoci)

freq2<-1-freq1
freq1qtl<-freq1
freq2qtl<-freq2
Snp_BreedB<-bin_snp(qtlLoci)

# alpha=abs(add_eff_1)+((freq2-freq1)*dom_eff)
alpha=add_eff_1+((freq2-freq1)*dom_eff)

#Breeding Value
xprogeny<-Snp_BreedB
xprogeny<-as.matrix(xprogeny)
q1<-xprogeny
for (i in 1:length(xprogeny[,1])){
ti<-xprogeny[i,]
two<-which(ti==2)
one<-which(ti==1)
zero<-which(ti==0)
q1[i,two]<-2*freq2[two]*alpha[two]
q1[i,one]<-(freq2[one]-freq1[one])*alpha[one]
q1[i,zero]<--2*freq1[zero]*alpha[zero]
}
xprogeny<-q1
tbv<- rowSums(xprogeny)
var(tbv)


#Dominance Deviation
xprogeny<-Snp_BreedB
xprogeny<-as.matrix(xprogeny)
q1<-xprogeny
for (i in 1:length(xprogeny[,1])){
ti<-xprogeny[i,]
two<-which(ti==2)
one<-which(ti==1)
zero<-which(ti==0)
q1[i,two]<-(freq2[two]**2)*-2*dom_eff[two]  #A1A1  -2q^2d   A1A2 2pqd   A2A2 -2p^2d
q1[i,one]<-2*freq1[one]*freq2[one]*dom_eff[one]
q1[i,zero]<--2*(freq1[zero]**2)*dom_eff[zero]
}
xprogeny<-q1
dom_dev<- rowSums(xprogeny)
var(dom_dev)

# scaling
ratio_add<-sqrt(var_add)/sqrt(var(tbv))
if(addTra==FALSE){
ratio_dom<-sqrt(var_dom)/sqrt(var(dom_dev))
} else {
ratio_dom<-0
}
dom_eff_new<-dom_eff*ratio_dom
dom_eff<-dom_eff_new
alpha_new<-alpha*ratio_add
alpha<-alpha_new
add_eff_1<-alpha-((freq2-freq1)*dom_eff)
add_eff_2<-add_eff_1*(-1)

#Breeding Value
# alpha=abs(add_eff_1)+((freq2-freq1)*dom_eff)
alpha=add_eff_1+((freq2-freq1)*dom_eff)
xprogeny<-Snp_BreedB
xprogeny<-as.matrix(xprogeny)
q1<-xprogeny
for (i in 1:length(xprogeny[,1])){
ti<-xprogeny[i,]
two<-which(ti==2)
one<-which(ti==1)
zero<-which(ti==0)
q1[i,two]<-2*freq2[two]*alpha[two]
q1[i,one]<-(freq2[one]-freq1[one])*alpha[one]
q1[i,zero]<--2*freq1[zero]*alpha[zero]
}
xprogeny<-q1
tbv<- rowSums(xprogeny)
var(tbv)

#Dominance Deviation
xprogeny<-Snp_BreedB
xprogeny<-as.matrix(xprogeny)
q1<-xprogeny
for (i in 1:length(xprogeny[,1])){
ti<-xprogeny[i,]
two<-which(ti==2)
one<-which(ti==1)
zero<-which(ti==0)
q1[i,two]<-(freq2[two]**2)*-2*dom_eff[two]  #A1A1  -2q^2d   A1A2 2pqd   A2A2 -2p^2d
q1[i,one]<-2*freq1[one]*freq2[one]*dom_eff[one]
q1[i,zero]<--2*(freq1[zero]**2)*dom_eff[zero]
}
xprogeny<-q1
dom_dev<- rowSums(xprogeny)
var(dom_dev)

variance_add_falc<-sum(2*freq1*freq2*(alpha**2))



vars<-c(var(tbv),var(dom_dev),variance_add_falc)
scaled_effects<-add_eff_1



# output

back<-list()
back[[1]]<-vars
back[[2]]<-data.frame(scaled_effects,alpha)
back[[3]]<-tbv
back[[4]]<-data.frame(freq1,freq2)


return(back)
}


# Function for 2-parent mating

mating<-function(individuals,noff,link_map,chrom_leng){


# This function is for mating among two parents. Each parent, which is
# one row in arg individuals, produces two strains (in total 4) and they are sampled
# to make a new offspring.

# Input: - individuals [This is a 2 row matrix where genotypes are as 11, ... 
#        where same allele is from sire and the other from dam]
#        - noff [number of offspring]
#        - link_map [This is output of functions in xbreed, three colomn data frame] 

# Output: Produced offspring in a Matrix with nrows=noff*2 and ncol=length(individuals)/2


offspring<-matrix(nrow=noff*2,ncol=(dim(individuals)[2]/2))


seqof1<-seq(1,noff*2,2)
seqof2<-seq(2,noff*2,2)
	for (T in 1:noff){
	ali1<-rec(individuals[1,],link_map,chrom_leng)
	ali2<-rec(individuals[2,],link_map,chrom_leng)
	ali<-rbind(ali1,ali2)
		
		res1<-sample(seq(1,4,2),1)
		res2<-sample(seq(2,4,2),1)
		while(res2==res1+1){
		res2<-sample(seq(2,4,2),1)
		}

		offspring[seqof1[T],]<-ali[res1,]
		offspring[seqof2[T],]<-ali[res2,]
	
	}
	
return(offspring)
}


# Function to create F1 (Internal Func for other functions)
Create_F1_DLF<-function(v1,v2,F1_No,sires,dams,qtl_effects,
noff,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var){


internal_F2<-list()

	indi<-as.matrix(rbind(v1,v2))
	kpo<-mating(individuals=indi,
	noff=noff,
	link_map=link_map,
	chrom_leng=chrom_length)
	dim(kpo)
	
	offspring_2str<-kpo
	dim(offspring_2str)
	
offspring_2str<-as.matrix(offspring_2str)	
offspring<-matrix(0,ncol=dim(offspring_2str)[2]*2,nrow=dim(offspring_2str)[1]/2)
dim(offspring)

sei1<-seq(1,length(offspring_2str[,1]),2)	
sei2<-seq(2,length(offspring_2str[,1]),2)	

seqr1<-seq(1,length(offspring[1,]),2)	
seqr2<-seq(2,length(offspring[1,]),2)	

offspring[,seqr1]<-offspring_2str[sei1,]
offspring[,seqr2]<-offspring_2str[sei2,]

dim(offspring)

offspring_2str[1:2,1:5]
offspring[1:1,1:10]

id<-1:dim(offspring)[1]
sire<-rep(sires,each=noff)
dam <-rep(dams,each=noff)
sex <-sample(c('F','M'),length(id),replace=T)
F1_id<-rep(F1_No,length(id))

#Extract qtl_genome for each offspring from total sequence based on positions
	nqtl_allele<-dim(linkage_map_qtl)[1]*2
	qtlsel1 <- matrix(ncol=nqtl_allele,nrow = length(id))

	pos<-link_map[,3]
	posqtl<-linkage_map_qtl[,3]

	offspringasli_1<-rbind(pos,offspring)
	index<-posqtl
	index_qtls <-which(offspringasli_1[1,]%in%index)

	qtlsel1<-offspringasli_1[,index_qtls]
	qtlsel1<-qtlsel1[-1,]
	qtlsel1<-as.matrix(qtlsel1)
	
	#Extract mrk_genome for each offspring from total sequence based on positions
	nmrk_allele<-dim(linkage_map_mrk)[1]*2
	mrksel1 <- matrix(ncol=nmrk_allele,nrow = length(id))

	pos<-link_map[,3]
	posmrk<-linkage_map_mrk[,3]

	offspringasli_1<-rbind(pos,offspring)
	index<-posmrk
	index_mrks <-which(offspringasli_1[1,]%in%index)

	mrksel1<-offspringasli_1[,index_mrks]
	mrksel1<-mrksel1[-1,]
	mrksel1<-as.matrix(mrksel1)
			

# # delete row below if doesnt work
  # offspringasli_1<-as.matrix(offspringasli_1)
  
qtlLoci<-qtlsel1
qtlLoci<-as.matrix(qtlLoci)

tbv<-matrix(ncol=4,nrow=length(id))
env<-matrix(ncol=4,nrow=length(id))
phen<-matrix(ncol=4,nrow=length(id))

add_eff_1<-qtl_effects
dom_eff<-rep(0,length(add_eff_1[,1]))
var_dom<-0

	for (j in 1:4){
		# TBV 
		tbv[,j]<-calc_TBV2(qtlLoci,add_eff_1[,j])
		
		# ENV 
		vv<-phen_var[j]-(var_add[j]+var_dom)
		env[,j]<-rnorm(length(tbv[,j]),mean=0,sd=sqrt(vv))
			
		# Phen
		phen[,j]<-tbv[,j]+env[,j]
	}

data<-data.frame(id,sire,dam,sex,F1_id,phen,tbv,env)
names(data)<-c('id','sire','dam','sex','F1_No','phen_1','phen_2','phen_3','phen_4',
'tbv_1','tbv_2','tbv_3','tbv_4',
'env_1','env_2','env_3','env_4')

sequ<-offspring

internal_F2[[1]]<-data
internal_F2[[2]]<-qtlsel1
internal_F2[[3]]<-mrksel1
internal_F2[[4]]<-sequ

	return(internal_F2)
}



# Function to create F2 (Internal Func for other functions)
Create_F2<-function(v1,v2,F2_No,sires,dams,qtl_effects,
noff,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var){


internal_F2<-list()
offspring_2str<-c()

	for (k in 1:dim(v1)[1]){
	indi<-as.matrix(rbind(v1[k,],v2[k,]))
	
	kpo<-mating(individuals=indi,
	noff=noff,
	link_map=link_map,
	chrom_leng=chrom_length)
	
	offspring_2str<-rbind(offspring_2str,kpo)	
	}

	
offspring_2str<-as.matrix(offspring_2str)	
offspring<-matrix(0,ncol=dim(offspring_2str)[2]*2,nrow=dim(offspring_2str)[1]/2)


sei1<-seq(1,length(offspring_2str[,1]),2)	
sei2<-seq(2,length(offspring_2str[,1]),2)	

seqr1<-seq(1,length(offspring[1,]),2)	
seqr2<-seq(2,length(offspring[1,]),2)	

offspring[,seqr1]<-offspring_2str[sei1,]
offspring[,seqr2]<-offspring_2str[sei2,]



id<-1:dim(offspring)[1]
sire<-rep(sires,each=noff)
dam <-rep(dams,each=noff)
sex <-sample(c('F','M'),length(id),replace=T)
F2_id<-rep(F2_No,length(id))


#Extract qtl_genome for each offspring from total sequence based on positions
	nqtl_allele<-dim(linkage_map_qtl)[1]*2
	qtlsel1 <- matrix(ncol=nqtl_allele,nrow = length(id))

	pos<-link_map[,3]
	posqtl<-linkage_map_qtl[,3]


	offspringasli_1<-rbind(pos,offspring)
	index<-posqtl
	index_qtls <-which(offspringasli_1[1,]%in%index)

	qtlsel1<-offspringasli_1[,index_qtls]
	qtlsel1<-qtlsel1[-1,]
	qtlsel1<-as.matrix(qtlsel1)

	
	#Extract mrk_genome for each offspring from total sequence based on positions
	nmrk_allele<-dim(linkage_map_mrk)[1]*2
	mrksel1 <- matrix(ncol=nmrk_allele,nrow = length(id))

	pos<-link_map[,3]
	posmrk<-linkage_map_mrk[,3]

	offspringasli_1<-rbind(pos,offspring)
	index<-posmrk
	index_mrks <-which(offspringasli_1[1,]%in%index)

	mrksel1<-offspringasli_1[,index_mrks]
	mrksel1<-mrksel1[-1,]
	mrksel1<-as.matrix(mrksel1)

		

# # delete row below if doesnt work
  # offspringasli_1<-as.matrix(offspringasli_1)
  
qtlLoci<-qtlsel1
qtlLoci<-as.matrix(qtlLoci)

tbv<-matrix(ncol=4,nrow=length(id))
env<-matrix(ncol=4,nrow=length(id))
phen<-matrix(ncol=4,nrow=length(id))

add_eff_1<-qtl_effects
dom_eff<-rep(0,length(add_eff_1[,1]))
var_dom<-0

	for (j in 1:4){
		# TBV 
		# tbv[,j]<-calc_TBV(qtlLoci,add_eff_1[,j],dom_eff)
		tbv[,j]<-calc_TBV2(qtlLoci,add_eff_1[,j])
		
		# ENV 
		vv<-phen_var[j]-(var_add[j]+var_dom)
		env[,j]<-rnorm(length(tbv[,j]),mean=0,sd=sqrt(vv))
			
		# Phen
		# tgv<-calc_TGV(qtlLoci,add_eff_1[,j],dom_eff)
		# phen[,j]<-tgv+env[,j]
		phen[,j]<-tbv[,j]+env[,j]
	}


data<-data.frame(id,sire,dam,sex,F2_id,
phen,tbv,env)
names(data)<-c('id','sire','dam','sex','F2_id',
'phen_1','phen_2','phen_3','phen_4',
'tbv_1','tbv_2','tbv_3','tbv_4',
'env_1','env_2','env_3','env_4')

sequ<-offspring

internal_F2[[1]]<-data
internal_F2[[2]]<-qtlsel1
internal_F2[[3]]<-mrksel1
internal_F2[[4]]<-sequ

	return(internal_F2)
}


# Function to create SYN 1 (Internal Func for other functions)
Create_Syn1<-function(size_polycross,SYN_0,Syn_id,qtl_effects,
noff,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var){


indi<-c()
indi_id<-c()
for (sycounter in 1:size_polycross){
U<-SYN_0[[sycounter]][[4]] #sequence
indi<-rbind(indi,U)
pes<-SYN_0[[sycounter]][[1]][1]
indi_id[sycounter]<-as.numeric(pes)
}
indi<-as.matrix(indi)
dim(indi)
indi_id


internal_syn<-list()


kpo<-X_parent_mating(size_polycross=size_polycross,individuals=indi,indi_id,noff=noff,
link_map=link_map,chrom_length)

	offspring_2str<-kpo[[1]]
	dim(offspring_2str)
	sires<-kpo[[2]][,1]
	dams<-kpo[[2]][,2]


	
offspring_2str<-as.matrix(offspring_2str)	
offspring<-matrix(0,ncol=dim(offspring_2str)[2]*2,nrow=dim(offspring_2str)[1]/2)
dim(offspring)

sei1<-seq(1,length(offspring_2str[,1]),2)	
sei2<-seq(2,length(offspring_2str[,1]),2)	

seqr1<-seq(1,length(offspring[1,]),2)	
seqr2<-seq(2,length(offspring[1,]),2)	

offspring[,seqr1]<-offspring_2str[sei1,]
offspring[,seqr2]<-offspring_2str[sei2,]



id<-1:dim(offspring)[1]
sire<-sires
dam <-dams
sex <-sample(c('F','M'),length(id),replace=T)
F2_id<-rep(Syn_id,length(id))

#Extract qtl_genome for each offspring from total sequence based on positions
	nqtl_allele<-dim(linkage_map_qtl)[1]*2
	qtlsel1 <- matrix(ncol=nqtl_allele,nrow = length(id))

	pos<-link_map[,3]
	posqtl<-linkage_map_qtl[,3]


	offspringasli_1<-rbind(pos,offspring)
	index<-posqtl
	index_qtls <-which(offspringasli_1[1,]%in%index)

	qtlsel1<-offspringasli_1[,index_qtls]
	qtlsel1<-qtlsel1[-1,]
	qtlsel1<-as.matrix(qtlsel1)

	
	#Extract mrk_genome for each offspring from total sequence based on positions
	nmrk_allele<-dim(linkage_map_mrk)[1]*2
	mrksel1 <- matrix(ncol=nmrk_allele,nrow = length(id))

	pos<-link_map[,3]
	posmrk<-linkage_map_mrk[,3]

	offspringasli_1<-rbind(pos,offspring)
	index<-posmrk
	index_mrks <-which(offspringasli_1[1,]%in%index)

	mrksel1<-offspringasli_1[,index_mrks]
	mrksel1<-mrksel1[-1,]
	mrksel1<-as.matrix(mrksel1)
		

# # delete row below if doesnt work
  # offspringasli_1<-as.matrix(offspringasli_1)
  
qtlLoci<-qtlsel1
qtlLoci<-as.matrix(qtlLoci)

tbv<-matrix(ncol=4,nrow=length(id))
env<-matrix(ncol=4,nrow=length(id))
phen<-matrix(ncol=4,nrow=length(id))

add_eff_1<-qtl_effects
dom_eff<-rep(0,length(add_eff_1[,1]))
var_dom<-0

	for (j in 1:4){
		# TBV 
		# tbv[,j]<-calc_TBV(qtlLoci,add_eff_1[,j],dom_eff)
		tbv[,j]<-calc_TBV2(qtlLoci,add_eff_1[,j])
		
		# ENV 
		vv<-phen_var[j]-(var_add[j]+var_dom)
		env[,j]<-rnorm(length(tbv[,j]),mean=0,sd=sqrt(vv))
			
		# Phen
		# tgv<-calc_TGV(qtlLoci,add_eff_1[,j],dom_eff)
		# phen[,j]<-tgv+env[,j]
		phen[,j]<-tbv[,j]+env[,j]
	}


data<-data.frame(id,sire,dam,sex,F2_id,
phen,tbv,env)
names(data)<-c('id','sire','dam','sex','Syn_id',
'phen_1','phen_2','phen_3','phen_4',
'tbv_1','tbv_2','tbv_3','tbv_4',
'env_1','env_2','env_3','env_4')

sequ<-offspring

internal_syn[[1]]<-data
internal_syn[[2]]<-qtlsel1
internal_syn[[3]]<-mrksel1
internal_syn[[4]]<-sequ

	return(internal_syn)
}



# Function to create SYN 2 (Internal Func for other functions)
Create_Syn2<-function(v1,v2,Syn2_No,sires,dams,qtl_effects,
noff,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var){

# Information should be explained below
# v1
# v2

internal_Syn2<-list()
offspring_2str<-c()

	for (k in 1:dim(v1)[1]){
	indi<-as.matrix(rbind(v1[k,],v2[k,]))
	
	kpo<-mating(individuals=indi,
	noff=noff,
	link_map=link_map,
	chrom_leng=chrom_length)
	
	offspring_2str<-rbind(offspring_2str,kpo)	
	}

	
offspring_2str<-as.matrix(offspring_2str)	
offspring<-matrix(0,ncol=dim(offspring_2str)[2]*2,nrow=dim(offspring_2str)[1]/2)


sei1<-seq(1,length(offspring_2str[,1]),2)	
sei2<-seq(2,length(offspring_2str[,1]),2)	

seqr1<-seq(1,length(offspring[1,]),2)	
seqr2<-seq(2,length(offspring[1,]),2)	

offspring[,seqr1]<-offspring_2str[sei1,]
offspring[,seqr2]<-offspring_2str[sei2,]



id<-1:dim(offspring)[1]
sire<-rep(sires,each=noff)
dam <-rep(dams,each=noff)
sex <-sample(c('F','M'),length(id),replace=T)
F2_id<-rep(Syn2_No,length(id))


#Extract qtl_genome for each offspring from total sequence based on positions
	nqtl_allele<-dim(linkage_map_qtl)[1]*2
	qtlsel1 <- matrix(ncol=nqtl_allele,nrow = length(id))

	pos<-link_map[,3]
	posqtl<-linkage_map_qtl[,3]


	offspringasli_1<-rbind(pos,offspring)
	index<-posqtl
	index_qtls <-which(offspringasli_1[1,]%in%index)

	qtlsel1<-offspringasli_1[,index_qtls]
	qtlsel1<-qtlsel1[-1,]
	qtlsel1<-as.matrix(qtlsel1)

	
	#Extract mrk_genome for each offspring from total sequence based on positions
	nmrk_allele<-dim(linkage_map_mrk)[1]*2
	mrksel1 <- matrix(ncol=nmrk_allele,nrow = length(id))

	pos<-link_map[,3]
	posmrk<-linkage_map_mrk[,3]

	offspringasli_1<-rbind(pos,offspring)
	index<-posmrk
	index_mrks <-which(offspringasli_1[1,]%in%index)

	mrksel1<-offspringasli_1[,index_mrks]
	mrksel1<-mrksel1[-1,]
	mrksel1<-as.matrix(mrksel1)

		

# # delete row below if doesnt work
  # offspringasli_1<-as.matrix(offspringasli_1)
  
qtlLoci<-qtlsel1
qtlLoci<-as.matrix(qtlLoci)

tbv<-matrix(ncol=4,nrow=length(id))
env<-matrix(ncol=4,nrow=length(id))
phen<-matrix(ncol=4,nrow=length(id))

add_eff_1<-qtl_effects
dom_eff<-rep(0,length(add_eff_1[,1]))
var_dom<-0

	for (j in 1:4){
		# TBV 
		# tbv[,j]<-calc_TBV(qtlLoci,add_eff_1[,j],dom_eff)
		tbv[,j]<-calc_TBV2(qtlLoci,add_eff_1[,j])
		
		# ENV 
		vv<-phen_var[j]-(var_add[j]+var_dom)
		env[,j]<-rnorm(length(tbv[,j]),mean=0,sd=sqrt(vv))
			
		# Phen
		# tgv<-calc_TGV(qtlLoci,add_eff_1[,j],dom_eff)
		# phen[,j]<-tgv+env[,j]
		phen[,j]<-tbv[,j]+env[,j]
	}


data<-data.frame(id,sire,dam,sex,F2_id,
phen,tbv,env)
names(data)<-c('id','sire','dam','sex','Syn2_No',
'phen_1','phen_2','phen_3','phen_4',
'tbv_1','tbv_2','tbv_3','tbv_4',
'env_1','env_2','env_3','env_4')

sequ<-offspring

#adding id as first colomn to qtl,mrk and seque
qtlsel1<-cbind(id,qtlsel1)
qtlsel1<-as.matrix(qtlsel1)

mrksel1<-cbind(id,mrksel1)
mrksel1<-as.matrix(mrksel1)

sequ<-cbind(id,sequ)
sequ<-as.matrix(sequ)

internal_Syn2[[1]]<-data
internal_Syn2[[2]]<-qtlsel1
internal_Syn2[[3]]<-mrksel1
internal_Syn2[[4]]<-sequ

	return(internal_Syn2)
}


# Function to create SYN 3 (Internal Func for other functions)
Create_Syn3<-function(v1,v2,Syn3_No,sires,dams,qtl_effects,
noff,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var){

# Information should be explained below
# v1
# v2

internal_Syn3<-list()
offspring_2str<-c()

	for (k in 1:dim(v1)[1]){
	indi<-as.matrix(rbind(v1[k,],v2[k,]))
	
	kpo<-mating(individuals=indi,
	noff=noff,
	link_map=link_map,
	chrom_leng=chrom_length)
	
	offspring_2str<-rbind(offspring_2str,kpo)	
	}

	
offspring_2str<-as.matrix(offspring_2str)	
offspring<-matrix(0,ncol=dim(offspring_2str)[2]*2,nrow=dim(offspring_2str)[1]/2)


sei1<-seq(1,length(offspring_2str[,1]),2)	
sei2<-seq(2,length(offspring_2str[,1]),2)	

seqr1<-seq(1,length(offspring[1,]),2)	
seqr2<-seq(2,length(offspring[1,]),2)	

offspring[,seqr1]<-offspring_2str[sei1,]
offspring[,seqr2]<-offspring_2str[sei2,]



id<-1:dim(offspring)[1]
sire<-rep(sires,each=noff)
dam <-rep(dams,each=noff)
sex <-sample(c('F','M'),length(id),replace=T)
F2_id<-rep(Syn3_No,length(id))


#Extract qtl_genome for each offspring from total sequence based on positions
	nqtl_allele<-dim(linkage_map_qtl)[1]*2
	qtlsel1 <- matrix(ncol=nqtl_allele,nrow = length(id))

	pos<-link_map[,3]
	posqtl<-linkage_map_qtl[,3]


	offspringasli_1<-rbind(pos,offspring)
	index<-posqtl
	index_qtls <-which(offspringasli_1[1,]%in%index)

	qtlsel1<-offspringasli_1[,index_qtls]
	qtlsel1<-qtlsel1[-1,]
	qtlsel1<-as.matrix(qtlsel1)

	
	#Extract mrk_genome for each offspring from total sequence based on positions
	nmrk_allele<-dim(linkage_map_mrk)[1]*2
	mrksel1 <- matrix(ncol=nmrk_allele,nrow = length(id))

	pos<-link_map[,3]
	posmrk<-linkage_map_mrk[,3]

	offspringasli_1<-rbind(pos,offspring)
	index<-posmrk
	index_mrks <-which(offspringasli_1[1,]%in%index)

	mrksel1<-offspringasli_1[,index_mrks]
	mrksel1<-mrksel1[-1,]
	mrksel1<-as.matrix(mrksel1)

# # delete row below if doesnt work
  # offspringasli_1<-as.matrix(offspringasli_1)
  
qtlLoci<-qtlsel1
qtlLoci<-as.matrix(qtlLoci)

tbv<-matrix(ncol=4,nrow=length(id))
env<-matrix(ncol=4,nrow=length(id))
phen<-matrix(ncol=4,nrow=length(id))

add_eff_1<-qtl_effects
dom_eff<-rep(0,length(add_eff_1[,1]))
var_dom<-0

	for (j in 1:4){
		# TBV 
		tbv[,j]<-calc_TBV2(qtlLoci,add_eff_1[,j])
		
		# ENV 
		vv<-phen_var[j]-(var_add[j]+var_dom)
		env[,j]<-rnorm(length(tbv[,j]),mean=0,sd=sqrt(vv))
			
		# Phen
		phen[,j]<-tbv[,j]+env[,j]
	}


data<-data.frame(id,sire,dam,sex,F2_id,
phen,tbv,env)
names(data)<-c('id','sire','dam','sex','Syn3_No',
'phen_1','phen_2','phen_3','phen_4',
'tbv_1','tbv_2','tbv_3','tbv_4',
'env_1','env_2','env_3','env_4')

sequ<-offspring

#adding id as first colomn to qtl,mrk and seque
qtlsel1<-cbind(id,qtlsel1)
qtlsel1<-as.matrix(qtlsel1)

mrksel1<-cbind(id,mrksel1)
mrksel1<-as.matrix(mrksel1)

sequ<-cbind(id,sequ)
sequ<-as.matrix(sequ)

internal_Syn3[[1]]<-data
internal_Syn3[[2]]<-qtlsel1
internal_Syn3[[3]]<-mrksel1
internal_Syn3[[4]]<-sequ

	return(internal_Syn3)
}


# Function to use BGLR to estimate effects	
Estimate_effects<-function(snp_reference,phen_reference){
	
y<-phen_reference
X<-snp_reference
y<-as.matrix(y)
X<-as.matrix(X)

X1<-snp_reference


#2  Setting the linear predictor
ETA<-list( 
list(X=X1, model='BRR')
)

#3  Fitting the model
answers<-BGLR(y=y,ETA=ETA, nIter=20000, burnIn=5000,verbose=FALSE)
G11<-answers$ETA[[1]]$b

ghat1<-G11
length(ghat1)
dhat<-rep(0,length(ghat1))
ghat2<-ghat1*-1
return(data.frame(ghat1,ghat2))

	}
	
	
	
Variant_Global<-function(n_all_var,nvar){

		# We want to start this set with the following varities:	
		index<-sample(1:n_all_var,nvar,rep=FALSE)
		file_V_data_set<-file_V_data[index]
		file_V_qtl_set<-file_V_qtl[index]
		file_V_mrk_set<-file_V_mrk[index]
		
		list_fun<-function() {
		list(list(), list(),list(), list())
		}

		Varit_M<-replicate(nvar, list_fun(), simplify=FALSE)

			
		for (i in 1:nvar){

			# Data
			a<-read.table(file_V_data_set[i],skip=1,
			colClasses = c(rep("integer", 3),'character',rep("NULL", 9)))
			
			a_1<-cbind(a,rep(i,length(a[,1])))
			names(a_1)<-c('id','sire','dam','sex','Varit_No')

			# QTL
			b<-read.table(file_V_qtl_set[i],skip=1)
			QTLs<-b[,-1]
			b<-as.matrix(b)
			QTLs<-as.matrix(QTLs)

			# MRK
			M<-read.table(file_V_mrk_set[i],skip=1)
			MRKs<-M[,-1]
			M<-as.matrix(M)
			MRKs<-as.matrix(MRKs)
			dim(MRKs)

			#Sequence
			No_mrk_allel<- dim(M)[2]-1
			No_qtl_allel<- dim(b)[2]-1
			Seq_mat<-matrix(0,ncol=(No_mrk_allel+No_qtl_allel),nrow=dim(M)[1])
			dim(Seq_mat)

			x1<-as.character(link_map[,1])
			x1<-rep(x1,each=2)
			x2<-substr(x1,1,1)
			index_qtl<-which(x2=='Q')
			index_mrk<-which(x2=='M')

			Seq_mat[,index_qtl]<-QTLs
			Seq_mat[,index_mrk]<-MRKs
			Seq_mat<-cbind(b[,1],Seq_mat)
			Seq_mat<-as.matrix(Seq_mat)

			Varit_M[[i]][[1]]<-a_1
			Varit_M[[i]][[2]]<-b
			Varit_M[[i]][[3]]<-M
			Varit_M[[i]][[4]]<-Seq_mat
			
		 cat('Reading Variety',i,"is done",fill=TRUE)
			
		}



		for (g_index in 1:length(Varit_M)){
		names(Varit_M[[g_index]])<-c('data','qtl','mrk','sequ')
		}


		Varit<-Varit_M
		return(Varit)

}



F1_Global<-function(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var){

F1<-list()

for (counter in 1:No_F1_Family){

a1<-sample(1:nvar,1,replace=TRUE)
a2<-sample(1:nvar,1,replace=TRUE)
while(a1==a2){
a2<-sample(1:nvar,1,replace=TRUE)
}

# Variety x
	id_sel_x<-sample(Varit[[a1]]$data[,1],
	nselected_each_variti,replace=FALSE)
	 
	# sequence selected individuals
	index<-match(id_sel_x,Varit[[a1]]$sequ[,1])
	v1<-Varit[[a1]]$sequ[index,]
	v1<-v1[-1] #remove 1st colomn that is ID

# Variety Y
	id_sel_y<-sample(Varit[[a2]]$data[,1],
	nselected_each_variti,replace=FALSE)
	 
	# sequence selected individuals
	index<-match(id_sel_y,Varit[[a2]]$sequ[,1])
	v2<-Varit[[a2]]$sequ[index,]
	v2<-v2[-1]

sires<-id_sel_x
dams<-id_sel_y

noff<-noff_m_v

F1[[counter]]<-Create_F1_DLF(v1,v2,F1_No=counter,sires,dams,qtl_effects,
noff,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
cat('Family F1:',counter,'is done',fill=TRUE)



}

for (g_index in 1:length(F1)){
names(F1[[g_index]])<-c('data','qtl','mrk','sequ')
}

 return(F1)

}



F2_Global<-function(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var){


F2<-list()
for (counter in 1:given_L){

x1<-F1[[counter]][[1]][,1]
s1<-sample(x1,size_F2,replace=TRUE)
s2<-sample(x1,size_F2,replace=TRUE)
testi<-s1==s2
amount<-which(testi==TRUE)
amount

while (length(amount)>0) {
s1<-sample(x1,size_F2,replace=TRUE)
s2<-sample(x1,size_F2,replace=TRUE)
testi<-s1==s2
amount<-which(testi==TRUE)
}

testi<-s1==s2
table(testi)

	# sequence selected individuals
	v1<-F1[[counter]][[4]][s1,]
    dim(v1)

	# sequence selected individuals
	v2<-F1[[counter]][[4]][s2,]
    dim(v2)

sires<-s1
dams<-s2

F2[[counter]]<-Create_F2(v1,v2,F2_No=counter,sires,dams,qtl_effects,
noff=noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

cat('Family F2:',counter,'is done',fill=TRUE)

}

for (g_index in 1:length(F2)){
names(F2[[g_index]])<-c('data','qtl','mrk','sequ')
}


return(F2)



}

	
F2_SP_Global<-function(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var){

#######################################
# CReation of F2_SP
######################################
	F2_SP<-list()
	for (counter in 1:given_L){

	x1<-F1_SP[[counter]][[1]][,1]
	s1<-sample(x1,size_F2_SP,replace=TRUE)
	s2<-sample(x1,size_F2_SP,replace=TRUE)
	testi<-s1==s2
	amount<-which(testi==TRUE)
	amount

	while (length(amount)>0) {
	s1<-sample(x1,size_F2_SP,replace=TRUE)
	s2<-sample(x1,size_F2_SP,replace=TRUE)
	testi<-s1==s2
	amount<-which(testi==TRUE)
	}

	testi<-s1==s2
	table(testi)

		# sequence selected individuals
		v1<-F1_SP[[counter]][[4]][s1,]
		dim(v1)

		# sequence selected individuals
		v2<-F1_SP[[counter]][[4]][s2,]
		dim(v2)

	sires<-s1
	dams<-s2

	F2_SP[[counter]]<-Create_F2(v1,v2,F2_No=counter,sires,dams,qtl_effects,
	noff=noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
	,chrom_length,var_add,phen_var)

	cat('Family F2_SP:',counter,'is done',fill=TRUE)

	}

# Sets[[Cyl]][[3]]<-F2_SP

	F2_SP[[1]][1]
	for (g_index in 1:length(F2_SP)){
	names(F2_SP[[g_index]])<-c('data','qtl','mrk','sequ')
	}

	return(F2_SP)

}



Global_Syn1<-function(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var){
# get mean of F2s for yeild
F2s_info<-matrix(nrow=length(F2),ncol=5)
phen_tr1<-c()
phen_tr3<-c()
phen_tr4<-c()

for (i in 1:length(F2)){
phen_tr1<-mean(F2[[i]][[1]]$phen_1)
phen_tr3<-mean(F2[[i]][[1]]$phen_3)
phen_tr4<-mean(F2[[i]][[1]]$phen_4)
shakhes<-(b_tr1*phen_tr1)+(b_tr3*phen_tr3)+(b_tr4*phen_tr4)

F2s_info[i,1]<-i
F2s_info[i,2]<-phen_tr1
F2s_info[i,3]<-phen_tr3
F2s_info[i,4]<-phen_tr4
F2s_info[i,5]<-shakhes
}

F2s_sorted<-F2s_info[order(-F2s_info[,5]),]
F2s_selected<-F2s_sorted[1:No_F2_Family_selected,]


# Select (size_polycross) plants randomly from each selected F2_SP for polycross

# # 
list_fun<-function() {
    list(list(), list(),list(), list())
}

 clonal_row<-replicate(nClonalrow, list_fun(), simplify=FALSE)
 

# version 2
as<-c()
as_list<-list()
counter<-1
for (j in 1:length(F2s_selected[,1])){

	z<-F2_SP[[F2s_selected[j,1]]][[1]]
	as<-sample(z[,1],size_polycross,replace=FALSE)
	as_list[[j]]<-as

	for (p in 1:length(as)){
			clonal_row[[counter]][[1]]<-F2_SP[[F2s_selected[j,1]]][[1]][as[p],] #data
			clonal_row[[counter]][[2]]<-F2_SP[[F2s_selected[j,1]]][[2]][as[p],] #qtl
			clonal_row[[counter]][[3]]<-F2_SP[[F2s_selected[j,1]]][[3]][as[p],] #mrk
			clonal_row[[counter]][[4]]<-F2_SP[[F2s_selected[j,1]]][[4]][as[p],] #sequ
			counter<-counter+1	
	}
		

}



#----------------------
# renaming id of Clonal Rows Start
#----------------------
for (i in 1:length(clonal_row)){
clonal_row[[i]][[1]]$id<-i
}
#----------------------
# renaming id of Clonal Rows Finish
#----------------------

#----------------------------------------
# GROUPING CLONAL ROWS (Selected SP from F2_SP)
#----------------------------------------
sorted_CR<-c()
for (i in 1:length(clonal_row)){
a<-clonal_row[[i]][[1]]
sorted_CR<-rbind(sorted_CR,a)
}

# 7 means phen_2 that is FT
sorted_CR<-sorted_CR[order(-sorted_CR[,7]),]
# sorted_CR[,1]

ngroups<-length(clonal_row)/size_polycross
sh1<-floor(ngroups/3) #3 is constant means FT low, Medi and high
sh2<-sh1
sh3<-ngroups-(sh1+sh2)
sh<-c(sh1,sh3,sh2)
sh
sum(sh)


abi<-sh*size_polycross
abi2<-cumsum(abi)
se1<-c(1,abi2[-3]+1)
se2<-abi2
se1
se2

sorted_CR1<-sorted_CR[se1[1]:se2[1],]
sorted_CR2<-sorted_CR[se1[2]:se2[2],]
sorted_CR3<-sorted_CR[se1[3]:se2[3],]

sorted_CR1_Sh <- sorted_CR1[sample(nrow(sorted_CR1)),]
sorted_CR2_Sh <- sorted_CR2[sample(nrow(sorted_CR2)),]
sorted_CR3_Sh <- sorted_CR3[sample(nrow(sorted_CR3)),]

# from each CR we sample No.=sh and size of each = size_polycross

# Low

se<-seq(1,length(sorted_CR1_Sh[,1]),size_polycross)
se1<-se
se2<-se-1
se2<-se2[-1]
se2<-c(se2,length(sorted_CR1_Sh[,1]))
se1
se2

group_ids_Low<-list()

	for (i in 1:sh[1]){
	group_ids_Low[[i]]<-sorted_CR1_Sh[se1[i]:se2[i],1]
	}
	
	
# Medium
se<-seq(1,length(sorted_CR2_Sh[,1]),size_polycross)
se1<-se
se2<-se-1
se2<-se2[-1]
se2<-c(se2,length(sorted_CR2_Sh[,1]))
se1
se2

group_ids_Med<-list()

	for (i in 1:sh[2]){
	group_ids_Med[[i]]<-sorted_CR2_Sh[se1[i]:se2[i],1]
	}
	
	# high
se<-seq(1,length(sorted_CR3_Sh[,1]),size_polycross)
se1<-se
se2<-se-1
se2<-se2[-1]
se2<-c(se2,length(sorted_CR3_Sh[,1]))
se1
se2

group_ids_hig<-list()

	for (i in 1:sh[3]){
	group_ids_hig[[i]]<-sorted_CR3_Sh[se1[i]:se2[i],1]
	}
	
	
	
output1<-matrix(unlist(group_ids_Low), ncol = 8, byrow = TRUE)
output2<-matrix(unlist(group_ids_Med), ncol = 8, byrow = TRUE)
output3<-matrix(unlist(group_ids_hig), ncol = 8, byrow = TRUE)

output_CR<-rbind(output1,output2,output3)	

ngroups<-length(clonal_row)/size_polycross

group_ids<-list()
for (i in 1:ngroups){
group_ids[[i]]<-output_CR[i,]
}

a<-sapply(clonal_row, "[[", c(1))
idis<-unlist(a[1,]) 


# Grouping clonal rows to (size_polycross)-parents based on similarity of FT
Four_parent_groups<-list()
for (i in 1:ngroups){
index<-match(group_ids[[i]],idis)
Four_parent_groups[[i]]<-clonal_row[index]
}


###################################
# CREATING MINI PLOTS (SYN1)
###################################
nplots<-length(Four_parent_groups)
mini_plots<-list()


	for (counter in 1:nplots){

	SYN_0<-Four_parent_groups[[counter]]
	Syn_id<-counter

	khor<-Create_Syn1(size_polycross,SYN_0,Syn_id,qtl_effects,
	noff=noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
	,chrom_length,var_add,phen_var)

	mini_plots[[counter]]<-khor

	 cat('SYN1 Number:',counter,'is done',fill=TRUE)
	}



for (g_index in 1:length(mini_plots)){
names(mini_plots[[g_index]])<-c('data','qtl','mrk','sequ')
}

return(mini_plots)

}


Global_Syn1_G<-function(F2_SP,nClonalrow,size_polycross,
eff_tr1,eff_tr3,eff_tr4,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var){



#################################
# Calc GEBV in Sp
#################################
	# Calc GEBV for tr1

mat_tr1<-c()
for (zp in 1:length(F2_SP)){
a<-F2_SP[[zp]][[3]]
mat_tr1<-rbind(mat_tr1,a)
}

ghat1<-eff_tr1[,1]
ghat2<-eff_tr1[,2]
dhat<-rep(0,length(ghat1))

gebv_tr1<-calc_GEBV2(mat_tr1,ghat1)

	# TRAIT 3
ghat1<-eff_tr3[,1]
ghat2<-eff_tr3[,2]

	# Calc GEBV for tr3
mat_tr3<-c()
for (zp in 1:length(F2_SP)){
a<-F2_SP[[zp]][[3]]
mat_tr3<-rbind(mat_tr3,a)
}
gebv_tr3<-calc_GEBV2(mat_tr3,ghat1)


	# TRAIT 4
ghat1<-eff_tr4[,1]
ghat2<-eff_tr4[,2]

	# Calc GEBV for tr4
mat_tr4<-c()
for (zp in 1:length(F2_SP)){
a<-F2_SP[[zp]][[3]]
mat_tr4<-rbind(mat_tr4,a)
}
gebv_tr4<-calc_GEBV2(mat_tr4,ghat1)
gebv_tr1<-as.numeric(gebv_tr1)
gebv_tr3<-as.numeric(gebv_tr3)
gebv_tr4<-as.numeric(gebv_tr4)



GS_data<-matrix(0,nrow=	length(gebv_tr1),ncol=6)
GS_data[,1]<-rep(1:dim(F2_SP[[1]][[1]])[1],length(F2_SP))
GS_data[,2]<-rep(1:length(F2_SP),each=dim(F2_SP[[1]][[1]])[1])
GS_data[,3]<-gebv_tr1
GS_data[,4]<-gebv_tr3
GS_data[,5]<-gebv_tr4
GS_data[,6]<-(b_tr1*gebv_tr1)+(b_tr3*gebv_tr3)+(b_tr4*gebv_tr4)

GS_data2<-GS_data[order(-GS_data[,6]),]
GS_data2_Selected<-GS_data2[1:nClonalrow,]


# # 
list_fun<-function() {
    list(list(), list(),list(), list())
}

 clonal_row<-replicate(nClonalrow, list_fun(), simplify=FALSE)
 

for (j in 1:length(clonal_row)){
clonal_row[[j]][[1]]<-F2_SP[[GS_data2_Selected[j,2]]][[1]][GS_data2_Selected[j,1],] #data
clonal_row[[j]][[2]]<-F2_SP[[GS_data2_Selected[j,2]]][[2]][GS_data2_Selected[j,1],] #qtl
clonal_row[[j]][[3]]<-F2_SP[[GS_data2_Selected[j,2]]][[3]][GS_data2_Selected[j,1],] #mrk
clonal_row[[j]][[4]]<-F2_SP[[GS_data2_Selected[j,2]]][[4]][GS_data2_Selected[j,1],] #sequ
}




#----------------------
# renaming id of Clonal Rows Start
#----------------------
for (i in 1:length(clonal_row)){
clonal_row[[i]][[1]]$id<-i
}

#----------------------
# renaming id of Clonal Rows Finish
#----------------------


#----------------------------------------
# GROUPING CLONAL ROWS (Selected SP from F2_SP)
#----------------------------------------
sorted_CR<-c()
for (i in 1:length(clonal_row)){
a<-clonal_row[[i]][[1]]
sorted_CR<-rbind(sorted_CR,a)
}

varTr1<-var(sorted_CR[,10])
varTr2<-var(sorted_CR[,11])
varTr3<-var(sorted_CR[,12])
varTr4<-var(sorted_CR[,13])
variances<-c(varTr1,varTr2,varTr3,varTr4)



# 7 means phen_2 that is FT
sorted_CR<-sorted_CR[order(-sorted_CR[,7]),]

ngroups<-length(clonal_row)/size_polycross
sh1<-floor(ngroups/3) #3 is constant means FT low, Medi and high
sh2<-sh1
sh3<-ngroups-(sh1+sh2)
sh<-c(sh1,sh3,sh2)
sh
sum(sh)


abi<-sh*size_polycross
abi2<-cumsum(abi)
se1<-c(1,abi2[-3]+1)
se2<-abi2
se1
se2

sorted_CR1<-sorted_CR[se1[1]:se2[1],]
sorted_CR2<-sorted_CR[se1[2]:se2[2],]
sorted_CR3<-sorted_CR[se1[3]:se2[3],]

sorted_CR1_Sh <- sorted_CR1[sample(nrow(sorted_CR1)),]
sorted_CR2_Sh <- sorted_CR2[sample(nrow(sorted_CR2)),]
sorted_CR3_Sh <- sorted_CR3[sample(nrow(sorted_CR3)),]

# from each CR we sample No.=sh and size of each = size_polycross
# Low

se<-seq(1,length(sorted_CR1_Sh[,1]),size_polycross)
se1<-se
se2<-se-1
se2<-se2[-1]
se2<-c(se2,length(sorted_CR1_Sh[,1]))
se1
se2

group_ids_Low<-list()

	for (i in 1:sh[1]){
	group_ids_Low[[i]]<-sorted_CR1_Sh[se1[i]:se2[i],1]
	}
	
	
# Medium
se<-seq(1,length(sorted_CR2_Sh[,1]),size_polycross)
se1<-se
se2<-se-1
se2<-se2[-1]
se2<-c(se2,length(sorted_CR2_Sh[,1]))
se1
se2

group_ids_Med<-list()

	for (i in 1:sh[2]){
	group_ids_Med[[i]]<-sorted_CR2_Sh[se1[i]:se2[i],1]
	}
	
	# high
se<-seq(1,length(sorted_CR3_Sh[,1]),size_polycross)
se1<-se
se2<-se-1
se2<-se2[-1]
se2<-c(se2,length(sorted_CR3_Sh[,1]))
se1
se2

group_ids_hig<-list()

	for (i in 1:sh[3]){
	group_ids_hig[[i]]<-sorted_CR3_Sh[se1[i]:se2[i],1]
	}
	
	
	
output1<-matrix(unlist(group_ids_Low), ncol = 8, byrow = TRUE)
output2<-matrix(unlist(group_ids_Med), ncol = 8, byrow = TRUE)
output3<-matrix(unlist(group_ids_hig), ncol = 8, byrow = TRUE)

output_CR<-rbind(output1,output2,output3)	

ngroups<-length(clonal_row)/size_polycross

group_ids<-list()
for (i in 1:ngroups){
group_ids[[i]]<-output_CR[i,]
}


##############################
# CREATING MINI PLOTS (SYN0)
##############################

a<-sapply(clonal_row, "[[", c(1))
idis<-unlist(a[1,]) 


# Grouping clonal rows to (size_polycross)-parents based on similarity of FT
Four_parent_groups<-list()
for (i in 1:ngroups){
index<-match(group_ids[[i]],idis)
Four_parent_groups[[i]]<-clonal_row[index]
}


######################################
# CREATING MINI PLOTS (SYN1)
######################################
nplots<-length(Four_parent_groups)
mini_plots<-list()


for (counter in 1:nplots){

SYN_0<-Four_parent_groups[[counter]]
Syn_id<-counter

khor<-Create_Syn1(size_polycross,SYN_0,Syn_id,qtl_effects,
noff=noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

mini_plots[[counter]]<-khor
 cat('SYN1 Number:',counter,'is done',fill=TRUE)
}


for (g_index in 1:length(mini_plots)){
names(mini_plots[[g_index]])<-c('data','qtl','mrk','sequ')
}

out_list<-list(mini_plots,variances)
return(out_list)
	
} #End of function



Global_Syn2<-function(mini_plots,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var){


######################################
# CREATING MINI PLOTS (SYN2)
######################################

Syn_2<-list()

for (counter in 1:length(mini_plots)){

x1<-mini_plots[[counter]][[1]][,1]
s1<-sample(x1,size_Syn1,replace=TRUE)
s2<-sample(x1,size_Syn1,replace=TRUE)
testi<-s1==s2
amount<-which(testi==TRUE)
amount

while (length(amount)>0) {
s1<-sample(x1,size_Syn1,replace=TRUE)
s2<-sample(x1,size_Syn1,replace=TRUE)
testi<-s1==s2
amount<-which(testi==TRUE)
}

testi<-s1==s2
table(testi)

	# sequence selected individuals
	v1<-mini_plots[[counter]][[4]][s1,]
    dim(v1)

	# sequence selected individuals
	v2<-mini_plots[[counter]][[4]][s2,]
    dim(v2)

sires<-s1
dams<-s2

Syn_2[[counter]]<-Create_Syn2(v1,v2,Syn2_No=counter,sires,dams,qtl_effects,
noff=noff_per_each_mating_in_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

cat('Syn2 number:',counter,'is done',fill=TRUE)

}


for (g_index in 1:length(Syn_2)){
names(Syn_2[[g_index]])<-c('data','qtl','mrk','sequ')
}


return(Syn_2)


}


Global_Syn3<-function(Syn_2,size_Syn3,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var){


SI<-c()
selected_syn2_id<-matrix(0,ncol=2,nrow=length(Syn_2))

	for (i in 1:length(Syn_2)){
		# tr1
		mean_tr1<-mean(Syn_2[[i]]$data$phen_1)
		
		# tr3
		mean_tr3<-mean(Syn_2[[i]]$data$phen_3)
		# tr4
		mean_tr4<-mean(Syn_2[[i]]$data$phen_4)

		SI[i]<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
		selected_syn2_id[i,1]<-i
		selected_syn2_id[i,2]<-SI[i]
	}



selected_syn2_id<-selected_syn2_id[order(-selected_syn2_id[,2]),]
selected_syn2_id<-selected_syn2_id[1:No_Selected_Syn2,1]
Syn_2_selected<-Syn_2[selected_syn2_id]

#--------------------------------
# CREATING MINI PLOTS (SYN3)
#--------------------------------

Syn_3<-list()

for (counter in 1:length(Syn_2_selected)){

x1<-Syn_2_selected[[counter]][[1]][,1]
s1<-sample(x1,size_Syn3,replace=TRUE)
s2<-sample(x1,size_Syn3,replace=TRUE)
testi<-s1==s2
amount<-which(testi==TRUE)
amount

while (length(amount)>0) {
s1<-sample(x1,size_Syn3,replace=TRUE)
s2<-sample(x1,size_Syn3,replace=TRUE)
testi<-s1==s2
amount<-which(testi==TRUE)
}

testi<-s1==s2


	# sequence selected individuals
	v1<-Syn_2_selected[[counter]][[4]][s1,]
	# sequence selected individuals
	v2<-Syn_2_selected[[counter]][[4]][s2,]


sires<-s1
dams<-s2

Syn_3[[counter]]<-Create_Syn3(v1,v2,Syn3_No=counter,sires,dams,qtl_effects,
noff=noff_per_each_mating_in_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

cat('Syn3 number:',counter,'is done',fill=TRUE)

}


for (g_index in 1:length(Syn_3)){
names(Syn_3[[g_index]])<-c('data','qtl','mrk','sequ')
}

return(Syn_3)


}


give_gebv_f2sp<-function(MIGIRAM){


taxi<-list()

for (upi in 1:length(MIGIRAM)){

F2_SP<-MIGIRAM[[upi]]

	mat_tr1<-c()
	for (zp in 1:length(F2_SP)){
	a<-F2_SP[[zp]][[3]]
	mat_tr1<-rbind(mat_tr1,a)
	}
	dim(mat_tr1)

	ghat1<-eff_tr1[,1]
	ghat2<-eff_tr1[,2]
	dhat<-rep(0,length(ghat1))

	gebv_tr1<-calc_GEBV2(mat_tr1,ghat1)
	summary(gebv_tr1)

		# TRAIT 3
	ghat1<-eff_tr3[,1]
	ghat2<-eff_tr3[,2]

		# Calc GEBV for tr3
	mat_tr3<-c()
	for (zp in 1:length(F2_SP)){
	a<-F2_SP[[zp]][[3]]
	mat_tr3<-rbind(mat_tr3,a)
	}
	dim(mat_tr3)
	gebv_tr3<-calc_GEBV2(mat_tr3,ghat1)


		# TRAIT 4
	ghat1<-eff_tr4[,1]
	ghat2<-eff_tr4[,2]

		# Calc GEBV for tr4
	mat_tr4<-c()
	for (zp in 1:length(F2_SP)){
	a<-F2_SP[[zp]][[3]]
	mat_tr4<-rbind(mat_tr4,a)
	}
		dim(mat_tr4)
	gebv_tr4<-calc_GEBV2(mat_tr4,ghat1)

	# GEBV
	gebv_tr1<-as.numeric(gebv_tr1)
	gebv_tr3<-as.numeric(gebv_tr3)
	gebv_tr4<-as.numeric(gebv_tr4)

	taxi[[upi]]<-list(gebv_tr1,gebv_tr3,gebv_tr4)
	cat('Analyses No:',upi,'Done',fill=TRUE)
}


return(taxi)

}


Give_gebv_F2<-function(MIGIRAM,nmrk,N_S_V,eff_tr1,eff_tr3,eff_tr4){

F2<-MIGIRAM

# # Sampling plants from each F2 to calculate GEBV 
Ref<-matrix(ncol=nmrk,nrow=length(F2))

	for (counter in 1:length(F2)){
	index<-sample(F2[[counter]]$data[,1],N_S_V,replace=FALSE)
	Gen_F2<-F2[[counter]]$mrk[index,]

		freq1<-Calc_freq(Gen_F2)
		freq2<-1-freq1
		alle_dosage<-freq1*2
		Ref[counter,]<-alle_dosage
	}

	
# Calc GEBV for tr1
	ghat1<-eff_tr1[,1]
	ghat2<-eff_tr1[,2]
	gebv_tr1<-Ref%*%ghat1

# Calc GEBV for tr3
	ghat1<-eff_tr3[,1]
	ghat2<-eff_tr3[,2]
	gebv_tr3<-Ref%*%ghat1
	
# Calc GEBV for tr4
	ghat1<-eff_tr4[,1]
	ghat2<-eff_tr4[,2]
	gebv_tr4<-Ref%*%ghat1
	
gebv_tr1<-as.numeric(gebv_tr1)
gebv_tr3<-as.numeric(gebv_tr3)
gebv_tr4<-as.numeric(gebv_tr4)

	taxi<-list(gebv_tr1,gebv_tr3,gebv_tr4)
return(taxi)

}



Family_Maker<-function(MIGIRAM,nvar,GS_data2_Selected){


list_fun<-function() {
    list(list(), list(),list(), list())
}

 tempi_family<-replicate(nvar, list_fun(), simplify=FALSE)
 
 
GS_da<-cbind(rep(1:nvar,each=10),GS_data2_Selected)

 
for (pedal in 1:nvar){

tily<-MIGIRAM[[1]]
DQ<-matrix(nrow=10,ncol=dim(tily[[1]][[2]])[2])
DM<-matrix(nrow=10,ncol=dim(tily[[1]][[3]])[2])
DS<-matrix(nrow=10,ncol=dim(tily[[1]][[4]])[2])
dim(DS)

aban<-subset(GS_da,GS_da[,1]==pedal)
aban<-aban[,-1]

Set_index<-aban[,1]
fam_index<-aban[,3]
indi_index<-aban[,2]

D2<-c()
for (j in 1:dim(aban)[1]){
F2_SP<-MIGIRAM[[Set_index[j]]]
eshih<-F2_SP[[fam_index[j]]][[1]][indi_index[j],] #data
D2<-rbind(D2,eshih)
DQ[j,]<-F2_SP[[fam_index[j]]][[2]][indi_index[j],] #qtl
DM[j,]<-F2_SP[[fam_index[j]]][[3]][indi_index[j],] #mrk
DS[j,]<-F2_SP[[fam_index[j]]][[4]][indi_index[j],] #sequ
}


tempi_family[[pedal]][[1]]<-D2
tempi_family[[pedal]][[2]]<-DQ
tempi_family[[pedal]][[3]]<-DM
tempi_family[[pedal]][[4]]<-DS
}



for (g_index in 1:length(tempi_family)){
names(tempi_family[[g_index]])<-c('data','qtl','mrk','sequ')
}

	# Make tempi_family similar output of Syn2 or 3
	for (u in 1:length(tempi_family)){
			
			localID<-tempi_family[[u]]$data[,1]

			# QTL
			x_qtl<-tempi_family[[u]][[2]]
			x_qtl<-cbind(localID,x_qtl)
			tempi_family[[u]][[2]]<-x_qtl

			# MRK
			x_mrk<-tempi_family[[u]][[3]]
			x_mrk<-cbind(localID,x_mrk)
			tempi_family[[u]][[3]]<-x_mrk

			# SEQ
			x_seq<-tempi_family[[u]][[4]]
			x_seq<-cbind(localID,x_seq)
			tempi_family[[u]][[4]]<-x_seq
	}
	
	return(tempi_family)

	
} #End of function


give_gebv_Syn2<-function(DAHANDEH,nmrk,N_S_V,eff_tr1,eff_tr3,eff_tr4){

taxi<-list()

for (upi in 1:length(DAHANDEH)){

Syn_2<-DAHANDEH[[upi]]
# # Sampling plants from each Syn2 to calculate GEBV for tr3 &4

Ref<-matrix(ncol=nmrk,nrow=length(Syn_2))
dim(Ref)

	for (counter in 1:length(Syn_2)){
	index<-sample(Syn_2[[counter]]$data[,1],N_S_V,replace=FALSE)
	index2<-match(index,Syn_2[[counter]]$mrk[,1])
	Gen_syn2<-Syn_2[[counter]]$mrk[index2,]
	Gen_syn2<-Gen_syn2[,-1]

		freq1<-Calc_freq(Gen_syn2)
		freq2<-1-freq1
		alle_dosage<-freq1*2
		Ref[counter,]<-alle_dosage
	}
	
	
# Calc GEBV for tr1
	ghat1<-eff_tr1[,1]
	ghat2<-eff_tr1[,2]
	gebv_tr1<-Ref%*%ghat1

# Calc GEBV for tr3
	ghat1<-eff_tr3[,1]
	ghat2<-eff_tr3[,2]
	gebv_tr3<-Ref%*%ghat1
	
# Calc GEBV for tr4
	ghat1<-eff_tr4[,1]
	ghat2<-eff_tr4[,2]
	gebv_tr4<-Ref%*%ghat1

	# # GEBV
gebv_tr1<-as.numeric(gebv_tr1)
gebv_tr3<-as.numeric(gebv_tr3)
gebv_tr4<-as.numeric(gebv_tr4)


	 taxi[[upi]]<-list(gebv_tr1,gebv_tr3,gebv_tr4)
	 cat('Analyses No:',upi,'Syn2 Done',fill=TRUE)
}


return(taxi)

}



Family_Extract_Selected_Syn2<-function(DAHANDEH,nvar,GS_data2_Selected){

		list_fun<-function() {
		list(list(), list(),list(), list())
		}

Update_data<-replicate(nvar, list_fun(), simplify=FALSE)
evalu<-GS_data2_Selected




			for(upd in 1:dim(evalu)[1]){
			set_number<-evalu[upd,1]
			Syn_number<-evalu[upd,2]
			Update_data[[upd]][[1]]<-DAHANDEH[[set_number]][[Syn_number]][[1]]
			Update_data[[upd]][[2]]<-DAHANDEH[[set_number]][[Syn_number]][[2]]
			Update_data[[upd]][[3]]<-DAHANDEH[[set_number]][[Syn_number]][[3]]
			Update_data[[upd]][[4]]<-DAHANDEH[[set_number]][[Syn_number]][[4]]
			}
		
	        tempi_family<-Update_data	
		
		for (g_index in 1:length(tempi_family)){
		names(tempi_family[[g_index]])<-c('data','qtl','mrk','sequ')
		}
		
	return(tempi_family)	
}


Give_gebv_Syn2_Single_Set<-function(DAHANDEH,nmrk,N_S_V,eff_tr1,eff_tr3,eff_tr4){

Syn_2<-DAHANDEH


Ref<-matrix(ncol=nmrk,nrow=length(Syn_2))
dim(Ref)

	for (counter in 1:length(Syn_2)){
	index<-sample(Syn_2[[counter]]$data[,1],N_S_V,replace=FALSE)
	index2<-match(index,Syn_2[[counter]]$mrk[,1])
	Gen_syn2<-Syn_2[[counter]]$mrk[index2,]
	Gen_syn2<-Gen_syn2[,-1]

		freq1<-Calc_freq(Gen_syn2)
		freq2<-1-freq1
		alle_dosage<-freq1*2
		Ref[counter,]<-alle_dosage
	}
	
	
# Calc GEBV for tr1
	ghat1<-eff_tr1[,1]
	ghat2<-eff_tr1[,2]
	gebv_tr1<-Ref%*%ghat1

# Calc GEBV for tr3
	ghat1<-eff_tr3[,1]
	ghat2<-eff_tr3[,2]
	gebv_tr3<-Ref%*%ghat1
	
# Calc GEBV for tr4
	ghat1<-eff_tr4[,1]
	ghat2<-eff_tr4[,2]
	gebv_tr4<-Ref%*%ghat1

	# # GEBV
gebv_tr1<-as.numeric(gebv_tr1)
gebv_tr3<-as.numeric(gebv_tr3)
gebv_tr4<-as.numeric(gebv_tr4)


	taxi<-list(gebv_tr1,gebv_tr3,gebv_tr4)
return(taxi)

}



Global_Syn3_G<-function(Syn_2,size_Syn3,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var){


Syn_2_selected<-Syn_2

Syn_3<-list()

for (counter in 1:length(Syn_2_selected)){

x1<-Syn_2_selected[[counter]][[1]][,1]
s1<-sample(x1,size_Syn3,replace=TRUE)
s2<-sample(x1,size_Syn3,replace=TRUE)
testi<-s1==s2
amount<-which(testi==TRUE)
amount

while (length(amount)>0) {
s1<-sample(x1,size_Syn3,replace=TRUE)
s2<-sample(x1,size_Syn3,replace=TRUE)
testi<-s1==s2
amount<-which(testi==TRUE)
}

testi<-s1==s2


	# sequence selected individuals
	v1<-Syn_2_selected[[counter]][[4]][s1,]
	# sequence selected individuals
	v2<-Syn_2_selected[[counter]][[4]][s2,]


sires<-s1
dams<-s2

Syn_3[[counter]]<-Create_Syn3(v1,v2,Syn3_No=counter,sires,dams,qtl_effects,
noff=noff_per_each_mating_in_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

cat('Syn3 number:',counter,'is done',fill=TRUE)

}


for (g_index in 1:length(Syn_3)){
names(Syn_3[[g_index]])<-c('data','qtl','mrk','sequ')
}

return(Syn_3)

	
} #End of function



Syn2_SP_Global<-function(Syn_SP,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var){


######################################
# CREATING MINI PLOTS (SYN2)
######################################

Syn_2_SP<-list()

for (counter in 1:length(Syn_SP)){

x1<-Syn_SP[[counter]][[1]][,1]
s1<-sample(x1,size_Syn1,replace=TRUE)
s2<-sample(x1,size_Syn1,replace=TRUE)
testi<-s1==s2
amount<-which(testi==TRUE)
amount

while (length(amount)>0) {
s1<-sample(x1,size_Syn1,replace=TRUE)
s2<-sample(x1,size_Syn1,replace=TRUE)
testi<-s1==s2
amount<-which(testi==TRUE)
}

testi<-s1==s2
table(testi)

	# sequence selected individuals
	v1<-Syn_SP[[counter]][[4]][s1,]
    dim(v1)

	# sequence selected individuals
	v2<-Syn_SP[[counter]][[4]][s2,]
    dim(v2)

sires<-s1
dams<-s2

Syn_2_SP[[counter]]<-Create_Syn2(v1,v2,Syn2_No=counter,sires,dams,qtl_effects,
noff=noff_per_each_mating_in_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

cat('Syn2 SP number:',counter,'is done',fill=TRUE)

}


for (g_index in 1:length(Syn_2_SP)){
names(Syn_2_SP[[g_index]])<-c('data','qtl','mrk','sequ')
}


return(Syn_2_SP)


}



give_gebv_Syn2sp<-function(MIGIRAM,eff_tr1,eff_tr3,eff_tr4){


taxi<-list()

for (upi in 1:length(MIGIRAM)){

F2_SP<-MIGIRAM[[upi]]

	mat_tr1<-c()
	for (zp in 1:length(F2_SP)){
	a<-F2_SP[[zp]][[3]]
	mat_tr1<-rbind(mat_tr1,a)
	}
	dim(mat_tr1)
	mat_tr1<-mat_tr1[,-1]

	ghat1<-eff_tr1[,1]
	ghat2<-eff_tr1[,2]
	dhat<-rep(0,length(ghat1))

	gebv_tr1<-calc_GEBV2(mat_tr1,ghat1)
	summary(gebv_tr1)

		# TRAIT 3
	ghat1<-eff_tr3[,1]
	ghat2<-eff_tr3[,2]

		# Calc GEBV for tr3
	mat_tr3<-c()
	for (zp in 1:length(F2_SP)){
	a<-F2_SP[[zp]][[3]]
	mat_tr3<-rbind(mat_tr3,a)
	}
	dim(mat_tr3)
		mat_tr3<-mat_tr3[,-1]
	gebv_tr3<-calc_GEBV2(mat_tr3,ghat1)


		# TRAIT 4
	ghat1<-eff_tr4[,1]
	ghat2<-eff_tr4[,2]

		# Calc GEBV for tr4
	mat_tr4<-c()
	for (zp in 1:length(F2_SP)){
	a<-F2_SP[[zp]][[3]]
	mat_tr4<-rbind(mat_tr4,a)
	}
		dim(mat_tr4)
	mat_tr4<-mat_tr4[,-1]
	gebv_tr4<-calc_GEBV2(mat_tr4,ghat1)

	# GEBV
	gebv_tr1<-as.numeric(gebv_tr1)
	gebv_tr3<-as.numeric(gebv_tr3)
	gebv_tr4<-as.numeric(gebv_tr4)

	taxi[[upi]]<-list(gebv_tr1,gebv_tr3,gebv_tr4)
	cat('Analyses No:',upi,'Done',fill=TRUE)
}


return(taxi)

}



Family_Maker_SynSP<-function(MIGIRAM,nvar,GS_data2_Selected){


list_fun<-function() {
    list(list(), list(),list(), list())
}

 tempi_family<-replicate(nvar, list_fun(), simplify=FALSE)
 
 
GS_da<-cbind(rep(1:nvar,each=10),GS_data2_Selected)

 
for (pedal in 1:nvar){

tily<-MIGIRAM[[1]]
DQ<-matrix(nrow=10,ncol=dim(tily[[1]][[2]])[2])
DM<-matrix(nrow=10,ncol=dim(tily[[1]][[3]])[2])
DS<-matrix(nrow=10,ncol=dim(tily[[1]][[4]])[2])
dim(DS)

aban<-subset(GS_da,GS_da[,1]==pedal)
aban<-aban[,-1]

Set_index<-aban[,1]
fam_index<-aban[,3]
indi_index<-aban[,2]

D2<-c()
for (j in 1:dim(aban)[1]){
F2_SP<-MIGIRAM[[Set_index[j]]]
eshih<-F2_SP[[fam_index[j]]][[1]][indi_index[j],] #data
D2<-rbind(D2,eshih)
DQ[j,]<-F2_SP[[fam_index[j]]][[2]][indi_index[j],] #qtl
DM[j,]<-F2_SP[[fam_index[j]]][[3]][indi_index[j],] #mrk
DS[j,]<-F2_SP[[fam_index[j]]][[4]][indi_index[j],] #sequ
}


tempi_family[[pedal]][[1]]<-D2
tempi_family[[pedal]][[2]]<-DQ
tempi_family[[pedal]][[3]]<-DM
tempi_family[[pedal]][[4]]<-DS
}



for (g_index in 1:length(tempi_family)){
names(tempi_family[[g_index]])<-c('data','qtl','mrk','sequ')
}


return(tempi_family)

	
} #End of function



# ---------------------------------
# Functions------Finish-------------
# ---------------------------------