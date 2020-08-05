

# Load needed packages
library(MASS)
library(BGLR)
library(data.table)
source('AllFunctions.r')



#-----------------------------
# GENERAL INPUTS -------------
#-----------------------------


n_all_var<-20
nvar<-20  #Number of Varities to be sampled from hp
nCycle<-36 #Number of Cycles (Sets)


nRep<-50

# Parameters for traits
# Trait 1:
h2_tr1<-0.03     #Rust
phen_var_tr1<-10

# Trait 2:
h2_tr2<-0.05     #FT
phen_var_tr2<-10

# Trait 3:
h2_tr3<-0.03     #DM
phen_var_tr3<-10

# Trait 4:
h2_tr4<-0.02
phen_var_tr4<-10 #GM

rg_tr34<-0.7  # GEnetic correlation between trait 3&4
#-----------------------------
# Parameters for HP
#-----------------------------
chrom_length<-100

#-----------------------------
# Parameters for Varieties
#-----------------------------
noff_m_v<-10       #20          # No. of offspring per mating 


#-----------------------------
# Parameters for F1
#-----------------------------
No_F1_Family<-250 #250

#-----------------------------
# Parameters for F2
#-----------------------------
noff_F2<-1    # No. of offspring per mating in F1 to produce F2 
size_F2<-40  
size_F2_SP<-40 
No_F2_Family_selected<-No_F1_Family/5 

#-----------------------------
# Parameters for Polycross 
#-----------------------------
size_polycross<-8
nClonalrow<-No_F2_Family_selected*size_polycross
Cloning_size<-1

#-----------------------------
# Parameters for Syn0
#-----------------------------
noff_per_each_syn<-2 # No. of off for each x-parent mating in Syn0 groups
                      # for the moment it should be even number
					  
#-----------------------------
# Parameters for Syn1
#-----------------------------
noff_per_each_mating_in_syn<-1 # No. of off for each mating in Syn1 groups to produce Syn2
size_Syn1<-40  #40 #600 original

#-----------------------------
# Parameters for Syn2
#-----------------------------
b_tr1<-1/3
b_tr3<-1/3  # Weight to trait 3 in SI
b_tr4<-1/3  # Weight to trait 4 in SI
N_S_V<-20   #No. of plants sampled per variety to calculate allele dosage per plot
# Note -- N_S_V should be equal or 
          # smaller than size_Syn1 and size_F2
#-----------------------------
# Parameters for Syn3
#-----------------------------
No_Selected_Syn2<-nvar
size_Syn3<-40  #40 #600 original

# ---------------------------------
# Files to be written
# ---------------------------------
Ref_F2_Set<-paste('Ref_F2_Set',1:nCycle,sep='_')
file_Acc1<-paste('Acc_tr134_F2SP_rep',1:nRep,sep='_')

# F1
file_TBV_F1<-paste('TBV_AV_F1_rep',1:nRep,sep='_')
file_VAR_F1<-paste('TBV_VAR_F1_rep',1:nRep,sep='_')


# syn2
file4_syn2<-paste('TBV_AV_Syn2_rep',1:nRep,sep='_')
file5_syn2<-paste('Phen_ACC_Syn2_rep',1:nRep,sep='_')
file6_syn2<-paste('TBV_VAR_Syn2_rep',1:nRep,sep='_')

# syn3
file4_syn3<-paste('TBV_AV_Syn3_rep',1:nRep,sep='_')
file5_syn3<-paste('Phen_ACC_Syn3_rep',1:nRep,sep='_')
file6_syn3<-paste('TBV_VAR_Syn3_rep',1:nRep,sep='_')

# Selected
file_Selected<-paste('Selected_Cyl',1:nCycle,sep='_')

# Variance within set
File_Var_within_Tr1<-paste('Var_within_Tr1_rep',1:nRep,sep='_')
File_Var_within_Tr2<-paste('Var_within_Tr2_rep',1:nRep,sep='_')
File_Var_within_Tr3<-paste('Var_within_Tr3_rep',1:nRep,sep='_')
File_Var_within_Tr4<-paste('Var_within_Tr4_rep',1:nRep,sep='_')


# ---------------------------------
# Files to be read
# ---------------------------------
# Files of QMSIM to be read
file_V_data<-paste('Var',1:n_all_var,sep='')
file_V_data<-paste(file_V_data,'_data_001.txt',sep='')
file_V_data


file_V_qtl<-paste('Var',1:n_all_var,sep='')
file_V_qtl<-paste(file_V_qtl,'_qtl_001.txt',sep='')
file_V_qtl

file_V_mrk<-paste('Var',1:n_all_var,sep='')
file_V_mrk<-paste(file_V_mrk,'_mrk_001.txt',sep='')
file_V_mrk


# ---------------------------------
# Needed Matrixes for saving results
# ---------------------------------
# ACC_tr34_plot Trait 1 and 3 and 4
ACC_tr134_F2_SP<-matrix(0,nrow=(nCycle),ncol=4)
ACC_tr134_F2_SP_Selected<-matrix(0,nrow=(nCycle),ncol=4)


	# FOR SAVING GENETIC GAIN DATA
	
# ALL F1
	TBV_F1_AV<-matrix(0,nrow=nCycle,ncol=4)
	TBV_F1_VAR<-matrix(0,nrow=nCycle,ncol=4)	

# ALL Syn 2
	TBV_Syn2_AV<-matrix(0,nrow=nCycle,ncol=4)
	TBV_Syn2_sd<-matrix(0,nrow=nCycle,ncol=4)
    TBV_Syn2_ACC<-matrix(0,nrow=nCycle,ncol=4)
	TBV_Syn2_VAR<-matrix(0,nrow=nCycle,ncol=4)
	
# ALL Syn 3
	TBV_Syn3_AV<-matrix(0,nrow=nCycle,ncol=4)
	TBV_Syn3_sd<-matrix(0,nrow=nCycle,ncol=4)
	TBV_Syn3_ACC<-matrix(0,nrow=nCycle,ncol=4)
	TBV_Syn3_VAR<-matrix(0,nrow=nCycle,ncol=4)


# --------------------------------------
# MAIN PROG----------------START--------
# --------------------------------------


# READ VARITIES DATA

	linkage_map_mrk<-read.table('lm_mrk_001.txt',skip=1)
	linkage_map_qtl<-read.table('lm_qtl_001.txt',skip=1)
	link_map_main<-rbind(linkage_map_mrk,linkage_map_qtl)
	lmap<-c()
	li<-length(unique(link_map_main[,2]))
	for (i in 1:li){
	link<-subset(link_map_main,link_map_main[,2]==i)
	link<-link[order(link[,3]),]
	lmap<-rbind(lmap,link)
	}
	link_map<-lmap


	# START REPLICATES
	for (irep in 1:nRep){
	
	list_fun<-function() {
    list(list(), list(),list(), list())
}

Varit_M<-replicate(n_all_var, list_fun(), simplify=FALSE)


	
	# For creating traits
	file_V_data_set<-file_V_data
	file_V_qtl_set<-file_V_qtl
    file_V_mrk_set<-file_V_mrk
	
for (i in 1:n_all_var){

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

dim(Varit_M[[i]][[1]])
dim(Varit_M[[i]][[4]])


Varit<-Varit_M
	
	
# Extracting QTL from all varities
all_qtl<-c()
for (i in 1:length(Varit)){
a<-Varit[[i]]$qtl
# remove 1st colomn as it is id 
a<-a[,-c(1)]
all_qtl<-rbind(all_qtl,a)
}

all_qtl<-as.matrix(all_qtl)
nqtl<-dim(all_qtl)[2]/2

# trait simulation
var_add_tr1<-h2_tr1*phen_var_tr1
var_add_tr2<-h2_tr2*phen_var_tr2
var_add_tr3<-h2_tr3*phen_var_tr3
var_add_tr4<-h2_tr4*phen_var_tr4

cov_tr34<-rg_tr34*(sqrt(var_add_tr3)*sqrt(var_add_tr4))

Sigma <- matrix(0,nrow=4,ncol=4)
diag(Sigma)<-c(var_add_tr1,var_add_tr2,var_add_tr3,var_add_tr4)
Sigma[3,4]<-cov_tr34
Sigma[4,3]<-cov_tr34

add_effects<-mvrnorm(n = nqtl, rep(0, 4), Sigma,empirical = TRUE)

# Trait 1
after_tr1<-scale_trait_effect(h2=h2_tr1,phen_var=phen_var_tr1
,nqtl=nqtl,sampled_effects=add_effects[,1],qtlMatrix=all_qtl)
add_eff_1_tr1<-after_tr1[[2]][,1]
tbv_tr1<-after_tr1[[3]]

# Trait 2
after_tr2<-scale_trait_effect(h2=h2_tr2,phen_var=phen_var_tr2
,nqtl=nqtl,sampled_effects=add_effects[,2],qtlMatrix=all_qtl)
add_eff_1_tr2<-after_tr2[[2]][,1]
tbv_tr2<-after_tr2[[3]]

# Trait 3
after_tr3<-scale_trait_effect(h2=h2_tr3,phen_var=phen_var_tr3
,nqtl=nqtl,sampled_effects=add_effects[,3],qtlMatrix=all_qtl)
add_eff_1_tr3<-after_tr3[[2]][,1]
tbv_tr3<-after_tr3[[3]]

# Trait 4
after_tr4<-scale_trait_effect(h2=h2_tr4,phen_var=phen_var_tr4
,nqtl=nqtl,sampled_effects=add_effects[,4],qtlMatrix=all_qtl)
add_eff_1_tr4<-after_tr4[[2]][,1]
tbv_tr4<-after_tr4[[3]]

freq1<-after_tr3[[4]][,1]
freq2<-after_tr4[[4]][,2]
alpha_tr3<-after_tr3[[2]][,2]
alpha_tr4<-after_tr4[[2]][,2]

# some test
	rr<-sum(2*freq1*freq2*alpha_tr3*alpha_tr4)/
	(sqrt(after_tr3[[1]][3])*sqrt(after_tr4[[1]][3]))
	rr
	cor(add_effects[,3],add_effects[,4])
	cor(add_eff_1_tr3,add_eff_1_tr4)
	cor(tbv_tr1,tbv_tr2)
	cor(tbv_tr1,tbv_tr3)
	cor(tbv_tr1,tbv_tr4)
	cor(tbv_tr2,tbv_tr3)
	cor(tbv_tr2,tbv_tr4)
	cor(tbv_tr3,tbv_tr4)


# Create phen, env and tbv for each trait
add_eff_1<-data.frame(add_eff_1_tr1,add_eff_1_tr2,add_eff_1_tr3,add_eff_1_tr4)
add_eff_1<-as.matrix(add_eff_1)

phen_var<-c(phen_var_tr1,
phen_var_tr2,
phen_var_tr3,
phen_var_tr4)

var_add<-c(var_add_tr1,
var_add_tr2,
var_add_tr3,
var_add_tr4)
var_add


for (i in 1:length(Varit)){
qtlLoci<-Varit[[i]]$qtl
qtlLoci<-qtlLoci[,-1] #1st colomn is ID

dim(Varit[[i]]$data)[1]
tbv<-matrix(ncol=4,nrow=dim(Varit[[i]]$data)[1])
env<-matrix(ncol=4,nrow=dim(Varit[[i]]$data)[1])
phen<-matrix(ncol=4,nrow=dim(Varit[[i]]$data)[1])

dom_eff<-rep(0,length(add_eff_1[,1]))
	for (j in 1:4){

		tbv[,j]<-calc_TBV2(qtlLoci,add_eff_1[,j])
			
		# ENV 
		var_dom<-0
		vv<-phen_var[j]-(var_add[j]+var_dom)
		env[,j]<-rnorm(length(tbv[,j]),mean=0,sd=sqrt(vv))
			

		phen[,j]<-tbv[,j]+env[,j]
	}

	# Returning simulated traits to individuals in each vareity
	a<-Varit[[i]]$data
	new_data<-cbind(a,phen,tbv,env)
	dim(new_data)
	names(new_data)[6:17] <- c('Phen_1','Phen_2','Phen_3','Phen_4',
	'tbv_1','tbv_2','tbv_3','tbv_4',
	'env_1','env_1','env_3','env_4')
	
Varit[[i]]$data<-new_data	

cat('variety number:',i,'is done',fill=TRUE)

}

dim(Varit[[1]]$data)

		
list_fun<-function() {
    list(list(), list(), list(), list(),list(), list(), list(), list()
	, list(), list(), list(), list())
}

Sets<-replicate(nCycle, list_fun(), simplify=FALSE)


for (g_index in 1:length(Sets)){
names(Sets[[g_index]])<-c('Y1','Y2','Y3','Y4','Y5','Y6',
'Y7','Y8','Y9','Y10','Y11','Y12')
}



step_names<-c('F1','F2_Y','F2_W','F2_G','F2_SP','S_Cross',
'Syn1','Syn2_Y','Syn2_W','Syn2_G','Mult','Syn3')

 plan_XL<-matrix(0,ncol=nCycle+length(step_names),nrow=nCycle)
 dim(plan_XL)
  
  for(i in 1:nCycle){
  plan_XL[i,i:(i+length(step_names)-1)]<-step_names
   }
   
   plan_XL[,1:8]



# Cycle should start here
for (Cyl in 1:nCycle){

	qtl_effects<-data.frame(add_eff_1_tr1,add_eff_1_tr2,add_eff_1_tr3,add_eff_1_tr4)
	qtl_effects<-as.matrix(qtl_effects)

	phen_var<-c(phen_var_tr1,
	phen_var_tr2,
	phen_var_tr3,
	phen_var_tr4)

	var_add<-c(var_add_tr1,
	var_add_tr2,
	var_add_tr3,
	var_add_tr4)
	

  if(Cyl==1){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}


Varit<-Variant_Global(n_all_var,nvar)
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)

Sets[[1]][[Cyl]]<-F1_out

cat('All steps of year:',Cyl,'Finished',fill=TRUE)

}  

	
	
  if(Cyl==2){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}


Varit<-Variant_Global(n_all_var,nvar)
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)

given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Sets[[2]][[Cyl]]<-F1_out
Sets[[1]][[Cyl]]<-F2_Yellow

cat('All steps of year:',Cyl,'Finished',fill=TRUE)
}  



  if(Cyl==3){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}


Varit<-Variant_Global(n_all_var,nvar)
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)

given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)

Sets[[3]][[Cyl]]<-F1_out
Sets[[2]][[Cyl]]<-F2_Yellow
Sets[[1]][[Cyl]]<-F2_W

cat('All steps of year:',Cyl,'Finished',fill=TRUE)
}  



  if(Cyl==4){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}


Varit<-Variant_Global(n_all_var,nvar)
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)

given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)

F2_Green<-Sets[[1]][[Cyl-2]]

Sets[[4]][[Cyl]]<-F1_out
Sets[[3]][[Cyl]]<-F2_Yellow
Sets[[2]][[Cyl]]<-F2_W
Sets[[1]][[Cyl]]<-F2_Green

cat('All steps of year:',Cyl,'Finished',fill=TRUE)
}  




  if(Cyl==5){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}


Varit<-Variant_Global(n_all_var,nvar)
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)

given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[2]][[Cyl-2]]

given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Sets[[5]][[Cyl]]<-F1_out
Sets[[4]][[Cyl]]<-F2_Yellow
Sets[[3]][[Cyl]]<-F2_W
Sets[[2]][[Cyl]]<-F2_Green
Sets[[1]][[Cyl]]<-F2_SP_out

cat('All steps of year:',Cyl,'Finished',fill=TRUE)

}  




  if(Cyl==6){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}


Varit<-Variant_Global(n_all_var,nvar)
# Select initial cross from initial material+good SF
sample_some_varit<-sample(1:nvar,10)
Varit_aval<-Varit[sample_some_varit]
Varit_aval[[1]][[1]]

# get mean for these
hanchi<-Sets[[1]][[Cyl-2]] #F2
SI<-c()
selected<-matrix(0,ncol=2,nrow=length(hanchi))

	for (i in 1:length(hanchi)){
		# tr1
		mean_tr1<-mean(hanchi[[i]]$data$phen_1)
		
		# tr3
		mean_tr3<-mean(hanchi[[i]]$data$phen_3)
		# tr4
		mean_tr4<-mean(hanchi[[i]]$data$phen_4)

		SI[i]<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
		selected[i,1]<-i
		selected[i,2]<-SI[i]
	}
selected
selected<-selected[order(-selected[,2]),]

selected_ID_F2<-selected[1:10,1]

Varit_dov<-Sets[[1]][[Cyl-1]]
# length(Varit_dov)
# lengths(Varit_dov)

Selected_F2_SP<-Varit_dov[selected_ID_F2]
# length(Selected_F2_SP)
# lengths(Selected_F2_SP)
# Selected_F2_SP[[1]][[1]]



# Make F2 selected to similar output of Syn2 or 3
for (u in 1:length(Selected_F2_SP)){
localID<-Selected_F2_SP[[u]]$data[,1]

# QTL
x_qtl<-Selected_F2_SP[[u]][[2]]
x_qtl<-cbind(localID,x_qtl)
Selected_F2_SP[[u]][[2]]<-x_qtl

# MRK
x_mrk<-Selected_F2_SP[[u]][[3]]
x_mrk<-cbind(localID,x_mrk)
Selected_F2_SP[[u]][[3]]<-x_mrk

# SEQ
x_seq<-Selected_F2_SP[[u]][[4]]
x_seq<-cbind(localID,x_seq)
Selected_F2_SP[[u]][[4]]<-x_seq
}
# test
Selected_F2_SP[[1]][[2]][,1:10]



Varit_all<-c(Varit_aval,Selected_F2_SP)
# length(Varit_all)
# lengths(Varit_all)
# Varit_all[[1]][[1]]
# Varit_all[[2]][[1]]
# Varit_all[[15]][[1]]
# Varit_all[[20]][[1]]

F1_out<-F1_Global(Varit_all,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)

given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[3]][[Cyl-2]]

given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)

Sets[[6]][[Cyl]]<-F1_out
Sets[[5]][[Cyl]]<-F2_Yellow
Sets[[4]][[Cyl]]<-F2_W
Sets[[3]][[Cyl]]<-F2_Green
Sets[[2]][[Cyl]]<-F2_SP_out
Sets[[1]][[Cyl]]<-S_Cross

cat('All steps of year:',Cyl,'Finished',fill=TRUE)
}  





  if(Cyl==7){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}


# Here we select F2_SP set 1 and 2

# Set 1
	# get mean for these
	hanchi<-Sets[[1]][[Cyl-3]] #F2
	SI<-c()
	selected<-matrix(0,ncol=2,nrow=length(hanchi))

		for (i in 1:length(hanchi)){
			# tr1
			mean_tr1<-mean(hanchi[[i]]$data$phen_1)
			
			# tr3
			mean_tr3<-mean(hanchi[[i]]$data$phen_3)
			# tr4
			mean_tr4<-mean(hanchi[[i]]$data$phen_4)

			SI[i]<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
			selected[i,1]<-i
			selected[i,2]<-SI[i]
		}
	selected
	selected<-selected[order(-selected[,2]),]
	selected_ID_F2_1<-selected[1:10,1]

	Varit_inter<-Sets[[1]][[Cyl-2]]
	# length(Varit_aval)
	# lengths(Varit_aval)

	Selected_F2_SP_AVAL<-Varit_inter[selected_ID_F2_1]
	# length(Selected_F2_SP_AVAL)
	# lengths(Selected_F2_SP_AVAL)
	# Selected_F2_SP_AVAL[[1]][[1]]


	# Make F2 selected to similar output of Syn2 or 3
	for (u in 1:length(Selected_F2_SP_AVAL)){
	localID<-Selected_F2_SP_AVAL[[u]]$data[,1]

	# QTL
	x_qtl<-Selected_F2_SP_AVAL[[u]][[2]]
	x_qtl<-cbind(localID,x_qtl)
	Selected_F2_SP_AVAL[[u]][[2]]<-x_qtl

	# MRK
	x_mrk<-Selected_F2_SP_AVAL[[u]][[3]]
	x_mrk<-cbind(localID,x_mrk)
	Selected_F2_SP_AVAL[[u]][[3]]<-x_mrk

	# SEQ
	x_seq<-Selected_F2_SP_AVAL[[u]][[4]]
	x_seq<-cbind(localID,x_seq)
	Selected_F2_SP_AVAL[[u]][[4]]<-x_seq
	}
	# test
	Selected_F2_SP_AVAL[[1]][[2]][,1:10]

# SEt 2
	# get mean for these
	hanchi<-Sets[[2]][[Cyl-2]] #F2
	SI<-c()
	selected<-matrix(0,ncol=2,nrow=length(hanchi))

		for (i in 1:length(hanchi)){
			# tr1
			mean_tr1<-mean(hanchi[[i]]$data$phen_1)
			
			# tr3
			mean_tr3<-mean(hanchi[[i]]$data$phen_3)
			# tr4
			mean_tr4<-mean(hanchi[[i]]$data$phen_4)

			SI[i]<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
			selected[i,1]<-i
			selected[i,2]<-SI[i]
		}
	selected
	selected<-selected[order(-selected[,2]),]
	selected_ID_F2_2<-selected[1:10,1]

	Varit_inter<-Sets[[2]][[Cyl-1]] #F2_SP
	Selected_F2_SP_DOV<-Varit_inter[selected_ID_F2_2]


	# Make F2 selected to similar output of Syn2 or 3
	for (u in 1:length(Selected_F2_SP_DOV)){
	localID<-Selected_F2_SP_DOV[[u]]$data[,1]

	# QTL
	x_qtl<-Selected_F2_SP_DOV[[u]][[2]]
	x_qtl<-cbind(localID,x_qtl)
	Selected_F2_SP_DOV[[u]][[2]]<-x_qtl

	# MRK
	x_mrk<-Selected_F2_SP_DOV[[u]][[3]]
	x_mrk<-cbind(localID,x_mrk)
	Selected_F2_SP_DOV[[u]][[3]]<-x_mrk

	# SEQ
	x_seq<-Selected_F2_SP_DOV[[u]][[4]]
	x_seq<-cbind(localID,x_seq)
	Selected_F2_SP_DOV[[u]][[4]]<-x_seq
	}
	# test
	Selected_F2_SP_DOV[[1]][[2]][,1:10]

	

Varit_all<-c(Selected_F2_SP_AVAL,Selected_F2_SP_DOV)
# length(Varit_all)
# lengths(Varit_all)


F1_out<-F1_Global(Varit_all,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)

given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[Cyl-3]][[Cyl-2]]

given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)


# Syn1
F2<-Sets[[Cyl-6]][[Cyl-3]]
F2_SP<-Sets[[Cyl-6]][[Cyl-2]]

all_Syn1_out<-Global_Syn1(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
Syn1_out<-all_Syn1_out[[1]]
S_Cross<-all_Syn1_out[[2]]



Sets[[7]][[Cyl]]<-F1_out
Sets[[6]][[Cyl]]<-F2_Yellow
Sets[[5]][[Cyl]]<-F2_W
Sets[[4]][[Cyl]]<-F2_Green
Sets[[3]][[Cyl]]<-F2_SP_out
Sets[[2]][[Cyl]]<-S_Cross
Sets[[1]][[Cyl]]<-Syn1_out

Sets[[1]][[Cyl-1]]<-S_Cross

cat('All steps of year:',Cyl,'Finished',fill=TRUE)

}  





  if(Cyl==8){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}


which_sets<-c(1,2,3)
which_step<-c(4,5,6)
which_sets
which_step
		
     		evalu_lst<-list()
		
			for(K_1 in 1:length(which_sets)){
			dakheli1<-Sets[[which_sets[K_1]]][[which_step[K_1]]]
				evalu_mat<-matrix(ncol=4,nrow=length(dakheli1))
				for (K_2 in 1:length(dakheli1)){
				# tr1
				mean_tr1<-mean(dakheli1[[K_2]]$data$phen_1)
				# tr3
				mean_tr3<-mean(dakheli1[[K_2]]$data$phen_3)
				# tr4
				mean_tr4<-mean(dakheli1[[K_2]]$data$phen_4)
				SI<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
				evalu_mat[K_2,1]<-which_sets[K_1]
				evalu_mat[K_2,2]<-which_step[K_1]
				evalu_mat[K_2,3]<-K_2
				evalu_mat[K_2,4]<-SI
				}
				evalu_lst[[K_1]]<-evalu_mat
			}
			
			
		evalu<-do.call("rbind", evalu_lst)
		colnames(evalu)<-c("Set No.","Step","Family No.",'Index')
		evalu
		
		evalu<-evalu[order(-evalu[,4]),]
		evalu	
		evalu<-evalu[1:No_Selected_Syn2,]
		evalu[,2]<-evalu[,2]+1
		
		write.table(evalu,file=file_Selected[Cyl])
				
		list_fun<-function() {
		list(list(), list(),list(), list())
		}

Update_data<-replicate(No_Selected_Syn2, list_fun(), simplify=FALSE)


			for(upd in 1:dim(evalu)[1]){
			set_number<-evalu[upd,1]
			Step_code<-evalu[upd,2]
			Syn_number<-evalu[upd,3]
			Update_data[[upd]][[1]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[1]]
			Update_data[[upd]][[2]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[2]]
			Update_data[[upd]][[3]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[3]]
			Update_data[[upd]][[4]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[4]]
			}
		
	        Varit<-Update_data	
		
		for (g_index in 1:length(Varit)){
		names(Varit[[g_index]])<-c('data','qtl','mrk','sequ')
		}
		
		

	# Make selected to similar output of Syn2 or 3
	for (u in 1:length(Varit)){
	localID<-Varit[[u]]$data[,1]

	# QTL
	x_qtl<-Varit[[u]][[2]]
	x_qtl<-cbind(localID,x_qtl)
	Varit[[u]][[2]]<-x_qtl

	# MRK
	x_mrk<-Varit[[u]][[3]]
	x_mrk<-cbind(localID,x_mrk)
	Varit[[u]][[3]]<-x_mrk

	# SEQ
	x_seq<-Varit[[u]][[4]]
	x_seq<-cbind(localID,x_seq)
	Varit[[u]][[4]]<-x_seq
	}
	# test
	Varit[[1]][[2]][,1:10]


F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)

given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[Cyl-3]][[Cyl-2]]

given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)


# Syn1
F2<-Sets[[Cyl-6]][[Cyl-3]]
F2_SP<-Sets[[Cyl-6]][[Cyl-2]]

all_Syn1_out<-Global_Syn1(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
Syn1_out<-all_Syn1_out[[1]]
S_Cross<-all_Syn1_out[[2]]

mini_plots<-Sets[[Cyl-7]][[Cyl-1]]
Syn2_out<-Global_Syn2(mini_plots,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)




Sets[[8]][[Cyl]]<-F1_out
Sets[[7]][[Cyl]]<-F2_Yellow
Sets[[6]][[Cyl]]<-F2_W
Sets[[5]][[Cyl]]<-F2_Green
Sets[[4]][[Cyl]]<-F2_SP_out
Sets[[3]][[Cyl]]<-S_Cross
Sets[[2]][[Cyl]]<-Syn1_out
Sets[[1]][[Cyl]]<-Syn2_out

Sets[[2]][[Cyl-1]]<-S_Cross

cat('All steps of year:',Cyl,'Finished',fill=TRUE)

}  



  if(Cyl==9){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}


which_sets<-c(1,2,3,4)
which_step<-c(4,5,6,7)
which_sets
which_step
		
     		evalu_lst<-list()
		
			for(K_1 in 1:length(which_sets)){
			dakheli1<-Sets[[which_sets[K_1]]][[which_step[K_1]]]
				evalu_mat<-matrix(ncol=4,nrow=length(dakheli1))
				for (K_2 in 1:length(dakheli1)){
				# tr1
				mean_tr1<-mean(dakheli1[[K_2]]$data$phen_1)
				# tr3
				mean_tr3<-mean(dakheli1[[K_2]]$data$phen_3)
				# tr4
				mean_tr4<-mean(dakheli1[[K_2]]$data$phen_4)
				SI<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
				evalu_mat[K_2,1]<-which_sets[K_1]
				evalu_mat[K_2,2]<-which_step[K_1]
				evalu_mat[K_2,3]<-K_2
				evalu_mat[K_2,4]<-SI
				}
				evalu_lst[[K_1]]<-evalu_mat
			}
			
			
		evalu<-do.call("rbind", evalu_lst)
		colnames(evalu)<-c("Set No.","Step","Family No.",'Index')
		evalu
		
		evalu<-evalu[order(-evalu[,4]),]
		evalu	
		evalu<-evalu[1:No_Selected_Syn2,]
		evalu[,2]<-evalu[,2]+1
		
		write.table(evalu,file=file_Selected[Cyl])
				
		list_fun<-function() {
		list(list(), list(),list(), list())
		}

Update_data<-replicate(No_Selected_Syn2, list_fun(), simplify=FALSE)


			for(upd in 1:dim(evalu)[1]){
			set_number<-evalu[upd,1]
			Step_code<-evalu[upd,2]
			Syn_number<-evalu[upd,3]
			Update_data[[upd]][[1]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[1]]
			Update_data[[upd]][[2]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[2]]
			Update_data[[upd]][[3]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[3]]
			Update_data[[upd]][[4]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[4]]
			}
		
	        Varit<-Update_data	
		
		for (g_index in 1:length(Varit)){
		names(Varit[[g_index]])<-c('data','qtl','mrk','sequ')
		}
		
		

	# Make selected to similar output of Syn2 or 3
	for (u in 1:length(Varit)){
	localID<-Varit[[u]]$data[,1]

	# QTL
	x_qtl<-Varit[[u]][[2]]
	x_qtl<-cbind(localID,x_qtl)
	Varit[[u]][[2]]<-x_qtl

	# MRK
	x_mrk<-Varit[[u]][[3]]
	x_mrk<-cbind(localID,x_mrk)
	Varit[[u]][[3]]<-x_mrk

	# SEQ
	x_seq<-Varit[[u]][[4]]
	x_seq<-cbind(localID,x_seq)
	Varit[[u]][[4]]<-x_seq
	}
	# test
	Varit[[1]][[2]][,1:10]

#F1
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)


# F2
given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[Cyl-3]][[Cyl-2]]

# F2_SP
given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)


# Syn1
F2<-Sets[[Cyl-6]][[Cyl-3]]
F2_SP<-Sets[[Cyl-6]][[Cyl-2]]

all_Syn1_out<-Global_Syn1(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
Syn1_out<-all_Syn1_out[[1]]
S_Cross<-all_Syn1_out[[2]]

# Syn2
mini_plots<-Sets[[Cyl-7]][[Cyl-1]]
Syn2_out<-Global_Syn2(mini_plots,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Syn2_W<-rnorm(5,1,10)

Sets[[9]][[Cyl]]<-F1_out
Sets[[8]][[Cyl]]<-F2_Yellow
Sets[[7]][[Cyl]]<-F2_W
Sets[[6]][[Cyl]]<-F2_Green
Sets[[5]][[Cyl]]<-F2_SP_out
Sets[[4]][[Cyl]]<-S_Cross
Sets[[3]][[Cyl]]<-Syn1_out
Sets[[3]][[Cyl-1]]<-S_Cross
Sets[[2]][[Cyl]]<-Syn2_out
Sets[[1]][[Cyl]]<-Syn2_W

cat('All steps of year:',Cyl,'Finished',fill=TRUE)

}  




  if(Cyl==10){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}


which_sets<-c(1:(Cyl-5))
which_step<-c(4:(Cyl-2))
which_sets
which_step
		
     		evalu_lst<-list()
		
			for(K_1 in 1:length(which_sets)){
			dakheli1<-Sets[[which_sets[K_1]]][[which_step[K_1]]]
				evalu_mat<-matrix(ncol=4,nrow=length(dakheli1))
				for (K_2 in 1:length(dakheli1)){
				# tr1
				mean_tr1<-mean(dakheli1[[K_2]]$data$phen_1)
				# tr3
				mean_tr3<-mean(dakheli1[[K_2]]$data$phen_3)
				# tr4
				mean_tr4<-mean(dakheli1[[K_2]]$data$phen_4)
				SI<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
				evalu_mat[K_2,1]<-which_sets[K_1]
				evalu_mat[K_2,2]<-which_step[K_1]
				evalu_mat[K_2,3]<-K_2
				evalu_mat[K_2,4]<-SI
				}
				evalu_lst[[K_1]]<-evalu_mat
			}
			
			
		evalu<-do.call("rbind", evalu_lst)
		colnames(evalu)<-c("Set No.","Step","Family No.",'Index')
		evalu
		
		evalu<-evalu[order(-evalu[,4]),]
		evalu	
		evalu<-evalu[1:No_Selected_Syn2,]
		evalu[,2]<-evalu[,2]+1
		
		write.table(evalu,file=file_Selected[Cyl])
				
		list_fun<-function() {
		list(list(), list(),list(), list())
		}

Update_data<-replicate(No_Selected_Syn2, list_fun(), simplify=FALSE)


			for(upd in 1:dim(evalu)[1]){
			set_number<-evalu[upd,1]
			Step_code<-evalu[upd,2]
			Syn_number<-evalu[upd,3]
			Update_data[[upd]][[1]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[1]]
			Update_data[[upd]][[2]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[2]]
			Update_data[[upd]][[3]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[3]]
			Update_data[[upd]][[4]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[4]]
			}
		
	        Varit<-Update_data	
		
		for (g_index in 1:length(Varit)){
		names(Varit[[g_index]])<-c('data','qtl','mrk','sequ')
		}
		
		

	# Make selected to similar output of Syn2 or 3
	for (u in 1:length(Varit)){
	localID<-Varit[[u]]$data[,1]

	# QTL
	x_qtl<-Varit[[u]][[2]]
	x_qtl<-cbind(localID,x_qtl)
	Varit[[u]][[2]]<-x_qtl

	# MRK
	x_mrk<-Varit[[u]][[3]]
	x_mrk<-cbind(localID,x_mrk)
	Varit[[u]][[3]]<-x_mrk

	# SEQ
	x_seq<-Varit[[u]][[4]]
	x_seq<-cbind(localID,x_seq)
	Varit[[u]][[4]]<-x_seq
	}
	# test
	Varit[[1]][[2]][,1:10]

#F1
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)


# F2
given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[Cyl-3]][[Cyl-2]]

# F2_SP
given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)


# Syn1
F2<-Sets[[Cyl-6]][[Cyl-3]]
F2_SP<-Sets[[Cyl-6]][[Cyl-2]]

all_Syn1_out<-Global_Syn1(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
Syn1_out<-all_Syn1_out[[1]]
S_Cross<-all_Syn1_out[[2]]

# Syn2
mini_plots<-Sets[[Cyl-7]][[Cyl-1]]
Syn2_out<-Global_Syn2(mini_plots,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Syn2_W<-rnorm(5,1,10)
Syn2_Green<-Sets[[Cyl-9]][[Cyl-2]]


Sets[[10]][[Cyl]]<-F1_out
Sets[[9]][[Cyl]]<-F2_Yellow
Sets[[8]][[Cyl]]<-F2_W
Sets[[7]][[Cyl]]<-F2_Green
Sets[[6]][[Cyl]]<-F2_SP_out
Sets[[5]][[Cyl]]<-S_Cross
Sets[[4]][[Cyl]]<-Syn1_out
Sets[[4]][[Cyl-1]]<-S_Cross
Sets[[3]][[Cyl]]<-Syn2_out
Sets[[2]][[Cyl]]<-Syn2_W
Sets[[1]][[Cyl]]<-Syn2_Green

cat('All steps of year:',Cyl,'Finished',fill=TRUE)
}  




  if(Cyl==11){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}


which_sets<-c(1:(Cyl-5))
which_step<-c(4:(Cyl-2))
which_sets
which_step
		
     		evalu_lst<-list()
		
			for(K_1 in 1:length(which_sets)){
			dakheli1<-Sets[[which_sets[K_1]]][[which_step[K_1]]]
				evalu_mat<-matrix(ncol=4,nrow=length(dakheli1))
				for (K_2 in 1:length(dakheli1)){
				# tr1
				mean_tr1<-mean(dakheli1[[K_2]]$data$phen_1)
				# tr3
				mean_tr3<-mean(dakheli1[[K_2]]$data$phen_3)
				# tr4
				mean_tr4<-mean(dakheli1[[K_2]]$data$phen_4)
				SI<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
				evalu_mat[K_2,1]<-which_sets[K_1]
				evalu_mat[K_2,2]<-which_step[K_1]
				evalu_mat[K_2,3]<-K_2
				evalu_mat[K_2,4]<-SI
				}
				evalu_lst[[K_1]]<-evalu_mat
			}
			
			
		evalu<-do.call("rbind", evalu_lst)
		colnames(evalu)<-c("Set No.","Step","Family No.",'Index')
		evalu
		
		evalu<-evalu[order(-evalu[,4]),]
		evalu	
		evalu<-evalu[1:No_Selected_Syn2,]
		evalu[,2]<-evalu[,2]+1
		
		write.table(evalu,file=file_Selected[Cyl])
				
		list_fun<-function() {
		list(list(), list(),list(), list())
		}

Update_data<-replicate(No_Selected_Syn2, list_fun(), simplify=FALSE)


			for(upd in 1:dim(evalu)[1]){
			set_number<-evalu[upd,1]
			Step_code<-evalu[upd,2]
			Syn_number<-evalu[upd,3]
			Update_data[[upd]][[1]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[1]]
			Update_data[[upd]][[2]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[2]]
			Update_data[[upd]][[3]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[3]]
			Update_data[[upd]][[4]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[4]]
			}
		
	        Varit<-Update_data	
		
		for (g_index in 1:length(Varit)){
		names(Varit[[g_index]])<-c('data','qtl','mrk','sequ')
		}
		
		

	# Make selected to similar output of Syn2 or 3
	for (u in 1:length(Varit)){
	localID<-Varit[[u]]$data[,1]

	# QTL
	x_qtl<-Varit[[u]][[2]]
	x_qtl<-cbind(localID,x_qtl)
	Varit[[u]][[2]]<-x_qtl

	# MRK
	x_mrk<-Varit[[u]][[3]]
	x_mrk<-cbind(localID,x_mrk)
	Varit[[u]][[3]]<-x_mrk

	# SEQ
	x_seq<-Varit[[u]][[4]]
	x_seq<-cbind(localID,x_seq)
	Varit[[u]][[4]]<-x_seq
	}
	# test
	Varit[[1]][[2]][,1:10]

#F1
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)


# F2
given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[Cyl-3]][[Cyl-2]]

# F2_SP
given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)


# Syn1
F2<-Sets[[Cyl-6]][[Cyl-3]]
F2_SP<-Sets[[Cyl-6]][[Cyl-2]]

all_Syn1_out<-Global_Syn1(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
Syn1_out<-all_Syn1_out[[1]]
S_Cross<-all_Syn1_out[[2]]

# Syn2
mini_plots<-Sets[[Cyl-7]][[Cyl-1]]
Syn2_out<-Global_Syn2(mini_plots,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Syn2_W<-rnorm(5,1,10)
Syn2_Green<-Sets[[Cyl-9]][[Cyl-2]]

mult<-rnorm(5,1,10)

Sets[[11]][[Cyl]]<-F1_out
Sets[[10]][[Cyl]]<-F2_Yellow
Sets[[9]][[Cyl]]<-F2_W
Sets[[8]][[Cyl]]<-F2_Green
Sets[[7]][[Cyl]]<-F2_SP_out
Sets[[6]][[Cyl]]<-S_Cross
Sets[[5]][[Cyl]]<-Syn1_out
Sets[[5]][[Cyl-1]]<-S_Cross
Sets[[4]][[Cyl]]<-Syn2_out
Sets[[3]][[Cyl]]<-Syn2_W
Sets[[2]][[Cyl]]<-Syn2_Green
Sets[[1]][[Cyl]]<-mult

cat('All steps of year:',Cyl,'Finished',fill=TRUE)

}  



  if(Cyl==12){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}

# For F2_SP
which_sets<-c(1:(Cyl-5))
which_step<-c(4:(Cyl-2))
which_sets
which_step

# For Syn2
which_sets2<-c(1:(Cyl-11))
which_step2<-c(10:(Cyl-2))
which_sets2
which_step2

which_sets<-c(which_sets2,which_sets)
which_step<-c(which_step2,which_step)
which_sets
which_step
		
     		evalu_lst<-list()
		
			for(K_1 in 1:length(which_sets)){
			dakheli1<-Sets[[which_sets[K_1]]][[which_step[K_1]]]
				evalu_mat<-matrix(ncol=4,nrow=length(dakheli1))
				for (K_2 in 1:length(dakheli1)){
				# tr1
				mean_tr1<-mean(dakheli1[[K_2]]$data$phen_1)
				# tr3
				mean_tr3<-mean(dakheli1[[K_2]]$data$phen_3)
				# tr4
				mean_tr4<-mean(dakheli1[[K_2]]$data$phen_4)
				SI<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
				evalu_mat[K_2,1]<-which_sets[K_1]
				evalu_mat[K_2,2]<-which_step[K_1]
				evalu_mat[K_2,3]<-K_2
				evalu_mat[K_2,4]<-SI
				}
				evalu_lst[[K_1]]<-evalu_mat
			}
			
			
				evalu<-do.call("rbind", evalu_lst)
		colnames(evalu)<-c("Set No.","Step","Family No.",'Index')
		evalu
		
		evalu<-evalu[order(-evalu[,4]),]
		evalu	
		evalu<-evalu[1:No_Selected_Syn2,]
		evalu

		jedi<-evalu[,2]-evalu[,1]
		x_f<-which(jedi==9)
		
		if(length(x_f)>0){
		evalu[-x_f,2]<-evalu[-x_f,2]+1		
		}
		
		if(length(x_f)==0){
		evalu[,2]<-evalu[,2]+1
		}
		

		
		write.table(evalu,file=file_Selected[Cyl])
				
		list_fun<-function() {
		list(list(), list(),list(), list())
		}

Update_data<-replicate(No_Selected_Syn2, list_fun(), simplify=FALSE)


			for(upd in 1:dim(evalu)[1]){
			set_number<-evalu[upd,1]
			Step_code<-evalu[upd,2]
			Syn_number<-evalu[upd,3]
			Update_data[[upd]][[1]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[1]]
			Update_data[[upd]][[2]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[2]]
			Update_data[[upd]][[3]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[3]]
			Update_data[[upd]][[4]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[4]]
			}
		
	        Varit<-Update_data	
		
		for (g_index in 1:length(Varit)){
		names(Varit[[g_index]])<-c('data','qtl','mrk','sequ')
		}
		
		

	# Make selected to similar output of Syn2 or 3
	for (u in 1:length(Varit)){
	
	if(u %in% x_f){
	Varit[[u]]<-Varit[[u]]
		} else {
			
			localID<-Varit[[u]]$data[,1]

			# QTL
			x_qtl<-Varit[[u]][[2]]
			x_qtl<-cbind(localID,x_qtl)
			Varit[[u]][[2]]<-x_qtl

			# MRK
			x_mrk<-Varit[[u]][[3]]
			x_mrk<-cbind(localID,x_mrk)
			Varit[[u]][[3]]<-x_mrk

			# SEQ
			x_seq<-Varit[[u]][[4]]
			x_seq<-cbind(localID,x_seq)
			Varit[[u]][[4]]<-x_seq
		}
	}
	
	# test
	# Varit[[9]][[2]][,1:10]

#F1
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)


# F2
given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[Cyl-3]][[Cyl-2]]

# F2_SP
given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)


# Syn1
F2<-Sets[[Cyl-6]][[Cyl-3]]
F2_SP<-Sets[[Cyl-6]][[Cyl-2]]

all_Syn1_out<-Global_Syn1(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
Syn1_out<-all_Syn1_out[[1]]
S_Cross<-all_Syn1_out[[2]]

# Syn2
mini_plots<-Sets[[Cyl-7]][[Cyl-1]]
Syn2_out<-Global_Syn2(mini_plots,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Syn2_W<-rnorm(5,1,10)
Syn2_Green<-Sets[[Cyl-9]][[Cyl-2]]

mult<-rnorm(5,1,10)

# Syn3
Syn_2<-Sets[[Cyl-11]][[Cyl-2]]
Syn3_out<-Global_Syn3(Syn_2,size_Syn3,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)


Sets[[12]][[Cyl]]<-F1_out
Sets[[11]][[Cyl]]<-F2_Yellow
Sets[[10]][[Cyl]]<-F2_W
Sets[[9]][[Cyl]]<-F2_Green
Sets[[8]][[Cyl]]<-F2_SP_out
Sets[[7]][[Cyl]]<-S_Cross
Sets[[6]][[Cyl]]<-Syn1_out
Sets[[6]][[Cyl-1]]<-S_Cross
Sets[[5]][[Cyl]]<-Syn2_out
Sets[[4]][[Cyl]]<-Syn2_W
Sets[[3]][[Cyl]]<-Syn2_Green
Sets[[2]][[Cyl]]<-mult
Sets[[1]][[Cyl]]<-Syn3_out

cat('All steps of year:',Cyl,'Finished',fill=TRUE)
}  




  if(Cyl==13){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}

# For F2_SP
which_sets<-c(1:(Cyl-5))
which_step<-c(4:(Cyl-2))
which_sets
which_step

# For Syn2
which_sets2<-c(1:(Cyl-11))
which_step2<-c(10:(Cyl-2))
which_sets2
which_step2

which_sets<-c(which_sets2,which_sets)
which_step<-c(which_step2,which_step)
which_sets
which_step
		
     		evalu_lst<-list()
		
			for(K_1 in 1:length(which_sets)){
			dakheli1<-Sets[[which_sets[K_1]]][[which_step[K_1]]]
				evalu_mat<-matrix(ncol=4,nrow=length(dakheli1))
				for (K_2 in 1:length(dakheli1)){
				# tr1
				mean_tr1<-mean(dakheli1[[K_2]]$data$phen_1)
				# tr3
				mean_tr3<-mean(dakheli1[[K_2]]$data$phen_3)
				# tr4
				mean_tr4<-mean(dakheli1[[K_2]]$data$phen_4)
				SI<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
				evalu_mat[K_2,1]<-which_sets[K_1]
				evalu_mat[K_2,2]<-which_step[K_1]
				evalu_mat[K_2,3]<-K_2
				evalu_mat[K_2,4]<-SI
				}
				evalu_lst[[K_1]]<-evalu_mat
			}
			
			
				evalu<-do.call("rbind", evalu_lst)
		colnames(evalu)<-c("Set No.","Step","Family No.",'Index')
		evalu
		
		evalu<-evalu[order(-evalu[,4]),]
		evalu	
		evalu<-evalu[1:No_Selected_Syn2,]
		evalu

		jedi<-evalu[,2]-evalu[,1]
		x_f<-which(jedi==9)
		
		if(length(x_f)>0){
		evalu[-x_f,2]<-evalu[-x_f,2]+1		
		}
		
		if(length(x_f)==0){
		evalu[,2]<-evalu[,2]+1
		}
		

		
		write.table(evalu,file=file_Selected[Cyl])
				
		list_fun<-function() {
		list(list(), list(),list(), list())
		}

Update_data<-replicate(No_Selected_Syn2, list_fun(), simplify=FALSE)


			for(upd in 1:dim(evalu)[1]){
			set_number<-evalu[upd,1]
			Step_code<-evalu[upd,2]
			Syn_number<-evalu[upd,3]
			Update_data[[upd]][[1]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[1]]
			Update_data[[upd]][[2]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[2]]
			Update_data[[upd]][[3]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[3]]
			Update_data[[upd]][[4]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[4]]
			}
		
	        Varit<-Update_data	
		
		for (g_index in 1:length(Varit)){
		names(Varit[[g_index]])<-c('data','qtl','mrk','sequ')
		}
		
		

	# Make selected to similar output of Syn2 or 3
	for (u in 1:length(Varit)){
	
	if(u %in% x_f){
	Varit[[u]]<-Varit[[u]]
		} else {
			
			localID<-Varit[[u]]$data[,1]

			# QTL
			x_qtl<-Varit[[u]][[2]]
			x_qtl<-cbind(localID,x_qtl)
			Varit[[u]][[2]]<-x_qtl

			# MRK
			x_mrk<-Varit[[u]][[3]]
			x_mrk<-cbind(localID,x_mrk)
			Varit[[u]][[3]]<-x_mrk

			# SEQ
			x_seq<-Varit[[u]][[4]]
			x_seq<-cbind(localID,x_seq)
			Varit[[u]][[4]]<-x_seq
		}
	}
	
	# test
	# Varit[[9]][[2]][,1:10]

#F1
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)


# F2
given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[Cyl-3]][[Cyl-2]]

# F2_SP
given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)


# Syn1
F2<-Sets[[Cyl-6]][[Cyl-3]]
F2_SP<-Sets[[Cyl-6]][[Cyl-2]]

all_Syn1_out<-Global_Syn1(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
Syn1_out<-all_Syn1_out[[1]]
S_Cross<-all_Syn1_out[[2]]

# Syn2
mini_plots<-Sets[[Cyl-7]][[Cyl-1]]
Syn2_out<-Global_Syn2(mini_plots,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Syn2_W<-rnorm(5,1,10)
Syn2_Green<-Sets[[Cyl-9]][[Cyl-2]]

mult<-rnorm(5,1,10)

# Syn3
Syn_2<-Sets[[Cyl-11]][[Cyl-2]]
Syn3_out<-Global_Syn3(Syn_2,size_Syn3,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)


Sets[[13]][[Cyl]]<-F1_out
Sets[[12]][[Cyl]]<-F2_Yellow
Sets[[11]][[Cyl]]<-F2_W
Sets[[10]][[Cyl]]<-F2_Green
Sets[[9]][[Cyl]]<-F2_SP_out
Sets[[8]][[Cyl]]<-S_Cross
Sets[[7]][[Cyl]]<-Syn1_out
Sets[[7]][[Cyl-1]]<-S_Cross
Sets[[6]][[Cyl]]<-Syn2_out
Sets[[5]][[Cyl]]<-Syn2_W
Sets[[4]][[Cyl]]<-Syn2_Green
Sets[[3]][[Cyl]]<-mult
Sets[[2]][[Cyl]]<-Syn3_out

cat('All steps of year:',Cyl,'Finished',fill=TRUE)
}  




  if(Cyl==14){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}

# For F2_SP
which_sets<-c(1:(Cyl-5))
which_step<-c(4:(Cyl-2))
which_sets
which_step

# For Syn2
which_sets2<-c(1:(Cyl-11))
which_step2<-c(10:(Cyl-2))
which_sets2
which_step2

which_sets<-c(which_sets2,which_sets)
which_step<-c(which_step2,which_step)
which_sets
which_step
		
     		evalu_lst<-list()
		
			for(K_1 in 1:length(which_sets)){
			dakheli1<-Sets[[which_sets[K_1]]][[which_step[K_1]]]
				evalu_mat<-matrix(ncol=4,nrow=length(dakheli1))
				for (K_2 in 1:length(dakheli1)){
				# tr1
				mean_tr1<-mean(dakheli1[[K_2]]$data$phen_1)
				# tr3
				mean_tr3<-mean(dakheli1[[K_2]]$data$phen_3)
				# tr4
				mean_tr4<-mean(dakheli1[[K_2]]$data$phen_4)
				SI<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
				evalu_mat[K_2,1]<-which_sets[K_1]
				evalu_mat[K_2,2]<-which_step[K_1]
				evalu_mat[K_2,3]<-K_2
				evalu_mat[K_2,4]<-SI
				}
				evalu_lst[[K_1]]<-evalu_mat
			}
			
			
				evalu<-do.call("rbind", evalu_lst)
		colnames(evalu)<-c("Set No.","Step","Family No.",'Index')
		evalu
		
		evalu<-evalu[order(-evalu[,4]),]
		evalu	
		evalu<-evalu[1:No_Selected_Syn2,]
		evalu

		jedi<-evalu[,2]-evalu[,1]
		x_f<-which(jedi==9)
		
		if(length(x_f)>0){
		evalu[-x_f,2]<-evalu[-x_f,2]+1		
		}
		
		if(length(x_f)==0){
		evalu[,2]<-evalu[,2]+1
		}
		

		
		write.table(evalu,file=file_Selected[Cyl])
				
		list_fun<-function() {
		list(list(), list(),list(), list())
		}

Update_data<-replicate(No_Selected_Syn2, list_fun(), simplify=FALSE)


			for(upd in 1:dim(evalu)[1]){
			set_number<-evalu[upd,1]
			Step_code<-evalu[upd,2]
			Syn_number<-evalu[upd,3]
			Update_data[[upd]][[1]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[1]]
			Update_data[[upd]][[2]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[2]]
			Update_data[[upd]][[3]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[3]]
			Update_data[[upd]][[4]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[4]]
			}
		
	        Varit<-Update_data	
		
		for (g_index in 1:length(Varit)){
		names(Varit[[g_index]])<-c('data','qtl','mrk','sequ')
		}
		
		

	# Make selected to similar output of Syn2 or 3
	for (u in 1:length(Varit)){
	
	if(u %in% x_f){
	Varit[[u]]<-Varit[[u]]
		} else {
			
			localID<-Varit[[u]]$data[,1]

			# QTL
			x_qtl<-Varit[[u]][[2]]
			x_qtl<-cbind(localID,x_qtl)
			Varit[[u]][[2]]<-x_qtl

			# MRK
			x_mrk<-Varit[[u]][[3]]
			x_mrk<-cbind(localID,x_mrk)
			Varit[[u]][[3]]<-x_mrk

			# SEQ
			x_seq<-Varit[[u]][[4]]
			x_seq<-cbind(localID,x_seq)
			Varit[[u]][[4]]<-x_seq
		}
	}
	
	# test
	# Varit[[9]][[2]][,1:10]

#F1
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)


# F2
given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[Cyl-3]][[Cyl-2]]

# F2_SP
given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)


# Syn1
F2<-Sets[[Cyl-6]][[Cyl-3]]
F2_SP<-Sets[[Cyl-6]][[Cyl-2]]

all_Syn1_out<-Global_Syn1(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
Syn1_out<-all_Syn1_out[[1]]
S_Cross<-all_Syn1_out[[2]]

# Syn2
mini_plots<-Sets[[Cyl-7]][[Cyl-1]]
Syn2_out<-Global_Syn2(mini_plots,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Syn2_W<-rnorm(5,1,10)
Syn2_Green<-Sets[[Cyl-9]][[Cyl-2]]

mult<-rnorm(5,1,10)

# Syn3
Syn_2<-Sets[[Cyl-11]][[Cyl-2]]
Syn3_out<-Global_Syn3(Syn_2,size_Syn3,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)


Sets[[14]][[Cyl]]<-F1_out
Sets[[13]][[Cyl]]<-F2_Yellow
Sets[[12]][[Cyl]]<-F2_W
Sets[[11]][[Cyl]]<-F2_Green
Sets[[10]][[Cyl]]<-F2_SP_out
Sets[[9]][[Cyl]]<-S_Cross
Sets[[8]][[Cyl]]<-Syn1_out
Sets[[8]][[Cyl-1]]<-S_Cross
Sets[[7]][[Cyl]]<-Syn2_out
Sets[[6]][[Cyl]]<-Syn2_W
Sets[[5]][[Cyl]]<-Syn2_Green
Sets[[4]][[Cyl]]<-mult
Sets[[3]][[Cyl]]<-Syn3_out


cat('All steps of year:',Cyl,'Finished',fill=TRUE)
}  



if(Cyl==15){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}

# For F2_SP
which_sets<-c(1:(Cyl-5))
which_step<-c(4:(Cyl-2))
which_sets
which_step

# For Syn2
which_sets2<-c(1:(Cyl-11))
which_step2<-c(10:(Cyl-2))
which_sets2
which_step2

which_sets<-c(which_sets2,which_sets)
which_step<-c(which_step2,which_step)
which_sets
which_step
		
     		evalu_lst<-list()
		
			for(K_1 in 1:length(which_sets)){
			dakheli1<-Sets[[which_sets[K_1]]][[which_step[K_1]]]
				evalu_mat<-matrix(ncol=4,nrow=length(dakheli1))
				for (K_2 in 1:length(dakheli1)){
				# tr1
				mean_tr1<-mean(dakheli1[[K_2]]$data$phen_1)
				# tr3
				mean_tr3<-mean(dakheli1[[K_2]]$data$phen_3)
				# tr4
				mean_tr4<-mean(dakheli1[[K_2]]$data$phen_4)
				SI<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
				evalu_mat[K_2,1]<-which_sets[K_1]
				evalu_mat[K_2,2]<-which_step[K_1]
				evalu_mat[K_2,3]<-K_2
				evalu_mat[K_2,4]<-SI
				}
				evalu_lst[[K_1]]<-evalu_mat
			}
			
			
		evalu<-do.call("rbind", evalu_lst)
		colnames(evalu)<-c("Set No.","Step","Family No.",'Index')
		evalu
		
		evalu<-evalu[order(-evalu[,4]),]
		evalu	
		evalu<-evalu[1:No_Selected_Syn2,]
		evalu

		jedi<-evalu[,2]-evalu[,1]
		x_f<-which(jedi==9)
		
		if(length(x_f)>0){
		evalu[-x_f,2]<-evalu[-x_f,2]+1		
		}
		
		if(length(x_f)==0){
		evalu[,2]<-evalu[,2]+1
		}
		

		
		write.table(evalu,file=file_Selected[Cyl])
				
		list_fun<-function() {
		list(list(), list(),list(), list())
		}

Update_data<-replicate(No_Selected_Syn2, list_fun(), simplify=FALSE)


			for(upd in 1:dim(evalu)[1]){
			set_number<-evalu[upd,1]
			Step_code<-evalu[upd,2]
			Syn_number<-evalu[upd,3]
			Update_data[[upd]][[1]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[1]]
			Update_data[[upd]][[2]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[2]]
			Update_data[[upd]][[3]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[3]]
			Update_data[[upd]][[4]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[4]]
			}
		
	        Varit<-Update_data	
		
		for (g_index in 1:length(Varit)){
		names(Varit[[g_index]])<-c('data','qtl','mrk','sequ')
		}
		
	# Make selected to similar output of Syn2 or 3
	for (u in 1:length(Varit)){
	
	if(u %in% x_f){
	Varit[[u]]<-Varit[[u]]
		} else {
			
			localID<-Varit[[u]]$data[,1]

			# QTL
			x_qtl<-Varit[[u]][[2]]
			x_qtl<-cbind(localID,x_qtl)
			Varit[[u]][[2]]<-x_qtl

			# MRK
			x_mrk<-Varit[[u]][[3]]
			x_mrk<-cbind(localID,x_mrk)
			Varit[[u]][[3]]<-x_mrk

			# SEQ
			x_seq<-Varit[[u]][[4]]
			x_seq<-cbind(localID,x_seq)
			Varit[[u]][[4]]<-x_seq
		}
	}

#F1
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)


# F2
given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[Cyl-3]][[Cyl-2]]

# F2_SP
given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)


# Syn1
F2<-Sets[[Cyl-6]][[Cyl-3]]
F2_SP<-Sets[[Cyl-6]][[Cyl-2]]

all_Syn1_out<-Global_Syn1(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
Syn1_out<-all_Syn1_out[[1]]
S_Cross<-all_Syn1_out[[2]]

# Syn2
mini_plots<-Sets[[Cyl-7]][[Cyl-1]]
Syn2_out<-Global_Syn2(mini_plots,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Syn2_W<-rnorm(5,1,10)
Syn2_Green<-Sets[[Cyl-9]][[Cyl-2]]

mult<-rnorm(5,1,10)

# Syn3
Syn_2<-Sets[[Cyl-11]][[Cyl-2]]
Syn3_out<-Global_Syn3(Syn_2,size_Syn3,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Sets[[Cyl]][[Cyl]]<-F1_out
Sets[[Cyl-1]][[Cyl]]<-F2_Yellow
Sets[[Cyl-2]][[Cyl]]<-F2_W
Sets[[Cyl-3]][[Cyl]]<-F2_Green
Sets[[Cyl-4]][[Cyl]]<-F2_SP_out
Sets[[Cyl-5]][[Cyl]]<-S_Cross
Sets[[Cyl-6]][[Cyl]]<-Syn1_out
Sets[[Cyl-6]][[Cyl-1]]<-S_Cross
Sets[[Cyl-7]][[Cyl]]<-Syn2_out
Sets[[Cyl-8]][[Cyl]]<-Syn2_W
Sets[[Cyl-9]][[Cyl]]<-Syn2_Green
Sets[[Cyl-10]][[Cyl]]<-mult
Sets[[Cyl-11]][[Cyl]]<-Syn3_out

cat('All steps of year:',Cyl,'Finished',fill=TRUE)
}  

if(Cyl==16){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}

# For F2_SP
which_sets<-c(1:(Cyl-5))
which_step<-c(4:(Cyl-2))
which_sets
which_step

# For Syn2
which_sets2<-c(1:(Cyl-11))
which_step2<-c(10:(Cyl-2))
which_sets2
which_step2

which_sets<-c(which_sets2,which_sets)
which_step<-c(which_step2,which_step)
which_sets
which_step
		
     		evalu_lst<-list()
		
			for(K_1 in 1:length(which_sets)){
			dakheli1<-Sets[[which_sets[K_1]]][[which_step[K_1]]]
				evalu_mat<-matrix(ncol=4,nrow=length(dakheli1))
				for (K_2 in 1:length(dakheli1)){
				# tr1
				mean_tr1<-mean(dakheli1[[K_2]]$data$phen_1)
				# tr3
				mean_tr3<-mean(dakheli1[[K_2]]$data$phen_3)
				# tr4
				mean_tr4<-mean(dakheli1[[K_2]]$data$phen_4)
				SI<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
				evalu_mat[K_2,1]<-which_sets[K_1]
				evalu_mat[K_2,2]<-which_step[K_1]
				evalu_mat[K_2,3]<-K_2
				evalu_mat[K_2,4]<-SI
				}
				evalu_lst[[K_1]]<-evalu_mat
			}
			
			
		evalu<-do.call("rbind", evalu_lst)
		colnames(evalu)<-c("Set No.","Step","Family No.",'Index')
		evalu
		
		evalu<-evalu[order(-evalu[,4]),]
		evalu	
		evalu<-evalu[1:No_Selected_Syn2,]
		evalu

		jedi<-evalu[,2]-evalu[,1]
		x_f<-which(jedi==9)
		
		if(length(x_f)>0){
		evalu[-x_f,2]<-evalu[-x_f,2]+1		
		}
		
		if(length(x_f)==0){
		evalu[,2]<-evalu[,2]+1
		}
		

		
		write.table(evalu,file=file_Selected[Cyl])
				
		list_fun<-function() {
		list(list(), list(),list(), list())
		}

Update_data<-replicate(No_Selected_Syn2, list_fun(), simplify=FALSE)


			for(upd in 1:dim(evalu)[1]){
			set_number<-evalu[upd,1]
			Step_code<-evalu[upd,2]
			Syn_number<-evalu[upd,3]
			Update_data[[upd]][[1]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[1]]
			Update_data[[upd]][[2]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[2]]
			Update_data[[upd]][[3]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[3]]
			Update_data[[upd]][[4]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[4]]
			}
		
	        Varit<-Update_data	
		
		for (g_index in 1:length(Varit)){
		names(Varit[[g_index]])<-c('data','qtl','mrk','sequ')
		}
		
	# Make selected to similar output of Syn2 or 3
	for (u in 1:length(Varit)){
	
	if(u %in% x_f){
	Varit[[u]]<-Varit[[u]]
		} else {
			
			localID<-Varit[[u]]$data[,1]

			# QTL
			x_qtl<-Varit[[u]][[2]]
			x_qtl<-cbind(localID,x_qtl)
			Varit[[u]][[2]]<-x_qtl

			# MRK
			x_mrk<-Varit[[u]][[3]]
			x_mrk<-cbind(localID,x_mrk)
			Varit[[u]][[3]]<-x_mrk

			# SEQ
			x_seq<-Varit[[u]][[4]]
			x_seq<-cbind(localID,x_seq)
			Varit[[u]][[4]]<-x_seq
		}
	}

#F1
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)


# F2
given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[Cyl-3]][[Cyl-2]]

# F2_SP
given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)


# Syn1
F2<-Sets[[Cyl-6]][[Cyl-3]]
F2_SP<-Sets[[Cyl-6]][[Cyl-2]]

all_Syn1_out<-Global_Syn1(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
Syn1_out<-all_Syn1_out[[1]]
S_Cross<-all_Syn1_out[[2]]

# Syn2
mini_plots<-Sets[[Cyl-7]][[Cyl-1]]
Syn2_out<-Global_Syn2(mini_plots,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Syn2_W<-rnorm(5,1,10)
Syn2_Green<-Sets[[Cyl-9]][[Cyl-2]]

mult<-rnorm(5,1,10)

# Syn3
Syn_2<-Sets[[Cyl-11]][[Cyl-2]]
Syn3_out<-Global_Syn3(Syn_2,size_Syn3,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Sets[[Cyl]][[Cyl]]<-F1_out
Sets[[Cyl-1]][[Cyl]]<-F2_Yellow
Sets[[Cyl-2]][[Cyl]]<-F2_W
Sets[[Cyl-3]][[Cyl]]<-F2_Green
Sets[[Cyl-4]][[Cyl]]<-F2_SP_out
Sets[[Cyl-5]][[Cyl]]<-S_Cross
Sets[[Cyl-6]][[Cyl]]<-Syn1_out
Sets[[Cyl-6]][[Cyl-1]]<-S_Cross
Sets[[Cyl-7]][[Cyl]]<-Syn2_out
Sets[[Cyl-8]][[Cyl]]<-Syn2_W
Sets[[Cyl-9]][[Cyl]]<-Syn2_Green
Sets[[Cyl-10]][[Cyl]]<-mult
Sets[[Cyl-11]][[Cyl]]<-Syn3_out

cat('All steps of year:',Cyl,'Finished',fill=TRUE)
}  

if(Cyl==17){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}

# For F2_SP
which_sets<-c(1:(Cyl-5))
which_step<-c(4:(Cyl-2))
which_sets
which_step

# For Syn2
which_sets2<-c(1:(Cyl-11))
which_step2<-c(10:(Cyl-2))
which_sets2
which_step2

which_sets<-c(which_sets2,which_sets)
which_step<-c(which_step2,which_step)
which_sets
which_step
		
     		evalu_lst<-list()
		
			for(K_1 in 1:length(which_sets)){
			dakheli1<-Sets[[which_sets[K_1]]][[which_step[K_1]]]
				evalu_mat<-matrix(ncol=4,nrow=length(dakheli1))
				for (K_2 in 1:length(dakheli1)){
				# tr1
				mean_tr1<-mean(dakheli1[[K_2]]$data$phen_1)
				# tr3
				mean_tr3<-mean(dakheli1[[K_2]]$data$phen_3)
				# tr4
				mean_tr4<-mean(dakheli1[[K_2]]$data$phen_4)
				SI<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
				evalu_mat[K_2,1]<-which_sets[K_1]
				evalu_mat[K_2,2]<-which_step[K_1]
				evalu_mat[K_2,3]<-K_2
				evalu_mat[K_2,4]<-SI
				}
				evalu_lst[[K_1]]<-evalu_mat
			}
			
			
		evalu<-do.call("rbind", evalu_lst)
		colnames(evalu)<-c("Set No.","Step","Family No.",'Index')
		evalu
		
		evalu<-evalu[order(-evalu[,4]),]
		evalu	
		evalu<-evalu[1:No_Selected_Syn2,]
		evalu

		jedi<-evalu[,2]-evalu[,1]
		x_f<-which(jedi==9)
		
		if(length(x_f)>0){
		evalu[-x_f,2]<-evalu[-x_f,2]+1		
		}
		
		if(length(x_f)==0){
		evalu[,2]<-evalu[,2]+1
		}
		

		
		write.table(evalu,file=file_Selected[Cyl])
				
		list_fun<-function() {
		list(list(), list(),list(), list())
		}

Update_data<-replicate(No_Selected_Syn2, list_fun(), simplify=FALSE)


			for(upd in 1:dim(evalu)[1]){
			set_number<-evalu[upd,1]
			Step_code<-evalu[upd,2]
			Syn_number<-evalu[upd,3]
			Update_data[[upd]][[1]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[1]]
			Update_data[[upd]][[2]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[2]]
			Update_data[[upd]][[3]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[3]]
			Update_data[[upd]][[4]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[4]]
			}
		
	        Varit<-Update_data	
		
		for (g_index in 1:length(Varit)){
		names(Varit[[g_index]])<-c('data','qtl','mrk','sequ')
		}
		
	# Make selected to similar output of Syn2 or 3
	for (u in 1:length(Varit)){
	
	if(u %in% x_f){
	Varit[[u]]<-Varit[[u]]
		} else {
			
			localID<-Varit[[u]]$data[,1]

			# QTL
			x_qtl<-Varit[[u]][[2]]
			x_qtl<-cbind(localID,x_qtl)
			Varit[[u]][[2]]<-x_qtl

			# MRK
			x_mrk<-Varit[[u]][[3]]
			x_mrk<-cbind(localID,x_mrk)
			Varit[[u]][[3]]<-x_mrk

			# SEQ
			x_seq<-Varit[[u]][[4]]
			x_seq<-cbind(localID,x_seq)
			Varit[[u]][[4]]<-x_seq
		}
	}

#F1
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)


# F2
given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[Cyl-3]][[Cyl-2]]

# F2_SP
given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)


# Syn1
F2<-Sets[[Cyl-6]][[Cyl-3]]
F2_SP<-Sets[[Cyl-6]][[Cyl-2]]

all_Syn1_out<-Global_Syn1(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
Syn1_out<-all_Syn1_out[[1]]
S_Cross<-all_Syn1_out[[2]]

# Syn2
mini_plots<-Sets[[Cyl-7]][[Cyl-1]]
Syn2_out<-Global_Syn2(mini_plots,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Syn2_W<-rnorm(5,1,10)
Syn2_Green<-Sets[[Cyl-9]][[Cyl-2]]

mult<-rnorm(5,1,10)

# Syn3
Syn_2<-Sets[[Cyl-11]][[Cyl-2]]
Syn3_out<-Global_Syn3(Syn_2,size_Syn3,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Sets[[Cyl]][[Cyl]]<-F1_out
Sets[[Cyl-1]][[Cyl]]<-F2_Yellow
Sets[[Cyl-2]][[Cyl]]<-F2_W
Sets[[Cyl-3]][[Cyl]]<-F2_Green
Sets[[Cyl-4]][[Cyl]]<-F2_SP_out
Sets[[Cyl-5]][[Cyl]]<-S_Cross
Sets[[Cyl-6]][[Cyl]]<-Syn1_out
Sets[[Cyl-6]][[Cyl-1]]<-S_Cross
Sets[[Cyl-7]][[Cyl]]<-Syn2_out
Sets[[Cyl-8]][[Cyl]]<-Syn2_W
Sets[[Cyl-9]][[Cyl]]<-Syn2_Green
Sets[[Cyl-10]][[Cyl]]<-mult
Sets[[Cyl-11]][[Cyl]]<-Syn3_out

cat('All steps of year:',Cyl,'Finished',fill=TRUE)
}  

if(Cyl==18){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}

# For F2_SP
which_sets<-c(1:(Cyl-5))
which_step<-c(4:(Cyl-2))
which_sets
which_step

# For Syn2
which_sets2<-c(1:(Cyl-11))
which_step2<-c(10:(Cyl-2))
which_sets2
which_step2

which_sets<-c(which_sets2,which_sets)
which_step<-c(which_step2,which_step)
which_sets
which_step
		
     		evalu_lst<-list()
		
			for(K_1 in 1:length(which_sets)){
			dakheli1<-Sets[[which_sets[K_1]]][[which_step[K_1]]]
				evalu_mat<-matrix(ncol=4,nrow=length(dakheli1))
				for (K_2 in 1:length(dakheli1)){
				# tr1
				mean_tr1<-mean(dakheli1[[K_2]]$data$phen_1)
				# tr3
				mean_tr3<-mean(dakheli1[[K_2]]$data$phen_3)
				# tr4
				mean_tr4<-mean(dakheli1[[K_2]]$data$phen_4)
				SI<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
				evalu_mat[K_2,1]<-which_sets[K_1]
				evalu_mat[K_2,2]<-which_step[K_1]
				evalu_mat[K_2,3]<-K_2
				evalu_mat[K_2,4]<-SI
				}
				evalu_lst[[K_1]]<-evalu_mat
			}
			
			
		evalu<-do.call("rbind", evalu_lst)
		colnames(evalu)<-c("Set No.","Step","Family No.",'Index')
		evalu
		
		evalu<-evalu[order(-evalu[,4]),]
		evalu	
		evalu<-evalu[1:No_Selected_Syn2,]
		evalu

		jedi<-evalu[,2]-evalu[,1]
		x_f<-which(jedi==9)
		
		if(length(x_f)>0){
		evalu[-x_f,2]<-evalu[-x_f,2]+1		
		}
		
		if(length(x_f)==0){
		evalu[,2]<-evalu[,2]+1
		}
		

		
		write.table(evalu,file=file_Selected[Cyl])
				
		list_fun<-function() {
		list(list(), list(),list(), list())
		}

Update_data<-replicate(No_Selected_Syn2, list_fun(), simplify=FALSE)


			for(upd in 1:dim(evalu)[1]){
			set_number<-evalu[upd,1]
			Step_code<-evalu[upd,2]
			Syn_number<-evalu[upd,3]
			Update_data[[upd]][[1]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[1]]
			Update_data[[upd]][[2]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[2]]
			Update_data[[upd]][[3]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[3]]
			Update_data[[upd]][[4]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[4]]
			}
		
	        Varit<-Update_data	
		
		for (g_index in 1:length(Varit)){
		names(Varit[[g_index]])<-c('data','qtl','mrk','sequ')
		}
		
	# Make selected to similar output of Syn2 or 3
	for (u in 1:length(Varit)){
	
	if(u %in% x_f){
	Varit[[u]]<-Varit[[u]]
		} else {
			
			localID<-Varit[[u]]$data[,1]

			# QTL
			x_qtl<-Varit[[u]][[2]]
			x_qtl<-cbind(localID,x_qtl)
			Varit[[u]][[2]]<-x_qtl

			# MRK
			x_mrk<-Varit[[u]][[3]]
			x_mrk<-cbind(localID,x_mrk)
			Varit[[u]][[3]]<-x_mrk

			# SEQ
			x_seq<-Varit[[u]][[4]]
			x_seq<-cbind(localID,x_seq)
			Varit[[u]][[4]]<-x_seq
		}
	}

#F1
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)


# F2
given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[Cyl-3]][[Cyl-2]]

# F2_SP
given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)


# Syn1
F2<-Sets[[Cyl-6]][[Cyl-3]]
F2_SP<-Sets[[Cyl-6]][[Cyl-2]]

all_Syn1_out<-Global_Syn1(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
Syn1_out<-all_Syn1_out[[1]]
S_Cross<-all_Syn1_out[[2]]

# Syn2
mini_plots<-Sets[[Cyl-7]][[Cyl-1]]
Syn2_out<-Global_Syn2(mini_plots,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Syn2_W<-rnorm(5,1,10)
Syn2_Green<-Sets[[Cyl-9]][[Cyl-2]]

mult<-rnorm(5,1,10)

# Syn3
Syn_2<-Sets[[Cyl-11]][[Cyl-2]]
Syn3_out<-Global_Syn3(Syn_2,size_Syn3,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Sets[[Cyl]][[Cyl]]<-F1_out
Sets[[Cyl-1]][[Cyl]]<-F2_Yellow
Sets[[Cyl-2]][[Cyl]]<-F2_W
Sets[[Cyl-3]][[Cyl]]<-F2_Green
Sets[[Cyl-4]][[Cyl]]<-F2_SP_out
Sets[[Cyl-5]][[Cyl]]<-S_Cross
Sets[[Cyl-6]][[Cyl]]<-Syn1_out
Sets[[Cyl-6]][[Cyl-1]]<-S_Cross
Sets[[Cyl-7]][[Cyl]]<-Syn2_out
Sets[[Cyl-8]][[Cyl]]<-Syn2_W
Sets[[Cyl-9]][[Cyl]]<-Syn2_Green
Sets[[Cyl-10]][[Cyl]]<-mult
Sets[[Cyl-11]][[Cyl]]<-Syn3_out

cat('All steps of year:',Cyl,'Finished',fill=TRUE)
}  

if(Cyl==19){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}

# For F2_SP
which_sets<-c(1:(Cyl-5))
which_step<-c(4:(Cyl-2))
which_sets
which_step

# For Syn2
which_sets2<-c(1:(Cyl-11))
which_step2<-c(10:(Cyl-2))
which_sets2
which_step2

which_sets<-c(which_sets2,which_sets)
which_step<-c(which_step2,which_step)
which_sets
which_step
		
     		evalu_lst<-list()
		
			for(K_1 in 1:length(which_sets)){
			dakheli1<-Sets[[which_sets[K_1]]][[which_step[K_1]]]
				evalu_mat<-matrix(ncol=4,nrow=length(dakheli1))
				for (K_2 in 1:length(dakheli1)){
				# tr1
				mean_tr1<-mean(dakheli1[[K_2]]$data$phen_1)
				# tr3
				mean_tr3<-mean(dakheli1[[K_2]]$data$phen_3)
				# tr4
				mean_tr4<-mean(dakheli1[[K_2]]$data$phen_4)
				SI<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
				evalu_mat[K_2,1]<-which_sets[K_1]
				evalu_mat[K_2,2]<-which_step[K_1]
				evalu_mat[K_2,3]<-K_2
				evalu_mat[K_2,4]<-SI
				}
				evalu_lst[[K_1]]<-evalu_mat
			}
			
			
		evalu<-do.call("rbind", evalu_lst)
		colnames(evalu)<-c("Set No.","Step","Family No.",'Index')
		evalu
		
		evalu<-evalu[order(-evalu[,4]),]
		evalu	
		evalu<-evalu[1:No_Selected_Syn2,]
		evalu

		jedi<-evalu[,2]-evalu[,1]
		x_f<-which(jedi==9)
		
		if(length(x_f)>0){
		evalu[-x_f,2]<-evalu[-x_f,2]+1		
		}
		
		if(length(x_f)==0){
		evalu[,2]<-evalu[,2]+1
		}
		

		
		write.table(evalu,file=file_Selected[Cyl])
				
		list_fun<-function() {
		list(list(), list(),list(), list())
		}

Update_data<-replicate(No_Selected_Syn2, list_fun(), simplify=FALSE)


			for(upd in 1:dim(evalu)[1]){
			set_number<-evalu[upd,1]
			Step_code<-evalu[upd,2]
			Syn_number<-evalu[upd,3]
			Update_data[[upd]][[1]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[1]]
			Update_data[[upd]][[2]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[2]]
			Update_data[[upd]][[3]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[3]]
			Update_data[[upd]][[4]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[4]]
			}
		
	        Varit<-Update_data	
		
		for (g_index in 1:length(Varit)){
		names(Varit[[g_index]])<-c('data','qtl','mrk','sequ')
		}
		
	# Make selected to similar output of Syn2 or 3
	for (u in 1:length(Varit)){
	
	if(u %in% x_f){
	Varit[[u]]<-Varit[[u]]
		} else {
			
			localID<-Varit[[u]]$data[,1]

			# QTL
			x_qtl<-Varit[[u]][[2]]
			x_qtl<-cbind(localID,x_qtl)
			Varit[[u]][[2]]<-x_qtl

			# MRK
			x_mrk<-Varit[[u]][[3]]
			x_mrk<-cbind(localID,x_mrk)
			Varit[[u]][[3]]<-x_mrk

			# SEQ
			x_seq<-Varit[[u]][[4]]
			x_seq<-cbind(localID,x_seq)
			Varit[[u]][[4]]<-x_seq
		}
	}

#F1
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)


# F2
given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[Cyl-3]][[Cyl-2]]

# F2_SP
given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)


# Syn1
F2<-Sets[[Cyl-6]][[Cyl-3]]
F2_SP<-Sets[[Cyl-6]][[Cyl-2]]

all_Syn1_out<-Global_Syn1(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
Syn1_out<-all_Syn1_out[[1]]
S_Cross<-all_Syn1_out[[2]]

# Syn2
mini_plots<-Sets[[Cyl-7]][[Cyl-1]]
Syn2_out<-Global_Syn2(mini_plots,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Syn2_W<-rnorm(5,1,10)
Syn2_Green<-Sets[[Cyl-9]][[Cyl-2]]

mult<-rnorm(5,1,10)

# Syn3
Syn_2<-Sets[[Cyl-11]][[Cyl-2]]
Syn3_out<-Global_Syn3(Syn_2,size_Syn3,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Sets[[Cyl]][[Cyl]]<-F1_out
Sets[[Cyl-1]][[Cyl]]<-F2_Yellow
Sets[[Cyl-2]][[Cyl]]<-F2_W
Sets[[Cyl-3]][[Cyl]]<-F2_Green
Sets[[Cyl-4]][[Cyl]]<-F2_SP_out
Sets[[Cyl-5]][[Cyl]]<-S_Cross
Sets[[Cyl-6]][[Cyl]]<-Syn1_out
Sets[[Cyl-6]][[Cyl-1]]<-S_Cross
Sets[[Cyl-7]][[Cyl]]<-Syn2_out
Sets[[Cyl-8]][[Cyl]]<-Syn2_W
Sets[[Cyl-9]][[Cyl]]<-Syn2_Green
Sets[[Cyl-10]][[Cyl]]<-mult
Sets[[Cyl-11]][[Cyl]]<-Syn3_out

cat('All steps of year:',Cyl,'Finished',fill=TRUE)
}  

if(Cyl==20){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}

# For F2_SP
which_sets<-c(1:(Cyl-5))
which_step<-c(4:(Cyl-2))
which_sets
which_step

# For Syn2
which_sets2<-c(1:(Cyl-11))
which_step2<-c(10:(Cyl-2))
which_sets2
which_step2

which_sets<-c(which_sets2,which_sets)
which_step<-c(which_step2,which_step)
which_sets
which_step
		
     		evalu_lst<-list()
		
			for(K_1 in 1:length(which_sets)){
			dakheli1<-Sets[[which_sets[K_1]]][[which_step[K_1]]]
				evalu_mat<-matrix(ncol=4,nrow=length(dakheli1))
				for (K_2 in 1:length(dakheli1)){
				# tr1
				mean_tr1<-mean(dakheli1[[K_2]]$data$phen_1)
				# tr3
				mean_tr3<-mean(dakheli1[[K_2]]$data$phen_3)
				# tr4
				mean_tr4<-mean(dakheli1[[K_2]]$data$phen_4)
				SI<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
				evalu_mat[K_2,1]<-which_sets[K_1]
				evalu_mat[K_2,2]<-which_step[K_1]
				evalu_mat[K_2,3]<-K_2
				evalu_mat[K_2,4]<-SI
				}
				evalu_lst[[K_1]]<-evalu_mat
			}
			
			
		evalu<-do.call("rbind", evalu_lst)
		colnames(evalu)<-c("Set No.","Step","Family No.",'Index')
		evalu
		
		evalu<-evalu[order(-evalu[,4]),]
		evalu	
		evalu<-evalu[1:No_Selected_Syn2,]
		evalu

		jedi<-evalu[,2]-evalu[,1]
		x_f<-which(jedi==9)
		
		if(length(x_f)>0){
		evalu[-x_f,2]<-evalu[-x_f,2]+1		
		}
		
		if(length(x_f)==0){
		evalu[,2]<-evalu[,2]+1
		}
		

		
		write.table(evalu,file=file_Selected[Cyl])
				
		list_fun<-function() {
		list(list(), list(),list(), list())
		}

Update_data<-replicate(No_Selected_Syn2, list_fun(), simplify=FALSE)


			for(upd in 1:dim(evalu)[1]){
			set_number<-evalu[upd,1]
			Step_code<-evalu[upd,2]
			Syn_number<-evalu[upd,3]
			Update_data[[upd]][[1]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[1]]
			Update_data[[upd]][[2]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[2]]
			Update_data[[upd]][[3]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[3]]
			Update_data[[upd]][[4]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[4]]
			}
		
	        Varit<-Update_data	
		
		for (g_index in 1:length(Varit)){
		names(Varit[[g_index]])<-c('data','qtl','mrk','sequ')
		}
		
	# Make selected to similar output of Syn2 or 3
	for (u in 1:length(Varit)){
	
	if(u %in% x_f){
	Varit[[u]]<-Varit[[u]]
		} else {
			
			localID<-Varit[[u]]$data[,1]

			# QTL
			x_qtl<-Varit[[u]][[2]]
			x_qtl<-cbind(localID,x_qtl)
			Varit[[u]][[2]]<-x_qtl

			# MRK
			x_mrk<-Varit[[u]][[3]]
			x_mrk<-cbind(localID,x_mrk)
			Varit[[u]][[3]]<-x_mrk

			# SEQ
			x_seq<-Varit[[u]][[4]]
			x_seq<-cbind(localID,x_seq)
			Varit[[u]][[4]]<-x_seq
		}
	}

#F1
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)


# F2
given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[Cyl-3]][[Cyl-2]]

# F2_SP
given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)


# Syn1
F2<-Sets[[Cyl-6]][[Cyl-3]]
F2_SP<-Sets[[Cyl-6]][[Cyl-2]]

all_Syn1_out<-Global_Syn1(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
Syn1_out<-all_Syn1_out[[1]]
S_Cross<-all_Syn1_out[[2]]

# Syn2
mini_plots<-Sets[[Cyl-7]][[Cyl-1]]
Syn2_out<-Global_Syn2(mini_plots,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Syn2_W<-rnorm(5,1,10)
Syn2_Green<-Sets[[Cyl-9]][[Cyl-2]]

mult<-rnorm(5,1,10)

# Syn3
Syn_2<-Sets[[Cyl-11]][[Cyl-2]]
Syn3_out<-Global_Syn3(Syn_2,size_Syn3,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Sets[[Cyl]][[Cyl]]<-F1_out
Sets[[Cyl-1]][[Cyl]]<-F2_Yellow
Sets[[Cyl-2]][[Cyl]]<-F2_W
Sets[[Cyl-3]][[Cyl]]<-F2_Green
Sets[[Cyl-4]][[Cyl]]<-F2_SP_out
Sets[[Cyl-5]][[Cyl]]<-S_Cross
Sets[[Cyl-6]][[Cyl]]<-Syn1_out
Sets[[Cyl-6]][[Cyl-1]]<-S_Cross
Sets[[Cyl-7]][[Cyl]]<-Syn2_out
Sets[[Cyl-8]][[Cyl]]<-Syn2_W
Sets[[Cyl-9]][[Cyl]]<-Syn2_Green
Sets[[Cyl-10]][[Cyl]]<-mult
Sets[[Cyl-11]][[Cyl]]<-Syn3_out

cat('All steps of year:',Cyl,'Finished',fill=TRUE)
}  

if(Cyl==21){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}

# For F2_SP
which_sets<-c(1:(Cyl-5))
which_step<-c(4:(Cyl-2))
which_sets
which_step

# For Syn2
which_sets2<-c(1:(Cyl-11))
which_step2<-c(10:(Cyl-2))
which_sets2
which_step2

which_sets<-c(which_sets2,which_sets)
which_step<-c(which_step2,which_step)
which_sets
which_step
		
     		evalu_lst<-list()
		
			for(K_1 in 1:length(which_sets)){
			dakheli1<-Sets[[which_sets[K_1]]][[which_step[K_1]]]
				evalu_mat<-matrix(ncol=4,nrow=length(dakheli1))
				for (K_2 in 1:length(dakheli1)){
				# tr1
				mean_tr1<-mean(dakheli1[[K_2]]$data$phen_1)
				# tr3
				mean_tr3<-mean(dakheli1[[K_2]]$data$phen_3)
				# tr4
				mean_tr4<-mean(dakheli1[[K_2]]$data$phen_4)
				SI<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
				evalu_mat[K_2,1]<-which_sets[K_1]
				evalu_mat[K_2,2]<-which_step[K_1]
				evalu_mat[K_2,3]<-K_2
				evalu_mat[K_2,4]<-SI
				}
				evalu_lst[[K_1]]<-evalu_mat
			}
			
			
		evalu<-do.call("rbind", evalu_lst)
		colnames(evalu)<-c("Set No.","Step","Family No.",'Index')
		evalu
		
		evalu<-evalu[order(-evalu[,4]),]
		evalu	
		evalu<-evalu[1:No_Selected_Syn2,]
		evalu

		jedi<-evalu[,2]-evalu[,1]
		x_f<-which(jedi==9)
		
		if(length(x_f)>0){
		evalu[-x_f,2]<-evalu[-x_f,2]+1		
		}
		
		if(length(x_f)==0){
		evalu[,2]<-evalu[,2]+1
		}
		

		
		write.table(evalu,file=file_Selected[Cyl])
				
		list_fun<-function() {
		list(list(), list(),list(), list())
		}

Update_data<-replicate(No_Selected_Syn2, list_fun(), simplify=FALSE)


			for(upd in 1:dim(evalu)[1]){
			set_number<-evalu[upd,1]
			Step_code<-evalu[upd,2]
			Syn_number<-evalu[upd,3]
			Update_data[[upd]][[1]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[1]]
			Update_data[[upd]][[2]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[2]]
			Update_data[[upd]][[3]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[3]]
			Update_data[[upd]][[4]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[4]]
			}
		
	        Varit<-Update_data	
		
		for (g_index in 1:length(Varit)){
		names(Varit[[g_index]])<-c('data','qtl','mrk','sequ')
		}
		
	# Make selected to similar output of Syn2 or 3
	for (u in 1:length(Varit)){
	
	if(u %in% x_f){
	Varit[[u]]<-Varit[[u]]
		} else {
			
			localID<-Varit[[u]]$data[,1]

			# QTL
			x_qtl<-Varit[[u]][[2]]
			x_qtl<-cbind(localID,x_qtl)
			Varit[[u]][[2]]<-x_qtl

			# MRK
			x_mrk<-Varit[[u]][[3]]
			x_mrk<-cbind(localID,x_mrk)
			Varit[[u]][[3]]<-x_mrk

			# SEQ
			x_seq<-Varit[[u]][[4]]
			x_seq<-cbind(localID,x_seq)
			Varit[[u]][[4]]<-x_seq
		}
	}

#F1
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)


# F2
given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[Cyl-3]][[Cyl-2]]

# F2_SP
given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)


# Syn1
F2<-Sets[[Cyl-6]][[Cyl-3]]
F2_SP<-Sets[[Cyl-6]][[Cyl-2]]

all_Syn1_out<-Global_Syn1(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
Syn1_out<-all_Syn1_out[[1]]
S_Cross<-all_Syn1_out[[2]]

# Syn2
mini_plots<-Sets[[Cyl-7]][[Cyl-1]]
Syn2_out<-Global_Syn2(mini_plots,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Syn2_W<-rnorm(5,1,10)
Syn2_Green<-Sets[[Cyl-9]][[Cyl-2]]

mult<-rnorm(5,1,10)

# Syn3
Syn_2<-Sets[[Cyl-11]][[Cyl-2]]
Syn3_out<-Global_Syn3(Syn_2,size_Syn3,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Sets[[Cyl]][[Cyl]]<-F1_out
Sets[[Cyl-1]][[Cyl]]<-F2_Yellow
Sets[[Cyl-2]][[Cyl]]<-F2_W
Sets[[Cyl-3]][[Cyl]]<-F2_Green
Sets[[Cyl-4]][[Cyl]]<-F2_SP_out
Sets[[Cyl-5]][[Cyl]]<-S_Cross
Sets[[Cyl-6]][[Cyl]]<-Syn1_out
Sets[[Cyl-6]][[Cyl-1]]<-S_Cross
Sets[[Cyl-7]][[Cyl]]<-Syn2_out
Sets[[Cyl-8]][[Cyl]]<-Syn2_W
Sets[[Cyl-9]][[Cyl]]<-Syn2_Green
Sets[[Cyl-10]][[Cyl]]<-mult
Sets[[Cyl-11]][[Cyl]]<-Syn3_out

cat('All steps of year:',Cyl,'Finished',fill=TRUE)
}  

if(Cyl==22){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}

# For F2_SP
which_sets<-c(1:(Cyl-5))
which_step<-c(4:(Cyl-2))
which_sets
which_step

# For Syn2
which_sets2<-c(1:(Cyl-11))
which_step2<-c(10:(Cyl-2))
which_sets2
which_step2

which_sets<-c(which_sets2,which_sets)
which_step<-c(which_step2,which_step)
which_sets
which_step
		
     		evalu_lst<-list()
		
			for(K_1 in 1:length(which_sets)){
			dakheli1<-Sets[[which_sets[K_1]]][[which_step[K_1]]]
				evalu_mat<-matrix(ncol=4,nrow=length(dakheli1))
				for (K_2 in 1:length(dakheli1)){
				# tr1
				mean_tr1<-mean(dakheli1[[K_2]]$data$phen_1)
				# tr3
				mean_tr3<-mean(dakheli1[[K_2]]$data$phen_3)
				# tr4
				mean_tr4<-mean(dakheli1[[K_2]]$data$phen_4)
				SI<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
				evalu_mat[K_2,1]<-which_sets[K_1]
				evalu_mat[K_2,2]<-which_step[K_1]
				evalu_mat[K_2,3]<-K_2
				evalu_mat[K_2,4]<-SI
				}
				evalu_lst[[K_1]]<-evalu_mat
			}
			
			
		evalu<-do.call("rbind", evalu_lst)
		colnames(evalu)<-c("Set No.","Step","Family No.",'Index')
		evalu
		
		evalu<-evalu[order(-evalu[,4]),]
		evalu	
		evalu<-evalu[1:No_Selected_Syn2,]
		evalu

		jedi<-evalu[,2]-evalu[,1]
		x_f<-which(jedi==9)
		
		if(length(x_f)>0){
		evalu[-x_f,2]<-evalu[-x_f,2]+1		
		}
		
		if(length(x_f)==0){
		evalu[,2]<-evalu[,2]+1
		}
		

		
		write.table(evalu,file=file_Selected[Cyl])
				
		list_fun<-function() {
		list(list(), list(),list(), list())
		}

Update_data<-replicate(No_Selected_Syn2, list_fun(), simplify=FALSE)


			for(upd in 1:dim(evalu)[1]){
			set_number<-evalu[upd,1]
			Step_code<-evalu[upd,2]
			Syn_number<-evalu[upd,3]
			Update_data[[upd]][[1]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[1]]
			Update_data[[upd]][[2]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[2]]
			Update_data[[upd]][[3]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[3]]
			Update_data[[upd]][[4]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[4]]
			}
		
	        Varit<-Update_data	
		
		for (g_index in 1:length(Varit)){
		names(Varit[[g_index]])<-c('data','qtl','mrk','sequ')
		}
		
	# Make selected to similar output of Syn2 or 3
	for (u in 1:length(Varit)){
	
	if(u %in% x_f){
	Varit[[u]]<-Varit[[u]]
		} else {
			
			localID<-Varit[[u]]$data[,1]

			# QTL
			x_qtl<-Varit[[u]][[2]]
			x_qtl<-cbind(localID,x_qtl)
			Varit[[u]][[2]]<-x_qtl

			# MRK
			x_mrk<-Varit[[u]][[3]]
			x_mrk<-cbind(localID,x_mrk)
			Varit[[u]][[3]]<-x_mrk

			# SEQ
			x_seq<-Varit[[u]][[4]]
			x_seq<-cbind(localID,x_seq)
			Varit[[u]][[4]]<-x_seq
		}
	}

#F1
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)


# F2
given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[Cyl-3]][[Cyl-2]]

# F2_SP
given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)


# Syn1
F2<-Sets[[Cyl-6]][[Cyl-3]]
F2_SP<-Sets[[Cyl-6]][[Cyl-2]]

all_Syn1_out<-Global_Syn1(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
Syn1_out<-all_Syn1_out[[1]]
S_Cross<-all_Syn1_out[[2]]

# Syn2
mini_plots<-Sets[[Cyl-7]][[Cyl-1]]
Syn2_out<-Global_Syn2(mini_plots,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Syn2_W<-rnorm(5,1,10)
Syn2_Green<-Sets[[Cyl-9]][[Cyl-2]]

mult<-rnorm(5,1,10)

# Syn3
Syn_2<-Sets[[Cyl-11]][[Cyl-2]]
Syn3_out<-Global_Syn3(Syn_2,size_Syn3,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Sets[[Cyl]][[Cyl]]<-F1_out
Sets[[Cyl-1]][[Cyl]]<-F2_Yellow
Sets[[Cyl-2]][[Cyl]]<-F2_W
Sets[[Cyl-3]][[Cyl]]<-F2_Green
Sets[[Cyl-4]][[Cyl]]<-F2_SP_out
Sets[[Cyl-5]][[Cyl]]<-S_Cross
Sets[[Cyl-6]][[Cyl]]<-Syn1_out
Sets[[Cyl-6]][[Cyl-1]]<-S_Cross
Sets[[Cyl-7]][[Cyl]]<-Syn2_out
Sets[[Cyl-8]][[Cyl]]<-Syn2_W
Sets[[Cyl-9]][[Cyl]]<-Syn2_Green
Sets[[Cyl-10]][[Cyl]]<-mult
Sets[[Cyl-11]][[Cyl]]<-Syn3_out

cat('All steps of year:',Cyl,'Finished',fill=TRUE)
}  

if(Cyl==23){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}

# For F2_SP
which_sets<-c(1:(Cyl-5))
which_step<-c(4:(Cyl-2))
which_sets
which_step

# For Syn2
which_sets2<-c(1:(Cyl-11))
which_step2<-c(10:(Cyl-2))
which_sets2
which_step2

which_sets<-c(which_sets2,which_sets)
which_step<-c(which_step2,which_step)
which_sets
which_step
		
     		evalu_lst<-list()
		
			for(K_1 in 1:length(which_sets)){
			dakheli1<-Sets[[which_sets[K_1]]][[which_step[K_1]]]
				evalu_mat<-matrix(ncol=4,nrow=length(dakheli1))
				for (K_2 in 1:length(dakheli1)){
				# tr1
				mean_tr1<-mean(dakheli1[[K_2]]$data$phen_1)
				# tr3
				mean_tr3<-mean(dakheli1[[K_2]]$data$phen_3)
				# tr4
				mean_tr4<-mean(dakheli1[[K_2]]$data$phen_4)
				SI<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
				evalu_mat[K_2,1]<-which_sets[K_1]
				evalu_mat[K_2,2]<-which_step[K_1]
				evalu_mat[K_2,3]<-K_2
				evalu_mat[K_2,4]<-SI
				}
				evalu_lst[[K_1]]<-evalu_mat
			}
			
			
		evalu<-do.call("rbind", evalu_lst)
		colnames(evalu)<-c("Set No.","Step","Family No.",'Index')
		evalu
		
		evalu<-evalu[order(-evalu[,4]),]
		evalu	
		evalu<-evalu[1:No_Selected_Syn2,]
		evalu

		jedi<-evalu[,2]-evalu[,1]
		x_f<-which(jedi==9)
		
		if(length(x_f)>0){
		evalu[-x_f,2]<-evalu[-x_f,2]+1		
		}
		
		if(length(x_f)==0){
		evalu[,2]<-evalu[,2]+1
		}
		

		
		write.table(evalu,file=file_Selected[Cyl])
				
		list_fun<-function() {
		list(list(), list(),list(), list())
		}

Update_data<-replicate(No_Selected_Syn2, list_fun(), simplify=FALSE)


			for(upd in 1:dim(evalu)[1]){
			set_number<-evalu[upd,1]
			Step_code<-evalu[upd,2]
			Syn_number<-evalu[upd,3]
			Update_data[[upd]][[1]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[1]]
			Update_data[[upd]][[2]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[2]]
			Update_data[[upd]][[3]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[3]]
			Update_data[[upd]][[4]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[4]]
			}
		
	        Varit<-Update_data	
		
		for (g_index in 1:length(Varit)){
		names(Varit[[g_index]])<-c('data','qtl','mrk','sequ')
		}
		
	# Make selected to similar output of Syn2 or 3
	for (u in 1:length(Varit)){
	
	if(u %in% x_f){
	Varit[[u]]<-Varit[[u]]
		} else {
			
			localID<-Varit[[u]]$data[,1]

			# QTL
			x_qtl<-Varit[[u]][[2]]
			x_qtl<-cbind(localID,x_qtl)
			Varit[[u]][[2]]<-x_qtl

			# MRK
			x_mrk<-Varit[[u]][[3]]
			x_mrk<-cbind(localID,x_mrk)
			Varit[[u]][[3]]<-x_mrk

			# SEQ
			x_seq<-Varit[[u]][[4]]
			x_seq<-cbind(localID,x_seq)
			Varit[[u]][[4]]<-x_seq
		}
	}

#F1
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)


# F2
given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[Cyl-3]][[Cyl-2]]

# F2_SP
given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)


# Syn1
F2<-Sets[[Cyl-6]][[Cyl-3]]
F2_SP<-Sets[[Cyl-6]][[Cyl-2]]

all_Syn1_out<-Global_Syn1(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
Syn1_out<-all_Syn1_out[[1]]
S_Cross<-all_Syn1_out[[2]]

# Syn2
mini_plots<-Sets[[Cyl-7]][[Cyl-1]]
Syn2_out<-Global_Syn2(mini_plots,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Syn2_W<-rnorm(5,1,10)
Syn2_Green<-Sets[[Cyl-9]][[Cyl-2]]

mult<-rnorm(5,1,10)

# Syn3
Syn_2<-Sets[[Cyl-11]][[Cyl-2]]
Syn3_out<-Global_Syn3(Syn_2,size_Syn3,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Sets[[Cyl]][[Cyl]]<-F1_out
Sets[[Cyl-1]][[Cyl]]<-F2_Yellow
Sets[[Cyl-2]][[Cyl]]<-F2_W
Sets[[Cyl-3]][[Cyl]]<-F2_Green
Sets[[Cyl-4]][[Cyl]]<-F2_SP_out
Sets[[Cyl-5]][[Cyl]]<-S_Cross
Sets[[Cyl-6]][[Cyl]]<-Syn1_out
Sets[[Cyl-6]][[Cyl-1]]<-S_Cross
Sets[[Cyl-7]][[Cyl]]<-Syn2_out
Sets[[Cyl-8]][[Cyl]]<-Syn2_W
Sets[[Cyl-9]][[Cyl]]<-Syn2_Green
Sets[[Cyl-10]][[Cyl]]<-mult
Sets[[Cyl-11]][[Cyl]]<-Syn3_out

cat('All steps of year:',Cyl,'Finished',fill=TRUE)
}  

if(Cyl==24){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}

# For F2_SP
which_sets<-c(1:(Cyl-5))
which_step<-c(4:(Cyl-2))
which_sets
which_step

# For Syn2
which_sets2<-c(1:(Cyl-11))
which_step2<-c(10:(Cyl-2))
which_sets2
which_step2

which_sets<-c(which_sets2,which_sets)
which_step<-c(which_step2,which_step)
which_sets
which_step
		
     		evalu_lst<-list()
		
			for(K_1 in 1:length(which_sets)){
			dakheli1<-Sets[[which_sets[K_1]]][[which_step[K_1]]]
				evalu_mat<-matrix(ncol=4,nrow=length(dakheli1))
				for (K_2 in 1:length(dakheli1)){
				# tr1
				mean_tr1<-mean(dakheli1[[K_2]]$data$phen_1)
				# tr3
				mean_tr3<-mean(dakheli1[[K_2]]$data$phen_3)
				# tr4
				mean_tr4<-mean(dakheli1[[K_2]]$data$phen_4)
				SI<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
				evalu_mat[K_2,1]<-which_sets[K_1]
				evalu_mat[K_2,2]<-which_step[K_1]
				evalu_mat[K_2,3]<-K_2
				evalu_mat[K_2,4]<-SI
				}
				evalu_lst[[K_1]]<-evalu_mat
			}
			
			
		evalu<-do.call("rbind", evalu_lst)
		colnames(evalu)<-c("Set No.","Step","Family No.",'Index')
		evalu
		
		evalu<-evalu[order(-evalu[,4]),]
		evalu	
		evalu<-evalu[1:No_Selected_Syn2,]
		evalu

		jedi<-evalu[,2]-evalu[,1]
		x_f<-which(jedi==9)
		
		if(length(x_f)>0){
		evalu[-x_f,2]<-evalu[-x_f,2]+1		
		}
		
		if(length(x_f)==0){
		evalu[,2]<-evalu[,2]+1
		}
		

		
		write.table(evalu,file=file_Selected[Cyl])
				
		list_fun<-function() {
		list(list(), list(),list(), list())
		}

Update_data<-replicate(No_Selected_Syn2, list_fun(), simplify=FALSE)


			for(upd in 1:dim(evalu)[1]){
			set_number<-evalu[upd,1]
			Step_code<-evalu[upd,2]
			Syn_number<-evalu[upd,3]
			Update_data[[upd]][[1]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[1]]
			Update_data[[upd]][[2]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[2]]
			Update_data[[upd]][[3]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[3]]
			Update_data[[upd]][[4]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[4]]
			}
		
	        Varit<-Update_data	
		
		for (g_index in 1:length(Varit)){
		names(Varit[[g_index]])<-c('data','qtl','mrk','sequ')
		}
		
	# Make selected to similar output of Syn2 or 3
	for (u in 1:length(Varit)){
	
	if(u %in% x_f){
	Varit[[u]]<-Varit[[u]]
		} else {
			
			localID<-Varit[[u]]$data[,1]

			# QTL
			x_qtl<-Varit[[u]][[2]]
			x_qtl<-cbind(localID,x_qtl)
			Varit[[u]][[2]]<-x_qtl

			# MRK
			x_mrk<-Varit[[u]][[3]]
			x_mrk<-cbind(localID,x_mrk)
			Varit[[u]][[3]]<-x_mrk

			# SEQ
			x_seq<-Varit[[u]][[4]]
			x_seq<-cbind(localID,x_seq)
			Varit[[u]][[4]]<-x_seq
		}
	}

#F1
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)


# F2
given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[Cyl-3]][[Cyl-2]]

# F2_SP
given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)


# Syn1
F2<-Sets[[Cyl-6]][[Cyl-3]]
F2_SP<-Sets[[Cyl-6]][[Cyl-2]]

all_Syn1_out<-Global_Syn1(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
Syn1_out<-all_Syn1_out[[1]]
S_Cross<-all_Syn1_out[[2]]

# Syn2
mini_plots<-Sets[[Cyl-7]][[Cyl-1]]
Syn2_out<-Global_Syn2(mini_plots,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Syn2_W<-rnorm(5,1,10)
Syn2_Green<-Sets[[Cyl-9]][[Cyl-2]]

mult<-rnorm(5,1,10)

# Syn3
Syn_2<-Sets[[Cyl-11]][[Cyl-2]]
Syn3_out<-Global_Syn3(Syn_2,size_Syn3,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Sets[[Cyl]][[Cyl]]<-F1_out
Sets[[Cyl-1]][[Cyl]]<-F2_Yellow
Sets[[Cyl-2]][[Cyl]]<-F2_W
Sets[[Cyl-3]][[Cyl]]<-F2_Green
Sets[[Cyl-4]][[Cyl]]<-F2_SP_out
Sets[[Cyl-5]][[Cyl]]<-S_Cross
Sets[[Cyl-6]][[Cyl]]<-Syn1_out
Sets[[Cyl-6]][[Cyl-1]]<-S_Cross
Sets[[Cyl-7]][[Cyl]]<-Syn2_out
Sets[[Cyl-8]][[Cyl]]<-Syn2_W
Sets[[Cyl-9]][[Cyl]]<-Syn2_Green
Sets[[Cyl-10]][[Cyl]]<-mult
Sets[[Cyl-11]][[Cyl]]<-Syn3_out

cat('All steps of year:',Cyl,'Finished',fill=TRUE)
}  

if(Cyl==25){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}

# For F2_SP
which_sets<-c(1:(Cyl-5))
which_step<-c(4:(Cyl-2))
which_sets
which_step

# For Syn2
which_sets2<-c(1:(Cyl-11))
which_step2<-c(10:(Cyl-2))
which_sets2
which_step2

which_sets<-c(which_sets2,which_sets)
which_step<-c(which_step2,which_step)
which_sets
which_step
		
     		evalu_lst<-list()
		
			for(K_1 in 1:length(which_sets)){
			dakheli1<-Sets[[which_sets[K_1]]][[which_step[K_1]]]
				evalu_mat<-matrix(ncol=4,nrow=length(dakheli1))
				for (K_2 in 1:length(dakheli1)){
				# tr1
				mean_tr1<-mean(dakheli1[[K_2]]$data$phen_1)
				# tr3
				mean_tr3<-mean(dakheli1[[K_2]]$data$phen_3)
				# tr4
				mean_tr4<-mean(dakheli1[[K_2]]$data$phen_4)
				SI<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
				evalu_mat[K_2,1]<-which_sets[K_1]
				evalu_mat[K_2,2]<-which_step[K_1]
				evalu_mat[K_2,3]<-K_2
				evalu_mat[K_2,4]<-SI
				}
				evalu_lst[[K_1]]<-evalu_mat
			}
			
			
		evalu<-do.call("rbind", evalu_lst)
		colnames(evalu)<-c("Set No.","Step","Family No.",'Index')
		evalu
		
		evalu<-evalu[order(-evalu[,4]),]
		evalu	
		evalu<-evalu[1:No_Selected_Syn2,]
		evalu

		jedi<-evalu[,2]-evalu[,1]
		x_f<-which(jedi==9)
		
		if(length(x_f)>0){
		evalu[-x_f,2]<-evalu[-x_f,2]+1		
		}
		
		if(length(x_f)==0){
		evalu[,2]<-evalu[,2]+1
		}
		

		
		write.table(evalu,file=file_Selected[Cyl])
				
		list_fun<-function() {
		list(list(), list(),list(), list())
		}

Update_data<-replicate(No_Selected_Syn2, list_fun(), simplify=FALSE)


			for(upd in 1:dim(evalu)[1]){
			set_number<-evalu[upd,1]
			Step_code<-evalu[upd,2]
			Syn_number<-evalu[upd,3]
			Update_data[[upd]][[1]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[1]]
			Update_data[[upd]][[2]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[2]]
			Update_data[[upd]][[3]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[3]]
			Update_data[[upd]][[4]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[4]]
			}
		
	        Varit<-Update_data	
		
		for (g_index in 1:length(Varit)){
		names(Varit[[g_index]])<-c('data','qtl','mrk','sequ')
		}
		
	# Make selected to similar output of Syn2 or 3
	for (u in 1:length(Varit)){
	
	if(u %in% x_f){
	Varit[[u]]<-Varit[[u]]
		} else {
			
			localID<-Varit[[u]]$data[,1]

			# QTL
			x_qtl<-Varit[[u]][[2]]
			x_qtl<-cbind(localID,x_qtl)
			Varit[[u]][[2]]<-x_qtl

			# MRK
			x_mrk<-Varit[[u]][[3]]
			x_mrk<-cbind(localID,x_mrk)
			Varit[[u]][[3]]<-x_mrk

			# SEQ
			x_seq<-Varit[[u]][[4]]
			x_seq<-cbind(localID,x_seq)
			Varit[[u]][[4]]<-x_seq
		}
	}

#F1
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)


# F2
given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[Cyl-3]][[Cyl-2]]

# F2_SP
given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)


# Syn1
F2<-Sets[[Cyl-6]][[Cyl-3]]
F2_SP<-Sets[[Cyl-6]][[Cyl-2]]

all_Syn1_out<-Global_Syn1(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
Syn1_out<-all_Syn1_out[[1]]
S_Cross<-all_Syn1_out[[2]]

# Syn2
mini_plots<-Sets[[Cyl-7]][[Cyl-1]]
Syn2_out<-Global_Syn2(mini_plots,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Syn2_W<-rnorm(5,1,10)
Syn2_Green<-Sets[[Cyl-9]][[Cyl-2]]

mult<-rnorm(5,1,10)

# Syn3
Syn_2<-Sets[[Cyl-11]][[Cyl-2]]
Syn3_out<-Global_Syn3(Syn_2,size_Syn3,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Sets[[Cyl]][[Cyl]]<-F1_out
Sets[[Cyl-1]][[Cyl]]<-F2_Yellow
Sets[[Cyl-2]][[Cyl]]<-F2_W
Sets[[Cyl-3]][[Cyl]]<-F2_Green
Sets[[Cyl-4]][[Cyl]]<-F2_SP_out
Sets[[Cyl-5]][[Cyl]]<-S_Cross
Sets[[Cyl-6]][[Cyl]]<-Syn1_out
Sets[[Cyl-6]][[Cyl-1]]<-S_Cross
Sets[[Cyl-7]][[Cyl]]<-Syn2_out
Sets[[Cyl-8]][[Cyl]]<-Syn2_W
Sets[[Cyl-9]][[Cyl]]<-Syn2_Green
Sets[[Cyl-10]][[Cyl]]<-mult
Sets[[Cyl-11]][[Cyl]]<-Syn3_out

cat('All steps of year:',Cyl,'Finished',fill=TRUE)
}  

if(Cyl==26){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}

# For F2_SP
which_sets<-c(1:(Cyl-5))
which_step<-c(4:(Cyl-2))
which_sets
which_step

# For Syn2
which_sets2<-c(1:(Cyl-11))
which_step2<-c(10:(Cyl-2))
which_sets2
which_step2

which_sets<-c(which_sets2,which_sets)
which_step<-c(which_step2,which_step)
which_sets
which_step
		
     		evalu_lst<-list()
		
			for(K_1 in 1:length(which_sets)){
			dakheli1<-Sets[[which_sets[K_1]]][[which_step[K_1]]]
				evalu_mat<-matrix(ncol=4,nrow=length(dakheli1))
				for (K_2 in 1:length(dakheli1)){
				# tr1
				mean_tr1<-mean(dakheli1[[K_2]]$data$phen_1)
				# tr3
				mean_tr3<-mean(dakheli1[[K_2]]$data$phen_3)
				# tr4
				mean_tr4<-mean(dakheli1[[K_2]]$data$phen_4)
				SI<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
				evalu_mat[K_2,1]<-which_sets[K_1]
				evalu_mat[K_2,2]<-which_step[K_1]
				evalu_mat[K_2,3]<-K_2
				evalu_mat[K_2,4]<-SI
				}
				evalu_lst[[K_1]]<-evalu_mat
			}
			
			
		evalu<-do.call("rbind", evalu_lst)
		colnames(evalu)<-c("Set No.","Step","Family No.",'Index')
		evalu
		
		evalu<-evalu[order(-evalu[,4]),]
		evalu	
		evalu<-evalu[1:No_Selected_Syn2,]
		evalu

		jedi<-evalu[,2]-evalu[,1]
		x_f<-which(jedi==9)
		
		if(length(x_f)>0){
		evalu[-x_f,2]<-evalu[-x_f,2]+1		
		}
		
		if(length(x_f)==0){
		evalu[,2]<-evalu[,2]+1
		}
		

		
		write.table(evalu,file=file_Selected[Cyl])
				
		list_fun<-function() {
		list(list(), list(),list(), list())
		}

Update_data<-replicate(No_Selected_Syn2, list_fun(), simplify=FALSE)


			for(upd in 1:dim(evalu)[1]){
			set_number<-evalu[upd,1]
			Step_code<-evalu[upd,2]
			Syn_number<-evalu[upd,3]
			Update_data[[upd]][[1]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[1]]
			Update_data[[upd]][[2]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[2]]
			Update_data[[upd]][[3]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[3]]
			Update_data[[upd]][[4]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[4]]
			}
		
	        Varit<-Update_data	
		
		for (g_index in 1:length(Varit)){
		names(Varit[[g_index]])<-c('data','qtl','mrk','sequ')
		}
		
	# Make selected to similar output of Syn2 or 3
	for (u in 1:length(Varit)){
	
	if(u %in% x_f){
	Varit[[u]]<-Varit[[u]]
		} else {
			
			localID<-Varit[[u]]$data[,1]

			# QTL
			x_qtl<-Varit[[u]][[2]]
			x_qtl<-cbind(localID,x_qtl)
			Varit[[u]][[2]]<-x_qtl

			# MRK
			x_mrk<-Varit[[u]][[3]]
			x_mrk<-cbind(localID,x_mrk)
			Varit[[u]][[3]]<-x_mrk

			# SEQ
			x_seq<-Varit[[u]][[4]]
			x_seq<-cbind(localID,x_seq)
			Varit[[u]][[4]]<-x_seq
		}
	}

#F1
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)


# F2
given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[Cyl-3]][[Cyl-2]]

# F2_SP
given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)


# Syn1
F2<-Sets[[Cyl-6]][[Cyl-3]]
F2_SP<-Sets[[Cyl-6]][[Cyl-2]]

all_Syn1_out<-Global_Syn1(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
Syn1_out<-all_Syn1_out[[1]]
S_Cross<-all_Syn1_out[[2]]

# Syn2
mini_plots<-Sets[[Cyl-7]][[Cyl-1]]
Syn2_out<-Global_Syn2(mini_plots,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Syn2_W<-rnorm(5,1,10)
Syn2_Green<-Sets[[Cyl-9]][[Cyl-2]]

mult<-rnorm(5,1,10)

# Syn3
Syn_2<-Sets[[Cyl-11]][[Cyl-2]]
Syn3_out<-Global_Syn3(Syn_2,size_Syn3,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Sets[[Cyl]][[Cyl]]<-F1_out
Sets[[Cyl-1]][[Cyl]]<-F2_Yellow
Sets[[Cyl-2]][[Cyl]]<-F2_W
Sets[[Cyl-3]][[Cyl]]<-F2_Green
Sets[[Cyl-4]][[Cyl]]<-F2_SP_out
Sets[[Cyl-5]][[Cyl]]<-S_Cross
Sets[[Cyl-6]][[Cyl]]<-Syn1_out
Sets[[Cyl-6]][[Cyl-1]]<-S_Cross
Sets[[Cyl-7]][[Cyl]]<-Syn2_out
Sets[[Cyl-8]][[Cyl]]<-Syn2_W
Sets[[Cyl-9]][[Cyl]]<-Syn2_Green
Sets[[Cyl-10]][[Cyl]]<-mult
Sets[[Cyl-11]][[Cyl]]<-Syn3_out

cat('All steps of year:',Cyl,'Finished',fill=TRUE)
}  

if(Cyl==27){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}

# For F2_SP
which_sets<-c(1:(Cyl-5))
which_step<-c(4:(Cyl-2))
which_sets
which_step

# For Syn2
which_sets2<-c(1:(Cyl-11))
which_step2<-c(10:(Cyl-2))
which_sets2
which_step2

which_sets<-c(which_sets2,which_sets)
which_step<-c(which_step2,which_step)
which_sets
which_step
		
     		evalu_lst<-list()
		
			for(K_1 in 1:length(which_sets)){
			dakheli1<-Sets[[which_sets[K_1]]][[which_step[K_1]]]
				evalu_mat<-matrix(ncol=4,nrow=length(dakheli1))
				for (K_2 in 1:length(dakheli1)){
				# tr1
				mean_tr1<-mean(dakheli1[[K_2]]$data$phen_1)
				# tr3
				mean_tr3<-mean(dakheli1[[K_2]]$data$phen_3)
				# tr4
				mean_tr4<-mean(dakheli1[[K_2]]$data$phen_4)
				SI<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
				evalu_mat[K_2,1]<-which_sets[K_1]
				evalu_mat[K_2,2]<-which_step[K_1]
				evalu_mat[K_2,3]<-K_2
				evalu_mat[K_2,4]<-SI
				}
				evalu_lst[[K_1]]<-evalu_mat
			}
			
			
		evalu<-do.call("rbind", evalu_lst)
		colnames(evalu)<-c("Set No.","Step","Family No.",'Index')
		evalu
		
		evalu<-evalu[order(-evalu[,4]),]
		evalu	
		evalu<-evalu[1:No_Selected_Syn2,]
		evalu

		jedi<-evalu[,2]-evalu[,1]
		x_f<-which(jedi==9)
		
		if(length(x_f)>0){
		evalu[-x_f,2]<-evalu[-x_f,2]+1		
		}
		
		if(length(x_f)==0){
		evalu[,2]<-evalu[,2]+1
		}
		

		
		write.table(evalu,file=file_Selected[Cyl])
				
		list_fun<-function() {
		list(list(), list(),list(), list())
		}

Update_data<-replicate(No_Selected_Syn2, list_fun(), simplify=FALSE)


			for(upd in 1:dim(evalu)[1]){
			set_number<-evalu[upd,1]
			Step_code<-evalu[upd,2]
			Syn_number<-evalu[upd,3]
			Update_data[[upd]][[1]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[1]]
			Update_data[[upd]][[2]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[2]]
			Update_data[[upd]][[3]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[3]]
			Update_data[[upd]][[4]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[4]]
			}
		
	        Varit<-Update_data	
		
		for (g_index in 1:length(Varit)){
		names(Varit[[g_index]])<-c('data','qtl','mrk','sequ')
		}
		
	# Make selected to similar output of Syn2 or 3
	for (u in 1:length(Varit)){
	
	if(u %in% x_f){
	Varit[[u]]<-Varit[[u]]
		} else {
			
			localID<-Varit[[u]]$data[,1]

			# QTL
			x_qtl<-Varit[[u]][[2]]
			x_qtl<-cbind(localID,x_qtl)
			Varit[[u]][[2]]<-x_qtl

			# MRK
			x_mrk<-Varit[[u]][[3]]
			x_mrk<-cbind(localID,x_mrk)
			Varit[[u]][[3]]<-x_mrk

			# SEQ
			x_seq<-Varit[[u]][[4]]
			x_seq<-cbind(localID,x_seq)
			Varit[[u]][[4]]<-x_seq
		}
	}

#F1
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)


# F2
given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[Cyl-3]][[Cyl-2]]

# F2_SP
given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)


# Syn1
F2<-Sets[[Cyl-6]][[Cyl-3]]
F2_SP<-Sets[[Cyl-6]][[Cyl-2]]

all_Syn1_out<-Global_Syn1(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
Syn1_out<-all_Syn1_out[[1]]
S_Cross<-all_Syn1_out[[2]]

# Syn2
mini_plots<-Sets[[Cyl-7]][[Cyl-1]]
Syn2_out<-Global_Syn2(mini_plots,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Syn2_W<-rnorm(5,1,10)
Syn2_Green<-Sets[[Cyl-9]][[Cyl-2]]

mult<-rnorm(5,1,10)

# Syn3
Syn_2<-Sets[[Cyl-11]][[Cyl-2]]
Syn3_out<-Global_Syn3(Syn_2,size_Syn3,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Sets[[Cyl]][[Cyl]]<-F1_out
Sets[[Cyl-1]][[Cyl]]<-F2_Yellow
Sets[[Cyl-2]][[Cyl]]<-F2_W
Sets[[Cyl-3]][[Cyl]]<-F2_Green
Sets[[Cyl-4]][[Cyl]]<-F2_SP_out
Sets[[Cyl-5]][[Cyl]]<-S_Cross
Sets[[Cyl-6]][[Cyl]]<-Syn1_out
Sets[[Cyl-6]][[Cyl-1]]<-S_Cross
Sets[[Cyl-7]][[Cyl]]<-Syn2_out
Sets[[Cyl-8]][[Cyl]]<-Syn2_W
Sets[[Cyl-9]][[Cyl]]<-Syn2_Green
Sets[[Cyl-10]][[Cyl]]<-mult
Sets[[Cyl-11]][[Cyl]]<-Syn3_out

cat('All steps of year:',Cyl,'Finished',fill=TRUE)
}  

if(Cyl==28){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}

# For F2_SP
which_sets<-c(1:(Cyl-5))
which_step<-c(4:(Cyl-2))
which_sets
which_step

# For Syn2
which_sets2<-c(1:(Cyl-11))
which_step2<-c(10:(Cyl-2))
which_sets2
which_step2

which_sets<-c(which_sets2,which_sets)
which_step<-c(which_step2,which_step)
which_sets
which_step
		
     		evalu_lst<-list()
		
			for(K_1 in 1:length(which_sets)){
			dakheli1<-Sets[[which_sets[K_1]]][[which_step[K_1]]]
				evalu_mat<-matrix(ncol=4,nrow=length(dakheli1))
				for (K_2 in 1:length(dakheli1)){
				# tr1
				mean_tr1<-mean(dakheli1[[K_2]]$data$phen_1)
				# tr3
				mean_tr3<-mean(dakheli1[[K_2]]$data$phen_3)
				# tr4
				mean_tr4<-mean(dakheli1[[K_2]]$data$phen_4)
				SI<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
				evalu_mat[K_2,1]<-which_sets[K_1]
				evalu_mat[K_2,2]<-which_step[K_1]
				evalu_mat[K_2,3]<-K_2
				evalu_mat[K_2,4]<-SI
				}
				evalu_lst[[K_1]]<-evalu_mat
			}
			
			
		evalu<-do.call("rbind", evalu_lst)
		colnames(evalu)<-c("Set No.","Step","Family No.",'Index')
		evalu
		
		evalu<-evalu[order(-evalu[,4]),]
		evalu	
		evalu<-evalu[1:No_Selected_Syn2,]
		evalu

		jedi<-evalu[,2]-evalu[,1]
		x_f<-which(jedi==9)
		
		if(length(x_f)>0){
		evalu[-x_f,2]<-evalu[-x_f,2]+1		
		}
		
		if(length(x_f)==0){
		evalu[,2]<-evalu[,2]+1
		}
		

		
		write.table(evalu,file=file_Selected[Cyl])
				
		list_fun<-function() {
		list(list(), list(),list(), list())
		}

Update_data<-replicate(No_Selected_Syn2, list_fun(), simplify=FALSE)


			for(upd in 1:dim(evalu)[1]){
			set_number<-evalu[upd,1]
			Step_code<-evalu[upd,2]
			Syn_number<-evalu[upd,3]
			Update_data[[upd]][[1]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[1]]
			Update_data[[upd]][[2]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[2]]
			Update_data[[upd]][[3]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[3]]
			Update_data[[upd]][[4]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[4]]
			}
		
	        Varit<-Update_data	
		
		for (g_index in 1:length(Varit)){
		names(Varit[[g_index]])<-c('data','qtl','mrk','sequ')
		}
		
	# Make selected to similar output of Syn2 or 3
	for (u in 1:length(Varit)){
	
	if(u %in% x_f){
	Varit[[u]]<-Varit[[u]]
		} else {
			
			localID<-Varit[[u]]$data[,1]

			# QTL
			x_qtl<-Varit[[u]][[2]]
			x_qtl<-cbind(localID,x_qtl)
			Varit[[u]][[2]]<-x_qtl

			# MRK
			x_mrk<-Varit[[u]][[3]]
			x_mrk<-cbind(localID,x_mrk)
			Varit[[u]][[3]]<-x_mrk

			# SEQ
			x_seq<-Varit[[u]][[4]]
			x_seq<-cbind(localID,x_seq)
			Varit[[u]][[4]]<-x_seq
		}
	}

#F1
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)


# F2
given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[Cyl-3]][[Cyl-2]]

# F2_SP
given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)


# Syn1
F2<-Sets[[Cyl-6]][[Cyl-3]]
F2_SP<-Sets[[Cyl-6]][[Cyl-2]]

all_Syn1_out<-Global_Syn1(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
Syn1_out<-all_Syn1_out[[1]]
S_Cross<-all_Syn1_out[[2]]

# Syn2
mini_plots<-Sets[[Cyl-7]][[Cyl-1]]
Syn2_out<-Global_Syn2(mini_plots,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Syn2_W<-rnorm(5,1,10)
Syn2_Green<-Sets[[Cyl-9]][[Cyl-2]]

mult<-rnorm(5,1,10)

# Syn3
Syn_2<-Sets[[Cyl-11]][[Cyl-2]]
Syn3_out<-Global_Syn3(Syn_2,size_Syn3,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Sets[[Cyl]][[Cyl]]<-F1_out
Sets[[Cyl-1]][[Cyl]]<-F2_Yellow
Sets[[Cyl-2]][[Cyl]]<-F2_W
Sets[[Cyl-3]][[Cyl]]<-F2_Green
Sets[[Cyl-4]][[Cyl]]<-F2_SP_out
Sets[[Cyl-5]][[Cyl]]<-S_Cross
Sets[[Cyl-6]][[Cyl]]<-Syn1_out
Sets[[Cyl-6]][[Cyl-1]]<-S_Cross
Sets[[Cyl-7]][[Cyl]]<-Syn2_out
Sets[[Cyl-8]][[Cyl]]<-Syn2_W
Sets[[Cyl-9]][[Cyl]]<-Syn2_Green
Sets[[Cyl-10]][[Cyl]]<-mult
Sets[[Cyl-11]][[Cyl]]<-Syn3_out

cat('All steps of year:',Cyl,'Finished',fill=TRUE)
}  

if(Cyl==29){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}

# For F2_SP
which_sets<-c(1:(Cyl-5))
which_step<-c(4:(Cyl-2))
which_sets
which_step

# For Syn2
which_sets2<-c(1:(Cyl-11))
which_step2<-c(10:(Cyl-2))
which_sets2
which_step2

which_sets<-c(which_sets2,which_sets)
which_step<-c(which_step2,which_step)
which_sets
which_step
		
     		evalu_lst<-list()
		
			for(K_1 in 1:length(which_sets)){
			dakheli1<-Sets[[which_sets[K_1]]][[which_step[K_1]]]
				evalu_mat<-matrix(ncol=4,nrow=length(dakheli1))
				for (K_2 in 1:length(dakheli1)){
				# tr1
				mean_tr1<-mean(dakheli1[[K_2]]$data$phen_1)
				# tr3
				mean_tr3<-mean(dakheli1[[K_2]]$data$phen_3)
				# tr4
				mean_tr4<-mean(dakheli1[[K_2]]$data$phen_4)
				SI<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
				evalu_mat[K_2,1]<-which_sets[K_1]
				evalu_mat[K_2,2]<-which_step[K_1]
				evalu_mat[K_2,3]<-K_2
				evalu_mat[K_2,4]<-SI
				}
				evalu_lst[[K_1]]<-evalu_mat
			}
			
			
		evalu<-do.call("rbind", evalu_lst)
		colnames(evalu)<-c("Set No.","Step","Family No.",'Index')
		evalu
		
		evalu<-evalu[order(-evalu[,4]),]
		evalu	
		evalu<-evalu[1:No_Selected_Syn2,]
		evalu

		jedi<-evalu[,2]-evalu[,1]
		x_f<-which(jedi==9)
		
		if(length(x_f)>0){
		evalu[-x_f,2]<-evalu[-x_f,2]+1		
		}
		
		if(length(x_f)==0){
		evalu[,2]<-evalu[,2]+1
		}
		

		
		write.table(evalu,file=file_Selected[Cyl])
				
		list_fun<-function() {
		list(list(), list(),list(), list())
		}

Update_data<-replicate(No_Selected_Syn2, list_fun(), simplify=FALSE)


			for(upd in 1:dim(evalu)[1]){
			set_number<-evalu[upd,1]
			Step_code<-evalu[upd,2]
			Syn_number<-evalu[upd,3]
			Update_data[[upd]][[1]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[1]]
			Update_data[[upd]][[2]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[2]]
			Update_data[[upd]][[3]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[3]]
			Update_data[[upd]][[4]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[4]]
			}
		
	        Varit<-Update_data	
		
		for (g_index in 1:length(Varit)){
		names(Varit[[g_index]])<-c('data','qtl','mrk','sequ')
		}
		
	# Make selected to similar output of Syn2 or 3
	for (u in 1:length(Varit)){
	
	if(u %in% x_f){
	Varit[[u]]<-Varit[[u]]
		} else {
			
			localID<-Varit[[u]]$data[,1]

			# QTL
			x_qtl<-Varit[[u]][[2]]
			x_qtl<-cbind(localID,x_qtl)
			Varit[[u]][[2]]<-x_qtl

			# MRK
			x_mrk<-Varit[[u]][[3]]
			x_mrk<-cbind(localID,x_mrk)
			Varit[[u]][[3]]<-x_mrk

			# SEQ
			x_seq<-Varit[[u]][[4]]
			x_seq<-cbind(localID,x_seq)
			Varit[[u]][[4]]<-x_seq
		}
	}

#F1
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)


# F2
given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[Cyl-3]][[Cyl-2]]

# F2_SP
given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)


# Syn1
F2<-Sets[[Cyl-6]][[Cyl-3]]
F2_SP<-Sets[[Cyl-6]][[Cyl-2]]

all_Syn1_out<-Global_Syn1(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
Syn1_out<-all_Syn1_out[[1]]
S_Cross<-all_Syn1_out[[2]]

# Syn2
mini_plots<-Sets[[Cyl-7]][[Cyl-1]]
Syn2_out<-Global_Syn2(mini_plots,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Syn2_W<-rnorm(5,1,10)
Syn2_Green<-Sets[[Cyl-9]][[Cyl-2]]

mult<-rnorm(5,1,10)

# Syn3
Syn_2<-Sets[[Cyl-11]][[Cyl-2]]
Syn3_out<-Global_Syn3(Syn_2,size_Syn3,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Sets[[Cyl]][[Cyl]]<-F1_out
Sets[[Cyl-1]][[Cyl]]<-F2_Yellow
Sets[[Cyl-2]][[Cyl]]<-F2_W
Sets[[Cyl-3]][[Cyl]]<-F2_Green
Sets[[Cyl-4]][[Cyl]]<-F2_SP_out
Sets[[Cyl-5]][[Cyl]]<-S_Cross
Sets[[Cyl-6]][[Cyl]]<-Syn1_out
Sets[[Cyl-6]][[Cyl-1]]<-S_Cross
Sets[[Cyl-7]][[Cyl]]<-Syn2_out
Sets[[Cyl-8]][[Cyl]]<-Syn2_W
Sets[[Cyl-9]][[Cyl]]<-Syn2_Green
Sets[[Cyl-10]][[Cyl]]<-mult
Sets[[Cyl-11]][[Cyl]]<-Syn3_out

cat('All steps of year:',Cyl,'Finished',fill=TRUE)
}  

if(Cyl==30){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}

# For F2_SP
which_sets<-c(1:(Cyl-5))
which_step<-c(4:(Cyl-2))
which_sets
which_step

# For Syn2
which_sets2<-c(1:(Cyl-11))
which_step2<-c(10:(Cyl-2))
which_sets2
which_step2

which_sets<-c(which_sets2,which_sets)
which_step<-c(which_step2,which_step)
which_sets
which_step
		
     		evalu_lst<-list()
		
			for(K_1 in 1:length(which_sets)){
			dakheli1<-Sets[[which_sets[K_1]]][[which_step[K_1]]]
				evalu_mat<-matrix(ncol=4,nrow=length(dakheli1))
				for (K_2 in 1:length(dakheli1)){
				# tr1
				mean_tr1<-mean(dakheli1[[K_2]]$data$phen_1)
				# tr3
				mean_tr3<-mean(dakheli1[[K_2]]$data$phen_3)
				# tr4
				mean_tr4<-mean(dakheli1[[K_2]]$data$phen_4)
				SI<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
				evalu_mat[K_2,1]<-which_sets[K_1]
				evalu_mat[K_2,2]<-which_step[K_1]
				evalu_mat[K_2,3]<-K_2
				evalu_mat[K_2,4]<-SI
				}
				evalu_lst[[K_1]]<-evalu_mat
			}
			
			
		evalu<-do.call("rbind", evalu_lst)
		colnames(evalu)<-c("Set No.","Step","Family No.",'Index')
		evalu
		
		evalu<-evalu[order(-evalu[,4]),]
		evalu	
		evalu<-evalu[1:No_Selected_Syn2,]
		evalu

		jedi<-evalu[,2]-evalu[,1]
		x_f<-which(jedi==9)
		
		if(length(x_f)>0){
		evalu[-x_f,2]<-evalu[-x_f,2]+1		
		}
		
		if(length(x_f)==0){
		evalu[,2]<-evalu[,2]+1
		}
		

		
		write.table(evalu,file=file_Selected[Cyl])
				
		list_fun<-function() {
		list(list(), list(),list(), list())
		}

Update_data<-replicate(No_Selected_Syn2, list_fun(), simplify=FALSE)


			for(upd in 1:dim(evalu)[1]){
			set_number<-evalu[upd,1]
			Step_code<-evalu[upd,2]
			Syn_number<-evalu[upd,3]
			Update_data[[upd]][[1]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[1]]
			Update_data[[upd]][[2]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[2]]
			Update_data[[upd]][[3]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[3]]
			Update_data[[upd]][[4]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[4]]
			}
		
	        Varit<-Update_data	
		
		for (g_index in 1:length(Varit)){
		names(Varit[[g_index]])<-c('data','qtl','mrk','sequ')
		}
		
	# Make selected to similar output of Syn2 or 3
	for (u in 1:length(Varit)){
	
	if(u %in% x_f){
	Varit[[u]]<-Varit[[u]]
		} else {
			
			localID<-Varit[[u]]$data[,1]

			# QTL
			x_qtl<-Varit[[u]][[2]]
			x_qtl<-cbind(localID,x_qtl)
			Varit[[u]][[2]]<-x_qtl

			# MRK
			x_mrk<-Varit[[u]][[3]]
			x_mrk<-cbind(localID,x_mrk)
			Varit[[u]][[3]]<-x_mrk

			# SEQ
			x_seq<-Varit[[u]][[4]]
			x_seq<-cbind(localID,x_seq)
			Varit[[u]][[4]]<-x_seq
		}
	}

#F1
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)


# F2
given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[Cyl-3]][[Cyl-2]]

# F2_SP
given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)


# Syn1
F2<-Sets[[Cyl-6]][[Cyl-3]]
F2_SP<-Sets[[Cyl-6]][[Cyl-2]]

all_Syn1_out<-Global_Syn1(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
Syn1_out<-all_Syn1_out[[1]]
S_Cross<-all_Syn1_out[[2]]

# Syn2
mini_plots<-Sets[[Cyl-7]][[Cyl-1]]
Syn2_out<-Global_Syn2(mini_plots,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Syn2_W<-rnorm(5,1,10)
Syn2_Green<-Sets[[Cyl-9]][[Cyl-2]]

mult<-rnorm(5,1,10)

# Syn3
Syn_2<-Sets[[Cyl-11]][[Cyl-2]]
Syn3_out<-Global_Syn3(Syn_2,size_Syn3,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Sets[[Cyl]][[Cyl]]<-F1_out
Sets[[Cyl-1]][[Cyl]]<-F2_Yellow
Sets[[Cyl-2]][[Cyl]]<-F2_W
Sets[[Cyl-3]][[Cyl]]<-F2_Green
Sets[[Cyl-4]][[Cyl]]<-F2_SP_out
Sets[[Cyl-5]][[Cyl]]<-S_Cross
Sets[[Cyl-6]][[Cyl]]<-Syn1_out
Sets[[Cyl-6]][[Cyl-1]]<-S_Cross
Sets[[Cyl-7]][[Cyl]]<-Syn2_out
Sets[[Cyl-8]][[Cyl]]<-Syn2_W
Sets[[Cyl-9]][[Cyl]]<-Syn2_Green
Sets[[Cyl-10]][[Cyl]]<-mult
Sets[[Cyl-11]][[Cyl]]<-Syn3_out

cat('All steps of year:',Cyl,'Finished',fill=TRUE)
}  

if(Cyl==31){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}

# For F2_SP
which_sets<-c(1:(Cyl-5))
which_step<-c(4:(Cyl-2))
which_sets
which_step

# For Syn2
which_sets2<-c(1:(Cyl-11))
which_step2<-c(10:(Cyl-2))
which_sets2
which_step2

which_sets<-c(which_sets2,which_sets)
which_step<-c(which_step2,which_step)
which_sets
which_step
		
     		evalu_lst<-list()
		
			for(K_1 in 1:length(which_sets)){
			dakheli1<-Sets[[which_sets[K_1]]][[which_step[K_1]]]
				evalu_mat<-matrix(ncol=4,nrow=length(dakheli1))
				for (K_2 in 1:length(dakheli1)){
				# tr1
				mean_tr1<-mean(dakheli1[[K_2]]$data$phen_1)
				# tr3
				mean_tr3<-mean(dakheli1[[K_2]]$data$phen_3)
				# tr4
				mean_tr4<-mean(dakheli1[[K_2]]$data$phen_4)
				SI<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
				evalu_mat[K_2,1]<-which_sets[K_1]
				evalu_mat[K_2,2]<-which_step[K_1]
				evalu_mat[K_2,3]<-K_2
				evalu_mat[K_2,4]<-SI
				}
				evalu_lst[[K_1]]<-evalu_mat
			}
			
			
		evalu<-do.call("rbind", evalu_lst)
		colnames(evalu)<-c("Set No.","Step","Family No.",'Index')
		evalu
		
		evalu<-evalu[order(-evalu[,4]),]
		evalu	
		evalu<-evalu[1:No_Selected_Syn2,]
		evalu

		jedi<-evalu[,2]-evalu[,1]
		x_f<-which(jedi==9)
		
		if(length(x_f)>0){
		evalu[-x_f,2]<-evalu[-x_f,2]+1		
		}
		
		if(length(x_f)==0){
		evalu[,2]<-evalu[,2]+1
		}
		

		
		write.table(evalu,file=file_Selected[Cyl])
				
		list_fun<-function() {
		list(list(), list(),list(), list())
		}

Update_data<-replicate(No_Selected_Syn2, list_fun(), simplify=FALSE)


			for(upd in 1:dim(evalu)[1]){
			set_number<-evalu[upd,1]
			Step_code<-evalu[upd,2]
			Syn_number<-evalu[upd,3]
			Update_data[[upd]][[1]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[1]]
			Update_data[[upd]][[2]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[2]]
			Update_data[[upd]][[3]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[3]]
			Update_data[[upd]][[4]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[4]]
			}
		
	        Varit<-Update_data	
		
		for (g_index in 1:length(Varit)){
		names(Varit[[g_index]])<-c('data','qtl','mrk','sequ')
		}
		
	# Make selected to similar output of Syn2 or 3
	for (u in 1:length(Varit)){
	
	if(u %in% x_f){
	Varit[[u]]<-Varit[[u]]
		} else {
			
			localID<-Varit[[u]]$data[,1]

			# QTL
			x_qtl<-Varit[[u]][[2]]
			x_qtl<-cbind(localID,x_qtl)
			Varit[[u]][[2]]<-x_qtl

			# MRK
			x_mrk<-Varit[[u]][[3]]
			x_mrk<-cbind(localID,x_mrk)
			Varit[[u]][[3]]<-x_mrk

			# SEQ
			x_seq<-Varit[[u]][[4]]
			x_seq<-cbind(localID,x_seq)
			Varit[[u]][[4]]<-x_seq
		}
	}

#F1
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)


# F2
given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[Cyl-3]][[Cyl-2]]

# F2_SP
given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)


# Syn1
F2<-Sets[[Cyl-6]][[Cyl-3]]
F2_SP<-Sets[[Cyl-6]][[Cyl-2]]

all_Syn1_out<-Global_Syn1(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
Syn1_out<-all_Syn1_out[[1]]
S_Cross<-all_Syn1_out[[2]]

# Syn2
mini_plots<-Sets[[Cyl-7]][[Cyl-1]]
Syn2_out<-Global_Syn2(mini_plots,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Syn2_W<-rnorm(5,1,10)
Syn2_Green<-Sets[[Cyl-9]][[Cyl-2]]

mult<-rnorm(5,1,10)

# Syn3
Syn_2<-Sets[[Cyl-11]][[Cyl-2]]
Syn3_out<-Global_Syn3(Syn_2,size_Syn3,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Sets[[Cyl]][[Cyl]]<-F1_out
Sets[[Cyl-1]][[Cyl]]<-F2_Yellow
Sets[[Cyl-2]][[Cyl]]<-F2_W
Sets[[Cyl-3]][[Cyl]]<-F2_Green
Sets[[Cyl-4]][[Cyl]]<-F2_SP_out
Sets[[Cyl-5]][[Cyl]]<-S_Cross
Sets[[Cyl-6]][[Cyl]]<-Syn1_out
Sets[[Cyl-6]][[Cyl-1]]<-S_Cross
Sets[[Cyl-7]][[Cyl]]<-Syn2_out
Sets[[Cyl-8]][[Cyl]]<-Syn2_W
Sets[[Cyl-9]][[Cyl]]<-Syn2_Green
Sets[[Cyl-10]][[Cyl]]<-mult
Sets[[Cyl-11]][[Cyl]]<-Syn3_out

cat('All steps of year:',Cyl,'Finished',fill=TRUE)
}  

if(Cyl==32){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}

# For F2_SP
which_sets<-c(1:(Cyl-5))
which_step<-c(4:(Cyl-2))
which_sets
which_step

# For Syn2
which_sets2<-c(1:(Cyl-11))
which_step2<-c(10:(Cyl-2))
which_sets2
which_step2

which_sets<-c(which_sets2,which_sets)
which_step<-c(which_step2,which_step)
which_sets
which_step
		
     		evalu_lst<-list()
		
			for(K_1 in 1:length(which_sets)){
			dakheli1<-Sets[[which_sets[K_1]]][[which_step[K_1]]]
				evalu_mat<-matrix(ncol=4,nrow=length(dakheli1))
				for (K_2 in 1:length(dakheli1)){
				# tr1
				mean_tr1<-mean(dakheli1[[K_2]]$data$phen_1)
				# tr3
				mean_tr3<-mean(dakheli1[[K_2]]$data$phen_3)
				# tr4
				mean_tr4<-mean(dakheli1[[K_2]]$data$phen_4)
				SI<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
				evalu_mat[K_2,1]<-which_sets[K_1]
				evalu_mat[K_2,2]<-which_step[K_1]
				evalu_mat[K_2,3]<-K_2
				evalu_mat[K_2,4]<-SI
				}
				evalu_lst[[K_1]]<-evalu_mat
			}
			
			
		evalu<-do.call("rbind", evalu_lst)
		colnames(evalu)<-c("Set No.","Step","Family No.",'Index')
		evalu
		
		evalu<-evalu[order(-evalu[,4]),]
		evalu	
		evalu<-evalu[1:No_Selected_Syn2,]
		evalu

		jedi<-evalu[,2]-evalu[,1]
		x_f<-which(jedi==9)
		
		if(length(x_f)>0){
		evalu[-x_f,2]<-evalu[-x_f,2]+1		
		}
		
		if(length(x_f)==0){
		evalu[,2]<-evalu[,2]+1
		}
		

		
		write.table(evalu,file=file_Selected[Cyl])
				
		list_fun<-function() {
		list(list(), list(),list(), list())
		}

Update_data<-replicate(No_Selected_Syn2, list_fun(), simplify=FALSE)


			for(upd in 1:dim(evalu)[1]){
			set_number<-evalu[upd,1]
			Step_code<-evalu[upd,2]
			Syn_number<-evalu[upd,3]
			Update_data[[upd]][[1]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[1]]
			Update_data[[upd]][[2]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[2]]
			Update_data[[upd]][[3]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[3]]
			Update_data[[upd]][[4]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[4]]
			}
		
	        Varit<-Update_data	
		
		for (g_index in 1:length(Varit)){
		names(Varit[[g_index]])<-c('data','qtl','mrk','sequ')
		}
		
	# Make selected to similar output of Syn2 or 3
	for (u in 1:length(Varit)){
	
	if(u %in% x_f){
	Varit[[u]]<-Varit[[u]]
		} else {
			
			localID<-Varit[[u]]$data[,1]

			# QTL
			x_qtl<-Varit[[u]][[2]]
			x_qtl<-cbind(localID,x_qtl)
			Varit[[u]][[2]]<-x_qtl

			# MRK
			x_mrk<-Varit[[u]][[3]]
			x_mrk<-cbind(localID,x_mrk)
			Varit[[u]][[3]]<-x_mrk

			# SEQ
			x_seq<-Varit[[u]][[4]]
			x_seq<-cbind(localID,x_seq)
			Varit[[u]][[4]]<-x_seq
		}
	}

#F1
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)


# F2
given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[Cyl-3]][[Cyl-2]]

# F2_SP
given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)


# Syn1
F2<-Sets[[Cyl-6]][[Cyl-3]]
F2_SP<-Sets[[Cyl-6]][[Cyl-2]]

all_Syn1_out<-Global_Syn1(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
Syn1_out<-all_Syn1_out[[1]]
S_Cross<-all_Syn1_out[[2]]

# Syn2
mini_plots<-Sets[[Cyl-7]][[Cyl-1]]
Syn2_out<-Global_Syn2(mini_plots,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Syn2_W<-rnorm(5,1,10)
Syn2_Green<-Sets[[Cyl-9]][[Cyl-2]]

mult<-rnorm(5,1,10)

# Syn3
Syn_2<-Sets[[Cyl-11]][[Cyl-2]]
Syn3_out<-Global_Syn3(Syn_2,size_Syn3,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Sets[[Cyl]][[Cyl]]<-F1_out
Sets[[Cyl-1]][[Cyl]]<-F2_Yellow
Sets[[Cyl-2]][[Cyl]]<-F2_W
Sets[[Cyl-3]][[Cyl]]<-F2_Green
Sets[[Cyl-4]][[Cyl]]<-F2_SP_out
Sets[[Cyl-5]][[Cyl]]<-S_Cross
Sets[[Cyl-6]][[Cyl]]<-Syn1_out
Sets[[Cyl-6]][[Cyl-1]]<-S_Cross
Sets[[Cyl-7]][[Cyl]]<-Syn2_out
Sets[[Cyl-8]][[Cyl]]<-Syn2_W
Sets[[Cyl-9]][[Cyl]]<-Syn2_Green
Sets[[Cyl-10]][[Cyl]]<-mult
Sets[[Cyl-11]][[Cyl]]<-Syn3_out

cat('All steps of year:',Cyl,'Finished',fill=TRUE)
}  

if(Cyl==33){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}

# For F2_SP
which_sets<-c(1:(Cyl-5))
which_step<-c(4:(Cyl-2))
which_sets
which_step

# For Syn2
which_sets2<-c(1:(Cyl-11))
which_step2<-c(10:(Cyl-2))
which_sets2
which_step2

which_sets<-c(which_sets2,which_sets)
which_step<-c(which_step2,which_step)
which_sets
which_step
		
     		evalu_lst<-list()
		
			for(K_1 in 1:length(which_sets)){
			dakheli1<-Sets[[which_sets[K_1]]][[which_step[K_1]]]
				evalu_mat<-matrix(ncol=4,nrow=length(dakheli1))
				for (K_2 in 1:length(dakheli1)){
				# tr1
				mean_tr1<-mean(dakheli1[[K_2]]$data$phen_1)
				# tr3
				mean_tr3<-mean(dakheli1[[K_2]]$data$phen_3)
				# tr4
				mean_tr4<-mean(dakheli1[[K_2]]$data$phen_4)
				SI<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
				evalu_mat[K_2,1]<-which_sets[K_1]
				evalu_mat[K_2,2]<-which_step[K_1]
				evalu_mat[K_2,3]<-K_2
				evalu_mat[K_2,4]<-SI
				}
				evalu_lst[[K_1]]<-evalu_mat
			}
			
			
		evalu<-do.call("rbind", evalu_lst)
		colnames(evalu)<-c("Set No.","Step","Family No.",'Index')
		evalu
		
		evalu<-evalu[order(-evalu[,4]),]
		evalu	
		evalu<-evalu[1:No_Selected_Syn2,]
		evalu

		jedi<-evalu[,2]-evalu[,1]
		x_f<-which(jedi==9)
		
		if(length(x_f)>0){
		evalu[-x_f,2]<-evalu[-x_f,2]+1		
		}
		
		if(length(x_f)==0){
		evalu[,2]<-evalu[,2]+1
		}
		

		
		write.table(evalu,file=file_Selected[Cyl])
				
		list_fun<-function() {
		list(list(), list(),list(), list())
		}

Update_data<-replicate(No_Selected_Syn2, list_fun(), simplify=FALSE)


			for(upd in 1:dim(evalu)[1]){
			set_number<-evalu[upd,1]
			Step_code<-evalu[upd,2]
			Syn_number<-evalu[upd,3]
			Update_data[[upd]][[1]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[1]]
			Update_data[[upd]][[2]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[2]]
			Update_data[[upd]][[3]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[3]]
			Update_data[[upd]][[4]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[4]]
			}
		
	        Varit<-Update_data	
		
		for (g_index in 1:length(Varit)){
		names(Varit[[g_index]])<-c('data','qtl','mrk','sequ')
		}
		
	# Make selected to similar output of Syn2 or 3
	for (u in 1:length(Varit)){
	
	if(u %in% x_f){
	Varit[[u]]<-Varit[[u]]
		} else {
			
			localID<-Varit[[u]]$data[,1]

			# QTL
			x_qtl<-Varit[[u]][[2]]
			x_qtl<-cbind(localID,x_qtl)
			Varit[[u]][[2]]<-x_qtl

			# MRK
			x_mrk<-Varit[[u]][[3]]
			x_mrk<-cbind(localID,x_mrk)
			Varit[[u]][[3]]<-x_mrk

			# SEQ
			x_seq<-Varit[[u]][[4]]
			x_seq<-cbind(localID,x_seq)
			Varit[[u]][[4]]<-x_seq
		}
	}

#F1
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)


# F2
given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[Cyl-3]][[Cyl-2]]

# F2_SP
given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)


# Syn1
F2<-Sets[[Cyl-6]][[Cyl-3]]
F2_SP<-Sets[[Cyl-6]][[Cyl-2]]

all_Syn1_out<-Global_Syn1(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
Syn1_out<-all_Syn1_out[[1]]
S_Cross<-all_Syn1_out[[2]]

# Syn2
mini_plots<-Sets[[Cyl-7]][[Cyl-1]]
Syn2_out<-Global_Syn2(mini_plots,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Syn2_W<-rnorm(5,1,10)
Syn2_Green<-Sets[[Cyl-9]][[Cyl-2]]

mult<-rnorm(5,1,10)

# Syn3
Syn_2<-Sets[[Cyl-11]][[Cyl-2]]
Syn3_out<-Global_Syn3(Syn_2,size_Syn3,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Sets[[Cyl]][[Cyl]]<-F1_out
Sets[[Cyl-1]][[Cyl]]<-F2_Yellow
Sets[[Cyl-2]][[Cyl]]<-F2_W
Sets[[Cyl-3]][[Cyl]]<-F2_Green
Sets[[Cyl-4]][[Cyl]]<-F2_SP_out
Sets[[Cyl-5]][[Cyl]]<-S_Cross
Sets[[Cyl-6]][[Cyl]]<-Syn1_out
Sets[[Cyl-6]][[Cyl-1]]<-S_Cross
Sets[[Cyl-7]][[Cyl]]<-Syn2_out
Sets[[Cyl-8]][[Cyl]]<-Syn2_W
Sets[[Cyl-9]][[Cyl]]<-Syn2_Green
Sets[[Cyl-10]][[Cyl]]<-mult
Sets[[Cyl-11]][[Cyl]]<-Syn3_out

cat('All steps of year:',Cyl,'Finished',fill=TRUE)
}  

if(Cyl==34){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}

# For F2_SP
which_sets<-c(1:(Cyl-5))
which_step<-c(4:(Cyl-2))
which_sets
which_step

# For Syn2
which_sets2<-c(1:(Cyl-11))
which_step2<-c(10:(Cyl-2))
which_sets2
which_step2

which_sets<-c(which_sets2,which_sets)
which_step<-c(which_step2,which_step)
which_sets
which_step
		
     		evalu_lst<-list()
		
			for(K_1 in 1:length(which_sets)){
			dakheli1<-Sets[[which_sets[K_1]]][[which_step[K_1]]]
				evalu_mat<-matrix(ncol=4,nrow=length(dakheli1))
				for (K_2 in 1:length(dakheli1)){
				# tr1
				mean_tr1<-mean(dakheli1[[K_2]]$data$phen_1)
				# tr3
				mean_tr3<-mean(dakheli1[[K_2]]$data$phen_3)
				# tr4
				mean_tr4<-mean(dakheli1[[K_2]]$data$phen_4)
				SI<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
				evalu_mat[K_2,1]<-which_sets[K_1]
				evalu_mat[K_2,2]<-which_step[K_1]
				evalu_mat[K_2,3]<-K_2
				evalu_mat[K_2,4]<-SI
				}
				evalu_lst[[K_1]]<-evalu_mat
			}
			
			
		evalu<-do.call("rbind", evalu_lst)
		colnames(evalu)<-c("Set No.","Step","Family No.",'Index')
		evalu
		
		evalu<-evalu[order(-evalu[,4]),]
		evalu	
		evalu<-evalu[1:No_Selected_Syn2,]
		evalu

		jedi<-evalu[,2]-evalu[,1]
		x_f<-which(jedi==9)
		
		if(length(x_f)>0){
		evalu[-x_f,2]<-evalu[-x_f,2]+1		
		}
		
		if(length(x_f)==0){
		evalu[,2]<-evalu[,2]+1
		}
		

		
		write.table(evalu,file=file_Selected[Cyl])
				
		list_fun<-function() {
		list(list(), list(),list(), list())
		}

Update_data<-replicate(No_Selected_Syn2, list_fun(), simplify=FALSE)


			for(upd in 1:dim(evalu)[1]){
			set_number<-evalu[upd,1]
			Step_code<-evalu[upd,2]
			Syn_number<-evalu[upd,3]
			Update_data[[upd]][[1]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[1]]
			Update_data[[upd]][[2]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[2]]
			Update_data[[upd]][[3]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[3]]
			Update_data[[upd]][[4]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[4]]
			}
		
	        Varit<-Update_data	
		
		for (g_index in 1:length(Varit)){
		names(Varit[[g_index]])<-c('data','qtl','mrk','sequ')
		}
		
	# Make selected to similar output of Syn2 or 3
	for (u in 1:length(Varit)){
	
	if(u %in% x_f){
	Varit[[u]]<-Varit[[u]]
		} else {
			
			localID<-Varit[[u]]$data[,1]

			# QTL
			x_qtl<-Varit[[u]][[2]]
			x_qtl<-cbind(localID,x_qtl)
			Varit[[u]][[2]]<-x_qtl

			# MRK
			x_mrk<-Varit[[u]][[3]]
			x_mrk<-cbind(localID,x_mrk)
			Varit[[u]][[3]]<-x_mrk

			# SEQ
			x_seq<-Varit[[u]][[4]]
			x_seq<-cbind(localID,x_seq)
			Varit[[u]][[4]]<-x_seq
		}
	}

#F1
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)


# F2
given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[Cyl-3]][[Cyl-2]]

# F2_SP
given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)


# Syn1
F2<-Sets[[Cyl-6]][[Cyl-3]]
F2_SP<-Sets[[Cyl-6]][[Cyl-2]]

all_Syn1_out<-Global_Syn1(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
Syn1_out<-all_Syn1_out[[1]]
S_Cross<-all_Syn1_out[[2]]

# Syn2
mini_plots<-Sets[[Cyl-7]][[Cyl-1]]
Syn2_out<-Global_Syn2(mini_plots,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Syn2_W<-rnorm(5,1,10)
Syn2_Green<-Sets[[Cyl-9]][[Cyl-2]]

mult<-rnorm(5,1,10)

# Syn3
Syn_2<-Sets[[Cyl-11]][[Cyl-2]]
Syn3_out<-Global_Syn3(Syn_2,size_Syn3,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Sets[[Cyl]][[Cyl]]<-F1_out
Sets[[Cyl-1]][[Cyl]]<-F2_Yellow
Sets[[Cyl-2]][[Cyl]]<-F2_W
Sets[[Cyl-3]][[Cyl]]<-F2_Green
Sets[[Cyl-4]][[Cyl]]<-F2_SP_out
Sets[[Cyl-5]][[Cyl]]<-S_Cross
Sets[[Cyl-6]][[Cyl]]<-Syn1_out
Sets[[Cyl-6]][[Cyl-1]]<-S_Cross
Sets[[Cyl-7]][[Cyl]]<-Syn2_out
Sets[[Cyl-8]][[Cyl]]<-Syn2_W
Sets[[Cyl-9]][[Cyl]]<-Syn2_Green
Sets[[Cyl-10]][[Cyl]]<-mult
Sets[[Cyl-11]][[Cyl]]<-Syn3_out

cat('All steps of year:',Cyl,'Finished',fill=TRUE)
}  

if(Cyl==35){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}

# For F2_SP
which_sets<-c(1:(Cyl-5))
which_step<-c(4:(Cyl-2))
which_sets
which_step

# For Syn2
which_sets2<-c(1:(Cyl-11))
which_step2<-c(10:(Cyl-2))
which_sets2
which_step2

which_sets<-c(which_sets2,which_sets)
which_step<-c(which_step2,which_step)
which_sets
which_step
		
     		evalu_lst<-list()
		
			for(K_1 in 1:length(which_sets)){
			dakheli1<-Sets[[which_sets[K_1]]][[which_step[K_1]]]
				evalu_mat<-matrix(ncol=4,nrow=length(dakheli1))
				for (K_2 in 1:length(dakheli1)){
				# tr1
				mean_tr1<-mean(dakheli1[[K_2]]$data$phen_1)
				# tr3
				mean_tr3<-mean(dakheli1[[K_2]]$data$phen_3)
				# tr4
				mean_tr4<-mean(dakheli1[[K_2]]$data$phen_4)
				SI<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
				evalu_mat[K_2,1]<-which_sets[K_1]
				evalu_mat[K_2,2]<-which_step[K_1]
				evalu_mat[K_2,3]<-K_2
				evalu_mat[K_2,4]<-SI
				}
				evalu_lst[[K_1]]<-evalu_mat
			}
			
			
		evalu<-do.call("rbind", evalu_lst)
		colnames(evalu)<-c("Set No.","Step","Family No.",'Index')
		evalu
		
		evalu<-evalu[order(-evalu[,4]),]
		evalu	
		evalu<-evalu[1:No_Selected_Syn2,]
		evalu

		jedi<-evalu[,2]-evalu[,1]
		x_f<-which(jedi==9)
		
		if(length(x_f)>0){
		evalu[-x_f,2]<-evalu[-x_f,2]+1		
		}
		
		if(length(x_f)==0){
		evalu[,2]<-evalu[,2]+1
		}
		

		
		write.table(evalu,file=file_Selected[Cyl])
				
		list_fun<-function() {
		list(list(), list(),list(), list())
		}

Update_data<-replicate(No_Selected_Syn2, list_fun(), simplify=FALSE)


			for(upd in 1:dim(evalu)[1]){
			set_number<-evalu[upd,1]
			Step_code<-evalu[upd,2]
			Syn_number<-evalu[upd,3]
			Update_data[[upd]][[1]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[1]]
			Update_data[[upd]][[2]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[2]]
			Update_data[[upd]][[3]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[3]]
			Update_data[[upd]][[4]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[4]]
			}
		
	        Varit<-Update_data	
		
		for (g_index in 1:length(Varit)){
		names(Varit[[g_index]])<-c('data','qtl','mrk','sequ')
		}
		
	# Make selected to similar output of Syn2 or 3
	for (u in 1:length(Varit)){
	
	if(u %in% x_f){
	Varit[[u]]<-Varit[[u]]
		} else {
			
			localID<-Varit[[u]]$data[,1]

			# QTL
			x_qtl<-Varit[[u]][[2]]
			x_qtl<-cbind(localID,x_qtl)
			Varit[[u]][[2]]<-x_qtl

			# MRK
			x_mrk<-Varit[[u]][[3]]
			x_mrk<-cbind(localID,x_mrk)
			Varit[[u]][[3]]<-x_mrk

			# SEQ
			x_seq<-Varit[[u]][[4]]
			x_seq<-cbind(localID,x_seq)
			Varit[[u]][[4]]<-x_seq
		}
	}

#F1
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)


# F2
given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[Cyl-3]][[Cyl-2]]

# F2_SP
given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)


# Syn1
F2<-Sets[[Cyl-6]][[Cyl-3]]
F2_SP<-Sets[[Cyl-6]][[Cyl-2]]

all_Syn1_out<-Global_Syn1(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
Syn1_out<-all_Syn1_out[[1]]
S_Cross<-all_Syn1_out[[2]]

# Syn2
mini_plots<-Sets[[Cyl-7]][[Cyl-1]]
Syn2_out<-Global_Syn2(mini_plots,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Syn2_W<-rnorm(5,1,10)
Syn2_Green<-Sets[[Cyl-9]][[Cyl-2]]

mult<-rnorm(5,1,10)

# Syn3
Syn_2<-Sets[[Cyl-11]][[Cyl-2]]
Syn3_out<-Global_Syn3(Syn_2,size_Syn3,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Sets[[Cyl]][[Cyl]]<-F1_out
Sets[[Cyl-1]][[Cyl]]<-F2_Yellow
Sets[[Cyl-2]][[Cyl]]<-F2_W
Sets[[Cyl-3]][[Cyl]]<-F2_Green
Sets[[Cyl-4]][[Cyl]]<-F2_SP_out
Sets[[Cyl-5]][[Cyl]]<-S_Cross
Sets[[Cyl-6]][[Cyl]]<-Syn1_out
Sets[[Cyl-6]][[Cyl-1]]<-S_Cross
Sets[[Cyl-7]][[Cyl]]<-Syn2_out
Sets[[Cyl-8]][[Cyl]]<-Syn2_W
Sets[[Cyl-9]][[Cyl]]<-Syn2_Green
Sets[[Cyl-10]][[Cyl]]<-mult
Sets[[Cyl-11]][[Cyl]]<-Syn3_out

cat('All steps of year:',Cyl,'Finished',fill=TRUE)
}  

if(Cyl==36){
   
params<-plan_XL[,Cyl]
inde_no<-which(params==0)
for (Kop in 1:length(inde_no)){
 Sets[[inde_no[Kop]]][[Cyl]]<-rnorm(10,1,20)
}

# For F2_SP
which_sets<-c(1:(Cyl-5))
which_step<-c(4:(Cyl-2))
which_sets
which_step

# For Syn2
which_sets2<-c(1:(Cyl-11))
which_step2<-c(10:(Cyl-2))
which_sets2
which_step2

which_sets<-c(which_sets2,which_sets)
which_step<-c(which_step2,which_step)
which_sets
which_step
		
     		evalu_lst<-list()
		
			for(K_1 in 1:length(which_sets)){
			dakheli1<-Sets[[which_sets[K_1]]][[which_step[K_1]]]
				evalu_mat<-matrix(ncol=4,nrow=length(dakheli1))
				for (K_2 in 1:length(dakheli1)){
				# tr1
				mean_tr1<-mean(dakheli1[[K_2]]$data$phen_1)
				# tr3
				mean_tr3<-mean(dakheli1[[K_2]]$data$phen_3)
				# tr4
				mean_tr4<-mean(dakheli1[[K_2]]$data$phen_4)
				SI<-(b_tr1*mean_tr1)+(b_tr3*mean_tr3)+(b_tr4*mean_tr4)
				evalu_mat[K_2,1]<-which_sets[K_1]
				evalu_mat[K_2,2]<-which_step[K_1]
				evalu_mat[K_2,3]<-K_2
				evalu_mat[K_2,4]<-SI
				}
				evalu_lst[[K_1]]<-evalu_mat
			}
			
			
		evalu<-do.call("rbind", evalu_lst)
		colnames(evalu)<-c("Set No.","Step","Family No.",'Index')
		evalu
		
		evalu<-evalu[order(-evalu[,4]),]
		evalu	
		evalu<-evalu[1:No_Selected_Syn2,]
		evalu

		jedi<-evalu[,2]-evalu[,1]
		x_f<-which(jedi==9)
		
		if(length(x_f)>0){
		evalu[-x_f,2]<-evalu[-x_f,2]+1		
		}
		
		if(length(x_f)==0){
		evalu[,2]<-evalu[,2]+1
		}
		

		
		write.table(evalu,file=file_Selected[Cyl])
				
		list_fun<-function() {
		list(list(), list(),list(), list())
		}

Update_data<-replicate(No_Selected_Syn2, list_fun(), simplify=FALSE)


			for(upd in 1:dim(evalu)[1]){
			set_number<-evalu[upd,1]
			Step_code<-evalu[upd,2]
			Syn_number<-evalu[upd,3]
			Update_data[[upd]][[1]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[1]]
			Update_data[[upd]][[2]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[2]]
			Update_data[[upd]][[3]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[3]]
			Update_data[[upd]][[4]]<-Sets[[set_number]][[Step_code]][[Syn_number]][[4]]
			}
		
	        Varit<-Update_data	
		
		for (g_index in 1:length(Varit)){
		names(Varit[[g_index]])<-c('data','qtl','mrk','sequ')
		}
		
	# Make selected to similar output of Syn2 or 3
	for (u in 1:length(Varit)){
	
	if(u %in% x_f){
	Varit[[u]]<-Varit[[u]]
		} else {
			
			localID<-Varit[[u]]$data[,1]

			# QTL
			x_qtl<-Varit[[u]][[2]]
			x_qtl<-cbind(localID,x_qtl)
			Varit[[u]][[2]]<-x_qtl

			# MRK
			x_mrk<-Varit[[u]][[3]]
			x_mrk<-cbind(localID,x_mrk)
			Varit[[u]][[3]]<-x_mrk

			# SEQ
			x_seq<-Varit[[u]][[4]]
			x_seq<-cbind(localID,x_seq)
			Varit[[u]][[4]]<-x_seq
		}
	}

#F1
F1_out<-F1_Global(Varit,No_F1_Family,nvar,noff_m_v,
link_map,linkage_map_qtl,linkage_map_mrk,qtl_effects
,chrom_length,var_add,phen_var)


# F2
given_L<-length(Sets[[Cyl-1]][[Cyl-1]])
F1<-Sets[[Cyl-1]][[Cyl-1]]
F2_Yellow<-F2_Global(given_L,F1,size_F2 ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

F2_W<-rnorm(5,1,10)
F2_Green<-Sets[[Cyl-3]][[Cyl-2]]

# F2_SP
given_L<-length(Sets[[Cyl-4]][[Cyl-4]])
F1_SP<-Sets[[Cyl-4]][[Cyl-4]]
F2_SP_out<-F2_SP_Global(given_L,F1_SP,size_F2_SP ,qtl_effects,
noff_F2,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

S_Cross<-rnorm(5,1,10)


# Syn1
F2<-Sets[[Cyl-6]][[Cyl-3]]
F2_SP<-Sets[[Cyl-6]][[Cyl-2]]

all_Syn1_out<-Global_Syn1(F2,F2_SP,No_F2_Family_selected,nClonalrow,size_polycross,
qtl_effects,noff_per_each_syn,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)
Syn1_out<-all_Syn1_out[[1]]
S_Cross<-all_Syn1_out[[2]]


# Syn2
mini_plots<-Sets[[Cyl-7]][[Cyl-1]]
Syn2_out<-Global_Syn2(mini_plots,size_Syn1,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Syn2_W<-rnorm(5,1,10)
Syn2_Green<-Sets[[Cyl-9]][[Cyl-2]]

mult<-rnorm(5,1,10)

# Syn3
Syn_2<-Sets[[Cyl-11]][[Cyl-2]]
Syn3_out<-Global_Syn3(Syn_2,size_Syn3,noff_per_each_mating_in_syn,
qtl_effects,link_map,linkage_map_qtl,linkage_map_mrk
,chrom_length,var_add,phen_var)

Sets[[Cyl]][[Cyl]]<-F1_out
Sets[[Cyl-1]][[Cyl]]<-F2_Yellow
Sets[[Cyl-2]][[Cyl]]<-F2_W
Sets[[Cyl-3]][[Cyl]]<-F2_Green
Sets[[Cyl-4]][[Cyl]]<-F2_SP_out
Sets[[Cyl-5]][[Cyl]]<-S_Cross
Sets[[Cyl-6]][[Cyl]]<-Syn1_out
Sets[[Cyl-6]][[Cyl-1]]<-S_Cross
Sets[[Cyl-7]][[Cyl]]<-Syn2_out
Sets[[Cyl-8]][[Cyl]]<-Syn2_W
Sets[[Cyl-9]][[Cyl]]<-Syn2_Green
Sets[[Cyl-10]][[Cyl]]<-mult
Sets[[Cyl-11]][[Cyl]]<-Syn3_out

cat('All steps of year:',Cyl,'Finished',fill=TRUE)
}  


} #End of Cycles



# -------------------
# All plots F1
# -------------------
tr1<-c()
tr2<-c()
tr3<-c()
tr4<-c()

Year<-1
Kol_set<-length(Sets)-11
for (Set_Count in 1:Kol_set){

F1<-Sets[[Set_Count]][[Year]]
length(F1)

	for (j in 1:length(F1)){
	# TBV
	tr1[j]<-mean(F1[[j]][[1]][,10])
	tr2[j]<-mean(F1[[j]][[1]][,11])
	tr3[j]<-mean(F1[[j]][[1]][,12])
	tr4[j]<-mean(F1[[j]][[1]][,13])

	}
	
	# Mean TBV
	TBV_F1_AV[Set_Count,1]<-mean(tr1)
	TBV_F1_AV[Set_Count,2]<-mean(tr2)
	TBV_F1_AV[Set_Count,3]<-mean(tr3)
	TBV_F1_AV[Set_Count,4]<-mean(tr4)
	

	# Var TBV
	TBV_F1_VAR[Set_Count,1]<-var(tr1)
	TBV_F1_VAR[Set_Count,2]<-var(tr2)
	TBV_F1_VAR[Set_Count,3]<-var(tr3)
	TBV_F1_VAR[Set_Count,4]<-var(tr4)

Year<-Year+1
}


# -------------------
# All plots Syn_2
# -------------------

tr1<-c()
tr2<-c()
tr3<-c()
tr4<-c()

cor_tr1<-c()
cor_tr2<-c()
cor_tr3<-c()
cor_tr4<-c()

Year<-10
Kol_set<-length(Sets)-11
for (Set_Count in 1:Kol_set){

Syn_2<-Sets[[Set_Count]][[Year]]
length(Syn_2)

	for (j in 1:length(Syn_2)){
	# TBV
	tr1[j]<-mean(Syn_2[[j]][[1]][,10])
	tr2[j]<-mean(Syn_2[[j]][[1]][,11])
	tr3[j]<-mean(Syn_2[[j]][[1]][,12])
	tr4[j]<-mean(Syn_2[[j]][[1]][,13])
	
	#Cor	
	cor_tr1[j]<-cor(Syn_2[[j]][[1]][,6],Syn_2[[j]][[1]][,10])
	cor_tr2[j]<-cor(Syn_2[[j]][[1]][,7],Syn_2[[j]][[1]][,11])
	cor_tr3[j]<-cor(Syn_2[[j]][[1]][,8],Syn_2[[j]][[1]][,12])
	cor_tr4[j]<-cor(Syn_2[[j]][[1]][,9],Syn_2[[j]][[1]][,13])
	}
	
	# Mean TBV
	TBV_Syn2_AV[Set_Count,1]<-mean(tr1)
	TBV_Syn2_AV[Set_Count,2]<-mean(tr2)
	TBV_Syn2_AV[Set_Count,3]<-mean(tr3)
	TBV_Syn2_AV[Set_Count,4]<-mean(tr4)
	
	# Mean ACC all traits
	TBV_Syn2_ACC[Set_Count,1]<-mean(cor_tr1)
	TBV_Syn2_ACC[Set_Count,2]<-mean(cor_tr2)
	TBV_Syn2_ACC[Set_Count,3]<-mean(cor_tr3)
	TBV_Syn2_ACC[Set_Count,4]<-mean(cor_tr4)
	
	# Var TBV
	TBV_Syn2_VAR[Set_Count,1]<-var(tr1)
	TBV_Syn2_VAR[Set_Count,2]<-var(tr2)
	TBV_Syn2_VAR[Set_Count,3]<-var(tr3)
	TBV_Syn2_VAR[Set_Count,4]<-var(tr4)


Year<-Year+1
}


# -------------------
# All plots Syn_3
# -------------------

tr1<-c()
tr2<-c()
tr3<-c()
tr4<-c()

cor_tr1<-c()
cor_tr2<-c()
cor_tr3<-c()
cor_tr4<-c()

Year<-12
 for (Set_Count in 1:Kol_set){
Syn_3<-Sets[[Set_Count]][[Year]]
length(Syn_3)

	for (j in 1:length(Syn_3)){
	# TBV
	tr1[j]<-mean(Syn_3[[j]][[1]][,10])
	tr2[j]<-mean(Syn_3[[j]][[1]][,11])
	tr3[j]<-mean(Syn_3[[j]][[1]][,12])
	tr4[j]<-mean(Syn_3[[j]][[1]][,13])
	
	#Cor	
	cor_tr1[j]<-cor(Syn_3[[j]][[1]][,6],Syn_3[[j]][[1]][,10])
	cor_tr2[j]<-cor(Syn_3[[j]][[1]][,7],Syn_3[[j]][[1]][,11])
	cor_tr3[j]<-cor(Syn_3[[j]][[1]][,8],Syn_3[[j]][[1]][,12])
	cor_tr4[j]<-cor(Syn_3[[j]][[1]][,9],Syn_3[[j]][[1]][,13])
	}
	
	# Mean TBV
	TBV_Syn3_AV[Set_Count,1]<-mean(tr1)
	TBV_Syn3_AV[Set_Count,2]<-mean(tr2)
	TBV_Syn3_AV[Set_Count,3]<-mean(tr3)
	TBV_Syn3_AV[Set_Count,4]<-mean(tr4)
	
	# Mean ACC all traits
	TBV_Syn3_ACC[Set_Count,1]<-mean(cor_tr1)
	TBV_Syn3_ACC[Set_Count,2]<-mean(cor_tr2)
	TBV_Syn3_ACC[Set_Count,3]<-mean(cor_tr3)
	TBV_Syn3_ACC[Set_Count,4]<-mean(cor_tr4)
	
	# Var TBV
	TBV_Syn3_VAR[Set_Count,1]<-var(tr1)
	TBV_Syn3_VAR[Set_Count,2]<-var(tr2)
	TBV_Syn3_VAR[Set_Count,3]<-var(tr3)
	TBV_Syn3_VAR[Set_Count,4]<-var(tr4)


Year<-Year+1
}




# -------------------
# Genetic variance within each set
# -------------------
Var_within_set_Tr1<-matrix(0,ncol=7,nrow=Kol_set)
Var_within_set_Tr2<-matrix(0,ncol=7,nrow=Kol_set)
Var_within_set_Tr3<-matrix(0,ncol=7,nrow=Kol_set)
Var_within_set_Tr4<-matrix(0,ncol=7,nrow=Kol_set)



tr1<-c()
tr2<-c()
tr3<-c()
tr4<-c()


 y_names<-c(1,2,5,7,8,12)
 
#  1: F1
#  2: F2
#  5: F2SP
#  7: Syn1
#  8: Syn2
# 12: Syn3
	  
 
 
 
 for (Set_Count in 1:Kol_set){
 

 for (pq in 1:length(y_names)){
 
idata<-Sets[[Set_Count]][[y_names[pq]]]
 
 	for (j in 1:length(idata)){
	# TBV
	tr1[j]<-mean(idata[[j]][[1]][,10])
	tr2[j]<-mean(idata[[j]][[1]][,11])
	tr3[j]<-mean(idata[[j]][[1]][,12])
	tr4[j]<-mean(idata[[j]][[1]][,13])
	}

Var_within_set_Tr1[Set_Count,pq]<-var(tr1)
Var_within_set_Tr2[Set_Count,pq]<-var(tr2)
Var_within_set_Tr3[Set_Count,pq]<-var(tr3)
Var_within_set_Tr4[Set_Count,pq]<-var(tr4)
    }
  

y_names<-y_names+1
}



	#  6: F2SP Selected
	
y_names<-6
 for (Set_Count in 1:Kol_set){
 
idata<-Sets[[Set_Count]][[y_names]]
 
Var_within_set_Tr1[Set_Count,7]<-idata[1]
Var_within_set_Tr2[Set_Count,7]<-idata[2]
Var_within_set_Tr3[Set_Count,7]<-idata[3]
Var_within_set_Tr4[Set_Count,7]<-idata[4]

y_names<-y_names+1
}


# Reorder as the steps within each set

# Original
#1  1: F1
#2  2: F2
#3  5: F2SP
#4  7: Syn1
#5  8: Syn2
#6 12: Syn3
#7  6: F2SP Selected

# Ordered
#1  1: F1
#2  2: F2
#3  5: F2SP
#7  6: F2SP Selected
#4  7: Syn1
#5  8: Syn2
#6 12: Syn3


ord<-c(1,2,3,7,4,5,6)
Var_within_set_Tr1<-Var_within_set_Tr1[,ord]
Var_within_set_Tr2<-Var_within_set_Tr2[,ord]
Var_within_set_Tr3<-Var_within_set_Tr3[,ord]
Var_within_set_Tr4<-Var_within_set_Tr4[,ord]

	
	
# ------------------------------
# Write some results to output
# -------------------------------

	# Results of replicates
# F1
	# TBV
	TBV_AV_F1<-TBV_F1_AV
	write.table(TBV_AV_F1,file=file_TBV_F1[irep])
		
	# Var genetic
	TBV_VAR_F1<-TBV_F1_VAR
	write.table(TBV_VAR_F1,file=file_VAR_F1[irep])	 


# Syn 2
	# TBV
	TBV_AV<-TBV_Syn2_AV
	write.table(TBV_AV,file=file4_syn2[irep])

	# ACC Phen 
	TBV_ACC<-TBV_Syn2_ACC
	write.table(TBV_ACC,file=file5_syn2[irep])
		
	# Var genetic
	TBV_VAR<-TBV_Syn2_VAR
	write.table(TBV_VAR,file=file6_syn2[irep])
	
	
	
# Syn 3
	# TBV
	TBV_AV<-TBV_Syn3_AV
	write.table(TBV_AV,file=file4_syn3[irep])

	# ACC Phen 
	TBV_ACC<-TBV_Syn3_ACC
	write.table(TBV_ACC,file=file5_syn3[irep])
	
	# Var genetic
	TBV_VAR<-TBV_Syn3_VAR
	write.table(TBV_VAR,file=file6_syn3[irep])
	
# Var within set	
	# Tr1
	write.table(Var_within_set_Tr1,file=File_Var_within_Tr1[irep])
	# Tr2
	write.table(Var_within_set_Tr2,file=File_Var_within_Tr2[irep])
	# Tr3
	write.table(Var_within_set_Tr3,file=File_Var_within_Tr3[irep])
	# Tr4
	write.table(Var_within_set_Tr4,file=File_Var_within_Tr4[irep])
	

cat('-------------------------------------',fill=TRUE)
cat('-------------------------------------',fill=TRUE)
cat('Replicate number:',irep,'is Finished',fill=TRUE)
cat('-------------------------------------',fill=TRUE)
cat('-------------------------------------',fill=TRUE)
	
	
	} #End of replicates
	

