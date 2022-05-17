##Reads interaction files obtained through the following data bases: APID, DOROTHEA and OMNIPATH for Homo sp.

apid <- read.delim("9606_noISI_Q2.txt",header=T, stringsAsFactors = F)

omni <- read.delim("omnipathdb.txt",header=T, stringsAsFactors = F)

omni <- omni[omni$is_directed==1,]

dorothea <- read.delim("dorothea_AB.txt",header=T, stringsAsFactors = F)

##Extracts lists of uniprot identifiers of proteins that interact with at least one of NF-kB sub units

# NFKBI (P105 or P50) - P19838
# RELA (P65) - Q04206

source("Auxiliar_functions.R")

subunits      <- c("P19838","Q04206")
apid_cln_pair <- c("UniprotID_A", "UniprotID_B")

filtered_apid   <- fd_objects_frm_cln_union(apid, cln = apid_cln_pair, objs = subunits, vice_versa = TRUE)
interactors_apid <- union(filtered_apid$vice, filtered_apid$versa)

##Extracts the same objects, but know, those must be NF-kB subunits' regulators

omni_n_dorothea_cln_pair   <- c("source", "target")
filtered_omni     <- fd_objects_frm_cln_union(omni,     cln = omni_n_dorothea_cln_pair, objs = subunits, vice_versa = TRUE)
filtered_dorothea <- fd_objects_frm_cln_union(dorothea, cln = omni_n_dorothea_cln_pair, objs = subunits, vice_versa = TRUE)

rm(subunits,omni_n_dorothea_cln_pair, apid_cln_pair)
##List with the union of targets that have as source either NFkB1 or RelA 

targets    <- union(filtered_omni$vice, filtered_dorothea$vice)

##Extraction of list of Uniprot indentifiers of proteins that regulate at least one of NF-kB's subunit
#So its a list of sources that have the subunits as targets

regulators <- union(filtered_omni$versa, filtered_dorothea$versa)

rm(filtered_apid, filtered_dorothea, filtered_omni, omni, dorothea, apid)
##Intersection between lists
  
targets_n_regulators     <- intersect(targets,   regulators)
targets_n_interactors    <- intersect(targets,   interactors_apid)
regulators_n_interactors <- intersect(regulators,interactors_apid) 

##Evaluates if the intersection lists are bigger then expected through FDR (false discovery rate) determinantion

#construction of variable universe
universe_apid     <- union(apid$UniprotID_A,apid$UniprotID_B)
universe_omni     <- union(omni$source,omni$target)
universe_dorothea <- union(dorothea$source,dorothea$target)
universe          <- union(universe_apid,union(universe_omni,universe_dorothea))

#sizes
t <- length(universe)
a <- length(targets)
r <- length(regulators)
i <- length(interactors_apid)
a_i <- length(targets_n_interactors)
a_r <- length(targets_n_regulators)
r_i <- length(regulators_n_interactors)

##p_values>>
#lower tail = F, so P[X > x]
#q = size - 1, 
p_a_i <- phyper(a_i-1,a,t-a,i,lower.tail=F)
p_a_r <- phyper(a_r-1,a,t-a,r,lower.tail=F)
p_r_i <- phyper(r_i-1,r,t-r,i,lower.tail=F) 

#It was found a large number of proteins that interact with NF-kB, larger than expected - there may be various retroaction cycles
#we obtained very low p_values, so our data is not random 

rm(a,r,i,t,a_i,a_r,r_i,p_a_i, p_a_r, p_r_i, universe_apid, universe_omni, universe_dorothea, universe)

#creation of a list that includes every candidate for NF-kB's retroaction cycles 
cyclesretro <- union(targets_n_interactors, targets_n_regulators) 

##Creation and organization of acquired RNAseq data of HeLa cancer cell line and various human tissues
library(org.Hs.eg.db)

#from: RNA-seq data of 675 commonly used human cancer cell lines, in Atlas expression
helaprot <- read.delim("HeLa_TPM_results.tsv",header=T, skip = 4 ,stringsAsFactors = F)
colnames(helaprot)[3] <- "HeLa"
helaprot <- helaprot[!is.na(helaprot$HeLa),]

#Usage of ENSEMBL identifiers, present in both experimental dataframes, instead of UNIPROT
cyclesretro_ensemble <- mapIds(org.Hs.eg.db, keys=cyclesretro,column="ENSEMBL", keytype = "UNIPROT")
rm(cyclesretro)
#from: RNA-seq data of 53 human tissues from paper GTEx, in Atlas expression
tissue_exp            <- read.delim("E-MTAB-5214-query-results.tsv",header=T, skip = 4 ,stringsAsFactors = FALSE)
intersection          <- intersect(cyclesretro_ensemble, intersect(helaprot$Gene.ID, tissue_exp$Gene.ID))

#Filtration of experimental data with intersection list
hela_cr        <- helaprot[index_gene_help(helaprot),]
tissue_cr      <- tissue_exp[index_gene_help(tissue_exp),]

rm(tissue_exp, helaprot)
tissue_filter  <- c("Gene.ID", "Gene.Name","liver", "stomach", "pancreas", "urinary.bladder",
                    "prostate.gland", "lung", "adrenal.gland" , "skeletal.muscle.tissue" , "amygdala", 
                    "breast", "uterus","ectocervix","endocervix")
tissue_cr      <- tissue_cr[,tissue_filter]

#Order data frames by ENSEMBL
tissue_cr <- tissue_cr[order(tissue_cr$Gene.ID),]
hela_cr <- hela_cr[order(hela_cr$Gene.ID),]

#Usage of a original function that does multiple correlation tests between HeLa RNAseq data and each chosen tissue

library(devtools)
library(easyGgplot2)

TEST <- Multiple_cor_tests(tissue_cr)