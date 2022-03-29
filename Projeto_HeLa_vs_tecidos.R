##Reads interaction files obtained through the following data bases: APID, DOROTHEA and OMNIPATH for Homo sp.

apid <- read.delim("9606_noISI_Q2.txt",header=T, stringsAsFactors = F)

omni <- read.delim("omnipathdb.txt",header=T, stringsAsFactors = F)

omni <- omni[omni$is_directed==1,]

dorothea <- read.delim("dorothea_AB.txt",header=T, stringsAsFactors = F)

##Extracts lists of uniprot identifiers of proteins that interact with at least one of NF-kB sub units

# NFKBI (P105 or P50) - P19838
# RELA (P65) - Q04206

source("R/Auxiliar_functions.R")

subunits      <- c("P19838","Q04206")
apid_cln_pair <- c("UniprotID_A", "UniprotID_B")

filtered_apid   <- fd_objects_frm_cln_union(apid, cln =  apid_cln_pair, objs = subunits, vice_versa = TRUE)
interators_apid <- union(filtered_apid$vice, filtered_apid$versa)

##Extracts the same objects, but know, those must be NF-kB subunits' regulators

omni_n_dorothea_cln_pair   <- c("source", "target")
filtered_omni     <- fd_objects_frm_cln_union(omni,     cln = omni_n_dorothea_cln_pair, objs = subunits, vice_versa = TRUE)
filtered_dorothea <- fd_objects_frm_cln_union(dorothea, cln = omni_n_dorothea_cln_pair, objs = subunits, vice_versa = TRUE)

##List with the union of targets that have as source either NFkB1 or RelA 

targets    <- union(filtered_omni$vice, filtered_dorothea$vice)

##Extraction of list of Uniprot indentifiers of proteins that regulate at least one of NF-kB's subunit
#So its a list of sources that have the subunits as targets

regulators <- union(filtered_omni$versa, filtered_dorothea$versa)

##Intersection between lists
  
targets_n_regulators     <- intersect(targets,   regulators)
targets_n_interactors    <- intersect(targets,   interactores_apid)
regulators_n_interactors <- intersect(regulators,interators_apid) 

<<ainda em PT>>
##Evaluates if the intersection lists are bigger then expected through FDR (false discovery rate) determinantion

#construction of variable universe
universe_apid     <- union(apid$UniprotID_A,apid$UniprotID_B)
universe_omni     <- union(omni$source,omni$target)
universe_dorothea <- union(dorothea$source,dorothea$target)
universe          <- union(universe_apid,union(universe_omni,universe_dorothea))

rm(universe_apid, universe_omni, universe_dorothea)

#tamanhos
t <- length(todos)
a <- length(targets)
r <- length(reguladores)
i <- length(interactores_apid)
a_i <- length(targets_n_regulators)
a_r <- length(alvos_e_reguladores)
r_i <- length(reguladores_e_interactores)

##p_values>>
#lower tail = F, logo  P[X > x]
#q = tamanho - 1, para incluir o inteiro de "listas duplas"
p_a_i <- phyper(a_i-1,a,t-a,i,lower.tail=F)
p_a_r <- phyper(a_r-1,a,t-a,r,lower.tail=F)
p_r_i <- phyper(r_i-1,r,t-r,i,lower.tail=F) #NÃO USAR - muitos dos reguladores 
#exercem a sua função por interação física

#Há um grande número de proteínas, maior do que o esperado, como interatores do NFkB - pode ter vários ciclos de retroação
#obtêm-se p_values muito baixos, diferenças significativas (não são dados aletórios)


##Retirar às três listas de genes/proteínas aquelas que aparecem em mais do que uma lista

ciclosretro <- union(alvos_e_interactores, alvos_e_reguladores) #Inclui todos os candidatos de compostos pertencentes a ciclos de retroação com o NF-kB


##Criar e intersetar a base de dados de interatores e reguladores com os dados de expressão por RNAseq para HeLa e tecidos humanos

library(org.Hs.eg.db)

#from: RNA-seq data of 675 commonly used human cancer cell lines, in Atlas expression
helaprot <- read.delim("HeLa_TPM_results.tsv",header=T, skip = 4 ,stringsAsFactors = F)
colnames(helaprot)[3] <- "HeLa"
helaprot <- helaprot[!is.na(helaprot$HeLa),]

#from: RNA-seq data of 53 human tissues from paper GTEx, in Atlas expression
tissue_exp           <- read.delim("E-MTAB-5214-query-results.tsv",header=T, skip = 4 ,stringsAsFactors = FALSE)
interseção           <- intersect(ciclosretro_ensemble, intersect(helaprot$Gene.ID, tissue_exp$Gene.ID))

#Uso de lista de indentificador ENSEMBL, presente nos dois data frames em estudo, em vez do UNIPROT
ciclosretro_ensemble <- mapIds(org.Hs.eg.db, keys=ciclosretro,column="ENSEMBL", keytype = "UNIPROT")

hela_cr        <- helaprot[index_gene_help(helaprot),]
tissue_cr      <- tissue_exp[index_gene_help(tissue_exp),]

lista          <- c("Gene.ID", "Gene.Name","liver", "pituitary.gland", "hypothalamus", "stomach", "pancreas", "urinary.bladder",
                    "prostate.gland", "adrenal.gland", "lung", "spleen", "amygdala", "vagina", "blood", 
                    "breast", "uterus", "ovary", "thyroid.gland","ectocervix","endocervix")
tissue_cr      <- tissue_cr[,lista]

#ordenação

tissue_cr <- tissue_cr[order(tissue_cr$Gene.ID),]
hela_cr <- hela_cr[order(hela_cr$Gene.ID),]

#Criação de função que faça testes múltiplos de correlação contra o data frame hela_cr

library(devtools)
library(easyGgplot2)


Multiple_cor_tests(tissue_cr)