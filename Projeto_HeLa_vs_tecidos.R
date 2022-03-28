##Reads interaction files obtained through the following data bases: APID, DOROTHEA and OMNIPATH for Homo sp.

apid <- read.delim("9606_noISI_Q2.txt",header=T, stringsAsFactors = F)

omni <- read.delim("omnipathdb.txt",header=T, stringsAsFactors = F)

omni <- omni[omni$is_directed==1,]

dorothea <- read.delim("dorothea_AB.txt",header=T, stringsAsFactors = F)

##Extracts lists of uniprot identifiers of proteins that interact with at least one of NF-kB sub units

# NFKBI (P105 or P50) - P19838
# RELA (P65) - Q04206

source("R/Auxiliar_functions.R")


linesAPID_N1 <- which(apid$UniprotID_A=="P19838")   
linesAPID_R1 <- which(apid$UniprotID_A=="Q04206")   
linesAPID_1 <- union(linesAPID_N1,linesAPID_R1)

interactors_apid_1 <- apid$UniprotID_B[linesAPID_1]

linesAPID_N2 <- which(apid$UniprotID_B=="P19838")   
linesAPID_R2 <- which(apid$UniprotID_B=="Q04206")   
linesAPID_2 <- union(linesAPID_N2,linesAPID_R2)  

interactors_apid_2 <- apid$UniprotID_A[linesAPID_2]

interactors_apid <- union(interactors_apid_1, interactors_apid_2)

rm(interactors_apid_1, interactors_apid_2, linesAPID_1, linesAPID_2, 
   linesAPID_N1, linesAPID_N2, linesAPID_R1, linesAPID_R2)

#We opted to remove, rm(), some unnecessary variables through out the programm

##Extracts the same objects, but know, those must be NF-kB subunits' regulators


lines_targets_omni_N <- which(omni$source=="P19838")   
lines_targets_omni_R <- which(omni$source=="Q04206")    

targets_omni <- union(omni$target[lines_targets_omni_N],omni$target[lines_targets_omni_R])

rm(lines_targets_omni_N,lines_targets_omni_R)

##List with the union of targets that have as source either NFkB1 or RelA

lines_targets_dorothea_N <- which(dorothea$source=="P19838")
lines_targets_dorothea_R <- which(dorothea$source=="Q04206")

targets_dorothea <- union(dorothea$target[lines_targets_dorothea_N],dorothea$target[lines_targets_dorothea_R])

targets <- union(targets_omni,targets_dorothea)
rm(llines_targets_dorothea_R,lines_targets_dorothea_N,targets_dorothea,targets_omni)

##Extraction of list of Uniprot indentifiers of proteins that regulate at least one of NF-kB's subunit
lines_regulators_omni_N <- which(omni$target=="P19838")
lines_regulators_omni_R <- which(omni$target=="Q04206")

regulators_omni <- union(omni$source[lines_regulators_omni_N],omni$source[lines_regulators_omni_R])


lines_regulators_dorothea_N <- which(dorothea$target=="P19838")
lines_regulators_dorothea_R <- which(dorothea$target=="Q04206")

regulators_dorothea <- union(dorothea$source[lines_regulators_dorothea_N],dorothea$source[linhas_reguladores_dorothea_R])

reguladores <- union(reguladores_omni,reguladores_dorothea)

rm(linhas_reguladores_omni_N,linhas_reguladores_omni_R,reguladores_omni, linhas_reguladores_dorothea_N, 
   linhas_reguladores_dorothea_R, reguladores_dorothea, reguladores_omni)

##Intersecção entre essas listas

alvos_e_reguladores <- intersect(alvos,reguladores)
alvos_e_interactores <- intersect(alvos,interactores_apid)
reguladores_e_interactores <- intersect(reguladores,interactores_apid)  

##Avaliar se as intersecções são maiores do que o esperado ao acaso com um teste hipergeométrico - FDR (false discovery rate)

#lista de todos os alvos
todos_apid <- union(apid$UniprotID_A,apid$UniprotID_B)
todos_omni <- union(omni$source,omni$target)
todos_dorothea <- union(dorothea$source,dorothea$target)
todos <- union(todos_apid,union(todos_omni,todos_dorothea))

rm(todos_apid, todos_omni, todos_dorothea)

#tamanhos
t <- length(todos)
a <- length(alvos)
r <- length(reguladores)
i <- length(interactores_apid)
a_i <- length(alvos_e_interactores)
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

rm(a_r, a_i, r_i, t, a, r, i)

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

#lista          <- colnames(tissue_exp[7:55])
lista          <- c("liver", "pituitary.gland", "hypothalamus", "stomach", "pancreas", "urinary.bladder", "prostate.gland", "adrenal.gland", "lung", "spleen", "amygdala", "vagina", "blood", "breast", "uterus", "ovary", "thyroid.gland","ectocervix","endocervix")
lista          <- append(lista, c("Gene.ID", 'Gene.Name'),0)  # 0 garante que adiciona no início da tabela
tissue_cr      <- tissue_cr[,lista]

#ordenação

tissue_cr <- tissue_cr[order(tissue_cr$Gene.ID),]
hela_cr <- hela_cr[order(hela_cr$Gene.ID),]

#Criação de função que faça testes múltiplos de correlação contra o data frame hela_cr

library(devtools)
library(easyGgplot2)


Multiple_cor_tests(tissue_cr)