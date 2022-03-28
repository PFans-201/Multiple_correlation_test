##Importar os ficheiros das interacções para o R

apid <- read.delim("9606_noISI_Q2.txt",header=T, stringsAsFactors = F)

omni <- read.delim("omnipathdb.txt",header=T, stringsAsFactors = F)

omni <- omni[omni$is_directed==1,]

dorothea <- read.delim("dorothea_AB.txt",header=T, stringsAsFactors = F)

##Extrair a lista de proteínas que interagem fisicamente com uma das subunidades do NFkB>>

# NFKBI (P105 -> P50) - P19838
# RELA (P65) - Q04206

linhasAPID_N1 <- which(apid$UniprotID_A=="P19838")   
linhasAPID_R1 <- which(apid$UniprotID_A=="Q04206")   
linhasAPID_1 <- union(linhasAPID_N1,linhasAPID_R1)

interactores_apid_1 <- apid$UniprotID_B[linhasAPID_1]

linhasAPID_N2 <- which(apid$UniprotID_B=="P19838")   
linhasAPID_R2 <- which(apid$UniprotID_B=="Q04206")   
linhasAPID_2 <- union(linhasAPID_N2,linhasAPID_R2)  

interactores_apid_2 <- apid$UniprotID_A[linhasAPID_2]

interactores_apid <- union(interactores_apid_1, interactores_apid_2)

rm(interactores_apid_1, interactores_apid_2, linhasAPID_1, linhasAPID_2, 
   linhasAPID_N1, linhasAPID_N2, linhasAPID_R1, linhasAPID_R2)

##Extrair a lista de genes que são regulados pelo NFkB>>

linhas_alvos_omni_N <- which(omni$source=="P19838")   
linhas_alvos_omni_R <- which(omni$source=="Q04206")    

alvos_omni <- union(omni$target[linhas_alvos_omni_N],omni$target[linhas_alvos_omni_R])

rm(linhas_alvos_omni_N,linhas_alvos_omni_R)

##lista com a união dos alvos que tem como fontes NFKB ou RELA>>

linhas_alvos_dorothea_N <- which(dorothea$source=="P19838")
linhas_alvos_dorothea_R <- which(dorothea$source=="Q04206")

alvos_dorothea <- union(dorothea$target[linhas_alvos_dorothea_N],dorothea$target[linhas_alvos_dorothea_R])

alvos <- union(alvos_omni,alvos_dorothea)
rm(linhas_alvos_dorothea_N,linhas_alvos_dorothea_R,alvos_dorothea,alvos_omni)

##Extrair a lista de genes que regulam o NFkB>>
linhas_reguladores_omni_N <- which(omni$target=="P19838")
linhas_reguladores_omni_R <- which(omni$target=="Q04206")

reguladores_omni <- union(omni$source[linhas_reguladores_omni_N],omni$source[linhas_reguladores_omni_R])


linhas_reguladores_dorothea_N <- which(dorothea$target=="P19838")
linhas_reguladores_dorothea_R <- which(dorothea$target=="Q04206")

reguladores_dorothea <- union(dorothea$source[linhas_reguladores_dorothea_N],dorothea$source[linhas_reguladores_dorothea_R])

reguladores <- union(reguladores_omni,reguladores_dorothea)

rm(linhas_reguladores_omni_N,linhas_reguladores_omni_R,reguladores_omni, linhas_reguladores_dorothea_N, 
   linhas_reguladores_dorothea_R, reguladores_dorothea, reguladores_omni)

##Fazer as intersecções entre essas listas>>

alvos_e_reguladores <- intersect(alvos,reguladores)
alvos_e_interactores <- intersect(alvos,interactores_apid)
reguladores_e_interactores <- intersect(reguladores,interactores_apid)  

##Avaliar se as intersecções são maiores do que o esperado ao acaso com um teste hipergeométrico>>

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

##Retirar às três listas de genes/proteínas aquelas que aparecem em mais do que uma lista>>

ciclosretro <- union(alvos_e_interactores, alvos_e_reguladores) #Inclui todos os candidatos de compostos pertencentes a ciclos de retroação com o NF-kB

##Criar e intersetar a base de dados de interatores e reguladores com os dados de expressão por RNAseq para HeLa e tecidos humanos

#library(clusterProfiler)
library(org.Hs.eg.db)

#from: RNA-seq of 675 commonly used human cancer cell lines, in Atlas expression
helaprot <- read.delim("HeLa_TPM_results.tsv",header=T, skip = 4 ,stringsAsFactors = F)
colnames(helaprot)[3] <- "HeLa"
helaprot <- helaprot[!is.na(helaprot$HeLa),]

ciclosretro_ensemble <- mapIds(org.Hs.eg.db, keys=ciclosretro,column="ENSEMBL", keytype = "UNIPROT")
tissue_exp           <- read.delim("E-MTAB-5214-query-results.tsv",header=T, skip = 4 ,stringsAsFactors = FALSE)
interseção           <- intersect(ciclosretro_ensemble, intersect(helaprot$Gene.ID, tissue_exp$Gene.ID))

index_gene_help <- function(tabela, amostra=interseção){
  linhas = c()
  for (i in 1:length(tabela$Gene.ID)){
    if (tabela$Gene.ID[i] %in% amostra){
     linhas <- append(linhas,i)
    }
  }
  return(linhas)
}

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
#library("ggpubr")

Multiple_cor_tests <- function(test_table){
  names <- colnames(test_table)
  df <- data.frame()
  not_possible <- c()
  for (i in 3:length(names)){
    correct = c()
    for (l in 1:length(test_table[,i])){
      if (is.na(test_table[l,i])==FALSE && test_table[l,i]>1) {
        correct <- append(correct,test_table[l,1])
      }
    }
    if (length(correct)<10){
      not_possible <- append(not_possible, names(i))
    }
    else{
      amostra_hl <- hela_cr[index_gene_help(hela_cr, amostra = correct),3]
      amostra_table_i <- test_table[index_gene_help(test_table, amostra = correct), i]
      
      teste <- cor.test(amostra_hl, amostra_table_i, method="pearson")
      df[i-2, 1] <- names[i]
      df[i-2, 2] <- teste$estimate
      df[i-2, 3] <- teste$p.value
      df[i-2, 4] <- length(correct)
    }
  }  
  
  if (length(not_possible)!= 0){
    print("The following tissues were not elegible for correlation test with the data for the cancerous cell line HeLa:", not_possible)
  }
  
  colnames(df) <- c("Tissues_to_compare","Correlation_coefficient","p_value", "Dimension")
  df <- df[order(df$Correlation_coefficient),]
  print(df)

  plot <- ggplot(df, aes(x = Correlation_coefficient, y=reorder(Tissues_to_compare, + Correlation_coefficient), fill = p_value)) +
          labs(x = "Coeficiente de correlação, cf. HeLa", y = "Tecidos estudados")+
          scale_fill_gradientn(trans = "log", breaks=c(1.0e-25, 1.0e-20, 1e-15, 1.0e-10, 1.0e-5, 1.0e-1), colors = colorspace::diverge_hcl(7))+
          geom_bar(stat="identity",color="black")+
          guides(fill = guide_colourbar(barwidth = 0.5, barheight = 10))+
          theme_classic()
  
  return(plot)
}

Multiple_cor_tests(tissue_cr)