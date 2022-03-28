
#Creation of function that finds objects, like Uniprot Identifiers, 
#that have a desired correspondence (objs)
#from (a) certain column(s) of a data frame (cln), and returns a list 
#of all the objects found (union)

fd_objects_frm_cln_union <- function(df,cln =  c(), objs = c()){
  list_of_lines <-  c()
  for (o in 1:length(objs)){
    lines <- c()
    for (c in 1:length(cln)){
      
    }
  }
}



#

index_gene_help <- function(tabela, amostra=interseção){
  linhas = c()
  for (i in 1:length(tabela$Gene.ID)){
    if (tabela$Gene.ID[i] %in% amostra){
      linhas <- append(linhas,i)
    }
  }
  return(linhas)
}




#

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
