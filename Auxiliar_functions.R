
#Creation of function that finds objects, like Uniprot Identifiers, 
#that have a desired correspondence (objs)
#from (a) certain column(s) of a data frame (cln), and returns a list 
#of all the objects found (union)

fd_objects_frm_cln_union <- function(df, cln =  c(search_, get_), objs = c(), vice_versa = TRUE){
  #ARGUMENTS:
  #df- dataframe/data; cln - colnumn names (str/char format) of df,
  #objs - vector that include the itemns to find,
  #vice_versa can be TRUE or FALSE, if TRUE (default), the function will 
  #search and get items from both columns indicated in cln
  
  list_of_lines  <-  c()
  list_of_lines_ <-  c()
 
  #NOTE: cln is a vetctor of paired columns, one to search information, and another to get information from
  for (o in objs){
    if (vice_versa == TRUE){
      list_of_lines  <- append(list_of_lines,  which(df[,cln[1]]==o))
      list_of_lines_ <- append(list_of_lines_, which(df[,cln[2]]==o))
    }
    else{
      list_of_lines  <- append(list_of_lines,  which(df[,cln[1]]==o))
    }
  }
  if (is.null(list_of_lines_) == FALSE){
    vice    <- unique(df[list_of_lines , cln[2]])
    versa   <- unique(df[list_of_lines_, cln[1]])
    
    #organizes vice and versa versions of cln in separate parts of a list called results
    results <- list("vice" = vice,"versa" = versa)
  }
  else {
    results <- unique(df[list_of_lines , cln[2]]) 
  }
  
  return(results)
}


  
#Function specialized in finding gene.id's from a given sample in a data frame
#It is a simplified version of the previous function and more directed for a specific problem 

index_gene_help <- function(df, sample=intersection){
  #ARGUMENTS:
  #df - data frame, sample - normally a vector of characters obtained by intersection of other data
  
  linhas = c()
  for (i in 1:length(df$Gene.ID)){
    if (df$Gene.ID[i] %in% sample){
      linhas <- append(linhas,i)
    }
  }
  return(linhas)
}



#Program useful to run multiple correlation tests between a main data, in this case, 
#RNAseq data for HeLa cell line, and various variables present in a data frame, corresponding here to
#human tissues from GTEx data, all acquired from Atlas Expression

Multiple_cor_tests <- function(test_table, base_table = hela_cr, xynames = c("Correlation coefficient, cf. HeLa", "Studied Tissues")){
  #ARGUMENTS:
  #test_table - table that contain multiple variables to test with base_table
  #base_table - hela_cr by default, xynames - vector that includes the names of x and y plot's axis labels, in this order 
  #(default values for the purpose of this whole program)
  
  
  names <- colnames(test_table)
  df <- data.frame()
  not_possible <- c()
  for (i in 3:length(names)){
    #in expression table, the first 2 columns concern gene id and gene name, that's why we start at 3
    corrected = c()
    for (l in 1:length(test_table[,i])){
      #we establish TPM inferior to 1.0 as not significant expression data for this test
      if (is.na(test_table[l,i])==FALSE && test_table[l,i]>1) {
        #the following character vector includes every ENSEMBL identifier which TPM followed the previous conditions 
        corrected <- append(corrected,test_table[l,1])
      }
    }
    if (length(corrected)<10){
      #we establish that if the dimension of the data for the test is inferior to 10, the data from the selected 
      #tissue cannot be used by this function, in other words, it cannot calculate the correlation coefficient  
      not_possible <- append(not_possible, names(i))
    }
    else{
      #acquisition of expression data for each test_table's column and the base's only column of interest
      sample_base <- base_table[index_gene_help(base_table, sample = corrected),3]
      sample_table_i <- test_table[index_gene_help(test_table, sample = corrected), i]
      
      #correlation test itself
      test  <- cor.test(sample_base, sample_table_i, method="spearman", exact= FALSE)
      test1 <- cor.test(sample_base, sample_table_i, method="kendall", exact= FALSE)
      
      #normality test
      test2 <- shapiro.test(sample_table_i) 
      
      #organization of information in the data frame, df, initially empty, in order to visualize it in the end
      df[i-2, 1] <- names[i]
      df[i-2, 2] <- test$estimate
      df[i-2, 3] <- test$p.value
      df[i-2, 4] <- test1$estimate
      df[i-2, 5] <- test2$p.value
      df[i-2, 6] <- length(corrected)
      df[i-2, 7] <- test2$p.value
    }
  }  
  
  if (length(not_possible)!= 0){
    print("The following tissues were not elegible for correlation test with the data for the cancerous cell line HeLa:", not_possible)
  }
  
  colnames(df) <- c("Tissues_to_compare","Spearman_Cor_coef","S_p_value", "Kendall_Cor_coef","K_p_value","Dimension", "SW Normality_test")
  df <- df[order(df$Spearman_Cor_coef),]
  
  x_name <- xynames[1]
  y_name <- xynames[2]
  
  #Function explained in ggplot 2 forums and main web page
  plot <- ggplot(df, aes(x = Spearman_Cor_coef, y=reorder(Tissues_to_compare, + Spearman_Cor_coef), fill = S_p_value)) +
    labs(x = x_name, y = y_name, fill = "p value")+
    scale_fill_gradientn(trans = "log",breaks=c(0.1,0.05,0.01, 0.001,0.0001), colors = colorspace::diverge_hcl(7))+
    geom_bar(stat="identity",color="black")+
    guides(fill = guide_colourbar(barwidth = 0.5, barheight = 10))+
    theme_classic()
   
  #organization of the plot and data frame in a list, giving the possibility to see both as a return for this function
  results <- list("plot"=plot, "dataframe" = df)
  return(results)
}
