# Multiple_correlation_test
<p>Program useful to run multiple correlation tests between a main data, in this case, RNAseq data for HeLa cell line, and various variables present in a data frame, in this case human tissues from GTEx data, all acquired from Atlas Expression.</p>

<p><b>Set working directory:</b> <br>Keep every .txt files in the same folder as the programms ('Auxiliar_funcions.R' and 'Projeto_HeLa_vs_Tecidos.R') then set the folder as the 'working directory' in the RStudio 'New Project' option.</p>

<p><b>R language version used:</b>  4.1.3</p>
<p><b>Libraries used:</b> </p>
<ul>
  <li>org.Hs.eg.db</li>
  <li>devtools</li>
  <li>easyGgplot2 or ggplot2</li>
</ul>

<p><b>Packages needed - install in console:</b><br>
<code>if (!requireNamespace("BiocManager", quietly = TRUE))</code><br>
<code>       install.packages("BiocManager")</code><br>
<code>BiocManager::install("org.Hs.eg.db")</code><br><br>  
<code>install.packages("devtools")</code><br>
<code>library(devtools)</code><br>
<code>install_github("kassambara/easyGgplot2")</code><br><br>
  
> :memo: **NOTE:** You might need RTools (version dependent of R language version) to install packages from sources like github 
<p>The following code installs a package that includes the library ggplot2 (without github):<br>  
<code>install.packages("tidyverse")</code></p>

<b>Raw data acquisition and description:</b>
| Data base         | Description                                                           | Document to source                            |
|-------------------|-----------------------------------------------------------------------|-----------------------------------------------|
| Apid              | Physical interactions protein-protein, human interactome <sup>i</sup> | <a href="http://cicblade.dep.usal.es:8080/APID/init.action " target="_blank">9606_noISI_Q2.txt</a>         |
| Omnipath          | Gene regulation and signaling network                                 | <a href="http://omnipathdb.org/interactions/" target="_blank">omnipathdb.txt</a>            |
| Dorothea          | Transcription factors and targets annotations                         | <a href="https://omnipathdb.org/interactions/?datasets=tfregulons&tfregulons_levels=A,B" target="_blank"> dorothea_AB.txt </a>  |
| Atlas Expression  | Protein expression - proteomics and transcriptomics                   | <ul><li><a href="https://www.ebi.ac.uk/gxa/experiments/E-MTAB-5214/Results " target="_blank">E-MTAB-5214-query-results.tsv</a> <li><a href="https://www.ebi.ac.uk/gxa/experiments/E-MTAB-2706/Results " target="_blank">HeLa_TPM_results.tsv</a></li></ul>|


<sup>i</sup> Data verified with 2 or more experimental evidences
## Authors
<p>Pedro Fanica, Beatriz Simões, Marta Cruz e Diogo Quaresma</p>
Faculdade de Ciências da Universidade de Lisboa - Bioquímica Experimental IV, Licenciatura em Bioquímica
