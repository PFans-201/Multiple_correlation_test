# Multiple_correlation_test
<p>Program useful to run multiple correlation tests between a main data, in this case, RNAseq data for HeLa cell line, and various variables present in a data frame, corresponding here to human tissues from GTEx data, all acquired from Atlas Expression.</p>

<p><b>Set working directory:</b> <br>Keep every .txt file in the same folder as the programm then set the folder as the 'working directory' in the RStudio 'New Project' option.</p>

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
  
<b>NOTE:</b> You might need RTools, to install packages from sources like github<br><br> 
The following code installs a package that includes the library ggplot2 (without github):<br>  
<code>install.packages("tidyverse")</code><br>
</p>

## Authors
<p>Pedro Fanica, Beatriz Simões, Marta Cruz e Diogo Quaresma</p>
Faculdade de Ciências da Universidade de Lisboa - Bioquímica Experimental IV, Licenciatura em Bioquímica
