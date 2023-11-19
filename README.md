# **Enriched Single-Nucleus RNA-Sequencing reveals unique attributes of distal convoluted tubule cells**
__*Xiao-Tong Su, Jeremiah Reyes, Anne Lackey, Yujiro Maeoka, Ryan Cornelius, Jim McCormick, Chao-Ling Yang, Hyun Jung, Paul Welling, *Jonathan Nelson, *David H. Ellison__  
*Corresponding authors  

If you use any of the code or workflows in this repository please cite our manuscript in xxx (will update).

Single cell sequencing data generated for this manuscript and finalized R objects can be downloaded in GEO. </br>
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228367

To ensure data accessibility to non-bioinformaticians, we made the DCT snRNA-Seq data available for further exploration via an interactive web tool generated using ShinyCell at https://ellisonlab.shinyapps.io/dct_shinycell/

Welcome to our GitHub repository!  
Here you will find analysis scripts for our manuscript. Please contact any corresponding authors, Drs. Xiao-Tong Su, Jonathan W. Nelson, and David H. Ellison, with questions or comments. 


<br/>
Thanks,
<br/>
Xiao-Tong,
<br/><br/>

Visit the Ellison lab website:<br/>
https://www.ohsu.edu/school-of-medicine/ellison-lab
<br/>

Check out the DCT snRNA-Seq data interactive web tool generated using ShinyCell at:<br/>
https://ellisonlab.shinyapps.io/dct_shinycell/
<br/>

**How to use the snRNA data analysis scripts**
1. Pre-process the aligned data and aggregate three control samples or all six samples (three control and three experimental. <br/>
Pre-process_Step1.R <br/>
Pre-process_Step2.R <br/>

3. Analyze the aggreaged and annotated dataset <br/>
Analysis_1_AddModuleScore.R <br/>
Analysis_2_Subpopulation.R <br/>
Analysis_3_Find_DEGs.R <br/>
Analysis_4_Pseudotime.R <br/>
Analysis_5_Pathway_analysis.R <br/>

4. Make figures for the manuscript. <br/>
Fig.1.R <br/>
Fig.2.R <br/>
Fig.3.R <br/>
Fig.4.R <br/>
Fig.5.R <br/>
Fig.6.R <br/>
