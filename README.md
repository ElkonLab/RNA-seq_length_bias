# Length_bias
Imports: cqn, edgeR


1. mkdir length_bias\
2. cd length_bias\
3. git clone https://github.com/ElkonLab/RNA-seq_length_bias.git\
or directly download the folder

For running example see Running_example.html
or create it by yourself in different format using R markdown with the file Running_example.Rmd

There are 4 steps in the anslysis:

# Step 1 : Filter genes

# function: identify_genes_with_leninfo_and_cpm_above_thersh()

Create filtered counts file with genes that have at least 1.0 cpm in all replicates of at least one biological condition in the dataset, and have a length and GC content data.
The output is cnts_data_Pgenes file, this file is being use as a counts file input for the next steps.

# Step 2 : compare biological conditions
The count file is the filtered counts file from the previous step (cnts_data_Pgenes file) and it need to be located in the wd.
The output includes 2 files:
treatment_vs_control_analysis_results.txt
biological_conditions_analysis_plots.pdf

# function: length_bias_analysis_compare_biological_conditions()

 Examine relationship between FC and Tx length using the following processing methods:
1. RPKM (add 1 pseudo count)
2. RPKM.qnorm
3. edgeR TMM
4. edgeR RLE
5. edgeR uq
The function assumes that the dataset contains 2 biological conditions (control and treatment) with equal number of replicates, ordered in the input files such that 1st col is Sym, then cols 2:(ncols/2) are for cond1 and (ncols/2)+1:ncols for cond2.  

# Step 3 : compare replicate samples 
The count file is the filtered counts file from the step 1 (cnts_data_Pgenes file) and it need to be located in the wd.
The analysis is on the treatment replicates. For control replicates set Rep.first.col to 2.
The output includes 2 files:
replicate_samples_analysis_results.txt
replicate_samples_analysis_plots.pdf

# function: length_bias_analysis_compare_replicate_samples()

Examine relationship between FC and Tx length for replicate samples using the following processing methods:
1. RPKM (add 1 pseudo count)
2. RPKM.qnorm
3. edgeR TMM
4. edgeR RLE
5. edgeR uq
The 1st col is Sym, then cols 2:n are for cond1 and (n+1):ncol for cond2.  
Set Rep.first.col=2 to analyze the three replicates of "control" samples
Set Rep.first.col=N.replicates+2 to analyze the three replicates of "treatmnet" samples

# Step 4 : compare biological and replicate samples after cqn normalization
The count file is the filtered counts file from the step 1 (cnts_data_Pgenes file) and it need to be located in the wd.
The analysis is on the treatment replicates. For control replicates set Rep.first.col to 2.
The output includes 4 files:
cqn_replicate_samples_analysis_results.txt
CQN_replicate_samples_analysis_plots.pdf
cqn_treatment_vs_control_analysis_results.txt
CQN_biological_conditions_analysis_plots.pdf

# function: length_bias_analysis_cqn_normalization()
Examine relationship between FC and Tx length for replicate samples  and biological conditions using the following processing methods:
1. CQN Normalization 
The function assumes that the dataset contains 2 biological conditions (control and treatment)
and each with N replicates, ordered in the input files such that 1st col is Sym, then cols 2:N.replicates are for cond1 and N.replicates+1:ncol for cond2.  
Set Rep.first.col=2 to analyze the replicates of "control" samples
Set Rep.first.col=N.replicates+2 to analyze the replicates of "treatmnet" samples
