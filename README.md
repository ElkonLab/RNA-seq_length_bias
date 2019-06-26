# Length_bias
Imports: cqn, edgeR

For full running example go to the Running_example.md file.

There are 4 steps in the anslysis:

# Step 1:

# function: identify_genes_with_leninfo_and_cpm_above_thresh()
Keep genes:
(1) Genes with length info
(2) At least 1.0 cpm in all replicates of at least one biological condition in the dataset

# Step 2:

# function: length_bias_analysis_compare_biological_conditions()

 Examine relationship between FC and Tx length using the following processing methods:
1. RPKM (add 1 pseudo count)
2. RPKM.qnorm
3. edgeR TMM
4. edgeR RLE
5. edgeR uq
The function assumes that the dataset contains 2 biological conditions (control and treatment) with equal numbers of replicates, ordered in the input files such Ens IDs are row names, the 1st col is Sym, then cols 2:(ncols/2) are for cond1 and (ncols/2)+1:ncols for cond2.  

# Step 3:

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

# Step 4:

# function: length_bias_analysis_cqn_normalization()
Examine relationship between FC and Tx length for replicate samples  and biological conditions using the following processing methods:
1. edgeR TMM
2. CQN Normalization - only with GC content
3. CQN Normalization - with GC content and genes length
The function assumes that the dataset contains 2 biological conditions (control and treatment)
and each with N replicates, ordered in the input files such that 1st col is Sym, then cols 2:N.replicates are for cond1 and N.replicates+1:ncol for cond2.  
Set Rep.first.col=2 to analyze the replicates of "control" samples
Set Rep.first.col=N.replicates+2 to analyze the replicates of "treatmnet" samples
