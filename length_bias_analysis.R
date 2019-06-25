library(edgeR)
library(cqn)

# Step 1:
# --------
# function: identify_genes_with_leninfo_and_cpm_above_thersh()
#
# Filter genes:
# (1) Genes with length info
# (2) At least 1.0 cpm in all replicates of at least one biological condition in the dataset
#


identify_genes_with_leninfo_and_cpm_above_thresh <- function(countsFile, Conds,N.samples.cutoff=3, CPM.CUTOFF=1.0, Hs.Tx.LenFile, DatasetTitle="GSE64233_TNFa_p65kd_samples") { 
  
  # Load count data
  data.cnts <- read.file(countsFile)
  
  # Load Gene length file
  Hs.genes.len <- read.table(file=Hs.Tx.LenFile, sep="\t", header=1, row.names=1, stringsAsFactors=F); 
  
  # Consider subset of genes with length info
  cGenes <- intersect(rownames(data.cnts), rownames(Hs.genes.len))
  
  # Calculate cpm levels (considering only genes with length info)
  data.cpm <- cpm(data.cnts[cGenes, 2:ncol(data.cnts)])
  
  # Identify "robustly expressed genes" in the dataset
  I.Pgenes <- find_genes_above_count_thresh_in_Conds_subset(counts.mat=data.cpm, Conds, N.cutoff=CPM.CUTOFF, N.samples.cutoff=N.samples.cutoff) 
  
  # The FILTERED count data (should be used as input to the 2nd step) 
  data.cnts.Pgenes <- data.cnts[rownames(data.cpm[I.Pgenes,]), ]
  
  outFile <- sprintf("%s_cnts_data_Pgenes.txt", DatasetTitle)
  str <- sprintf("\nWriting Pgenes count file to %s\n\n", outFile)
  cat(str, "\n")
  write.table(data.cnts.Pgenes, file=outFile, sep="\t", quote=F, col.names=NA)
}
#-----------------------------------------------------

#
# Step 2:
#---------
# function: length_bias_analysis_compare_biological_conditions()
#
# Examine relationship between FC and Tx length using the following processing methods:
# 1. RPKM (add 1 pseudo count)
# 2. RPKM.qnorm
# 3. edgeR TMM
# 4. edgeR RLE
# 5. edgeR uq
#
#
# the function assumes that the dataset contains 2 biological conditions (control and treatment)
# with equal number of replicates, ordered in the input files such that ,have Ens id as row names, 1st col is Sym, then cols 2:(ncols/2) are for cond1 and (ncols/2)+1:ncols for cond2.  
#

length_bias_analysis_compare_biological_conditions <- function (countsFile.Pgenes, Hs.Tx.LenFile, DatasetTitle="TNFa.p65kd", X.lim=c(8, 15), Y.lim=c(-2,2)) {
  
  #pdfFile <- sprintf("%s_biological_conditions_analysis_plots.pdf", DatasetTitle);
  #pdf(pdfFile)
  par(mfrow=c(2,3))
  P.bio <- array(dim=c(1, 6))
  Spearman.r.bio <- array(dim=c(1, 6))
  
  colnames(P.bio)          <- c("RPKM", "RPKM.qnorm", "edgeR.TMM", "edgeR.RLE", "RLE.RPKM", "UQ.RPKM")
  colnames(Spearman.r.bio) <- c("RPKM", "RPKM.qnorm", "edgeR.TMM", "edgeR.RLE", "RLE.RPKM", "UQ.RPKM")
  
  
  # Load Pgenes count data
  data.cnts.Pgenes <- read.file(countsFile.Pgenes)
  #get the number of replicates
  n_replicates <- (ncol(data.cnts.Pgenes)-1) / 2
  # Load Gene length file
  Hs.genes.len <- read.table(file=Hs.Tx.LenFile, sep="\t", header=1, row.names=1, stringsAsFactors=F); 
  
  # Consider subset of genes with length info
  cGenes <- intersect(rownames(data.cnts.Pgenes), rownames(Hs.genes.len))
  
  # Create DGEList list object: it contains the 3 componenets: (1) counts, (2) group = condition set and (3) genes = with gene annots
  #groupf <- factor(c(rep("c",3), rep("tnfa",3), rep("c.p65kd",3), rep("tnfa.p65kd",3)))
  groupf <- factor(c(rep("control",n_replicates), rep("treatment",n_replicates)))
  dge <- DGEList(counts=data.cnts.Pgenes[cGenes, 2:ncol(data.cnts.Pgenes)], group=groupf, genes=Hs.genes.len[cGenes, ])
  
  
  # 
  # RPKM analysis analysis
  #
  RPKM <- RPKM_norm(dge)
 
  Conds <- list("control"=1:n_replicates, "treatment"=(n_replicates+1):(n_replicates*2))
  RPKM.logFC <- calcLogFC(RPKM,Conds)
  
  main.title <- sprintf("%s: RPKM", DatasetTitle);
  r <- plotFCvsLen(logTxLen=log(dge$genes[, "Tx.len"],2), logFC=RPKM.logFC, main.title,xLab="Length (log2)",X.lim, Y.lim)
  P.bio[1, 1]          <- r$p.value; 
  Spearman.r.bio[1, 1] <- r$estimate; 
  
  # RPKM+qnorm
  RPKM.qnorm       <- normalizeQuantiles(RPKM)
  RPKM.qnorm.logFC <- calcLogFC(RPKM.qnorm,Conds)
  
  main.title <- sprintf("%s: RPKM-qnorm", DatasetTitle);
  r <- plotFCvsLen(logTxLen=log(dge$genes[, "Tx.len"],2), logFC=RPKM.qnorm.logFC, main.title,xLab="Length (log2)", X.lim, Y.lim)
  P.bio[1, 2]          <- r$p.value; 
  Spearman.r.bio[1, 2] <- r$estimate; 
  
  # 
  # edgeR-TMM analysis
  #
  # Design matrix
  design <- model.matrix(~groupf)
  
  #TMM norm with Estimation of dispersion
  dge.tmm <-norm_dge(dge,"TMM", design)
  
  # Fit NB model and perform LR test (here group factor has 2 levels, hence the fit has 2 coef. First is the Intercept. 
  lrt <- fitNBandLRtest(dge.tmm,design)
  
  # Results including log2FC are in lrt$table
  
  edgeR.TMM.logFC <- lrt$table[,1];  
  main.title <- sprintf("%s: edgeR-TMM", DatasetTitle);
  r <- plotFCvsLen(logTxLen=log(dge$genes[, "Tx.len"],2), logFC=edgeR.TMM.logFC, main.title,xLab="Length (log2)", X.lim, Y.lim)
  P.bio[1, 3]          <- r$p.value; 
  Spearman.r.bio[1, 3] <- r$estimate; 
  
  # edgeR+RLE
  dge.rle  <- norm_dge(dge,"RLE", design)
  
  # Fit NB model and perform LR test (here group factor has 2 levels, hence the fit has 2 coef. First is the Intercept. 
  lrt <- fitNBandLRtest(dge.rle,design)
  
  # Results including log2FC are in lrt$table
 
  edgeR.RLE.logFC <- lrt$table[,1];  
  
  main.title <- sprintf("%s: edgeR-RLE", DatasetTitle);
  r <- plotFCvsLen(logTxLen=log(dge$genes[, "Tx.len"],2), logFC=edgeR.RLE.logFC, main.title,xLab="Length (log2)", X.lim, Y.lim)
  P.bio[1, 4]          <- r$p.value; 
  Spearman.r.bio[1, 4] <- r$estimate; 
  
  
  # RPKM+RLE
  RPKM.rle <- rpkm(dge.rle,  gene.length=dge.rle$genes[, "Tx.len"])
  RPKM.rle.logFC <- calcLogFC(RPKM.rle,Conds)
  
  main.title <- sprintf("%s: RLE-RPKM", DatasetTitle);
  r <- plotFCvsLen(logTxLen=log(dge$genes[, "Tx.len"],2), logFC=RPKM.rle.logFC, main.title,xLab="Length (log2)", X.lim, Y.lim)
  P.bio[1, 5]          <- r$p.value; 
  Spearman.r.bio[1, 5] <- r$estimate; 
  
  
  # RPKM+UQ
  dge.uq  <- norm_dge(dge,"upperquartile",estimateDispersion = FALSE)
  RPKM.uq <- rpkm(dge.uq,  gene.length=dge.uq$genes[, "Tx.len"])
  RPKM.uq.logFC <- calcLogFC(RPKM.uq,Conds)
  
  main.title <- sprintf("%s: UQ-RPKM", DatasetTitle);
  r <- plotFCvsLen(logTxLen=log(dge$genes[, "Tx.len"],2), logFC=RPKM.uq.logFC, main.title,xLab="Length (log2)", X.lim, Y.lim)
  P.bio[1, 6]          <- r$p.value; 
  Spearman.r.bio[1, 6] <- r$estimate; 
  
  res <- rbind(Spearman.r.bio, P.bio)
  rownames(res) <- c("Spearman.r", "pval")
  
  outFile <- sprintf("%s_treatment_vs_control_analysis_results.txt", DatasetTitle);
  write.table(res, file=outFile, sep="\t", quote=F, col.names=NA)
  str <- sprintf("\nResults recorded in %s\n", outFile)
  cat(str, "\n")

  #invisible(dev.off())
  
  res
  
}
#---------------------------------------------------------------------------
# 
# Step 3
#---------
# function: length_bias_analysis_compare_replicate_samples()
#
# Examine relationship between FC and Tx length for replicate samples using the following processing methods:
# 1. RPKM (add 1 pseudo count)
# 2. RPKM.qnorm
# 3. edgeR TMM
# 4. edgeR RLE
# 5. edgeR uq
#
#
# The function assumes that the dataset contains 2 biological conditions (control and treatment)
# and each with n replicates, ordered in the input files such that Ens Ids are row names,1st col is Sym, then cols 2:n are for cond1 and (n+1):ncol for cond2.  
#
# Set Rep.first.col=2 to analyze the three replicates of "control" samples
# Set Rep.first.col=N.replicates+2 to analyze the three replicates of "treatmnet" samples
#
#------------------------------------------------------------

length_bias_analysis_compare_replicate_samples <- function(countsFile.Pgenes, Hs.Tx.LenFile,N.replicates=3, DatasetTitle="TNFa.p65kd", Rep.first.col=5, X.lim=c(8, 15), Y.lim=c(-3,3)) {
  
  #pdfFile <- sprintf("%s_replicate_samples_analysis_plots.pdf", DatasetTitle);
  n_pairs <- (N.replicates * (N.replicates-1))/2
 # pdf(pdfFile)
  par(mfrow=c(3,3))
  P.rep <- array(dim=c(5, n_pairs))
  Spearman.r.rep <- array(dim=c(5, n_pairs))
  rep_names <- c()
  for(i in 1:N.replicates){for(j in 1:i){if (i!=j) {rep_names<-c(rep_names,sprintf("rep%dvs%d",i,j))}}}
  colnames(P.rep)          <- rep_names
  colnames(Spearman.r.rep) <- rep_names
  
  rownames(P.rep)          <- c("edgeR.TMM", "RPKM", "RPKM.qnorm", "RLE.RPKM", "UQ.RPKM")
  rownames(Spearman.r.rep) <- c("edgeR.TMM", "RPKM", "RPKM.qnorm", "RLE.RPKM", "UQ.RPKM")
  
  Rep.last.col  <- Rep.first.col + N.replicates-1; 
  
  # Load Pgenes count data
  data.cnts.Pgenes <- read.file(countsFile.Pgenes)
  
  # Load Gene length file
  Hs.genes.len <- read.table(file=Hs.Tx.LenFile, sep="\t", header=1, row.names=1, stringsAsFactors=F); 
  
  # Consider subset of genes with length info
  cGenes <- intersect(rownames(data.cnts.Pgenes), rownames(Hs.genes.len))
  
  # For dataset with N replicates
  groupf <- factor(1:N.replicates);
  
  # Create DGEList list object: it contains the 3 componenets: (1) counts, (2) group = condition set and (3) genes = with gene annots
  dge <- DGEList(counts=data.cnts.Pgenes[cGenes, Rep.first.col:Rep.last.col], group=groupf, genes=Hs.genes.len[cGenes, ]) 
  log2.Tx.Len <- log(dge$genes[, "Tx.len"],2) #log2 of genes length
  
  # 
  # edgeR analysis
  #
  
  dge.tmm <- norm_dge(dge,"TMM",estimateDispersion = FALSE)
  # Design matrix
  design <- model.matrix(~groupf)
  
  # fit
  fit <- glmFit(dge.tmm, dispersion=0) 
  edgeR.logFC <- fit$coef[, 2]
  
  # Tests: comparisons of all of the replicates to replicate 1
  lrt <- glmLRT(fit, coef=2:N.replicates)
  # Results including log2FC are in lrt$table (comparisons to replicate 1)
  FC.data <- list()
  m <- N.replicates-1
  ind <- 1
  for(i in 1:m){ #comaprisons to replicate 1
    rep.p <- cor.test(log2.Tx.Len, lrt$table[,i], method="spearman")
    P.rep[1, i] <- rep.p$p.value; 
    Spearman.r.rep[1, i] <- rep.p$estimate; 
    FC.data[[ind]] <- lrt$table[,i]
    ind <- ind + 1
  }
  if(N.replicates > 2){
    for(i in 2:N.replicates){
      for(j in 2:i){
        if (i!=j) {
          x <-  FCBetweenReps(fit,N.replicates,i,j) #comaprisons between the differnt pairs of replicates
          rep.p <- cor.test(log2.Tx.Len, x$table[,1], method="spearman")
          P.rep[1,ind] <- rep.p$p.value; 
          Spearman.r.rep[1, ind] <- rep.p$estimate; 
          FC.data[[ind]] <- x$table[,1]
          ind <- ind + 1
        }
      }
    }
  }
  
  edgeR_TMM_plots <- plot_replicate_samples(log2.Tx.Len,FC.data,"edgeR-TMM",DatasetTitle,rep_names,xLab="Length (log2)", X.lim, Y.lim)
  
  
  
  # 
  # RPKM analysis
  #
  RPKM <- RPKM_norm(dge)
  FC.data <- FC_replicate_samples(RPKM)
  RPKM_plots <- plot_replicate_samples(log2.Tx.Len,FC.data,"RPKM",DatasetTitle,rep_names,xLab="Length (log2)", X.lim, Y.lim)
  
  
  for(i in seq(length(FC.data))){
    rep.p <- cor.test(log2.Tx.Len, FC.data[[i]], method="spearman")
    P.rep[2, i] <- rep.p$p.value; 
    Spearman.r.rep[2, i] <- rep.p$estimate; 
  }
  
  
  
  # 
  # RPKM.qnorm analysis
  #
  RPKM.qnorm <- normalizeQuantiles(RPKM)
  FC.data <- FC_replicate_samples(RPKM.qnorm)
  RPKM.qnorm_plots <- plot_replicate_samples(log2.Tx.Len,FC.data,"RPKM.qnorm",DatasetTitle,rep_names,xLab="Length (log2)", X.lim, Y.lim)
  
  for(i in seq(length(FC.data))){
    rep.p <- cor.test(log2.Tx.Len, FC.data[[i]], method="spearman")
    P.rep[3, i] <- rep.p$p.value; 
    Spearman.r.rep[3, i] <- rep.p$estimate; 
  }
  
  
  # 
  # RPKM.RLE analysis
  #
  dge.rle  <- norm_dge(dge,"RLE",design,estimateDispersion = FALSE)
  RPKM.rle <- rpkm(dge.rle,  gene.length=dge.rle$genes[, "Tx.len"])
  RPKM.rle <- RPKM.rle+1; 
  FC.data <- FC_replicate_samples(RPKM.rle)
  
  for(i in seq(length(FC.data))){
    rep.p <- cor.test(log2.Tx.Len, FC.data[[i]], method="spearman")
    P.rep[4, i] <- rep.p$p.value; 
    Spearman.r.rep[4, i] <- rep.p$estimate; 
  }
  
  
  # 
  # RPKM.UQ analysis
  #
  dge.uq  <- norm_dge(dge,"upperquartile",design,estimateDispersion = FALSE)
  RPKM.uq <- rpkm(dge.uq,  gene.length=dge.uq$genes[, "Tx.len"])
  RPKM.uq <- RPKM.uq+1; 
  FC.data <- FC_replicate_samples(RPKM.uq)
  
  for(i in seq(length(FC.data))){
    rep.p <- cor.test(log2.Tx.Len, FC.data[[i]], method="spearman")
    P.rep[5, i] <- rep.p$p.value; 
    Spearman.r.rep[5, i] <- rep.p$estimate; 
  }
  
  
  res <- data.frame("r"=Spearman.r.rep, "p"=P.rep)
  outFile <- sprintf("%s_replicate_samples_analysis_results.txt", DatasetTitle);
  write.table(res, file=outFile, sep="\t", quote=F, col.names=NA)
  
  str <- sprintf("\nResults recorded in %s\n", outFile)
  cat(str, "\n")
  
  #invisible(dev.off())
  
  res
}

#--------------------------------------------------------------
#
# Step 4:
#---------
# function: length_bias_analysis_cqn_normalization()
#
# Examine relationship between FC and Tx length for replicate samples  and biological conditions using the following processing methods:
# 1. CQN Normalization 
#
# The function assumes that the dataset contains 2 biological conditions (control and treatment)
# and each with n replicates, ordered in the input files such that Ens Ids are row names,1st col is Sym, then cols 2:n are for cond1 and (n+1):ncol for cond2.   
#
# Set Rep.first.col=2 to analyze the replicates of "control" samples
# Set Rep.first.col=N.replicates+2 to analyze the replicates of "treatmnet" samples
#


length_bias_analysis_cqn_normalization <- function (countsFile.Pgenes, Hs.Tx.GC.Content.And.LenFile,N.replicates=3, DatasetTitle="TNFa.p65kd", Rep.first.col=5, X.lim=c(8, 15), Y.lim=c(-4,4)) {
  
 # pdfFile <- sprintf("%s_CQN_replicate_samples_analysis_plots.pdf", DatasetTitle);
  n_pairs <- (N.replicates * (N.replicates-1))/2
 # pdf(pdfFile)
  par(mfrow=c(3,n_pairs))
  P.rep <- array(dim=c(3, n_pairs))
  Spearman.r.rep <- array(dim=c(3, n_pairs))
  rep_names <- c()
  for(i in 1:N.replicates){for(j in 1:i){if (i!=j) {rep_names<-c(rep_names,sprintf("rep%dvs%d",i,j))}}}
  colnames(P.rep)          <- rep_names
  colnames(Spearman.r.rep) <- rep_names
  
  rownames(P.rep)          <- c("edgeR.TMM", "CQN_GC", "CQN_GC_Length")
  rownames(Spearman.r.rep) <- c("edgeR.TMM", "CQN_GC", "CQN_GC_Length") 
  
  Rep.last.col  <- Rep.first.col + N.replicates - 1; 
  
  # Load Pgenes count data
  data.cnts.Pgenes <- read.file(countsFile.Pgenes)
  # Load Gene length file
  Hs.genes.len.gc <- read.table(file=Hs.Tx.GC.Content.And.LenFile, sep="\t", header=1, row.names=1, stringsAsFactors=F); 
  Hs.genes.len.gc <- Hs.genes.len.gc[,3:4] 
  # Consider subset of genes with length info
  cGenes <- intersect(rownames(data.cnts.Pgenes), rownames(Hs.genes.len.gc))
  
  # For dataset with N replicates
  groupf <- factor(1:N.replicates);
  
  # Create DGEList list object: it contains the 3 componenets: (1) counts, (2) group = condition set and (3) genes = with gene annots
  dge <- DGEList(counts=data.cnts.Pgenes[cGenes, Rep.first.col:Rep.last.col], group=groupf, genes=Hs.genes.len.gc[cGenes, ]) 
  log2.Tx.Len <- log(dge$genes[, "Tx.len"],2)
  # Design matrix
  design <- model.matrix(~groupf)
  design
  # 
  # edgeR analysis
  #
  
  dge.tmm <- norm_dge(dge,"TMM",design,estimateDispersion = FALSE)
  
  # fit
  fit <- glmFit(dge.tmm, dispersion=0) 
  
  # Tests: comparisons of all replicates to replicate 1
  lrt <- glmLRT(fit, coef=2:N.replicates)
  
  # Results including log2FC are in lrt$table (comparisons to replicate 1)
  FC.data <- list()
  m <- N.replicates-1
  ind <- 1
  for(i in 1:m){ # cor test of all replicates FC vs. rep1 
    rep.p <- cor.test(log2.Tx.Len, lrt$table[,i], method="spearman")
    P.rep[1, i] <- rep.p$p.value; 
    Spearman.r.rep[1, i] <- rep.p$estimate; 
    FC.data[[ind]] <- lrt$table[,i]
    ind <- ind + 1
  }
  if(N.replicates > 2){ # cor test for all pairs of replicates
    for(i in 2:N.replicates){
      for(j in 2:i){
        if (i!=j) {
          x <- FCBetweenReps(fit,N.replicates,i,j)
          rep.p <- cor.test(log2.Tx.Len, x$table[,1], method="spearman")
          P.rep[1,ind] <- rep.p$p.value; 
          Spearman.r.rep[1, ind] <- rep.p$estimate; 
          FC.data[[ind]] <- x$table[,1]
          ind <- ind + 1
        }
      }
    }
  }
  
  edgeR_TMM_plots <- plot_replicate_samples(log2.Tx.Len,FC.data,"edgeR-TMM",DatasetTitle,rep_names,xLab="Length (log2)", X.lim, Y.lim)
  
  colnames(Hs.genes.len.gc) <- c("length","GC.content")
  # Consider subset of genes with length info
  counts <- data.cnts.Pgenes[cGenes,]
  counts.matrix <- as.matrix(counts[,2:ncol(counts)])
  counts <- merge(counts,Hs.genes.len.gc,by="row.names")
  
  cqn_res <- cqn(counts = counts.matrix,x = counts$GC.content,lengths = counts$length,lengthMethod = "fixed") # cqn normalization only GC
  
  CQN_norm <- cqn_res$y + cqn_res$offset # values are in log2 
  offset_res <- cqn_res$offset
  
  # Examine relationship between FC and Tx length for replicate samples after CQN normalization - only GC
  Rep.first.col <- Rep.first.col - 1 
  Rep.last.col <- Rep.first.col + N.replicates -1
  CQN_norm_replicates <- CQN_norm[,Rep.first.col:Rep.last.col]
  
  FC.data <- FC_replicate_samples((2^CQN_norm_replicates) + 1) #add 1 pseudo count 
  CQN.qnorm_plots <- plot_replicate_samples(log(counts$length,2),FC.data,"CQN_GC",DatasetTitle,rep_names,xLab="Length (log2)", X.lim, Y.lim)
  for(i in seq(length(FC.data))){
    rep.p <- cor.test(log2.Tx.Len, FC.data[[i]], method="spearman")
    P.rep[2, i] <- rep.p$p.value; 
    Spearman.r.rep[2, i] <- rep.p$estimate; 
  }
  
  # Examine relationship between FC and Tx length for replicate samples after CQN normalization -  GC and Tx.Len
  
  cqn_res <- cqn(counts = counts.matrix,x = counts$GC.content,lengths = counts$length) # cqn normalization  
  
  CQN_norm <- cqn_res$y + cqn_res$offset # values are in log2 
  offset_res <- cqn_res$offset
 
  # Examine relationship between FC and Tx length for replicate samples after CQN normalization 
  CQN_norm_replicates <- CQN_norm[,Rep.first.col:Rep.last.col]
  
  FC.data <- FC_replicate_samples((2^CQN_norm_replicates) + 1) # add 1 for pseudo count
  plot_replicate_samples(log(counts$length,2),FC.data,"CQN_GC_Length",DatasetTitle,rep_names,xLab="Length (log2)", X.lim, Y.lim)
  #invisible(dev.off())
  
  for(i in seq(length(FC.data))){
    rep.p <- cor.test(log2.Tx.Len, FC.data[[i]], method="spearman")
    P.rep[3, i] <- rep.p$p.value; 
    Spearman.r.rep[3, i] <- rep.p$estimate; 
  }
  
  res_replicates <- data.frame("r"=Spearman.r.rep, "p"=P.rep)
  outFile <- sprintf("%s_cqn_replicate_samples_analysis_results.txt", DatasetTitle);
  write.table(res_replicates, file=outFile, sep="\t", quote=F, col.names=NA)
  str <- sprintf("\nResults recorded in %s\n", outFile)
  cat(str, "\n")
  
  #pdfFile <- sprintf("%s_CQN_biological_conditions_analysis_plots.pdf", DatasetTitle);
  #pdf(pdfFile)
  par(mfrow=c(3,3))
  P.bio <- array(dim=c(1, 3))
  Spearman.r.bio <- array(dim=c(1, 3))
  
  colnames(P.bio)          <- c("edgeR.TMM", "CQN_GC", "CQN_GC_Length")
  colnames(Spearman.r.bio) <- c("edgeR.TMM", "CQN_GC", "CQN_GC_Length")
  
  
  Conds <- list("control"=1:N.replicates, "treatment"=(N.replicates+1):(N.replicates*2))
  # edgeR-TMM analysis
  
  groupf <- factor(c(rep("control",N.replicates), rep("treatment",N.replicates)))  # Load Gene length file
  Hs.genes.len.gc <- read.table(file=Hs.Tx.GC.Content.And.LenFile, sep="\t", header=1, row.names=1, stringsAsFactors=F); 
  Hs.genes.len.gc <- Hs.genes.len.gc[,3:4] 
  dge <- DGEList(counts=data.cnts.Pgenes[cGenes, 2:ncol(data.cnts.Pgenes)], group=groupf, genes=Hs.genes.len.gc[cGenes, ])
  
  # Design matrix
  design <- model.matrix(~groupf)
  
  # TMM norm +Estimate dispersion
  dge.tmm <-norm_dge(dge,"TMM",design,estimateDispersion = TRUE)
  
  # Fit NB model and perform LR test (here group factor has 2 levels, hence the fit has 2 coef. First is the Intercept. 
  lrt <- fitNBandLRtest(dge.tmm,design)
  
  # Results including log2FC are in lrt$table
  head(lrt$table)
  
  edgeR.TMM.logFC <- lrt$table[,1];  
  main.title <- sprintf("%s: edgeR-TMM", DatasetTitle);
  r <- plotFCvsLen(logTxLen=log2.Tx.Len, logFC=edgeR.TMM.logFC, main.title,xLab = "Length (log2)", X.lim, Y.lim)
  P.bio[1, 1]          <- r$p.value; 
  Spearman.r.bio[1, 1] <- r$estimate; 
  
  colnames(Hs.genes.len.gc) <-  c("length","GC.content")
  # Examine relationship between FC and Tx length for biological conditions after CQN normalization with GC 
  cqn_res <- cqn(counts = counts.matrix,x = counts$GC.content,lengths = counts$length, lengthMethod = "fixed") # cqn normalization - fixed length  
  CQN_norm <- cqn_res$y + cqn_res$offset # values are in log2 
  CQN.logFC <- calcLogFC((2^CQN_norm)+1,Conds) #add 1 pseudo count
  main.title <- sprintf("%s: CQN_GC", DatasetTitle);
  r <- plotFCvsLen(logTxLen = log(counts$length,2),logFC = CQN.logFC, main.title = main.title ,xLab = "Length (log2)", X.lim = X.lim, Y.lim = Y.lim)
  P.bio[1, 2]          <- r$p.value; 
  Spearman.r.bio[1, 2] <- r$estimate; 
  
  # Examine relationship between FC and Tx length for biological conditions after CQN normalization with GC and Length
  cqn_res <- cqn(counts = counts.matrix,x = counts$GC.content,lengths = counts$length) # cqn normalization  
  CQN_norm <- cqn_res$y + cqn_res$offset # values are in log2 
  CQN.logFC <- calcLogFC((2^CQN_norm)+1,Conds) #add 1 pseudo count
  main.title <- sprintf("%s: CQN_GC_Length", DatasetTitle);
  r <- plotFCvsLen(logTxLen = log(counts$length,2),logFC = CQN.logFC, main.title = main.title ,xLab = "Length (log2)", X.lim = X.lim, Y.lim = Y.lim)
  #invisible(dev.off())
  
  P.bio[1, 3]          <- r$p.value; 
  Spearman.r.bio[1, 3] <- r$estimate; 
  
  res_bio <- rbind(Spearman.r.bio, P.bio)
  rownames(res_bio) <- c("Spearman.r", "pval")
  
  outFile <- sprintf("%s_cqn_treatment_vs_control_analysis_results.txt", DatasetTitle)
  write.table(res_bio, file=outFile, sep="\t", quote=F, col.names=NA)
  str <- sprintf("\nResults recorded in %s\n", outFile)
  cat(str, "\n")
  cat("cqn treatment vs control analysis results:\n")
  print(res_bio)
  cat("\ncqn replicates results:\n")
  print(res_replicates)
}


# Auxiliary functions:
#
#-----------------------------------------------------
#
# function: read.file
# ==========================
#
# Input: 
#	1. tab delimited file with header line; Row.names in row number 1.  
#	
#

read.file <- function (file.name) {
  read.table(file=file.name, sep="\t", header=1, row.names=1, stringsAsFactors=F)	
}

#-----------------------------------------------------
#
# Input:
#
#      1. expresion matrix (only counts/cpm data, no sym cols etc)
#      2. List defining condition sets
# 
#       Conds <- list("control"=1:N.replicates, "treatment"=(N.replicates+1):(N.replicates*2))
#
# OutPut:
#    
#      1. Geneids with atleast N.cutoff reads in atleast N.samples.cutoff of a certain cond set. 
#

find_genes_above_count_thresh_in_Conds_subset <- function(counts.mat, Conds, N.cutoff=20, N.samples.cutoff=2) {
  
  m <- (counts.mat > N.cutoff)
  
  Cs <- Conds[[1]]
  A <- apply(m[, Cs], 1, sum)
  
  for ( i in 2:length(Conds) ) {
    Cs <- Conds[[i]]
    A <- cbind( A, apply(m[, Cs], 1, sum) )
  }
  
  Is <- which( apply(A, 1, max)>=N.samples.cutoff  )
  
  
  cat("\nNum of genes passing filter: ", length(Is), "\n\n"); 
  
  rownames(counts.mat)[Is]
}

#-----------------------------------------------------
#function: plotFCvsLen()
# Input:
#
#      1. x axis var (logTxLen or GC content)
#      2.logFC values
#      3. plot title
#      4. x axis label
#      5.X.lim
#      6.Y.lim
#      7.cex 
# 
#  OutPut:
#    
#      1. Produces a plot of FC vs X axis var, return spearman cor.test results
#
plotFCvsLen <- function(logTxLen, logFC, main.title,xLab="Length (log2)", X.lim, Y.lim, cex.Main=1.0) {
  
  plot(logTxLen, logFC, pch=19, cex=0.6, main=main.title, xlab=xLab, ylab="FC (log2)", xlim=X.lim, ylim=Y.lim, cex.main=cex.Main)
  lm <- tryCatch(lm( logFC ~  logTxLen),error=function(e) NA) 
  if(!is.na(lm)){
    abline(lm, col="red", lwd=1)  
  }
  r <- cor.test(logTxLen, logFC, method="spearman")
  str <- sprintf("\nr=%.2f\np=%.2e", r$estimate, r$p.value)
  text(x=X.lim[1], y=Y.lim[2], str, adj=0, col="black", cex=0.8)
  
  return(r)
  
}
#-----------------------------------------------------

#function: fitNBandLRtest()
#Fit NB model and perform LR test (here group factor has 2 levels, hence the fit has 2 coef. First is the Intercept.
#
fitNBandLRtest <- function(dge,design){
  
  fit <- glmFit(dge, design)
  lrt <- glmLRT(fit, coef=2) # group2 vs group1
  return(lrt)
}
#-----------------------------------------------------
# calculates FC between two conditions
calcLogFC <- function(exprData,Conds){
  logFC <- log( apply(exprData[, Conds[[2]]], 1, mean) / apply(exprData[, Conds[[1]]], 1, mean), 2)
  return(logFC)
}
#-----------------------------------------------------
RPKM_norm <- function(dge){
  dge.none <- calcNormFactors(dge,  method="none") 
  RPKM <- rpkm(dge.none, gene.length=dge$genes[, "Tx.len"])
  RPKM <- RPKM + 1; # Pseudo count
  return(RPKM)
}
#-----------------------------------------------------
#
# function: norm_dge()
# apply normalization method on DGE list object
# return normalized (with/without estimation of dispersion) DGE object
#
norm_dge <- function(dge,normMethod,design,estimateDispersion=TRUE){
  # normaliztion
  dge.norm <- calcNormFactors(dge,  method=normMethod) 
  dge.norm$samples
  
  # Estimate dispersion
  if(estimateDispersion==TRUE){
    dge.norm <- estimateDisp(dge.norm, design)
  }
  return(dge.norm)
}

#-----------------------------------------------------
#
# function: FC_replicate_samples()
#
# Calculates the fold change between replicate samples. 
# returns list of vectors of the FC values between all pairs of replicate sample 
#
FC_replicate_samples <- function(expr.data){
  res <- list()
  rep_names <-c()
  ind <- 1
  for(i in 1:ncol(expr.data)){
    for(j in 1:i){
      if (i!=j) {
        expr.data.logFC <- log( expr.data[, i]/expr.data[, j], 2)
        res[[ind]] <- expr.data.logFC
        rep_names <- c(rep_names,sprintf("rep%dvs%d",i,j))
        ind = ind +1
      }
    }
  }
  names(res) <- rep_names
  return(res)
}

#------------------------------------------------------------

plot_replicate_samples <- function(log2.Tx.Len,FC.data,normalization.method, DatasetTitle,rep_names, xLab="Length (log2)", X.lim=c(8, 15), Y.lim=c(-4,4)){
  res <- c()
  for(i in seq(length(FC.data))){
    if (i == 10){ plot.new()}
    expr.data.logFC <- FC.data[[i]]
    rep_name <- rep_names[i]
    main.title <- sprintf("%s %s\n%s", DatasetTitle,rep_name, normalization.method);
    r <- plotFCvsLen(logTxLen=log2.Tx.Len, logFC=expr.data.logFC, main.title,xLab, X.lim, Y.lim, cex.Main=1.0)
    p1 <- recordPlot()
    res <- c(res,p1)
  }
  return(res)
}

#------------------------------------------------------------
FCBetweenReps <- function(fit,N.replicates,i,j){
  contrast_vals <- rep_len(0,N.replicates)
  contrast_vals[i] <- 1
  contrast_vals[j] <- -1
  x <- glmLRT(fit, contrast=contrast_vals)
  return(x)
}
#------------------------------------------------------------
