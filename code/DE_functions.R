require(edgeR)
require(DESeq2)
require(MAST)

DE_process <- function(gcount, ncount, filtcpm=10, filtpercent=0.2, perm=F){
  # perm true of false, if true will perform a permuted version.
  coldata <- data.frame(row.names = c(colnames(gcount),colnames(ncount)),
                        condition=c(rep('G',dim(gcount)[2]),rep('N',dim(ncount)[2])))
  countall <- cbind(gcount,ncount)
  if (perm==T){
    coldata$condition <- permute(coldata$condition)
  }
  y <- DGEList(counts= countall,group=coldata$condition)
  keep <- rowSums(edgeR::cpm(y)[,y$samples$group=="N"]>filtcpm) >= dim(ncount)[2] * filtpercent
  y <- y[keep, keep.lib.sizes=FALSE]
  return(y)
}

run_edgeR_qlf <- function(y,plotit=T) {
  y <- calcNormFactors(y)
  group= y$samples[,"group"]
  design <- model.matrix(~group)
  y <- estimateDisp(y,design)
  
  fitqlf <- glmQLFit(y,design)
  qlf <- glmQLFTest(fitqlf,coef=2)
  out <- topTags(qlf, n=Inf, adjust.method = "BH")
  
  if (plotit==T) {
    summ_pvalues(qlf$table$PValue)
    outsig <- subset(out$table,FDR <0.1)
    print(paste0("There are ",dim(outsig)[1], " genes passed FDR <0.1 cutoff"))
    print(knitr::kable(signif(as.matrix(head(out$table[order(out$table$PValue),])),digit=2)))
  }
  
  return(out)
}

run_ttest <- function(y) {
  countfl <- y$counts
  group= y$samples[,"group"]
  gcount <- countfl[,group=="G"]
  ncount <- countfl[,group=="N"]
  pv <- rep(1,dim(countfl)[1])
  for (i in 1:dim(countfl)[1]){
    a <- t.test(countfl[i,1:dim(gcount)[2]], y=countfl[i,1:dim(ncount)[2]])
    pv[i] <- a$p.value
  }
  pv
}

# run_deseq2 <- function(y) {
  # coldata <- data.frame(row.names = y$counts,
  #                       condition=y$samples[,"group"])
  # countfl <- y$counts
run_deseq2 <- function(countfl,coldata) {
  dds = DESeqDataSetFromMatrix(countData = countfl, colData = coldata,design = ~condition)
  dds = estimateSizeFactors(dds)
  ddsWARD = DESeq(dds)
  resWARD = results(ddsWARD)
  return(resWARD)
}

run_MASTcpmDetRate <- function(y) {
  countfl <- y$counts
  grp <-y$samples[,"group"]
  names(grp) <- rownames(y$samples)
  cdr <- scale(colSums(countfl > 0)/dim(countfl)[1])
  dge <- DGEList(counts = countfl)
  dge <- edgeR::calcNormFactors(dge)
  cpms <- edgeR::cpm(dge)
  sca <- FromMatrix(exprsArray = log2(cpms + 1),
                    cData = data.frame(wellKey = names(grp),
                                       grp = grp, cdr = cdr))
  zlmdata <- zlm.SingleCellAssay(~cdr + grp, sca)
  mast <- lrTest(zlmdata, "grp")
  
  # hist(mast[, "hurdle", "Pr(>Chisq)"], 50)
  df = data.frame(pval = mast[, "hurdle", "Pr(>Chisq)"],
                  row.names = names(mast[, "hurdle", "Pr(>Chisq)"]))
  df$fdr <- p.adjust(df$pval, method="BH")
  b <- getLogFC(zlmdata)
  df <- cbind(df, b[b$contrast=="grpN",c("logFC", "varLogFC","z"), with=F])
  df <- df[order(df$pval),]
  return(df)
}