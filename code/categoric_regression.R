categoric <- function(mtx,tg.cells,neg.cells){
  pval <- rep(NA,nrow(mtx))
  beta <- rep(NA,nrow(mtx))
  se <- rep(NA,nrow(mtx))
  subset <- mtx[,c(tg.cells,neg.cells)]
  for (i in 1:nrow(subset)){
    df <- data.frame(row.names = c(tg.cells,neg.cells),
                     residual = as.numeric(subset[i,]),
                     condition=c(rep(1,length(tg.cells)),
                                 rep(0,length(neg.cells))))
    cat.lm <- summary(lm(residual~factor(condition), data = df))
    beta[i] <- cat.lm$coefficients[2,1]
    se[i] <- cat.lm$coefficients[2,2]
    pval[i] <- cat.lm$coefficients[2,4]
  }
  names(beta) <- rownames(mtx)
  names(se) <- rownames(mtx)
  names(pval) <- rownames(mtx)
  return(list(beta=beta,se=se,pval=pval))
}