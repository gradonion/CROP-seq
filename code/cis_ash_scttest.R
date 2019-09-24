ash_on_scttest <- function(summary_per_gene,cisgenes.df,option){
  beta.mtx <- sapply(summary_per_gene, function(x){x$beta})
  se.mtx <- sapply(summary_per_gene, function(x){x$se})
  pval.mtx <- sapply(summary_per_gene, function(x){x$pval})
  for (i in 1:nrow(cisgenes.df)) {
    if (option=='enhancer'){
      tmp.column <- cisgenes.df$locus[i]
    } else if (option=='gRNA') {
      tmp.column <- cisgenes.df$gRNA[i]
    } else {
      print('Please provide a valid \'option\' type: enhancer or gRNA.')
      return(1)
    }
    tmp.gene <- cisgenes.df$cisGene[i]
    cisgenes.df$ttest.beta[i] <- beta.mtx[tmp.gene,tmp.column]
    cisgenes.df$ttest.se[i] <- se.mtx[tmp.gene,tmp.column]
    cisgenes.df$ttest.pval[i] <- pval.mtx[tmp.gene,tmp.column]
  }
  cisgenes.df <- na.omit(cisgenes.df)
  cat(nrow(cisgenes.df),'data points:')
  cisgenes.df$zscore <- cisgenes.df$ttest.beta/cisgenes.df$ttest.se
  plot1 <- ggplot(cisgenes.df,aes(zscore)) + geom_histogram(aes(y=..density..),bins = 20) +
    labs(x = 'z-score from t test',y='frequency') + 
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1),col=2)
  
  ash.ttest <- ash(cisgenes.df$ttest.beta,cisgenes.df$ttest.se,mixcompdist = '-uniform')
  cisgenes.df$ash.beta <- ash.ttest$result$PosteriorMean
  cisgenes.df$lfsr <- get_lfsr(ash.ttest)
  # cisgenes.df$lfdr <- get_lfdr(ash.ttest)
  log10_lfsr <- log10(get_lfsr(ash.ttest))
  plot2 <- qplot(cisgenes.df$ttest.beta, cisgenes.df$ash.beta,
                 color=log10_lfsr, main = 'ash shrinkage on t test result',
                 xlab = 'log fold change', ylab = 'beta after shrinkage') + 
    geom_abline(intercept = 0,slope = 1,linetype = "dotted")
  grid.arrange(plot1,plot2,ncol=2)
  cisgenes.df <- cisgenes.df[order(cisgenes.df$lfsr),]
  cisgenes.df <- format(cisgenes.df,digits=3)
  print(kable(cisgenes.df,row.names = F) %>% kable_styling() %>%
          scroll_box(width = "100%", height = "300px"))
  ash.stats <- data.frame(pi=ash.ttest$fitted_g$pi, 
                          interval_a=ash.ttest$fitted_g$a, interval_b=ash.ttest$fitted_g$b)
  cat('ash estimated mixture proportions and corresponding intervals:')
  print(kable(ash.stats[ash.stats$pi!=0,], row.names = F, digits = 3) %>% 
          kable_styling(position = 'center',full_width = T) )
  return(ash.ttest)
}
