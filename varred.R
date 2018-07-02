var.red <- function(distmatrix,df_vars,npermut=1000,plevel=0.05,pcrit=pcrit_rda)
{ require(vegan)
  if (plevel>0.05) plevel <- 0.05
  #if (pcrit >0.1 ) pcrit  <- 0.1
  rows <- rownames(df_vars)
  df_vars <- df_vars[,colSums(df_vars)!=0]
  spe.rda.1 <- dbrda(distmatrix ~ 1, data = df_vars, add=T)
  spe.rda <- dbrda(distmatrix ~ ., data = df_vars, add=T)
  a <- anova(spe.rda)
  pval <- a[1,4]
  cat("\ndbrda pvalue preciding ordiR2step: ",pval,"\n")
  if (pval > pcrit) 
    { cat(paste("warning: dbrda not significant, pvalue > ",pcrit,", no vars added\n",sep="")) 
      df_vars <- as.data.frame(df_vars[,0])
      rownames(df_vars) <- rows
      return(df_vars)
    }
  if (pval > plevel) 
    cat(paste("warning: dbrda significance anova: ",plevel," < pvalue <= ",pcrit,"\n",sep=""))
  # if the analysis is "significant", compute the adjusted R2
  spe.fwd <- ordiR2step(spe.rda.1, scope = formula(spe.rda),trace=F, Pin=plevel, 
             R2permutations=npermut, permutations=how(nperm=npermut))
  #significant vars
  fwd_res <- spe.fwd$terminfo$terms
  fwd_chr <- as.character(fwd_res)
  fwd_chr <- fwd_chr[2]
  fwd_chr.split <- unlist(strsplit(fwd_chr,"[+ ]"))
  fwd_chr.split <- subset(fwd_chr.split,fwd_chr.split!="")
  cat("sign. vars   following ordiR2step:\n")
  print(fwd_chr.split)
  indices <- match(fwd_chr.split,colnames(df_vars),nomatch=0)
  df_vars <- as.data.frame(df_vars[,indices])
  if (ncol(df_vars)==1) 
  {colnames(df_vars)[1] <- fwd_chr.split; rownames(df_vars) <- rows }
  spe.rda <- dbrda(distmatrix ~ ., data = df_vars, add=T)
  a <- anova(spe.rda)
  pval <- a[1,4]
  cat("dbrda pvalue following ordiR2step: ",pval,"\n")
  return(df_vars)
}
