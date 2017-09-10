library(RColorBrewer)
library(ggplot2)
library(data.table)
library(dtplyr)
library(dplyr)
library(DESeq2)

"%.%" <- function(f,g) function(x) f(g(x))
rowScale <- t %.% scale %.% t

dftabler <- function(x,rn="ens_gene") data.table(x,keep.rownames=rn)

# countfile: from featureCounts (column names must be sample-names)
# tablefile: 1st column=sample-names, 2nd=condition, 3rd=time, ...
getdds <- function(countfile,tbl,model=~time+cond,reduced=NULL,...)
{
    mat <- read.table(countfile,row.names=1,header=TRUE)[,rownames(tbl)]
    if(is.null(tbl$rep)) tbl$rep <- "Rep1"
    if(is.null(tbl$clone)) tbl$clone <- with(tbl,paste(cond,rep,sep="_"))
    rownames(mat) <- sub("\\..*$","",rownames(mat))
    dds <- DESeqDataSetFromMatrix(mat,tbl,model)
    if(is.null(reduced)) dds <- DESeq(dds,test="Wald")
    else dds <- DESeq(dds,test="LRT",reduced=reduced,...)
    return(dds)
}

# Preserve the orders of the factors
getdef <- function(deftable)
  read.table(deftable,row.names=NULL,head=TRUE) %>%
  lapply(function(x) factor(x,unique(x))) %>%
  data.frame(row.names=1)

timecoursePlots <- function(res,rld,k=4,alpha=0.1,showsd=TRUE,bpal="Set2")
{
  tbl <- colData(rld)
  pal <- brewer.pal(k,bpal)
  # Find time-course expression patterns
  A <- rowScale(assay(rld)[na.omit(rownames(res)[res$padj<alpha]),])
  km <- kmeans(A,k,algorithm="Lloyd",nstart=200,iter.max=1000)
  cl <- factor(paste0("c",km$cluster))
  names(cl) <- names(km$cluster)
  m <- reshape2::melt(A)
  m <- data.frame(m,tbl[m$Var2,],cluster=cl[m$Var1])
  m <- m %>% group_by(cond,time,rep,clone,cluster) %>%
    summarise(size=n(),mean=mean(value),sd=sd(value))
  pd <- position_dodge(0.2)
  clabel <- sprintf("%s (%s)", levels(cl),prettyNum(table(cl),","))
  names(clabel) <- levels(cl)
  p <- ggplot(m,aes(time,mean,shape=cond,colour=cluster,group=clone,linetype=cond)) +
    geom_point(size=3,position=pd) + geom_line() +
    geom_hline(yintercept=0,linetype=2,colour="grey") +
    scale_color_brewer(palette=bpal) + theme_bw() +
    xlab("Time") + ylab("Relative expression level") +
    facet_wrap(~cluster, labeller=labeller(cluster=clabel))
  if(showsd) p <- p + geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=.4,position=pd)
  print(p)
  # LDA: inspcet clusters
  ld <- MASS::lda(A,km$cluster)
  op <- par(mfrow=c(1,2))
  plot(A%*%ld$scaling[,1:2],pch=16,col=paste0(pal[km$cluster],"66"),asp=1) 
  abline(h=0,v=0,lty=2,col="grey")
  #legend("topleft",legend=levels(cl),col=pal,pch=16)
  barplot(ld$svd,names.arg=colnames(ld$scaling),las=2,ylab="Singular value")
  par(op)
  ldsc <- reshape2::melt(ld$scaling)
  colnames(ldsc) <- c('sample','LD','loading')
  ldsc <- data.frame(tbl[as.character(ldsc$sample),],ldsc)
  p <- ggplot(ldsc,aes(time,loading,shape=cond,linetype=cond,group=clone,colour=LD)) +
    geom_line() + geom_point(size=3) +
    geom_hline(yintercept=0,linetype=2,colour="grey") +
    facet_wrap(~LD) + theme_bw()
  print(p)
  data.table(ens_gene=rownames(res),tbl_dt(res)) %>%
    filter(padj < alpha) %>%
    mutate(cluster=cl[ens_gene])
}

#if(!exists("e2g")){
#  library(biomaRt)
#  ensembl <- useMart("ENSEMBL_MART_ENSEMBL",host="asia.ensembl.org")
#  mart <- useDataset("mmusculus_gene_ensembl",mart=ensembl)
#  e2g <- getBM(attributes = c("ensembl_gene_id","external_gene_name","gene_biotype"), mart = mart)
#  e2g <- dplyr::rename(e2g, ens_gene = ensembl_gene_id, ext_gene = external_gene_name, biotype = gene_biotype)
#}

#matplot(cbind(x,x-0.3,x-0.3+c(0,-0.2,-0.4)),type='b',pch=1,axes=FALSE,lwd=2,lty=c(1,2,1),col=c(1,2,2),ylab="");arrows(c(1,2,3),c(x[1],x[2]-0.3,x[3]-0.3),c(1,2,3),c(x[1]-0.3,x[2]-0.5,x[3]-0.7),lwd = 2,length = .2,col="blue");abline(h=x[1],lty=2)