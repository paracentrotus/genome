require(RColorBrewer)
require(pheatmap)
require(tidyverse)
require(wesanderson)

plotSynt <- function(synt,labx,laby){
  synt  %>% select(s2gi,s2chr) %>%  group_by(s2chr) %>%  summarise(max=max(s2gi))%>%
    arrange(max) %>% mutate(mids=lag(max,default=0)+0.5*(max-lag(max,default=0)))-> s2chLim
  synt  %>% select(s1gi,s1chr) %>%  group_by(s1chr) %>%  summarise(max=max(s1gi))%>%
    arrange(max) %>% mutate(mids=lag(max,default=0)+0.5*(max-lag(max,default=0)))-> s1chLim
  #palette(colorRampPalette(brewer.pal(6,'BuPu'))(nrow(s2chLim)))
  #palette(colorRampPalette(brewer.pal(12,'Paired'))(nrow(s1chLim)))
  plab=paste("(",nrow(synt)," orthologues)",sep="")
  labxp=paste(labx,plab,sep=" ")
  #plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=synt$s2chr,xlab=labxp,ylab=laby)
  plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=synt$s1chr,xlab=labx,ylab=laby,xlim=c(0,max(s1chLim$max)))
  #plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=ifelse(synt$scol=='S',as.factor(synt$s1chr),'lightgrey'),xlab=labx,ylab=laby,xlim=c(0,max(s1chLim$max)))
  plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=ifelse(synt$scol=='S','#810f7c','lightgrey'),xlab=labx,ylab=laby,xlim=c(0,max(s1chLim$max)))
  
  abline(v=0,lty=3,lwd=0.5)
  abline(h=0,lty=3,lwd=0.5)
  
  abline(v=s1chLim$max,lty=3,lwd=0.5)
  abline(h=s2chLim$max,lty=3,lwd=0.5)
  
  axis(3,at=s1chLim$mids,labels=s1chLim$s1chr,las=2,cex.axis=0.6)
  axis(4,at=s2chLim$mids,labels=s2chLim$s2chr,las=1,cex.axis=0.6)
  
}

plotSyntr <- function(synt,labx,laby){
  synt  %>% select(s2gi,s2chr) %>%  group_by(s2chr) %>%  summarise(max=max(s2gi)) %>%
    arrange(max) %>% mutate(mids=lag(max,default=0)+0.5*(max-lag(max,default=0)))-> s2chLim
  synt  %>% select(s1gi,s1chr) %>%  group_by(s1chr) %>%  summarise(max=max(s1gi)) %>%
    arrange(max) %>% mutate(mids=lag(max,default=0)+0.5*(max-lag(max,default=0)))-> s1chLim
  print(nrow(s2chLim))
  #palette(colorRampPalette(brewer.pal(12,'Paired'))(nrow(s2chLim)))
  palette(colorRampPalette(brewer.pal(12,'Paired'))(nrow(s2chLim)))
  #plab=paste("(",nrow(synt)," orthologues)",sep="")
  #labxp=paste(labx,plab,sep=" ")
  #plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=synt$s2chr,xlab=labxp,ylab=laby)
  #plot(synt$s2gi,synt$s1gi,pch=20,cex=0.4,col=as.factor(synt$s2chr),xlab=labx,ylab=laby)
  plot(synt$s2gi,synt$s1gi,pch=20,cex=0.4,col=ifelse(synt$scol=='S',as.factor(synt$s2chr),'lightgrey'),xlab=labx,ylab=laby)

  abline(v=0,lty=3,lwd=0.5)
  abline(h=0,lty=3,lwd=0.5)
  
  abline(h=s1chLim$max,lty=3,lwd=0.5)
  abline(v=s2chLim$max,lty=3,lwd=0.5)
  
  axis(4,at=s1chLim$mids,labels=s1chLim$s1chr,las=1,cex.axis=0.6)
  axis(3,at=s2chLim$mids,labels=s2chLim$s2chr,las=2,cex.axis=0.6)
  
}

plotSyntCLG <- function(synt,labx,laby,c){
  synt  %>% select(s2gi,s2chr) %>%  group_by(s2chr) %>%  summarise(max=max(s2gi))%>%
    arrange(max) %>% mutate(mids=lag(max,default=0)+0.5*(max-lag(max,default=0)))-> s2chLim
  synt  %>% select(s1gi,s1chr) %>%  group_by(s1chr) %>%  summarise(max=max(s1gi))%>%
    arrange(max) %>% mutate(mids=lag(max,default=0)+0.5*(max-lag(max,default=0)))-> s1chLim
  #palette(c(brewer.pal(10,'Paired'),brewer.pal(7,'Dark2')))
  #palette(colorRampPalette(brewer.pal(12,'Paired'))(length(levels(synt$clg))))
  #palette(colorRampPalette(brewer.pal(12,'Paired'))(nrow(s1chLim)))
  
  plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=synt$clg,xlab=labx,ylab=laby)
  #plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=synt$s1chr,xlab=labx,ylab=laby)
  
  abline(v=0,lty=3,lwd=0.5)
  abline(h=0,lty=3,lwd=0.5)
  
  abline(v=s1chLim$max,lty=3,lwd=0.5)
  abline(h=s2chLim$max,lty=3,lwd=0.5)
  
  axis(3,at=s1chLim$mids,labels=s1chLim$s1chr,las=2,cex.axis=0.8)
  axis(4,at=s2chLim$mids,labels=s2chLim$s2chr,las=1,cex.axis=0.4)
  
}

plotSyntCLGoc <- function(synt,labx,laby,c){
  synt  %>% select(s2gi,s2chr) %>%  group_by(s2chr) %>%  summarise(max=max(s2gi))%>%
    arrange(max) %>% mutate(mids=lag(max,default=0)+0.5*(max-lag(max,default=0)))-> s2chLim
  synt  %>% select(s1gi,s1chr) %>%  group_by(s1chr) %>%  summarise(max=max(s1gi))%>%
    arrange(max) %>% mutate(mids=lag(max,default=0)+0.5*(max-lag(max,default=0)))-> s1chLim
  #palette(c(brewer.pal(10,'Paired'),brewer.pal(7,'Dark2')))
  #palette(colorRampPalette(brewer.pal(12,'Paired'))(length(levels(synt$clg))))
  #palette(colorRampPalette(brewer.pal(12,'Paired'))(nrow(s1chLim)))
  
  plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=ifelse(synt$scol=='S',synt$clg,'lightgrey'),xlab=labx,ylab=laby)
  #plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=synt$s1chr,xlab=labx,ylab=laby)
  
  abline(v=0,lty=3,lwd=0.5)
  abline(h=0,lty=3,lwd=0.5)
  
  abline(v=s1chLim$max,lty=3,lwd=0.5)
  abline(h=s2chLim$max,lty=3,lwd=0.5)
  
  axis(3,at=s1chLim$mids,labels=s1chLim$s1chr,las=2,cex.axis=0.7)
  axis(4,at=s2chLim$mids,labels=s2chLim$s2chr,las=1,cex.axis=0.7)
  
}

plotSyntCLGr <- function(synt,labx,laby){
  synt  %>% select(s2gi,s2chr) %>%  group_by(s2chr) %>%  summarise(max=max(s2gi))%>%
    arrange(max) %>% mutate(mids=lag(max,default=0)+0.5*(max-lag(max,default=0)))-> s2chLim
  synt  %>% select(s1gi,s1chr) %>%  group_by(s1chr) %>%  summarise(max=max(s1gi))%>%
    arrange(max) %>% mutate(mids=lag(max,default=0)+0.5*(max-lag(max,default=0)))-> s1chLim
  palette(c(brewer.pal(10,'Paired'),brewer.pal(7,'Dark2')))
    #palette(colorRampPalette(brewer.pal(12,'Paired'))(length(levels(synt$clg))))
  #palette(colorRampPalette(brewer.pal(12,'Paired'))(nrow(s1chLim)))
  
  #plot(synt$s1gi,synt$s2gi,pch=20,cex=0.4,col=synt$clg,xlab=labx,ylab=laby)
  #plot(synt$s2gi,synt$s1gi,pch=20,cex=0.4,col=synt$clg,xlab=labx,ylab=laby)
  plot(synt$s2gi,synt$s1gi,pch=20,cex=0.4,col=ifelse(synt$scol=='S',synt$clg,'lightgrey'),xlab=labx,ylab=laby)
  
  abline(v=0,lty=3,lwd=0.5)
  abline(h=0,lty=3,lwd=0.5)
  
  abline(h=s1chLim$max,lty=3,lwd=0.5)
  abline(v=s2chLim$max,lty=3,lwd=0.5)
  
  axis(4,at=s1chLim$mids,labels=s1chLim$s1chr,las=1,cex.axis=0.6)
  axis(3,at=s2chLim$mids,labels=s2chLim$s2chr,las=2,cex.axis=0.6)
  
}

ordSynt <- function(synt){
  synt %>% count(s1chr,s2chr) %>% spread(s2chr,n,fill=0) %>% column_to_rownames(var="s1chr")-> chr_mat
  hmp<-pheatmap(as.matrix(chr_mat))
  y_chrom<-rownames(as.matrix(chr_mat))[hmp$tree_row$order]
  x_chrom<-colnames(as.matrix(chr_mat))[hmp$tree_col$order]
  #print(x_chrom)
  synt %>% mutate(s1chrO=match(s1chr,y_chrom)) %>% 
    arrange(s1chrO,s1gp) %>% mutate(s1gi=row_number(s1chrO)) %>%
    mutate(s2chrO=match(s2chr,x_chrom)) %>% 
    arrange(s2chrO,s2gp) %>%  mutate(s2gi=row_number(s2chrO)) -> synt.ord
  return(synt.ord)
  #return(c(x_chrom,y_chrom))
}

ordSyntJO <- function(synt){
  synt %>% count(s1chr,s2chr) %>% spread(s2chr,n,fill=0) %>% column_to_rownames(var="s1chr")-> chr_mat
  hmp<-pheatmap(as.matrix(chr_mat))
  y_chrom<-rownames(as.matrix(chr_mat))[hmp$tree_row$order]
  x_chrom<-colnames(as.matrix(chr_mat))[hmp$tree_col$order]
  #print(x_chrom)
  synt %>% mutate(s1chrO=match(s1chr,y_chrom)) %>% 
    arrange(s1chrO,s1gp) %>% mutate(s1gi=row_number(s1chrO)) %>%
    mutate(s2chrO=match(s2chr,x_chrom)) %>% 
    arrange(s2chrO,s2gp) %>%  mutate(s2gi=row_number(s2chrO)) -> synt.ord
  #return(synt.ord)
  return(c(x_chrom,y_chrom))
}

ordSyntMan <- function(synt,x_chrom,y_chrom){
    synt %>% mutate(s1chrO=match(s1chr,y_chrom)) %>% 
      arrange(s1chrO,s1gp) %>% mutate(s1gi=row_number(s1chrO)) %>%
      mutate(s2chrO=match(s2chr,x_chrom)) %>% 
      arrange(s2chrO,s2gp) %>%  mutate(s2gi=row_number(s2chrO)) -> synt.ord
    return(synt.ord)
}

#### TRY TO MAKE SOME FUNCTIONS
testEnrichCLG <- function(synt){
  synt %>% 
    group_by(s1chr) %>% mutate(s1tot=n()) %>%
    group_by(s2chr) %>% mutate(s2tot=n())%>%
    group_by(s1chr,s2chr,s1tot,s2tot) %>% tally() -> synt.chrn
  synt.no=sum(synt.chrn$n)
  synt.chrn %>% rowwise() %>% 
    mutate(fishp=fisher.test(matrix(c(n,s1tot-n,s2tot-n,synt.no),ncol=2),alternative='greater')$p.value) %>%
    ungroup() %>% mutate(padj=p.adjust(fishp,method="bonferroni")) %>% mutate(scol=ifelse(padj<.05,'S','NS')) -> synt.test
  return(synt.test)
}

testEnrich <- function(synt){
  synt %>% 
    group_by(s1chr) %>% mutate(s1tot=n()) %>%
    group_by(s2chr) %>% mutate(s2tot=n())%>%
    group_by(s1chr,s2chr,s1tot,s2tot) %>% tally() -> synt.chrn
  synt.no=sum(synt.chrn$n)
  synt.chrn %>% rowwise() %>% 
    mutate(fishp=fisher.test(matrix(c(n,s1tot-n,s2tot-n,synt.no),ncol=2),alternative='greater')$p.value) %>%
    ungroup() %>% mutate(padj=p.adjust(fishp,method="bonferroni")) %>% mutate(scol=ifelse(padj<.05,'S','NS')) -> synt.test
  return(synt.test)
} 

gchkExp <- function(synt.exp){
  synt.exp %>% arrange(s1chr,s1gi) %>% #filter(s1chr=="Scaffold_2100") %>% 
    group_by(s1chr) %>% mutate(chunk=cut_width(s1gi,width=20)) %>% ungroup() %>% filter(scol=="S") %>% 
    select(s2chr,chunk,s1chr,s1gp) %>% 
    group_by(s1chr,chunk) %>% mutate(xmin=min(s1gp),xmax=max(s1gp)) %>%
    group_by(s2chr,s1chr,chunk,xmin,xmax) %>% tally() %>% 
    group_by(chunk) %>% mutate(ymin=(cumsum(n)-n)/sum(n),ymax=cumsum(n)/(sum(n))) -> synt.gchk
  return(synt.gchk)
}

palette(c(brewer.pal(12,'Paired'),brewer.pal(7,'Dark2'),brewer.pal(3,'Set3')))
gcol=c(brewer.pal(12,'Paired'),brewer.pal(7,'Dark2'),brewer.pal(3,'Set3'))
palette(c(brewer.pal(12,'Paired'),brewer.pal(7,'Dark2')))
gcol=c(brewer.pal(12,'Paired'),brewer.pal(7,'Dark2'))
c(brewer.pal(12,'Paired'),brewer.pal(7,'Dark2'))

bflo_col=c("#1F78B4","#1B9E77","#33A02C","#B2DF8A","#FB9A99","#E31A1C","#FFFF99","#B15928","#D95F02","#7570B3","#E7298A","#66A61E","#CAB2D6","#6A3D9A","#E6AB02","#A6761D","#8DD3C7","#666666","#FF7F00","#FDBF6F","#BEBADA","#ffd92f")
palette(bflo_col)

pmax_col=c("#A6CEE3","#1F78B4","#1B9E77","#33A02C","#B2DF8A","#FB9A99","#FFFF99","#B15928","#D95F02","#7570B3","#E7298A","#66A61E","#6A3D9A","#E6AB02","#A6761D","#666666","#FF7F00","#BEBADA")
palette(pmax_col)

par(resetPar())     ## reset
par(mar=c(5,5,6,6))

gcol=c(brewer.pal(10,'Paired'),brewer.pal(7,'Dark2'))
setwd("~/Dropbox/Genomes/Synteny/synt-forge/Urchin/")
#
hsa_mmu_synt=read.table("Hsa_Mmu_synt.txt",h=T,stringsAsFactors=T)
hsa_mmu_synt.test=testEnrich(hsa_mmu_synt)
hsa_mmu_synt %>% inner_join(hsa_mmu_synt.test) -> hsa_mmu_synt.exp

plotSynt(hsa_mmu_synt.exp,"Homo","Mus")
#amphioxus
pliv_bflo_synt=read.table("Pliv_Bflo_synt.txt",h=T,stringsAsFactors=T)
plotSyntCLG(pliv_bflo_synt,'Paracentrotus','Branchiostoma')
pliv_bflo_synt.test=testEnrich(pliv_bflo_synt)
pliv_bflo_synt %>% inner_join(pliv_bflo_synt.test) -> pliv_bflo_synt.exp
x_chrom=c('BFL_1','BFL_5','BFL_16','BFL_3','BFL_6','BFL_13','BFL_11','BFL_2','BFL_7','BFL_12','BFL_19','BFL_4','BFL_8','BFL_15','BFL_9','BFL_14','BFL_17','BFL_18','BFL_10')
y_chrom=c('chrom1','chrom2','chrom4','chrom18','chrom7','chrom8','chrom11','chrom9','chrom12','chrom15','chrom3','chrom14','chrom17','chrom6','chrom10','chrom16','chrom5','chrom13')
#y_chrom=c('Scaffold_174','Scaffold_2100','Scaffold_3433','Scaffold_3430','Scaffold_3434','Scaffold_3425','Scaffold_3429','Scaffold_3428','Scaffold_3426','Scaffold_3432','Scaffold_649','Scaffold_1721','Scaffold_3435','Scaffold_218','Scaffold_674','Scaffold_3431','Scaffold_2974','Scaffold_641')
plotSyntCLGoc(ordSynt(pliv_bflo_synt.exp),'Paracentrotus','Branchiostoma')
plotSyntCLGoc(ordSyntMan(pliv_bflo_synt.exp,x_chrom,y_chrom),'Paracentrotus','Branchiostoma')
plotSyntCLG(ordSyntMan(pliv_bflo_synt.exp,x_chrom,y_chrom),'Paracentrotus','Branchiostoma')

legend('bottomright',levels(spur_bflo_synt.exp$clg),fill=bflo_col)
levels(pliv_bflo_synt.exp$clg)
cols=data.frame(row.names=levels(pliv_bflo_synt.exp$clg),col=bflo_col)

#pecten
pliv_pmax_synt=read.table("Pliv_Pmax_synt.txt",h=T,stringsAsFactors=T)
pliv_pmax_synt.test=testEnrich(pliv_pmax_synt)
pliv_pmax_synt %>% inner_join(pliv_pmax_synt.test) -> pliv_pmax_synt.exp
plotSyntCLGoc(pliv_pmax_synt.exp,'Paracentrotus','Pecten')
plotSyntCLGoc(ordSynt(pliv_pmax_synt.exp),'Paracentrotus','Pecten')
x_chrom=c('scaffold_4','scaffold_14','scaffold_2','scaffold_7','scaffold_13','scaffold_15','scaffold_18','scaffold_3','scaffold_10','scaffold_11','scaffold_17','scaffold_16','scaffold_1','scaffold_19','scaffold_5','scaffold_6','scaffold_8','scaffold_9','scaffold_12')
#y_chrom=c('Scaffold_174','Scaffold_2100','Scaffold_649','Scaffold_3430','Scaffold_3433','Scaffold_3434','Scaffold_3428','Scaffold_3425','Scaffold_1721','Scaffold_3429','Scaffold_3432','Scaffold_218','Scaffold_3426','Scaffold_3435','Scaffold_674','Scaffold_2974','Scaffold_3431','Scaffold_641')
y_chrom=c('chrom1','chrom2','chrom18','chrom13','chrom4','chrom7','chrom8','chrom12','chrom11','chrom6','chrom9','chrom15','chrom10','chrom3','chrom16','chrom5','chrom17','chrom14')
plotSyntCLG(ordSyntMan(pliv_pmax_synt.exp,x_chrom,y_chrom),'Paracentrotus','Pecten')

plotSyntCLGoc(ordSyntMan(pliv_pmax_synt.exp,x_chrom,y_chrom),'Paracentrotus','Pecten')
legend('bottomright',levels(pliv_pmax_synt.exp$clg),fill=pmax_col)

# Spurp
pliv_spur_synt=read.table("Pliv_Spur_synt.txt",h=T,stringsAsFactors=T)
pliv_spur_synt.test=testEnrich(pliv_spur_synt)
pliv_spur_synt %>% inner_join(pliv_spur_synt.test) -> pliv_spur_synt.exp
plotSynt(ordSynt(pliv_spur_synt.exp),'Paracentrotus','Strongylocentrotus')
write.table(pliv_spur_synt.exp,file='Pliv_Spur_synt_ext.txt',quote=F,sep='\t')
ordSyntJO(pliv_spur_synt.exp)
# Lvar
pliv_lvar_synt=read.table("Pliv_Lvar_synt_a.txt",h=T,stringsAsFactors=T)
pliv_lvar_synt.test=testEnrich(pliv_lvar_synt)
pliv_lvar_synt %>% inner_join(pliv_lvar_synt.test) -> pliv_lvar_synt.exp
plotSynt(ordSynt(pliv_lvar_synt.exp),'Paracentrotus','Lytechinus')
plotSynt((pliv_lvar_synt.exp),'Paracentrotus','Lytechinus')

pliv_pmax_synt.exp
pliv_lvar_synt.exp %>% select(-clg) %>% mutate(clg=pliv_pmax_synt.exp[match(pliv_lvar_synt.exp$s1g,pliv_pmax_synt.exp$s1g),]$clg) -> pliv_lvar_synt.ext

write.table(pliv_lvar_synt.ext,file='Pliv_Lvar_synt_extpm.txt',quote=F,sep='\t')

pliv_nvec_synt=read.table("Pliv_Nvec_synt.txt",h=T,stringsAsFactors=T)
pliv_nvec_synt.test=testEnrich(pliv_nvec_synt)
pliv_nvec_synt %>% inner_join(pliv_nvec_synt.test) -> pliv_nvec_synt.exp
plotSynt(ordSynt(pliv_nvec_synt.exp),'Paracentrotus','Nematostella')


# Spur-Lvar
spur_lvar_synt=read.table("Spur_Lvar_synt.txt",h=T,stringsAsFactors=T)
spur_lvar_synt.test=testEnrich(spur_lvar_synt)
spur_lvar_synt %>% inner_join(spur_lvar_synt.test) -> spur_lvar_synt.exp
plotSynt(ordSynt(spur_lvar_synt.exp),'Strongylocentrotus','Lytechinus')
#plotSynt((pliv_lvar_synt.exp),'Paracentrotus','Lytechinus')
write.table(spur_lvar_synt.exp,file='Spur_Lvar_synt_ext.txt',quote=F,sep='\t')

# Spur-Blo
spur_bflo_synt=read.table("Spur_Bflo_synt.txt",h=T,stringsAsFactors=T)
#plotSyntCLG(pliv_bflo_synt,'Paracentrotus','Branchiostoma')
spur_bflo_synt.test=testEnrich(spur_bflo_synt)
spur_bflo_synt %>% inner_join(spur_bflo_synt.test) -> spur_bflo_synt.exp
plotSyntCLGoc(ordSynt(spur_bflo_synt.exp),'Strongylocentrotus','Branchiostoma')
levels(as.factor(spur_bflo_synt.exp$s1chr))
spur_bflo_synt.exp
spur_bflo_synt[match(spur_lvar_synt.exp$s1g,spur_bflo_synt$s1g),]$clg

spur_lvar_synt.exp %>% mutate(clg=spur_bflo_synt[match(spur_lvar_synt.exp$s1g,spur_bflo_synt$s1g),]$clg) -> spur_lvar_synt.ext
write.table(spur_lvar_synt.ext,file='Spur_Lvar_synt_ext.txt',quote=F,sep='\t')

# Lvar-Blo
lvar_bflo_synt=read.table("Lvar_Bflo_synt.txt",h=T,stringsAsFactors=T)
#plotSyntCLG(pliv_bflo_synt,'Paracentrotus','Branchiostoma')
lvar_bflo_synt.test=testEnrich(lvar_bflo_synt)
lvar_bflo_synt %>% inner_join(lvar_bflo_synt.test) -> lvar_bflo_synt.exp
plotSyntCLGoc(ordSynt(lvar_bflo_synt.exp),'Lytechinus','Branchiostoma')

lvar_pmax_synt=read.table("Lvar_Pmax_synt.txt",h=T,stringsAsFactors=T)
lvar_pmax_synt.test=testEnrich(lvar_pmax_synt)
lvar_pmax_synt %>% inner_join(lvar_pmax_synt.test) -> lvar_pmax_synt.exp
plotSyntCLGoc(ordSynt(lvar_pmax_synt.exp),'Lytechinus','Pecten')

bflo_pmax_synt=read.table("Bflo_Pmax_synt.txt",h=T,stringsAsFactors=T)
bflo_pmax_synt.test=testEnrich(bflo_pmax_synt)
bflo_pmax_synt %>% inner_join(bflo_pmax_synt.test) -> bflo_pmax_synt.exp
plotSyntCLGoc(ordSynt(bflo_pmax_synt.exp),'Brachiostoma','Pecten')

