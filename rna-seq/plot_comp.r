
require(tidyverse)
require(pheatmap)
require(RColorBrewer)

setwd("~/Dropbox/Genomes/Urchin/Paracentrotus_genome/Analyses/Transcriptomes/mfuzz")
setwd("~/Downloads/")



#elts=read.table('blan_drer_phyper.tsv',h=T)
elts=read.table('Revision/pliv_spur_phyper.tsv',h=T)
elts=read.table('Revision/pliv_blan_phyper.tsv',h=T)

#elts %>% select(cl1,cl2,pv) %>% pivot_wider(names_from = cl1,id_cols=cl2,values_from=pv)%>%
#  na_if(0) %>% column_to_rownames(var='cl2') -> elts.m

#diag(elts.m)<-NA

pal= colorRampPalette(brewer.pal(8, "YlGn"))(100)
pal2= colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)
#pal= colorRampPalette(brewer.pal(8, "YlOrBr"))(100)

elts %>% mutate(paj=p.adjust(pv,method='BH')) %>% dplyr::select(cl1,cl2,paj) %>%
  pivot_wider(names_from = cl1,id_cols=cl2,values_from=paj)%>%
  column_to_rownames(var='cl2') -> elts.m

elts %>% mutate(paj=p.adjust(pv,method='BH')) %>% dplyr::select(cl1,cl2,paj) %>%
  group_by(cl1) %>% mutate(best=max(-log10(paj)),mv=ifelse(-log10(paj)==best,'*','')) %>%
  pivot_wider(names_from = cl1,id_cols=cl2,values_from=mv) %>%
  column_to_rownames(var='cl2') -> elts.l

plsp.m<-elts.m
plsp.l<-elts.l
plbl.m<-elts.m
plbl.l<-elts.l

table(-log10(elts.m)<0.01)    
#elts.d=dist(-log10(elts.m),method='ward.D2')
#ct=hclust(dist(-log10(elts.m)),method='ward.D2')     
plsp_ht<-pheatmap(-log10(plsp.m),col=pal,breaks = seq(0, 50, length.out = 100),display_numbers = plsp.l)
plbl_ht<-pheatmap(-log10(plbl.m),col=pal,breaks = seq(0, 20, length.out = 100),display_numbers = elts.l)

plsp_hto_sp<-rownames(plsp.m[plsp_ht$tree_row[['order']],])
plsp_hto_pl<-colnames(plsp.m[,plsp_ht$tree_col[['order']]])

plbl_hto_pl<-colnames(plbl.m[,plbl_ht$tree_col[['order']]])
plbl_hto_bl<-rownames(plbl.m[plbl_ht$tree_row[['order']],])


spur_ctr=read.table("Spur.mfuzz.centr.txt",h=T,row.names=1,sep='\t')
pheatmap(spur_ctr[plsp_hto_sp,],cluster_cols = F,col=pal2,cluster_rows=F)
spur_ctr[rownames(spur_ctr),]


pliv_ctr=read.table("Pliv.mfuzz.centr.txt",h=T,row.names=1,sep='\t')
pheatmap(pliv_ctr[rev(plsp_hto_pl),],cluster_cols = F,col=pal2,cluster_rows=F)
pheatmap(pliv_ctr[rev(plbl_hto_pl),],cluster_cols = F,col=pal2,cluster_rows=F)

spur_ctr[rownames(spur_ctr),]

blan_ctr=read.table("Blan.mfuzz.centr.txt",h=T,row.names=1,sep='\t')
pheatmap(blan_ctr[plbl_hto_bl,],cluster_cols = F,col=pal2,cluster_rows=F)
spur_ctr[rownames(spur_ctr),]

#### ENRICHMENT TFBS
pal3= colorRampPalette(brewer.pal(8, "BuPu"))(100)

tf_cl<-read.table("Revision/tf_enrich_mfuzz18.tsv",h=T)
tf_cl[,2:19]
tf_cl<-tf_cl[,2:19]
table(apply(tf_cl,1,min)<0.05)
table(apply(tf_cl,2,min)<0.05)
tf_cl.f<-tf_cl[apply(tf_cl,1,min)<0.05,apply(tf_cl,2,min)<0.05]
tf_cl.f1<-tf_cl[apply(tf_cl,1,min)<0.01,apply(tf_cl,2,min)<0.01]

pheatmap(-log10(tf_cl.f1),col=pal3,breaks = seq(0,5, length.out = 100))
?pheatmap
#-apply(modstrathyp.c,2,qnorm)
pheatmap(-apply(elts.m,2,qnorm),col=pal,breaks = seq(0, 10, length.out = 100),cluster_rows = ct,cluster_cols = ct))

pheatmap(as.matrix(-log10(elts.m)),col=pal)
pheatmap(as.matrix(-log10(elts.m)),breaks = seq(0, 24, length.out = 100),col=pal)
#pheatmap(as.matrix(-log10(elts.m)),breaks = seq(2, 24, length.out = 100),col=pal)
pal= colorRampPalette(brewer.pal(8, "YlOrRd"))(100)
pal= colorRampPalette(brewer.pal(8, "YlGn"))(100)




