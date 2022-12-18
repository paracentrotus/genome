require(DESeq2)
require(pheatmap)
require(RColorBrewer)
require(tidyverse)
require(Mfuzz)
require(wesanderson)



setwd("~/Dropbox/Genomes/Urchin/Paracentrotus_genome/Analyses/Transcriptomes/mfuzz/Revision/")

### RUNNING MFUZZ FPKM

Pliv.fpkm=read.table("Pliv_rename_fpkms.tsv",h=T)

Pliv.fpkm=Pliv.fpkm[apply(Pliv.fpkm,1,max)>2,]

Pliv.fpkm.red=as.matrix(Pliv.fpkm[,c(1,6,9,11,13,15,16,17)])
Pliv.fpkm.log=as.matrix(log10(Pliv.fpkm.red+1))

par(mfrow=c(1,1))
# I tried to run it with all ther data, but the signal is a bit fuzzy and it looks like it wants to make 
# too many clusters; I tried to pick select stages that make more sense (in the heatmap)

Pliv.set <- new("ExpressionSet",exprs=Pliv.fpkm.log)
Pliv.set <- filter.std(Pliv.set,min.std=0)
Pliv.set <- standardise(Pliv.set)
Pliv.m <- mestimate(Pliv.set)
#stg.R=stg[c(1,6,9,11,13,15,16,17)]
stg.Rn=c('Egg','32C','EBl','SBl','MBl','Ga','Pr','Pl')

# determine the nombre of clusters by computing the distance between centroids 
Pliv.dmin=Dmin(Pliv.set,Pliv.m,crange=seq(4,48,2),repeats=3,visu=TRUE)
plot(seq(4,48,2),Pliv.dmin,pch=19)
Pliv.clust <-mfuzz(Pliv.set,c=18,m=Pliv.m)
pdf("Pliv.mfuzz.plot.pdf",paper="a4")
mfuzz.plot(Pliv.set,cl=Pliv.clust,mfrow=c(3,3),colo=brewer.pal(9,'GnBu'),time.labels=stg.Rn,new.window=FALSE)
pheatmap(cor(t(Pliv.clust[[1]])))
pheatmap(cor(Pliv.clust$membership))
pal=rev(colorRampPalette(brewer.pal(9,'RdBu'))(100))

pheatmap(Pliv.clust$centers,cluster_cols = F,col=pal)
head(Pliv.clust$membership)
dev.off()
head(Pliv.clust$centers)
Pliv.over <- overlap(Pliv.clust)
overlap.plot(Pliv.clust ,over=Pliv.over,thres=0.05)

write.table(as.data.frame(Pliv.clust$centers),file='Pliv.mfuzz.centr.txt',quote=F,sep='\t')
write.table(as.data.frame(Pliv.clust$cluster),file='Pliv.mfuzz.clust.txt',quote=F,sep='\t')

###### Strongylocentrotus 

Spur.fpkm=read.table('Spur_rename_fpkms.tsv',h=T)

Spur.fpkm=Spur.fpkm[apply(Spur.fpkm,1,max)>2,]
Spur.fpkm.red=as.matrix(Spur.fpkm[,c(1:4,6:8,11)])
#colnames(Spur.fpkm.red)
#Spur.fpkm.log=as.matrix(log10(Spur.fpkm+1))
Spur.fpkm.log.red=as.matrix(log10(Spur.fpkm.red+1))
# I tried to run it with all ther data, but the signal is a bit fuzzy and it looks like it wants to make 
# too many clusters; I tried to pick select stages that make more sense (in the heatmap)
head(Spur.fpkm.log.red)
Spur.set <- new("ExpressionSet",exprs=Spur.fpkm.log.red)
Spur.set <- filter.std(Spur.set,min.std=0)
Spur.set <- standardise(Spur.set)
Spur.m <- mestimate(Spur.set)
stg.R=stg[c(1,6,9,11,13,15,16,17)]
Spur.stg.red=c('eggs','cleavage','swim_blast','mes_blast','mid_gastr','late_gastr','prism','larva')
colnames(Spur.fpkm.log.red)
# determine the nombre of clusters by computing the distance between centroids 
Spur.dmin=Dmin(Spur.set,Spur.m,crange=seq(4,48,2),repeats=3,visu=TRUE)
plot(seq(4,48,2),Spur.dmin,pch=19)

Spur.clust <-mfuzz(Spur.set,c=18,m=Spur.m)

pheatmap(cor(t(Spur.clust[[1]])))
pheatmap(cor(Spur.clust$membership))

pdf("Spur.mfuzz.plot.pdf",paper="a4")
mfuzz.plot(Spur.set,cl=Spur.clust,mfrow=c(3,3),time.labels=Spur.stg.red,new.window=FALSE,colo=brewer.pal(9,'GnBu'),)
dev.off()

Spur.over <- overlap(Spur.clust)
overlap.plot(Spur.clust ,over=Spur.over,thres=0.05)
write.table(as.data.frame(Spur.clust$cluster),file='Spur.mfuzz.clust.txt',quote=F,sep='\t')
write.table(as.data.frame(Spur.clust$centers),file='Spur.mfuzz.centr.txt',quote=F,sep='\t')

pal2=rev(colorRampPalette(brewer.pal(9,'RdBu'))(100))
colnames(Spur.clust$centers)<-Spur.stg.red
pheatmap(Spur.clust$centers,cluster_cols = F,col=pal2)


###### Branchiostoma 
Blan.fpkm=read.table('../../data/Bla_cRPKMs_sel.tsv')

colnames(Blan.fpkm)
#colnames(Blan.fpkm)=gsub("Bl_", "", names(Blan.fpkm))
#Blan.stg=gsub("\\_.", "", names(Blan.fpkm))

Blan.kstg=c('Eggs','32cells','8h','15h','21h','36h','60h','PreM')
Blan.fpkm %>% select(Bl_Eggs_A,Bl_32cells_A,Bl_8h_a,Bl_15h_A,Bl_21h_h,Bl_36h_d,Bl_60h_A,Bl_PreM_u) -> Blan.fpkm.red
colnames(Blan.fpkm.red)<-Blan.kstg
Blan.fpkm.red=Blan.fpkm.red[apply(Blan.fpkm.red,1,max)>2,]

Blan.fpkm.log.red=as.matrix(log10(Blan.fpkm.red+1))

par(mfrow=c(1,1))

Blan.set <- new("ExpressionSet",exprs=Blan.fpkm.log.red)
Blan.set <- filter.std(Blan.set,min.std=0)
Blan.set <- standardise(Blan.set)
Blan.m <- mestimate(Blan.set)

# determine the nombre of clusters by computing the distance between centroids 
Blan.dmin=Dmin(Blan.set,Blan.m,crange=seq(4,48,2),repeats=3,visu=TRUE)
plot(seq(4,48,2),Blan.dmin,pch=19)

Blan.clust <-mfuzz(Blan.set,c=18,m=Blan.m)
pdf("Blan.mfuzz.plot.pdf",paper="a4")
mfuzz.plot(Blan.set,cl=Blan.clust,mfrow=c(3,3),new.window=FALSE,colo=brewer.pal(9,'GnBu'))
dev.off()
write.table(as.data.frame(Blan.clust$cluster),file='Blan.mfuzz.clust.txt',quote=F,sep='\t')
write.table(as.data.frame(Blan.clust$centers),file='Blan.mfuzz.centr.txt',quote=F,sep='\t')

pal2=rev(colorRampPalette(brewer.pal(9,'RdBu'))(100))
colnames(Blan.clust$centers)<-Spur.stg.red
pheatmap(Blan.clust$centers,cluster_cols = F,col=pal2)

pheatmap(cor(t(Blan.clust[[1]])))
pheatmap(cor(Blan.clust$membership))


Blan.over <- overlap(Blan.clust)
overlap.plot(Blan.clust ,over=Blan.over,thres=0.05)

