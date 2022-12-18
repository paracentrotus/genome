require("DESeq2")
require('tidyverse')
library("RColorBrewer")
library("pheatmap")
library("ggplot2")
require("ComplexHeatmap")
require("wesanderson")
#require("viridis")

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

setwd("/Users/fmarletaz/Dropbox/Genomes/Urchin/ATAC-urchin/re-stages")

peaks=read.table("stage_peaks_anc.tsv",h=T,sep='\t',row.names=1)
head(peaks)
head(peaks)
peaks %>% #filter(Pattern %in% selPat) %>% 
  select(X16.cells,EarlyBlastula,LateBlastula,Gastrula,Prism,Pluteus)  -> ma
m = make_comb_mat(ma)
c_ord=c(16,12,6,13,10,15,11,14,3,8,5,9,7,2,4,1)
UpSet(m,comb_order = c_ord,set_order = c('EB','LB','G','PR'))

UpSet(m,set_order = colnames(m),comb_order = c(1:16))
UpSet(m,set_order = colnames(m),comb_order = c(),decreasing = T))
order(comb_size(m),decreasing = T)

sort(c_ord)
c(18,)
comb_size(m)
as.data.frame(comb_size(m))
m[c(40,)]
UpSet(m[comb_size(m) >= 450],
      comb_order = order(comb_size(m[comb_size(m) >= 450]),decreasing = T),
      set_order = colnames(m))
UpSet(m[filt],
      comb_order = order(comb_size(m[filt]),decreasing = T),
      set_order = colnames(m[filt]))

UpSet(m,comb_order = cord,set_order = colnames(m))

order(comb_size(m))
as.data.frame(comb_size(m)) %>% mutate(id = row_number())
c(18,12,17,11,16,10,15,9,14,8,13)
cord=c(13,8,14,9,15,10,16,11,17,12,18,6,5,7,4,3,2,1)
m
peaks %>% head()
abs(peaks$Dist)

as.character(peaks$Pattern)
peaks %>% select(Status,Pattern) %>% group_by(Pattern) %>% tally %>% arrange(desc(n)) %>% head(n=20)
peaks %>% select(Status,Pattern) %>% group_by(Status) %>% tally %>% arrange(desc(n)) %>% head(n=20)

selpat=c("YNNNNN","YYNNNN","NYNNNN","NYYNNN","NNYNNN",
         "NNYYNN","NNNYNN","NNNYYN","NNNNYN","NNNNYY",
         "NNNNNY")


peaks %>% select(Status,Pattern) %>% group_by(Pattern) %>% tally %>% arrange(desc(n)) %>% filter(n>450) %>%select(Pattern) -> selPat#-> locPerStg
selPat=unlist(selPat)
selPat2=c('NYYNNN','NYYYNN','NNNYNN','NNNYYN','NNNNYN',"NYYYYN","NNNNNY","NNNNYY","NNNYYY","YYYNNN",
          "NYYYYY","YYYYNN",'YYYYYN','YYYYYY')
selPat=c("NYNNNN","NYYYNN","NYYNNN","NNYNNN","NNYYNN","NNYYYN","NNNYNN","NNNYYN", "NNNNYN",
          "NYYYYN","NNNYYY","NNYYYY","YYYYNN","YYYYYN","NYYYYY","YYYYYY")
selPat=c("YNNNNN","YYNNNN","NYNNNN","NYYNNN","NNYNNN",
         "NNYYNN","NNNYNN","NNNYYN","NNNNYN","NNNNYY",
         "NNNNNY","NNNYYY","NNYYYN",
         "NNYYYY","NYYYYN","YYYYYN","NYYYYY","YYYYYY")
selPat=c("NYNNNN","NYYNNN","NNYNNN","NNYYNN","NNNYNN","NNNYYN", "NNNNYN","NYYYYN","NYYYYY","YYYYYY")

peaks %>% select(Status,Pattern) %>% filter(Pattern %in% selPat) %>% group_by(Status,Pattern)  -> locPerPat
locPerPat$Pattern=factor(locPerPat$Pattern,levels=selPat)
ggplot(locPerPat, aes(Pattern,fill=Status))+geom_bar()+theme_bw()+scale_fill_brewer(palette='GnBu',direction=-1)


peaks %>% filter(Status!='tss') %>% select(Dist,Pattern) %>% filter(Pattern %in% selPat) %>% mutate(dist_tss=abs(Dist))  -> distPerPat
distPerPat$Pattern=factor(distPerPat$Pattern,levels=selPat)
ggplot(distPerPat, aes(dist_tss,group=Pattern,color=Pattern))+geom_density()+  xlim(0, 50000)#+scale_x_continuous(trans='log10')
ggplot(distPerPat, aes(x=Pattern,y=dist_tss,fill=Pattern))+
  #geom_boxplot(outlier.size =0,outlier.color='lightgrey')+  ylim(0, 50000)+
  scale_fill_brewer(palette='Set3')+
  geom_violin()+ stat_summary(fun.data=data_summary)+
  scale_y_continuous(trans='log2')#+ geom_jitter(shape=16, position=position_jitter(0.2))


peaks %>% select(Status,Pattern) %>% filter(Pattern %in% selPat) %>% group_by(Status,Pattern) %>%
  tally %>% ungroup() %>% group_by(Pattern) %>% mutate(per = n / sum(n)*100) -> locPerPatPer
locPerPatPer$Pattern=factor(locPerPatPer$Pattern,levels=selPat)

ggplot(locPerPatPer,aes(x=Pattern,y=per,fill=Status))+geom_bar(stat="identity")+theme_bw()+scale_fill_brewer(palette='GnBu',direction=-1)

peaks %>% select(Status) %>% group_by(Status) %>% tally %>% mutate(sp='Pli') -> Plpk

ospLoc=read.table("loc_BlaDre.txt",h=T)
ospLoc=rbind(ospLoc,Plpk)
ospLoc %>% group_by(sp) %>% mutate(per = n / sum(n)*100) -> ospLoc
#ospLoc$Pattern=factor(ospLoc$Pattern,levels=selPat)
#  arrange(desc(n)) %>%  filter(n>500) %>% select(Pattern) %>% as_vector() -> selPat#-> locPerStg
ggplot(ospLoc,aes(x=sp,y=per,fill=Status))+geom_bar(stat="identity")+theme_bw()+scale_fill_brewer(palette='GnBu',direction=-1)+
  geom_text(aes(label=n),position = position_stack(vjust = 0.5))
?scale_fill_brewer
locPerPat
# Motifs

colsig=colorRampPalette(brewer.pal(9,'BuPu'))(100)

enrichPval=read.table("enrich_motifs.tsv",h=T)
enrichCt=read.table("mat_tf_pattern.txt",h=T)
enrichCt %>% pivot_longer(cols=NYYYYN:YYYYYY,names_to="Pattern",values_to="NbPk") -> enrichCtLg
enrichPval %>% pivot_longer(cols=NYYYYN:YYYYYY,names_to="Pattern",values_to="Pval") -> enrichPvalLg
enrichPvalLg$NbPk=enrichCtLg$NbPk

enrichPvalLg %>% filter(Pval<0.05) %>% filter(Pattern!='YYYYYY') %>% select(Motif) -> signif_motifs
enrichPvalLg %>% filter(Pval<0.05) %>%filter(Pattern!='YYYYYY') %>%
  group_by(Pattern) %>% top_n(20,Pval) %>% select(Motif) -> signif_motifs_red

enrichPvalLg %>% filter(Pval<0.05) %>%   group_by(Pattern) %>% top_n(30,Pval) %>%
  ungroup() %>% select(Motif) -> signif_motifs_red

length(signif_motifs_red$Motif)
enrichPvalLg %>% filter(Pval<0.05)  %>%  group_by(Pattern) %>% top_n(20) %>% tally
 
enrichPvalLg %>% filter(Pval<0.05) %>% filter(Pattern!='YYYYYY' & Pattern!='NYYYYY') %>% 
  group_by(Pattern) %>% top_n(20,Pval) -> signif_motifs_red2
  group_by(Pattern) %>% tally()

enrichPvalLg %>% filter(Pval<0.05) %>% filter(Motif %in% as.vector(signif_motifs_red2$Motif)) %>% 
  filter(Pattern %in% selPat) -> enrichPvalLg.filt2
signif_motifs_red=as.vector(signif_motifs_red$Motif)
signif_motifs=as.vector(signif_motifs$Motif)

enrichPvalLg.filt2

enrichPval%>% filter(Motif %in% as.vector(signif_motifs)) %>% column_to_rownames(var = "Motif") -> enrichPval.filt
enrichPval%>% filter(Motif %in% as.vector(signif_motifs_red)) %>% column_to_rownames(var = "Motif") -> enrichPval.red.filt
enrichPval%>% filter(Motif %in% as.vector(signif_motifs_red2$Motif)) %>% column_to_rownames(var = "Motif") -> enrichPval.red2.filt
selPat
pheatmap(as.matrix(-log10(enrichPval.filt+1)),color=cividis(100))
pheatmap(as.matrix(-log10(enrichPval.red.filt+1)),color=cividis(100))
clred2=pheatmap(as.matrix(-log10(enrichPval.red2.filt+1)),color=colsig)
clord2=rownames(enrichPval.red2.filt[clred2$tree_row$order,])
nmot=read.table("names_motifs.txt",row.names=1)
names_row_red=paste(rownames(enrichPval.red.filt),nmotl[rownames(enrichPval.red.filt)])
ggplot(enrichPvalLg.filt2,aes(x=Pattern, y=NMotif, size=NbPk, color=-log10(Pval+1))) +
  geom_point(alpha=0.5) +  scale_radius(range=c(1,16),name='Number of Hit') +
  scale_color_distiller(palette='BuPu',direction =1) + theme_bw() + theme(axis.text.x = element_text(angle = -45, hjust = 0))
enrichPvalLg.filt2$Motif=factor(enrichPvalLg.filt2$Motif,levels=rev(clord2))

names_row_red=paste(rownames(enrichPval.red.filt),nmotl[rownames(enrichPval.red.filt)])
enrichPvalLg.filt2$NMotif=nmotl[rownames(enrichPval.red.filt)]

enrichPvalLg.filt2$NMotif=paste(enrichPvalLg.filt2$Motif,nmotl[as.vector(enrichPvalLg.filt2$Motif)])


peaks %>% filter(Status!='tss' & !is.na(Dist)) %>% select(Scaf,Start,End,Pattern) %>%  
  unite("pos", Start:End, sep= "-",remove = FALSE) %>%
  unite("pkp", Scaf:pos, sep= ":",remove = FALSE) %>%
  select(-Scaf,-pos,-Start,-End) %>% rownames_to_column(var='peak') -> peaks_sel

write.table(peaks_sel,file='peaks_pattern.txt',quote=F,)

## CONSERVATION 

head(peaks)
#col=colorRampPalette(c('#ccece6','#99d8c9','#66c2a4','#41ae76','#238b45','#006d2c'))(11)
col <- wes_palette("Zissou1", 12, type = "continuous")
selPat=c("YYNNNN","NYNNNN","NYYNNN","NNYNNN","NNYYNN","NNNYNN","NNNYYN", "NNNNYN","NNNNNY","NYYYYN","NYYYYY","YYYYYY")

peaks %>% filter(Status!='tss' & !is.na(Dist) & Repeat==0) %>% select(Phast,Pattern,Size,Dist,Status) %>%
  filter(Pattern %in% selPat)-> consPerPat
consPerPat$Pattern=factor(consPerPat$Pattern,levels=selPat)
#ggplot(consPerPat, aes(cons,group=Pattern,color=Pattern))+geom_density()+  xlim(0, 50000)#+scale_x_continuous(trans='log10')
ggplot(consPerPat, aes(x=Pattern,y=Phast,fill=Pattern))+
  #geom_boxplot(outlier.size =0,outlier.color='lightgrey')+
  #geom_density(alpha=0.3)+
  geom_violin(scale='width')+ stat_summary(fun.data=data_summary)+  
  scale_fill_manual(values=col)+theme_light()
#len(set_component_height())
#+geom_violin()+ stat_summary(fun.data=data_summary)#+
scale_y_continuous(trans='log2')#+ geom_jitter(shape=16, position=position_jitter(0.2))

