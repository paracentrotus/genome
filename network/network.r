require(igraph)
require(tidyverse)
require(pheatmap)
library(viridis)
require(RColorBrewer)
setwd("/Users/fmarletaz/Dropbox (UCL)/Genomes/Urchin/Paracentrotus_genome/Analyses/Network")

procNet <- function(nwk) {
  qt<-quantile(cleav.nwk$prob,probs = c(0,.10,.25,.5,.75,.9,1))
  nwk %>% filter(prob>qt[6]) %>% separate(tf_target,c('tf1','tf2')) -> nwk.f
  nwk.g<-graph_from_data_frame(nwk.f)
  return(nwk.g)
}

cleav.nwk=read_tsv("ansr_16-cells_network.tsv.gz")
earbl.nwk=read_tsv("ansr_EarlyBlastula_network.tsv.gz")
latbl.nwk=read_tsv("ansr_LateBlastula_network.tsv.gz")
gastr.nwk=read_tsv("ansr_Gastrula_network.tsv.gz")
prism.nwk=read_tsv("ansr_Prism_network.tsv.gz")
plute.nwk=read_tsv("ansr_Pluteus_network.tsv.gz")

# filter network 
cleav.g<-procNet(cleav.nwk)
earbl.g<-procNet(earbl.nwk)
latbl.g<-procNet(latbl.nwk)
gastr.g<-procNet(gastr.nwk)
prism.g<-procNet(prism.nwk)
plute.g<-procNet(plute.nwk)

# compute degrees

cleav.dg<-degree(cleav.g, v = V(cleav.g),mode = c("out"),loops = T, normalized = T)
earbl.dg<-degree(earbl.g, v = V(earbl.g),mode = c("out"),loops = T, normalized = T)
latbl.dg<-degree(latbl.g, v = V(latbl.g),mode = c("out"),loops = T, normalized = T)
gastr.dg<-degree(gastr.g, v = V(gastr.g),mode = c("out"),loops = T, normalized = T)
prism.dg<-degree(prism.g, v = V(prism.g),mode = c("out"),loops = T, normalized = T)
plute.dg<-degree(plute.g, v = V(plute.g),mode = c("out"),loops = T, normalized = T)

# plot centrality by TF category
m2f<-read_tsv("Pliv_PqN3S.gimme.vertebrate.v5.0.motif2factors.txt")
m2f %>% separate(Motif,c('c1','c2','c3','fam_tf','tfid'),sep="[.]",remove=F) %>%
  select(-c1,-c2,-c3,-tfid,-Curated) %>% distinct(Factor,.keep_all = T) %>%  
  column_to_rownames(var='Factor') -> m2fc

as.data.frame(cleav.dg) %>% rename(deg=cleav.dg) %>% rownames_to_column(var="gid") %>% mutate(tf=m2fc[gid,]$fam_tf) %>% filter(!is.na(tf))-> cleav.dga
as.data.frame(earbl.dg) %>% rename(deg=earbl.dg) %>% rownames_to_column(var="gid") %>% mutate(tf=m2fc[gid,]$fam_tf) %>%  filter(!is.na(tf))-> earbl.dga
as.data.frame(latbl.dg) %>% rename(deg=latbl.dg) %>% rownames_to_column(var="gid") %>% mutate(tf=m2fc[gid,]$fam_tf) %>%  filter(!is.na(tf))-> latbl.dga
as.data.frame(gastr.dg) %>% rename(deg=gastr.dg) %>% rownames_to_column(var="gid") %>% mutate(tf=m2fc[gid,]$fam_tf) %>%  filter(!is.na(tf))-> gastr.dga
as.data.frame(prism.dg) %>% rename(deg=prism.dg) %>% rownames_to_column(var="gid") %>% mutate(tf=m2fc[gid,]$fam_tf) %>%  filter(!is.na(tf))-> prism.dga
as.data.frame(plute.dg) %>% rename(deg=plute.dg) %>% rownames_to_column(var="gid") %>% mutate(tf=m2fc[gid,]$fam_tf) %>%  filter(!is.na(tf))-> plute.dga

as.data.frame(gastr.dg) %>% rename(deg=gastr.dg) %>% rownames_to_column(var="gid") %>% mutate(tf=m2fc[gid,]$fam_tf) %>% mutate(tf=replace_na(tf,"non-TF")) -> gasta.dge

gasta.dge %>% mutate(status=ifelse(tf=='non-TF','non-TF','TF')) -> gasta.dge
ggplot(gasta.dge,aes(x=status,y=deg,fill=status)) + geom_violin() + theme_light()

cleav.dga %>% filter(deg>0) %>% group_by(tf) %>% summarise(ctr=median(deg),n_gen=n()) %>% mutate(stg='16-cells') -> cleav.med
earbl.dga %>% filter(deg>0) %>% group_by(tf) %>% summarise(ctr=median(deg),n_gen) %>% mutate(stg='EarlyBlastula') -> earbl.med
latbl.dga %>% filter(deg>0) %>% group_by(tf) %>% summarise(ctr=median(deg)) %>% mutate(stg='LateBlastula') -> latbl.med
gastr.dga %>% filter(deg>0) %>% group_by(tf) %>% summarise(ctr=median(deg)) %>% mutate(stg='Gastrula') -> gastr.med
prism.dga %>% filter(deg>0) %>% group_by(tf) %>% summarise(ctr=median(deg)) %>% mutate(stg='Prism') -> prism.med
plute.dga %>% filter(deg>0) %>% group_by(tf) %>% summarise(ctr=median(deg)) %>% mutate(stg='Pluteus') -> plute.med

stg.med <- rbind(cleav.med,earbl.med,latbl.med,gastr.med,prism.med,plute.med) 
stg.med %>% pivot_wider(id_cols=stg,names_from=tf,values_from=ctr,values_fill=0) %>% column_to_rownames(var='stg')-> stf.med.mat
scale(stf.med.mat) -> stf.med.mat.sc
pheatmap(t(stf.med.mat),cluster_cols = F,col=magma(100))

pheatmap(t(stf.med.mat.sc),cluster_cols = F,col=pal)
pal=rev(colorRampPalette(brewer.pal(9,'RdBu'))(100))


# plot top 10 factors for each network  

gname<-read.table("Pliv.gname.txt")
colnames(gname)=c('gnm','gid')
gname %>% column_to_rownames(var='gid') -> gname

topNet <- function(g.df,g) {
  as.data.frame(g.df) -> g.df
  g.df$gnm<-gname[rownames(g.df),]
  names(g.df)=c('ctr','gnm')
  g.df %>% arrange(desc(ctr)) %>% slice_max(ctr,n=10)-> g.top
  g.gtp <- induced.subgraph(graph=g,vids=rownames(g.top))
  g.gtp<-set.vertex.attribute(g.gtp, "name", value=g.top$gnm)
  return(g.gtp)
}

cleav.gtp<-topNet(cleav.dg,cleav.g)
earbl.gtp<-topNet(earbl.dg,earbl.g)
latbl.gtp<-topNet(latbl.dg,latbl.g)
gastr.gtp<-topNet(gastr.dg,gastr.g)
prism.gtp<-topNet(prism.dg,prism.g)
plute.gtp<-topNet(plute.dg,plute.g)

par(mfrow=c(3,2),mar=c(1.5,1.5,1.5,1.5))
plot(cleav.gtp,edge.arrow.size = .5,frame.color='red',main="16-cells")
plot(earbl.gtp,edge.arrow.size = .5,size=5,main="Early Blastula")
plot(latbl.gtp,edge.arrow.size = .5,main="Late Blastula")
plot(gastr.gtp,edge.arrow.size = .5,main="Gastrula")
plot(prism.gtp,edge.arrow.size = .5,main="Prism")
plot(plute.gtp,edge.arrow.size = .5,main="Pluteus")
?plot
#plot(gastr.gtp,layout=layout_with_kk(gastr.gtp),edge.arrow.size = .5)

cleav.gtp<-topNet(cleav.dg,cleav.g)
plot(cleav.gtp,edge.arrow.size = .5)

# plot endomesoderme et 

endomes<-read.table("Endomesoderm.txt")
names(endomes)=c('layer','name','gid','core')
endomes %>% filter(core==1)-> emd

selNet <- function(g,gl) {
  gl<-gl[gl$gid %in% V(latbl.g)$name,]
  g.msd<-induced.subgraph(graph=g,vids=gl$gid)
  g.msd.rn<-set.vertex.attribute(g.msd, "name", value=gl$name)
  return(g.msd.rn)
}

cleav.msd<-selNet(cleav.g,emd)
earbl.msd<-selNet(earbl.g,emd)
latbl.msd<-selNet(latbl.g,emd)
gastr.msd<-selNet(gastr.g,emd)
prism.msd<-selNet(prism.g,emd)
plute.msd<-selNet(plute.g,emd)

plot(cleav.msd,edge.arrow.size = .5,main='16-cells')#,layout=layout_with_kk(cleav.msd))
plot(earbl.msd,edge.arrow.size = .5,main='Early Blastula')#,layout=layout_with_kk(gastr.msd.rn))
plot(latbl.msd,edge.arrow.size = .5,main='Late Blastula')
plot(gastr.msd,edge.arrow.size = .5,main='Gastrula')
plot(prism.msd,edge.arrow.size = .5,main='Prism')
plot(plute.msd,edge.arrow.size = .5,main='Pluteus')
#gastr.msd<-induced.subgraph(graph=gastr.g,vids=emd$gid)
#gastr.msd.rn<-set.vertex.attribute(gastr.msd, "name", value=emd$name)

skelt<-read.table("skeletogenic.txt")
names(skelt)=c('name','gid')
cleav.skl<-selNet(cleav.g,skelt)
earbl.skl<-selNet(earbl.g,skelt)
latbl.skl<-selNet(latbl.g,skelt)
gastr.skl<-selNet(gastr.g,skelt)
prism.skl<-selNet(prism.g,skelt)
plute.skl<-selNet(plute.g,skelt)

plot(cleav.skl,edge.arrow.size = .5,main='16-cells')#,layout=layout_with_kk(cleav.skl))
plot(earbl.skl,edge.arrow.size = .5,main='Early Blastula')#,layout=layout_with_kk(gastr.skl.rn))
plot(latbl.skl,edge.arrow.size = .5,main='Late Blastula')
plot(gastr.skl,edge.arrow.size = .5,main='Gastrula')
plot(prism.skl,edge.arrow.size = .5,main='Prism')
plot(plute.skl,edge.arrow.size = .5,main='Pluteus')


gastr.skl<-induced.subgraph(graph=gastr.gp,vids=skelt$gid)
gastr.skl.rn<-set.vertex.attribute(gastr.skl, "name", value=skelt$name)
plot(gastr.skl.rn)

gastr.df %>% dplyr::arrange(desc(gastr.df)) %>% filter(grepl("Zic",gnm))
test<-subgraph.edges(graph=gastr.g, eids=which(E(gastr.g)$label=="PL05901"), delete.vertices = TRUE)
test<-subgraph(graph=gastr.g, eids=which(V(gastr.g)$name=="PL24291"))
'PL21602'
selegoV <- ego(gastr.gp, order=1, nodes = which(V(gastr.g)$name=="PL24291"), mode = "out", mindist = 0)
sox.i<-names(unlist(selegoV))[names(unlist(selegoV)) %in% as.vector(m2f$Factor)]
selegoV <- ego(gastr.gp, order=1, nodes = which(V(gastr.g)$name %in% sox.i), mode = "out", mindist = 0)

selegoG <- induced_subgraph(gastr.gp,unlist(selegoV))
plot(selegoG)
sox.i
ego_size(gastr.gp,order=1,mode="out",nodes = which(V(gastr.g)$name=="PL24291"))
sg1 <- decompose.graph(gastr.g,mode="weak")
g2 <- induced.subgraph(graph=gastr.gp,vids=as.vector(sox.i))
g2 <- induced.subgraph(graph=gastr.gp,vids=as.vector(gastr.top))
g2 <- induced.subgraph(graph=gastr.gp,vids=rownames(gastr.top))
gastr.top
plot(g2)
neighverts <- unique(unlist(sapply(sg1,FUN=function(s){if(any(V(s)$name %in% subv)) V(s)$name else NULL})))
g3 <- induced.subgraph(graph=g1,vids=neighverts)
plot(g3)
neighbors(graph.gp, ,
which(V(gastr.gp)$name=="PL03351")
test.g<-degree(test, v = V(test),mode = c("out"),loops = TRUE, normalized = FALSE)
plot(test)
#, vertex.size=test.g*6, vertex.color=rgb(0.1,0.7,0.8,0.5) )


