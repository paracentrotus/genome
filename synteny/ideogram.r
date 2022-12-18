require(RIdeogram)
require(tidyverse)
library("wesanderson")

setwd("~/Dropbox/Genomes/Synteny/synt-forge/Urchin/")

# EXAMPLE

data(karyotype_dual_comparison, package="RIdeogram")
head(karyotype_dual_comparison)
data(synteny_dual_comparison, package="RIdeogram")
head(synteny_dual_comparison)

# REST

lg=read.table("Pliv_Spur_synt_ext.txt",h=T)

head(lg)
lg %>%  group_by(s1chr) %>%   summarise(End=max(s1gp)) %>% arrange(-End) %>% rename(Chr=s1chr) %>% 
  mutate(Start=1,species='Paracentrotus',fill=969696,size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_pliv

lg  %>% group_by(s2chr) %>%   summarise(End=max(s2gp)) %>% arrange(-End) %>% rename(Chr=s2chr) %>% 
  mutate(Start=1,species='Strongylocentrotus',fill=969696,size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_spur


k_pliv_ord=c('chrom1','chrom2','chrom3','chrom4','chrom5','chrom6','chrom7','chrom8','chrom9','chrom10','chrom11','chrom12','chrom13','chrom14','chrom15','chrom16','chrom17','chrom18')
k_spur_ord=c('scaffold_5','scaffold_1','scaffold_6','scaffold_18','scaffold_17','scaffold_10','scaffold_20','scaffold_7',
  'scaffold_3','scaffold_14','scaffold_2','scaffold_13','scaffold_21','scaffold_9','scaffold_11',
  'scaffold_8','scaffold_19','scaffold_15','scaffold_12','scaffold_16','scaffold_4')

k_pliv=k_pliv[order(match(k_pliv$Chr,k_pliv_ord)),]
k_spur=k_spur[order(match(k_spur$Chr,k_spur_ord)),]

rbind(k_pliv,k_spur) -> chroms


lg$Species_1=match(lg$s1chr,k_pliv$Chr)
lg$Species_2=match(lg$s2chr,k_spur$Chr)

lg %>% filter(scol=='S') %>% select(Species_1,s1gp,Species_2,s2gp) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill='cccccc') %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) -> lgt


ideogram(karyotype = chroms, synteny = lgt)

# Lytechinus vs. Strongylocentrutus

lg=read.table("Spur_Lvar_synt_ext.txt",h=T)
head(lg)

lg %>%  group_by(s1chr) %>%   summarise(End=max(s1gp)) %>% arrange(-End) %>% rename(Chr=s1chr) %>% 
  mutate(Start=1,species='Strongylocentrotus',fill='Strongylocentrotus',size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_spur

lg  %>% group_by(s2chr) %>%   summarise(End=max(s2gp)) %>% arrange(-End) %>% rename(Chr=s2chr) %>% 
  mutate(Start=1,species='Lytechinus',fill='Lytechinus',size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_lvar

k_spur_ord=c('scaffold_5','scaffold_1','scaffold_6','scaffold_18','scaffold_17','scaffold_20','scaffold_10','scaffold_7',
             'scaffold_3','scaffold_14','scaffold_2','scaffold_13','scaffold_21','scaffold_9','scaffold_11',
             'scaffold_8','scaffold_19','scaffold_15','scaffold_12','scaffold_16','scaffold_4')
k_lvar_ord=c('chr7','chr8','chr3','chr9','chr19','chr6','chr18','chr1','chr11','chr10',
             'chr2','chr17','chr4','chr5','chr13','chr12','chr16','chr14','chr15')

k_spur=k_spur[order(match(k_spur$Chr,k_spur_ord)),]
k_lvar=k_lvar[order(match(k_lvar$Chr,k_lvar_ord)),]

rbind(k_spur,k_lvar) -> chroms

lg$Species_1=match(lg$s1chr,k_spur$Chr)
lg$Species_2=match(lg$s2chr,k_lvar$Chr)

lg %>% filter(scol=='S') %>% select(Species_1,s1gp,Species_2,s2gp) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill='cccccc') %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) -> lgt

ideogram(karyotype = chroms, synteny = lgt)

### Spur lyt by CLG

read.table('spur_alg_assign.txt',h=T) %>% column_to_rownames('clg') -> spu
cols=read.table('colors.txt',row.names=1)
lg=read.table("Spur_Lvar_synt_ext.txt",h=T)
spu 
head(lgr)
head(lgrf)
lg %>% filter(s1chr=='scaffold_20')
#lg <-lgrf

lg %>%  group_by(s1chr) %>%   summarise(End=max(s1gp)) %>% arrange(-End) %>% rename(Chr=s1chr) %>% 
  mutate(Start=1,species='Strongylocentrotus',fill='Strongylocentrotus',size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_spur

lg  %>% group_by(s2chr) %>%   summarise(End=max(s2gp)) %>% arrange(-End) %>% rename(Chr=s2chr) %>% 
  mutate(Start=1,species='Lytechinus',fill='Lytechinus',size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_lvar

k_spur_ord=c('scaffold_5','scaffold_1','scaffold_6','scaffold_18','scaffold_17','scaffold_20','scaffold_10','scaffold_7',
             'scaffold_3','scaffold_14','scaffold_2','scaffold_13','scaffold_21','scaffold_9','scaffold_11',
             'scaffold_8','scaffold_19','scaffold_15','scaffold_12','scaffold_16','scaffold_4')
k_lvar_ord=c('chr7','chr8','chr3','chr9','chr19','chr6','chr18','chr1','chr11','chr10',
             'chr2','chr17','chr4','chr5','chr13','chr12','chr16','chr14','chr15')

k_spur=k_spur[order(match(k_spur$Chr,k_spur_ord)),]
k_lvar=k_lvar[order(match(k_lvar$Chr,k_lvar_ord)),]


rbind(k_spur,k_lvar) -> chroms
lg$Species_1=match(lg$s1chr,k_spur$Chr)
lg$Species_2=match(lg$s2chr,k_lvar$Chr)
spu[(spu$chr=='scaffold_1'),]$clg
spu[(spu$clg==clg),]$chr

lg %>% filter(spu[clg,]==s1chr & !is.na(clg)) -> lg

lg %>% filter(scol=='S' & !is.na(clg)) %>% select(Species_1,s1gp,Species_2,s2gp,clg) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill=cols[match(clg,rownames(cols),),]) %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clg) -> lgt

spu[(spu$chr=='scaffold_1'),]$clg

#rcols=data.frame(row.names=rownames(cols),col=str_replace(cols$col,'#',''))
match(cols,lg$clg)
head(lgt)
ideogram(karyotype = chroms, synteny = lgt)
convertSVG("chromosome.svg", device = "png")




# Lytechinus 

read.table('pliv_alg_assign_lg.txt',h=T) %>% column_to_rownames('clg') -> pla
cols=read.table('colors.txt',row.names=1)
lgr=read.table("Pliv_Lvar_synt_ext.txt",h=T)
levels(as.factor(lgr$clg))
#lgr %>% filter(pla[clg,]==s1chr & !is.na(clg)) -> lgrf
lgr
lg <-lgr
  
lg %>%  group_by(s1chr) %>%   summarise(End=max(s1gp)) %>% arrange(-End) %>% rename(Chr=s1chr) %>% 
  mutate(Start=1,species='Paracentrotus',fill='Paracentrotus',size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_pliv

lg  %>% group_by(s2chr) %>%   summarise(End=max(s2gp)) %>% arrange(-End) %>% rename(Chr=s2chr) %>% 
  mutate(Start=1,species='Lytechinus',fill='Lytechinus',size=12,color=252525) %>% relocate(Start,.after=Chr)-> k_lvar
 
k_pliv_ord=c('chrom1','chrom2','chrom3','chrom4','chrom5','chrom6','chrom7','chrom8','chrom9','chrom10','chrom11','chrom12',
             'chrom13','chrom14','chrom15','chrom16','chrom17','chrom18')

k_lvar_ord=c('chr7','chr8','chr3','chr9','chr19','chr6','chr18','chr1','chr11','chr10',
             'chr2','chr17','chr4','chr5','chr13','chr12','chr16','chr14','chr15')

k_pliv=k_pliv[order(match(k_pliv$Chr,k_pliv_ord)),]
k_lvar=k_lvar[order(match(k_lvar$Chr,k_lvar_ord)),]


rbind(k_pliv,k_lvar) -> chroms
lg$Species_1=match(lg$s1chr,k_pliv$Chr)
lg$Species_2=match(lg$s2chr,k_lvar$Chr)

lg %>% filter(scol=='S' & !is.na(clg)) %>% select(Species_1,s1gp,Species_2,s2gp,clg) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill=cols[match(clg,rownames(cols),),]) %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clg) -> lgt

lg %>% filter(scol=='S' & !is.na(clg)) %>% select(Species_1,s1gp,Species_2,s2gp,clg) %>% rename(Start_1=s1gp,Start_2=s2gp) %>%
  mutate(End_1=Start_1+1000,End_2=Start_2+1000,fill='cccccc') %>%
  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2) %>% select(-clg) -> lgt

head(lgt)  
'cccccc'
#rcols=data.frame(row.names=rownames(cols),col=str_replace(cols$col,'#',''))

ideogram(karyotype = chroms, synteny = lgt)
convertSVG("chromosome.svg", device = "png")

### TEST DENSITY Gene

lg %>%  group_by(s1chr) %>%   summarise(End=max(s1gp)) %>% arrange(-End) %>% rename(Chr=s1chr) %>% 
  mutate(Start=1) %>% relocate(Start,.after=Chr)-> k_pliv

gene_density <- GFFex(input = "../gtfs/Pliv_aH2p.gtf", karyotype = 'Pliv_karyo.txt', feature = "gene", window = 1000000)

!is.na(lg$clg)


