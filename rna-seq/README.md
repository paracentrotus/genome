# Stage-specific expression modules 

Temporal clustering was performed using `mfuzz` using the code described in `mfuzz.r` using the gene expression data included in `Bla_fpkms_sel.tsv.gz`,`Pliv_rename_fpkms.tsv.gz`,`Spur_rename_fpkms.tsv.gz`. 
For each species, we plotted the minimum centroid distance (Dmin) for increasing number of clusters (`{sp}_Dmin.pdf`), the correlation of centroids (`{sp}_cor_centroids`) and the correlation of membership values (`{sp}_cor_membership.pdf`) representing the shared gene content of fuzzy clusters. 

The enrichment of genes belonging to the same gene families among pairs of cluster was calculated as described in `co-enrichment.ipynb` and subsequently plotted using ` plot_comp.r`. 

