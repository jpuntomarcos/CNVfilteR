# CNVfilteR
R package to remove false positives of CNV calling tools by using SNV calls


## Introduction

Many tools for germline copy number variant (CNV) detection 
from NGS data have been developed. Usually, these tools were 
designed for different input data like WGS, WES or
panel data, and their performance may depend on the CNV size. Available
benchmarks show that all these tools obtain false positives, sometimes 
reaching a very high number of them.

With the aim of reducing the number of false positives,
[CNVfilteR](http://bioconductor.org/packages/CNVfilteR) identifies those 
germline CNVs that can be discarded. This task is performed by using the 
germline single nucleotide variant (SNV) calls that are usually 
obtained in common NGS pipelines. As VCF field interpretation is key 
when working with these files, CNVfilteR specifically supports 
VCFs produced by VarScan2, Strelka/Strelka2, freeBayes, HaplotypeCaller, 
UnifiedGenotyper and Torrent Variant Caller. Additionally, results can be 
plotted using the functions provided by the R/Bioconductor packages
[karyoploteR](http://bioconductor.org/packages/karyoploteR/) and 
[CopyNumberPlots](http://bioconductor.org/packages/CopyNumberPlots/).




## How to use it

Documentation ([vignette](http://bioconductor.org/packages/devel/bioc/vignettes/CNVfilteR/inst/doc/CNVfilteR.html) and [user manual](http://bioconductor.org/packages/devel/bioc/manuals/CNVfilteR/man/CNVfilteR.pdf)) are available at the [CNVfilteR](http://bioconductor.org/packages/CNVfilteR) 
Bioconductor site, allow with instructions to install it.


## Citation

CNVfilteR was developed by José Marcos Moreno-Cabrera and Bernat Gel. For citations,
please refer to the [Bioinformatics paper](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btab356/6275259): 

José Marcos Moreno-Cabrera, Jesús del Valle, Elisabeth Castellanos, Lidia Feliubadaló, Marta Pineda, Eduard Serra, Gabriel Capellá, Conxi Lázaro, Bernat Gel, *CNVfilteR: an R/Bioconductor package to identify false positives produced by germline NGS CNV detection tools*, Bioinformatics, 2021;, btab356, [https://doi.org/10.1093/bioinformatics/btab356](https://doi.org/10.1093/bioinformatics/btab356)


