Thank you for your interest!! Here, I will briefly describe my analysis in the paper "Contribution of CRISPRable DNA on human complex traits" where we investigated the 21 Cas enriched genomic regions' contribution to human complex traits and diseases.

## Heritability enrichment analysis
Analysis of the human genome was based on the GRCh37 assembly, which can be found [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890711-GRCh37-hg19-b37-humanG1Kv37-Human-Reference-Discrepancies#humanG1Kv37). To simplify the process, we considered 22 autosomal chromosomes end to end as one single piece with a total of 2,881,033,286 bases. We cut the entire piece into 20,000 segments, resulting in about 144,052 bases per segment. We counted the number of PAM sequences in each segment, where the reverse strand was also considered. For example, when counting the number of the NGG sequence, we counted the number of both 5’-NGG-3’ and 5’-CCN-3’.
Briefly, running `fecth.py` can get the position of the PAM sequence. 
Taking **GG** for example, 
``` 
python fecth.py human_g1k_v37_bk.fasta GG > GG_hg19_pos.tsv
```
can get **GG** positions into **GG_hg19_pos.tsv** file.

Considering the reverse strand, we also run 
```
python fecth.py human_g1k_v37_bk.fasta CC > CC_hg19_pos.tsv
``` 
to get reversing **GG** positions.

Annotation of Cas enriched regions is based on the number of individual PAM within each segment. For Cas with more than one PAM sequence, we selected the top 2,000 segments that have the highest sum of all its PAMs, denoting Cas enriched regions. These regions were saved into the 'Top10.bed' file.

To investigate the magnitude of these Cas-enriched regions' contribution to human complex traits, we applied stratified linkage disequilibrium (LD) score regression (S-LDSC) to partition the heritability of each human complex trait.
Run the following codes in the LDSC environment to do the heritability enrichment analysis, you can modify it to analyze other Cas/PAM easily.
```
bash make_annot_gzip_Cas.sh
bash compute_l2_Cas.sh
bash Partition_heritability_Cas.sh
```

## Functional annotations
The functional annotations for the autosomes are available [here](https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_baselineLD_v2.2_ldscores.tgz), thanks to the efforts of the Alkes group from the Broad Institute. In our analysis, we test the enrichment of SNPs annotated for our Cas on the regulatory elements using Fisher's exact test. See `Cas_SNP_baseline_OR.R`.


For the X chromosome, we kept the segments in the same length (144 kb) as we did for the autosomes, resulting in 1070 segments in total. Top10% of the 1070 segments were selected as the Cas enriched genomic regions on the X chromosome. Same as the autosomes, SNPs present in the top10% regions were annotated as 1 and other SNPs were annotated as 0 for each Cas. We then queried 125,497 Hapmap3 SNPs on https://www.snp-nexus.org/v4/ for four gene annotations from UCSC, eight epigenetic markers from Roadmap, and five regulator elements from Ensembl. The odds ratio was also obtained from Fisher's exact test.


## References
1. Bulik-Sullivan, B., Loh, PR., Finucane, H. et al. LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. Nat Genet 47, 291–295 (2015). https://doi.org/10.1038/ng.3211
2. Finucane, H., Bulik-Sullivan, B., Gusev, A. et al. Partitioning heritability by functional annotation using genome-wide association summary statistics. Nat Genet 47, 1228–1235 (2015). https://doi.org/10.1038/ng.3404
3. Jorge Oscanoa, Lavanya Sivapalan, Emanuela Gadaleta, Abu Z Dayem Ullah, Nicholas R Lemoine, Claude Chelala, SNPnexus: a web server for functional annotation of human genome sequence variation (2020 update), Nucleic Acids Research, Volume 48, Issue W1, 02 July 2020, Pages W185–W192, https://doi.org/10.1093/nar/gkaa420

