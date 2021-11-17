

## Count PAMs on 22 chromosomes

We used the GRCh37 assembly of the human genome, which can be found [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890711-GRCh37-hg19-b37-humanG1Kv37-Human-Reference-Discrepancies#humanG1Kv37). 
Briefly, run `fecth.py` can get position of PAM sequence. 
Taking **NGG** for example, 
``` 
python fecth.py human_g1k_v37_bk.fasta GG > NGG_hg19_pos.tsv
```
can get **NGG** positions into **NGG_hg19_pos.tsv** file.

Considering the reverse strand, we also run 
```
python fecth.py human_g1k_v37_bk.fasta CC > CC_hg19_pos.tsv
``` 
to get reversing **NGG** positions.

Then, using `make_bed_PAM.R`, we split **NGG** positions into 22 files for 22 chromosomes.

## Cut the genome into 20,000 segments
We simply think the 22 chromosomes as one. Then we cut it into 20,000 segments, and count the number of PAM sequence in each segment. Next, we select segments contain the largest number of PAM sequence, save those locations in one `Top10.bed` file for LDSC. 
See `PAMs_top10.R`.


## Running LDSC
The code in `make_annot_with_bed.R` is modified from the [LD Score Estimation Tutorial](https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial) to compute LD scores for all 22 chromosomes.
**Partitioning heritability** of 28 human complex traits also followed the instructions from [Partitioned Heritability](https://github.com/bulik/ldsc/wiki/Partitioned-Heritability) of LDSC software.









