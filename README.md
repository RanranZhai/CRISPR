

## count PAMs on 22 chromosomes

We used the GRCh37 assembly of the uman genome, which can be found [here](). 
Briefly, run `fecth.py` can get position of PAM sequence. For example, 
``` 
python fecth.py human_g1k_v37_bk.fasta GG > GG_hg19_pos.tsv
```
can get **GG** positions into **GG_hg19_pos.tsv** file.

Considering the reverse strand, we also run `python fecth.py human_g1k_v37_bk.fasta CC > CC_hg19_pos.tsv` to get reversing **GG** positions.









