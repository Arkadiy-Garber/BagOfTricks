# BagOfTricks
A collection of python and bash scripts for various bioinformatics-related tasks

## FindFindMeHemes.py
#### Prediction of heme-binding motifs
    FindMeHemes.py -contigs contigs.fna -outDir ./ -mode metagenome

## Seqs.py
#### Print to stdout assembly stats
    Seqs.py -fasta assembly.fasta

## checkm-quality.sh
#### A wrapper for five CheckM commands to generate alignments, trees, and summary stats (completion/redundancy scores) for a directory of genomes or bins (metagenome-assembled genomes)
    checkm-quality.sh .fa path/to/bins/ 16

## coding-density.py
#### Takes as input a GFF file that represents a genome. Outputs to stdout a tab-delimited summary of the following information (in order)
file name - coding density - genome size (in Mb) - number of coding genes - number of pseudogenes (if included in annotation) - number of transposases - number of mutS genes - number of mutL genes - number of recA genes

    coding-density.py -gff genome.gff

## daltons.py
#### Calculate exact weights in daltons from protein sequences
    daltons.py -f proteins.faa

## masker.py
#### Mask alignment positions with too many gaps (user-defined fraction)
    masker.py -i alignment.fa -o alignment.masked.fa -m 0.5

## ribosome.py
#### Translates user-provided DNA sequences to proteins
    ribosome.py -i coding_genes.ffn -o protein_translations.faa -x y





