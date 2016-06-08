# Sparc: a sparsity-based consensus algorithm for long erroneous sequencing reads

https://peerj.com/articles/2016.pdf

test data:
https://sourceforge.net/projects/sparc-consensus/files/testdata/


Parameters:

b: backbone file.

m: the reads mapping files produced by blasr, using option -m 5. (A blasr example command: blasr -nproc 32 query.reads.fasta backbone.fasta -bestn 1 -m 5 -minMatch 19 -out backbone.mapped.m5)

k: k-mer size (suggested range: [1,2]).

c: coverage threshold (range: [1,5], suggest: 2).

t: adaptive threshold (suggested range[0.0,0.3]).

g: skip size, the larger the value, the more memory efficient the algorithm is (suggested range: [1,3]).

HQ_Prefix: Shared prefix of the high quality read names. (e.g. if the sec-gen sequences have names >Contig_xxx, then ‘Contig’ is a shared prefix of the high quality reads)

boost: boosting weight for the high quality reads (suggested range: [1,10]).

Example command: 

Using third-gen data only:

Sparc b Backbone.fa m backbone.mapped.m5 k 2 g 2 c 2 t 0.1 o ConsensusOutput

Using hybrid data:

Sparc b backbone.fasta m backbone.mapped.m5 k 2 g 2 c 2 t 0.1 HQ_Prefix Contig boost 5 o ConsensusOutput



