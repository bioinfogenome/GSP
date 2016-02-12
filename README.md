# GSP
A tool for Genome Specific Primers design in polyploid species

In a polyploid species, genes from the homeologous genomes exhibit a very high sequence similarity, especially in the coding regions. This makes it difficult to design genome-specific primers to amplify specific sequences from the individual genomes in the polyploid genome background. Development of genome-specific primers for important genes in polyploid species is very useful and critical not only for the study of sequence diversity and association mapping of genes in natural populations, but also for the development of gene-based functional markers for marker-assisted breeding. GSP is a powerful web-based platform for designing genome-specific primers that can distinguish the sequences among different genomes in a polypoid species.

Web server: http://probes.pw.usda.gov/GSP/

# Installation
1. Download and upzip.
2. change the path to upzip folder.
3. cmake ./src
4. make

# Usage
Usage: GSP
-r    blast table result path

-d    sequence database path (blast database path)
-a    fasta file path (for design specefic priemrs in multiple sequences only)
-b    bed path
-m    muscle path
-p    primer3 path
-t    primer3 parameters file path
-o    output path
-f    hit sequence flanking length (default: 200)
-s    product min size (default: 200)
-l    product max size (default: 1000)
-c    different site in primer (default: 2)
-e    different site in 3 end of primer (default: No)
