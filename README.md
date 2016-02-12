# GSP
A tool for Genome Specific Primers design in polyploid species

In a polyploid species, genes from the homeologous genomes exhibit a very high sequence similarity, especially in the coding regions. This makes it difficult to design genome-specific primers to amplify specific sequences from the individual genomes in the polyploid genome background. Development of genome-specific primers for important genes in polyploid species is very useful and critical not only for the study of sequence diversity and association mapping of genes in natural populations, but also for the development of gene-based functional markers for marker-assisted breeding. GSP is a powerful web-based platform for designing genome-specific primers that can distinguish the sequences among different genomes in a polypoid species.

Web server: http://probes.pw.usda.gov/GSP/

# Dependencies
There is number of additional dependencies not provided by GSP authors. Additional programs include:  
1. <b>bedtools:</b> a powerful toolset for genome arithmetic (http://bedtools.readthedocs.org/en/latest/).
2. <b>Muscle:</b> Multiple sequence alignment program (http://www.drive5.com/muscle/).
3. <b>Primer3:</b> program for designing PCR primers (http://primer3.sourceforge.net/).


# Installation
1. Download and upzip.
2. change the path to upzip folder.
3. cmake ./src
4. make

# Usage
<b>Usage:</b> GSP  
<b>-r</b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;blast table result path  
<b>-d</b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;sequence database path (blast database path)  
<b>-a</b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;fasta file path (for design specefic priemrs in multiple sequences only)  
<b>-b</b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;bed path  
<b>-m</b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;muscle path  
<b>-p</b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;primer3 path  
<b>-t</b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;primer3 parameters file path  
<b>-o</b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;output path  
<b>-f</b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;hit sequence flanking length (default: 200)  
<b>-s</b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;product min size (default: 200)  
<b>-l</b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;product max size (default: 1000)  
<b>-c</b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;different site in primer (default: 2)  
<b>-e</b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;different site in 3 end of primer (default: No)  
