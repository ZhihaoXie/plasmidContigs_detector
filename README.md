# plasmidContigs\_detector

bacteria plasmid detect.


# Installation

`git clone git@github.com:ZhihaoXie/plasmidContigs_detector.git`


# Requests

+ Python3
+ Biopython, numpy and matplotlib
+ Perl and Bioperl
+ cBar: http://csbl.bmb.uga.edu/~ffzhou/cBar/
+ abricate: https://github.com/tseemann/abricate
+ Interproscan: https://www.ebi.ac.uk/interpro/
+ Glimmer: http://ccb.jhu.edu/software/glimmer/index.shtml
+ Prodigal: http://prodigal.ornl.gov/


# Usage

```
python3 plasmidContigs_detector.py [options]

options:
  -f CONTIGS   contigs fasta
  -o OUTDIR    the output dir
  -x PREFIX    the output file prefix
  -p PFAM      the pfam of plasmid copy family, example plasmid_pfam.id
  -c CBAR      the cBar path, you can get it from
               http://csbl.bmb.uga.edu/~ffzhou/cBar/
  -a ABRICATE  the abricate path, you can get it from
               https://github.com/tseemann/abricate
  -i IPR       Interproscan path, get it from https://www.ebi.ac.uk/interpro/

```

# Copyright

Zhihao Xie,  xiezhihao1122@outlook.com

