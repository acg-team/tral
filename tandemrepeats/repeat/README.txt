# (C) 2012 Elke Schaper

Scoring module tandem repeats


This module provides the basic functions to calculate a range of scores, pvalues and divergences on tandem repeat muliple sequence alignments as described in 
1. Schaper,E., Kajava,A.V., Hauser,A. and Anisimova,M. (2012) Repeat or not repeat?--Statistical validation of tandem repeat prediction in genomic sequences. Nucleic Acids Res., 10.1093/nar/gks726.

#######################
Requirements:
#######################

Python 3

Modules:
bisect
collections
copy
csv
re
logging
math
numpy
scipy
scipy.linalg
scipy.special
scipy.stats

The requirement can be tested by running:
<path>/python3 <path>/test_dependencies.py


######################
USAGE
######################

1. Set paths.py
DATAROOT: Points to data-directory, containing score distribution files for pvalue calculation, and paml substitution matrix files.
RESULTROOT: Points to results-directory, where results should be saved.

2. Set repeat_main.py
All other parameters are directly set within the repeat_main.py file.
Details can be found within that file.

3. Execute repeat_main.py
<path>/python3 repeat_main.py


 





