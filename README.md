
[<img src="https://img.shields.io/pypi/v/tral.svg?branch=master">](https://pypi.python.org/pypi/tral)

# Tandem Repeat Annotation Library

TRAL is a highly modularized, flexible sequence tandem repeats annotation Python2/3 library.

- Large scale annotation of tandem repeats with *de novo* detectors, and sequence profile models  
- Statistical significance testing, overlap detection, and filtering of annotations  
- Refinement of tandem repeat annotations with circular profile hidden Markov models  
- User-defined output formats  

The source code is [documented on GitHub IO].

### Version

2.0

### Installation

It is recommended to use TRAL with the [docker image] or install it locally with the [easy_setup] system.

### License

GPL-2.0

### Dependencies

Some of TRAL's functions depend on external software ([Installation instructions for dependencies]). This includes creation of sequence profile hidden Markov models, alignment of tandem repeat units, and *de novo* repeat detection.

### Citation

Delucchi, M., Näf, P., Bliven, S., & Anisimova, M. (2021). [TRAL 2.0: Tandem Repeat Detection with Circular Profile Hidden Markov Models and Evolutionary Aligner](https://www.frontiersin.org/articles/10.3389/fbinf.2021.691865). *Frontiers in Bioinformatics*, DOI:  10.3389/fbinf.2021.691865

E Schaper, A Korsunsky, J Pecerska, A Messina, R Murri, H Stockinger, S Zoller, I Xenarios, and M Anisimova (2015). [TRAL: Tandem Repeat Annotation Library](http://bioinformatics.oxfordjournals.org/content/early/2015/05/17/bioinformatics.btv306.abstract). *Bioinformatics*. DOI:  10.1093/bioinformatics/btv306

[documented on GitHub IO]:https://acg-team.github.io/tral/
[docker image]:https://github.com/acg-team/tral/packages
[easy_setup]:https://github.com/acg-team/tral/tree/develop/easy_setup
[Installation instructions for dependencies]:https://acg-team.github.io/tral/install_external.html#install-external
[Pypi]:https://pypi.python.org/pypi
[pip]:https://pip.pypa.io/en/latest/
