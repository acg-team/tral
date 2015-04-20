[<img src="https://img.shields.io/pypi/v/tral.svg?branch=master">](https://pypi.python.org/pypi/tral)

# Tandem Repeat Annotation Library

TRAL is a highly modularized, flexible sequence tandem repeats annotation Python2/3 library.

  - Large scale annotation of tandem repeats with *de novo* detectors, and sequence profile models
  - Statistical significance testing, overlap detection, and filtering of annotations
  - Refinement of tandem repeat annotations with circular profile hidden Markov models
  - User-defined output formats

The source code is [documented on GitHub IO].

### Version
0.3.4


### Installation

TRAL is available on [Pypi] and can be installed with [pip] for Python>=3.2:

```sh
$ pip install tral
```

See also more extensive [Installation instructions].


### License

GPL


### Dependencies

Some of TRAL's functions depend on external software ([Installation instrutions for dependencies]). This includes creation of sequence profile hidden Markov models, alignment of tandem repeat units, and *de novo* repeat detection.



[documented on GitHub IO]:http://elkeschaper.github.io/tral/
[Installation instructions]:http://elkeschaper.github.io/tral/install.html#install
[Installation instructions for dependencies]:http://elkeschaper.github.io/tral/install_external.html#install-external
[Pypi]:https://pypi.python.org/pypi
[pip]:https://pip.pypa.io/en/latest/
