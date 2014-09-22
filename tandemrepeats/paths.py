import os

# Where are the executables stored? E.g. the tandem repeat detection algorithms?
EXECROOT = "/usr/local/bin"
# Where is the code located?
CODEROOT = "/I/have/stored/the/code/here/TandemRepeats"
# Where are the data-files, e.g. for pValue calculation stored?
DATAROOT = os.path.join(CODEROOT, "data")
# (Needed?) Where shall results be stored?
RESULTROOT = "/Save/my/results/here"
# (Needed?) What is the root of the file system?
ROOT = "/Root/of/my/filesystem"
