# test-sfm
An automated model test system for DNA sequence-function models.

Test-sfm uses a database of 16000 unique characterized genetic systems to run Python-wrapped, sequence-function models, quantify model accuracies, accept or reject proposed mechanistic hypotheses, and identify sources of model error. This package is easily modifiable to expand the genetic system database, calculate additional statistical test metrics, and test new and improved gene expression models implemented in nearly any language.

## Getting Started

### Dependencies
Python packages used are listed below. You can install the first three packages as a part of the [SciPy Stack](https://www.scipy.org/install.html).
* pandas - Database management
* scipy - Statistics calculations 
* numpy - General purpose numerical computing
* [scikit-learn](http://scikit-learn.org/stable/install.html) - Machine learning

[ViennaRNA](https://www.tbi.univie.ac.at/RNA/) - A C code library for the prediction and comparison of RNA secondary structures. ViennaRNA is wrapped with /models/PyVRNA.py for use in modeling and machine learning analysis.

### Installing
Install with the following command:
```
sudo python setup.py install
```
The model test system can then be imported by:
```
import testsfm
```

### Usage
If you would like to use the provided genetic system database, the best way is to clone the complete repository, navigate to /testsfm, and run the datbase initialization module (initdb.py). In Linux:
```
git clone https://github.com/reisalex/test-sfm
cd /testsfm
python initdb.py
```

See /examples for usage. 

Please cite:

Alexander C. Reis, and Howard M. Salis
An automated model test system for systematic development and
improvement of gene expression models, In Preparation (2017).

