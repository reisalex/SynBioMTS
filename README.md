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
Install with the following:
```
git clone https://github.com/reisalex/test-sfm
cd test-sfm
sudo python setup.py install
```
The model test system can then be imported by:
```
import testsfm
```

### Usage
If you would like to use the provided genetic system database, the best way is to navigate to /testsfm, and run the datbase initialization module (initdb.py):
```
cd /testsfm
python initdb.py
```

To use the model test system, you have to do the following:
1. Wrap the model with a Python function.
2. Create a models `Container` object and pass the wrapped functions with the `add` method.
3. Specify the functional form between the model predictor and the system function with the `setform` method.
4. Create a `ModelTest` object and pass the models Container.
5. Run model calculations and statistics with `run`.

```
import testsfm

# Wrap the model with a function
import RBS_Calculator_v2
def RBSCalcv2(sequence,temperature):
    rRNA = 'ACCTCCTTA'
    model = RBS_Calculator_v2.RBS_Calculator(mRNA=sequence,rRNA=rRNA)
    model.temp = temperature
    model.run()
    RBS = model.output() # simplified for the example
    
    # Results should be returned as a dictionary
    # The keys will become labels in the resulting pandas dataframe
    results = {
    'TIR': RBS.tir,
    'dG_total': RBS.dG_total,
    'dG_mRNA_rRNA': RBS.dG_mRNA_rRNA',  
    'dG_mRNA': RBS.dG_mRNA
    }
    return results

if __name__ == "__main__":
    
    # create models Container object
    models = testsfm.interface.Container()
    
    # add the model(s)
    models.add(RBSCalcv2)
    
    # specify the form of each model
    # RBS_Calculator is a thermodynamic model where: Protein ~ K*exp(-0.45*dG_total)
    models.setform(['RBSCalcv2'],x='dG_total',y='PROT.MEAN',yScale='ln',a1=-0.45)
    
    # create test system object
    testsystem = testsfm.analyze.ModelTest(models,'geneticsystems.db')
    
    # run model predictions and statistics calculations 
    testsystem.run()
    
    # if you want to shelve the model calculations
    # testsystem.run(calcsFilename='savedcalcs.db')
```

You can specify filters to run predictions on a subset of genetic systems with shared properties:
```
filters = { 'ORGANISM': ['Escherichia coli'],
            ''DATASET': ['Beck_PLoS_2016']
}
testsystem = testsfm.analyze.ModelTest(models,'geneticsystems.db',filters)
```

You can specify the number of processes to force single process model calculations. The model test system uses the number of available CPUs by default.
```
testsystem = testsfm.analyze.ModelTest(models,'geneticsystems.db',nprocesses=1)
```

If you simply want to run model predictions, you can use `predict` instead of `run` which includes statistics:
```
testsystem.predict()
```

See /examples for more detailed examples.

### Statistics
By default, the model test system will run statistics assuming the model predictor and the system function share a linear relationship. Specifically, `analyze.statistics()` calls a custom function `linear_complete` from the stats module to compute the following:
* Fitted slope and y-intercept with outliers removed (via MAD method)
* Relative model error (Apparent Value/Predicted Value)
* Pearson & Spearman correlation coefficients
* One-sided model error cummulative distribution function
* Kullback-Leibler divergence

If you want to run futher statistics, you can import the stats module:
```
from testsfm import stats
```
You can always add additional stats functions as needed.

### Exporting
Export to Excel is as simple as:
```
test.to_excel('filename')
```
By default pandas exports with the labels (columns) alphabetized. The model test system overrides the default export if you specify the labels. See /examples/labels for the ones I use:
```
test.run()

with open("labels/labels1.txt","r") as f:
    predictLabels = [x.strip('\n') for x in f.readlines()]

with open("labels/labels_stats.txt","r") as f:
    statsLabels = [x.strip('\n') for x in f.readlines()]

test.to_excel('Beck',predictLabels,statsLabels) 
```

## Acknowledgements

Thanks to Howard M Salis (Penn State), Iman Farasat (Merck), Amin Espah Borujeni (MIT), Tian Tian (JBEI), Daniel Goodman (Harvard), Sri Kosuri (UCLA), Robert Egbert (Berkeley), Mark Mimmee (MIT), and Heather Beck (Vienna) for publishing/providing high quality characterization data. A special thanks to Daniel Goodman for discussion on Flow-seq and for providing additional information on the 2013 Flow-seq datasets.

If you use test-sfm, please cite:

Alexander C. Reis, and Howard M. Salis An automated model test system for systematic development and improvement of gene expression models, In Preparation (2017).
