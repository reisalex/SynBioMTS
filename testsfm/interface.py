"""
Model interface to import models, write model wrappers, collect model output
and run statistics specific to model type

Copyright 2017 Alexander C. Reis, Howard M. Salis, all rights reserved.
  
"""

import stats
from functools import partial

class Models(dict):
    __name__ = "Models"

    def __init__(self):
        pass

    def register(self, alias, model, *args, **kargs):
        '''Register a model with the interface under the name *alias*. You
        may provide default arguments that will be passed automatically when
        calling the registered model. Fixed arguments can then be overriden
        at model call time.
        Inputs:
        alias (string) = The name the operator will take in the interface. If
                         the alias already exist it will overwrite the the
                         operator already present.
        model (object) = The model to which refer the alias.
        argument       = One or more argument (and keyword argument) to pass
                         automatically to the registered model when called,
                         optional.

        The following code block is an example of how the toolbox is used. ::
            >>> def RBS_Calculator_v2_0(sequence, organism, temp):
            ...     # run model code here
            ...     print seq, organism, temperature
            ...
            >>> models = Interface()
            >>> models.register("RBSCalc2",RBS_Calculator_v2_0,'ACTAGC',temp=37.0)
            >>> models.RBSCalc2('Escherichia coli')
            'ACTAGC' 'Escherichia coli' 37.0

        The registered model will be given the attributes :attr:`__name__`
        set to the alias and :attr:`__doc__` set to the original model's
        documentation. The :attr:`__dict__` attribute will also be updated
        with the original model's instance dictionary, if any. '''

        # Adapted from DEAP's base.Toolbox
        # https://github.com/DEAP/deap

        pmodel = partial(model, *args, **kargs)
        pmodel.__name__ = alias
        pmodel.__doc__ = model.__doc__

        if hasattr(model, "__dict__") and not isinstance(model, type):
            # Some functions don't have a dictionary, in these cases
            # simply don't copy it. Moreover, if the model is actually
            # a class, we do not want to copy the dictionary.
            pmodel.__dict__.update(model.__dict__.copy())

        if not pmodel.args:
            pmodel.co_varnames = model.__code__.co_varnames
        else:
            pmodel.co_varnames = pmodel.args

        # print pmodel.func
        # print model.__code__.co_varnames
        print pmodel.co_varnames
        # print pmodel.args
        print pmodel.keywords

        setattr(self, alias, pmodel)
        self[alias] = pmodel

    def unregister(self, alias):
        '''Unregister *alias* from the model interface.
        alias (string) = The name of the operator to remove from the interface.
        '''
        delattr(self, alias)
        self.pop(alias)

    def _listed(self):
        return sorted(self.keys())

    def __add__(self,another):
        assert another.__name__ == "Models", "Must add two Models to combine."
        for alias,pmodel in another.iteritems():
            self.register(alias,pmodel)
        return self

    def __sub__(self,another):
        assert another.__name__ == "Models", "Must subtract two Models to remove."
        for alias,pmodel in another.iteritems():
            self.unregister(alias,pmodel)
        return self


if __name__ == "__main__":
    pass