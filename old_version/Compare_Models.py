import sys,os,re
import scipy,xlrd,xlwt
from time import time
import numpy as np
from math import *
from scipy import stats,polyfit
from scipy.optimize import minimize
from itertools import product
import cPickle as pickle
import copy_reg
import types
import multiprocessing as mp

sys.path.append('../../')
from DNAc import *

from Model_Interface import _call_model


class Compare_Models(object):
    
    use_MPI = False
    use_reported_start = True
    recalculate = False
    
    def __init__(self,models,datasets=None):
        
        handle = open('mRNA_database/mRNA_database.p','rb')
        self.db = pickle.load(handle)
        handle.close()
        
        handle = open('model_calcs.p','rb')
        self.model_calcs = pickle.load(handle)
        handle.close()
        
        if datasets is None:
            datasets = self.db.keys()
        
        if not all( [s in self.db for s in datasets] ):
            raise Exception('''At least one input dataset is not in the database.
                               Double check & run _import.py if needed.''')
        
        self.datasets = datasets
        self.models = models
    
    def set_options(self, use_MPI=False, use_reported_start=True, recalculate=False):
        self.use_MPI = use_MPI
        self.use_reported_start = use_reported_start
        self.recalculate = recalculate
    
    def run(self):
        
        inputList = self._define_inputList()
        
        if self.use_MPI:
            
            # MPI via mpi4py
            #from mpi4py import MPI
            #from MPI_pool import Pool
            #pool = Pool(MPI.COMM_WORLD)
            #pool.start()
            #calcs = pool.map(_call_model, inputList)
            #pool.stop()

            # Python's multiprocessing
            ncpus = mp.cpu_count()
            pool = mp.Pool(processes=ncpus,maxtasksperchild=1)
            t0 = time()
            output = pool.map(_call_model,inputList)
            tf = time()
            pool.close()
            pool.join()
            dt = (tf-t0)/60.0
            print "TOTAL RUN TIME = {} minutes".format(dt)

        else:
            output = map(_call_model,inputList)
            
        #=== Decompress data into appropriate lists ===#
        
        # Set up datastructure if:
        # (1) not previously calculated or (2) if recalculating
        for paper in self.datasets:
            if paper not in self.model_calcs.keys():
                self.model_calcs[paper] = {}            
            for model in self.models:
                if model not in self.model_calcs[paper].keys() or self.recalculate:
                    self.model_calcs[paper][model] = {}
        
        # Add data!            
        for i in range(len(inputList)):
            
            paper = inputList[i]['paper']
            model = inputList[i]['model']
            n = inputList[i]['subset']
            results_dic = output[i]
            k = self.db[paper]['num_subsets']
            
            # Save values
            for key in results_dic.keys():
                if key not in self.model_calcs[paper][model].keys():
                    self.model_calcs[paper][model][key] = [[] for _ in range(k)]
                self.model_calcs[paper][model][key][n].append(results_dic[key])
                
        #=== Calculate statistics ===#
        # Always recalculate stats!
        
        for paper,model in product(self.datasets,self.models):
            for n in range(self.db[paper]['num_subsets']):
                
                if model in ['RBS Calculator v1.0','UTR Designer','RBS Calculator v1.1','RBS Calculator v2.0','RBS Calculator v2.1','RBS_Calc_v2_UTR_Designer']:
                    self._calc_thermodynamic_stats((paper,n,model))
                elif model in ['RBS Designer']:
                    self._calc_efficiency_stats((paper,n,model))
                elif model in ['EMOPEC']:
                    self._calc_general_stats((paper,n,model))
                else:
                    print "WARNING: Model, {}, has not been programmed for stats calcs...".format(model)
                    pass
        
        
    def _define_inputList(self):
        
        # Add subsets & indexes for each mRNA
        def _add(inputList,paper,model):
            for n in range(self.db[paper]['num_subsets']):
                for i in range(self.db[paper]['subset_size'][n]):
                    
                    data = {
                    'paper': paper,
                    'subset': n,
                    'indx': i,
                    'model': model,
                    'mRNA': self.db[paper]['mRNA'][n][i],
                    'start_pos': self.db[paper]['start_pos'][n][i],
                    'temperature': self.db[paper]['temp'][n][i],
                    'rRNA': self.db[paper]['rRNA'][n][i]}
                    
                    inputList.append(data)
            return inputList
        
        #=== Set up input list of models and datasets for map ===#
        inputList = []
        for paper,model in product(self.datasets,self.models):
            
            #recalculate everything
            if self.recalculate: 
                inputList = _add(inputList,paper,model)
            
            #recalculate only if the paper or model is new
            elif paper not in self.model_calcs or model not in self.model_calcs[paper]:
                inputList = _add(inputList,paper,model)
        
        return inputList
        
    def _calc_thermodynamic_stats(self,var,fixed_beta=True):
        
        #=== Useful functions ====#
        
        calc_dG_apparent = lambda K,beta,rate: (log(rate) - log(K))/(-beta)
        calc_TIR = lambda K,beta,dG: K*exp(-beta*dG)
        
        # Linear least squares (LSQ) linear regression b/w dG_total & fluorescence
        def _LSQ_linreg(dG_total,Exp):
            (beta,logK) = polyfit(dG_total,np.log(Exp),1)
            beta = -beta
            K = exp(logK)
            return beta,K
        
        def _fit_K(dG_total,Exp):
            beta = 0.45
            LSQ = lambda K: np.sum( ( np.log(Exp)+beta*dG_total-K[0] )**2.0 )
            res = minimize(LSQ,2500.0,bounds=[(0.0,None)])
            K = np.exp(res.x[0])
            return beta,K
        
        paper,n,model = var
        
        Expression = np.array(self.db[paper]['Expression'][n]).astype(float)
        dG_total = np.array(self.model_calcs[paper][model]['dG_total'][n]).astype(float)
        
        # Remove nonpredicted sequences (1x10^12)
        indxs = [i for i in range(len(dG_total)) if dG_total[i]!=1e12]
        Expression = Expression[indxs]
        dG_total = dG_total[indxs]
        
        #=== Determine outliers with initial fit ===#
        
        if fixed_beta: beta,K = _fit_K(dG_total,Expression)
        else:          beta,K = _LSQ_linreg(dG_total,Expression)        
        
        dG_apparent = np.array([calc_dG_apparent(K,beta,E) for E in Expression])
        ddG = dG_total - dG_apparent
        ddG_abs = np.absolute(ddG)
        
        outliers = self.find_outliers(ddG_abs)
        keepers = np.invert(outliers)
        
        Exp_trimmed = Expression[keepers]
        dG_total_trimmed = dG_total[keepers]
        
        #=== Reapply fit with trimmed dataset ===#
        
        if fixed_beta: beta,K = _fit_K(dG_total_trimmed,Exp_trimmed)
        else:          beta,K = _LSQ_linreg(dG_total_trimmed,Exp_trimmed)
        
        dG_apparent = np.array([calc_dG_apparent(K,beta,E) for E in Expression])
        ddG = dG_total - dG_apparent
        ddG_abs = np.absolute(ddG)
        
        #=== Calculate correlation statistics ===#
        
        (rho,spearman_p) = stats.spearmanr(dG_total,dG_apparent)
        (R,pearson_p) = stats.pearsonr(dG_total,dG_apparent)
        
        #=== Calculate TIRs and TIR-fold change errors (fce) ===#
        
        # Get all Expression again
        Expression = np.array(self.db[paper]['Expression'][n]).astype(float)
        dG_total = np.array(self.model_calcs[paper][model]['dG_total'][n]).astype(float)
        
        dG_apparent = np.array([calc_dG_apparent(K,beta,E) for E in Expression])
        ddG = dG_total - dG_apparent
        ddG_abs = np.absolute(ddG)        
        
        TIR = [calc_TIR(K,beta,dG) for dG in dG_total]
        TIR_fce = Expression/np.array(TIR)
        
        #=== Save information ===#
        
        print paper,model,n,beta,K,R**2,pearson_p
        
        keyList = ['beta','K','R^2','spearman_rho','spearman_p','pearson_r','pearson_p',
                   'TIR','fold_change_error','ddG','ddG_abs','dG_apparent']
        
        for k in keyList:
            if k not in self.model_calcs[paper][model]:
                self.model_calcs[paper][model][k] = [None]*self.db[paper]['num_subsets']
        
        self.model_calcs[paper][model]['beta'][n]         = beta
        self.model_calcs[paper][model]['K'][n]            = K
        self.model_calcs[paper][model]['R^2'][n]          = R**2
        self.model_calcs[paper][model]['spearman_rho'][n] = rho
        self.model_calcs[paper][model]['spearman_p'][n]   = spearman_p
        self.model_calcs[paper][model]['pearson_r'][n]    = R
        self.model_calcs[paper][model]['pearson_p'][n]    = pearson_p
        
        self.model_calcs[paper][model]['TIR'][n]               = list(TIR)
        self.model_calcs[paper][model]['fold_change_error'][n] = list(TIR_fce)
        self.model_calcs[paper][model]['ddG'][n]               = list(ddG)
        self.model_calcs[paper][model]['ddG_abs'][n]           = list(ddG_abs)
        self.model_calcs[paper][model]['dG_apparent'][n]       = list(dG_apparent)
        
    def _calc_general_stats(self,var):
        
        #=== Useful functions ===#
        
        def _fit_K(apparent,predicted):
            fxn = lambda K: np.sum( ( apparent-K[0]*predicted )**2.0 )
            res = minimize(fxn,max(apparent),bounds=[(0.0,None)])
            K = res.x[0]
            return K
        
        paper,n,model = var
        
        apparent  = np.array(self.db[paper]['Expression'][n])
        predicted = np.array(self.model_calcs[paper][model]['Expression'][n])
        predicted[predicted==0] = 1e-12 # remove zeroes
        
        #=== Determine outliers with initial fit ===#
        
        K = _fit_K(apparent,predicted)
        predicted_scaled = K*predicted
        error = np.absolute( apparent / predicted_scaled )
        
        outliers = self.find_outliers(error)
        keepers = np.invert(outliers)
        
        # remove outlier removal strategy on FlowSeq data, not working well
        if len(predicted) > 1000:
            keepers = np.array( [1]*len(predicted) )
        
        apparent_trimmed = apparent[keepers]
        predicted_trimmed = predicted[keepers]
        
        #=== Reapply fit with trimmed dataset ===#
        
        K = _fit_K(apparent_trimmed,predicted_trimmed)
        predicted_scaled = K*predicted
        error = np.absolute( apparent / predicted_scaled )
        
        #=== Calculate correlation statistics ===#
        
        (rho,spearman_p) = stats.spearmanr(np.log10(predicted),np.log10(apparent))
        (R,pearson_p) = stats.pearsonr(np.log10(predicted),np.log10(apparent))        

        #=== Save information ===#
        
        print paper,model,n,K,R**2
        
        keyList = ['K','R^2','spearman_rho','spearman_p','pearson_r','pearson_p',
                   'fold_change_error','Expression']
        
        for k in keyList:
            if k not in self.model_calcs[paper][model]:
                self.model_calcs[paper][model][k] = [None]*self.db[paper]['num_subsets']
        
        self.model_calcs[paper][model]['K'][n] = K
        self.model_calcs[paper][model]['R^2'][n]          = R**2
        self.model_calcs[paper][model]['spearman_rho'][n] = rho
        self.model_calcs[paper][model]['spearman_p'][n]   = spearman_p
        self.model_calcs[paper][model]['pearson_r'][n]    = R
        self.model_calcs[paper][model]['pearson_p'][n]    = pearson_p
        
        self.model_calcs[paper][model]['Expression'][n] = list(predicted_scaled)
        self.model_calcs[paper][model]['fold_change_error'][n] = list(error)

    def _calc_efficiency_stats(self,var):
        
        #=== Useful functions ===#
        
        def _fit_K(apparent,efficiency):
            
            fxn = lambda K: np.sum( ( apparent - K[0]*efficiency )**2.0 )
            res = minimize(fxn,max(apparent),bounds=[(0.0,None)])
            K = res.x[0]
            return K
        
        paper,n,model = var
        
        apparent  = np.array(self.db[paper]['Expression'][n])
        efficiency = np.array(self.model_calcs[paper][model]['translation efficiency'][n])
        efficiency[efficiency==0] = 1e-12 # remove zeroes
        
        # Remove no prediction sequences (required for RBS Designer)
        indxs = [i for i in range(len(efficiency)) if efficiency[i]!=-1]
        
        calculated = []
        for x in efficiency:
            if x!=-1: calculated.append(1)
            else: calculated.append(0)
        
        apparent = apparent[indxs].astype(float)
        efficiency = efficiency[indxs].astype(float)
        
        #=== Determine outliers with initial fit ===#
        
        K = _fit_K(apparent,efficiency)
        predicted = K*efficiency
        error = np.absolute( apparent / predicted )
                
        outliers = self.find_outliers(error)
        keepers = np.invert(outliers)
        
        apparent_trimmed = apparent[keepers]
        efficiency_trimmed = efficiency[keepers]
        
        #=== Reapply fit with trimmed dataset ===#
        
        K = _fit_K(apparent_trimmed,efficiency_trimmed)
        predicted = K*efficiency
        error = np.absolute( apparent / predicted )
        
        #=== Calculate correlation statistics ===#
        
        (rho,spearman_p) = stats.spearmanr(np.log10(predicted),np.log10(apparent))
        (R,pearson_p) = stats.pearsonr(np.log10(predicted),np.log10(apparent))        

        #=== Save information ===#
        
        print paper,model,n,K,R**2
        
        keyList = ['K','R^2','spearman_rho','spearman_p','pearson_r','pearson_p',
                   'fold_change_error','Expression','calculated']
        
        for k in keyList:
            if k not in self.model_calcs[paper][model]:
                self.model_calcs[paper][model][k] = [None]*self.db[paper]['num_subsets']
        
        self.model_calcs[paper][model]['K'][n] = K
        self.model_calcs[paper][model]['R^2'][n]          = R**2
        self.model_calcs[paper][model]['spearman_rho'][n] = rho
        self.model_calcs[paper][model]['spearman_p'][n]   = spearman_p
        self.model_calcs[paper][model]['pearson_r'][n]    = R
        self.model_calcs[paper][model]['pearson_p'][n]    = pearson_p
        
        efficiency = self.model_calcs[paper][model]['translation efficiency'][n]
        efficiency = [-1 if v==None else v for v in efficiency[:]]
        Expression = list(K*np.array(efficiency))
        Expression = [v if v>=0 else -1 for v in Expression[:]]
        
        all_error = []
        j = 0
        for i in calculated:
            if i==1:
                all_error.append(error[j])
                j+=1
            else:
                all_error.append(-1)
        
        self.model_calcs[paper][model]['Expression'][n]        = Expression
        self.model_calcs[paper][model]['fold_change_error'][n] = list(all_error)
        self.model_calcs[paper][model]['calculated'][n]     = list(calculated)
        
    def find_outliers(self, vals, threshold=2):
        ''' Uses a median-absolute-deviation (MAD) test to identify outliers
        Input is a Python list of values (ddG)
        Output: a ddG-list-length boolean list (TRUE = outlier)'''
        
        median = np.median(vals)
        diff = np.absolute(vals - median)
        MAD = np.median(diff)
        score = 0.6745 * diff / MAD
        outliers = score > threshold
        
        return outliers
        

def main(models,datasets):
    '''Run [models] on [datasets] using Compare_Models class'''
    
    CM = Compare_Models(models,datasets)
    CM.set_options(use_MPI=True,recalculate=True)
    
    # t0 = time()
    CM.run()
    # tf = time()
    # dt = (tf-t0)/60.0
    # print "TOTAL RUN TIME = {} minutes".format(dt)
    
    # Open existing pickle file and save data to it
    # pkl_file = open('model_calcs.p','rb')
    # my_dict = pickle.load(pkl_file)
    # pkl_file.close()
    
    # model_calcs = CM.model_calcs
    # for paper in datasets:
    #     for m in models:
    #         my_dict[paper][m] = model_calcs[paper][m]
    # handle = open('model_calcs.p','wb')
    # pickle.dump(my_dict,handle)
    # handle.close()
    
    # Write completely brand new pickle file
    # handle = open('fast_model_calcs.pkl','wb')
    # model_calcs = CM.model_calcs
    # pickle.dump(model_calcs,handle)
    # handle.close()
    
if __name__ == "__main__":
    
    # All models
    models = ['RBS Calculator v1.0',
              'RBS Calculator v1.1',
              'RBS Calculator v2.0',
              'RBS Calculator v2.1',
              'UTR Designer',
              'RBS Designer',
              'EMOPEC']
    
    #models = ['RBS_Calc_v2_UTR_Designer']
    models = ['RBS Calculator v2.1']
    
    # Temporary for model hypothesis testing
    # handle = open('hypotheses.pkl','wb')
    # model_calcs = {}
    # pickle.dump(model_calcs,handle)
    # handle.close()
    
    # Individually characterized sequences
    datasets = ['EspahBorujeni_NAR_2013',
                'EspahBorujeni_NAR_2015',
                'EspahBorujeni_JACS_2016',
                'EspahBorujeni_Footprint',
                'Salis_Nat_Biotech_2009',
                'Farasat_MSB_2014',
                'Tian_NAR_2015',
                'Mimee_Cell_Sys_2015',
                'Bonde_NatMethods_IC_2016']

    # Flow-seq datasets run and save separately
    #datasets = ['Kosuri_PNAS_2013',
                #'Goodman_Science_2013']
    
    main(models,datasets)
