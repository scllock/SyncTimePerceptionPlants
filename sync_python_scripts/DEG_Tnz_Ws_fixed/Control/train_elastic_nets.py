import numpy as np
import pandas as pd
from sklearn.linear_model import ElasticNetCV, LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import pickle

from multiprocessing import Pool
import functools

import argparse

def custom_loocv_single(
    sample_pos: int,
    sample_names: list[str],
    input_pheno: np.array,
    input_expr: pd.DataFrame,
    model_fn,
    verbose = True
):
    """Truncates the dataset based on inputs then returns the output of a single
    fit from model_fn.
    """
    
    sample = sample_names[sample_pos]
    if verbose: print(f'Performing fitting without sample {sample}')

    # Remove the chosen sample from the list, then create the dataset
    remaining_positions = list(range(len(sample_names)))
    remaining_positions.remove(sample_pos)
        
    # make truncated datasets
    trunc_pheno = input_pheno[remaining_positions]
    trunc_expr = input_expr.iloc[:, remaining_positions]

    y_values = trunc_pheno
    x_values = np.transpose(trunc_expr.values)

    # Check for NaN values in x_values and y_values
    if np.any(np.isnan(x_values)):
        raise ValueError(f"NaN values found in x_values for sample {sample}")
    if np.any(np.isnan(y_values)):
        raise ValueError(f"NaN values found in y_values for sample {sample}")


    # get the output of the fitting (eg ---.fit(x_values, y_values))
    return(model_fn(x_values, y_values))

def run_model_with_parallel_loocv(
    input_expr: pd.DataFrame,
    input_pheno: np.array,
    model_fn,
    verbose = True
) -> dict:
    """Returns a dictionary with the model results, where the keys are IDs of 
    each sample and values are the outputs of model_fn.

    Args:
        input_expr (pd.DataFrame): represents the gene expression values.
        input_pheno (np.array): represents the phenotypes to be predicted - eg floats
            for regression or bools for classification.
        model_fn: A function which will take an array of expression values and an
            array of phenotypes, and return an arbitrary modelling output.
        verbose (bool, optional): whether to print steps. Defaults to True.
    """

    # validate whether there are the same number of samples as phenotypes
    sample_names = input_expr.columns.values.tolist()
    if not (len(input_pheno) == len(sample_names)):
        raise(ValueError("Length of samples and phenotypes must be equal"))

    results_dict = dict.fromkeys(sample_names)
    range_of_samples = range(len(sample_names))

    # given all the sample indices, run the loocv parallelised
    pool = Pool()
    result_async = pool.map_async(
        functools.partial(  # use this to keep data inputs fixed
            custom_loocv_single,
            sample_names = sample_names,
            input_pheno = input_pheno,
            input_expr = input_expr,
            model_fn = model_fn,
            verbose = verbose
        ), 
        range_of_samples
    )  # this will take time
    result_async.wait() # until completion
    
    list_of_results = result_async.get()
    
    # Cast this into a results dictionary
    for new_key, new_val in zip(results_dict.keys(), list_of_results):
        results_dict[new_key] = new_val
    
    return results_dict


def elastic_net_model_fn(
    x_values: np.array,
    y_values: np.array,
    l1_list = [.005,.01,.015,.02, .03, .04, .05, .06, .07, .08, .09, .095, .099, .1, .15, .2, .25, .3, .35, .4, .45, .5, .55, .6, .65, .7, .75, .8, .85, .9, .95, .99, 1],
    max_iter = 1000  # Increase the number of iterations
) -> ElasticNetCV:

    regr = ElasticNetCV(
        l1_ratio = l1_list,
        cv = 5,  # automatically chooses set to 5-fold validation
        random_state = 123,  # ensures a repeatable outcome
        verbose = 1,
        n_jobs = 1,  # do NOT parallelise the 5-fold CV
        max_iter = max_iter  # Set max iterations here
    )


    # fit with CV
    return regr.fit(x_values, y_values)


# regularised logistic regressions with enet
def logistic_model_fn(
    x_values: np.array,
    y_values: np.array,  # dtype must be bool
    l1_list =  [.1, .5, .7, .9, .95, .99, 1]
) -> LogisticRegressionCV:

    regr = LogisticRegressionCV(
        Cs = 10,
        l1_ratios = l1_list,
        penalty = 'elasticnet',
        cv = 5,  # automatically chooses set to 5-fold validation
        random_state = 123,
        solver = 'saga',
        # scoring = 'neg_log_loss',  # use accuracy instead
        max_iter = 10**3,
        verbose = 0, # too many output lines to be useful if verbose = 1!
        n_jobs = 1 # do NOT parallelise the 5-fold CV
    )

    return regr.fit(x_values, y_values)


if __name__ == "__main__":

    # load the datasets
    expr_df = pd.read_csv("data/Zscore_DEG_Ws_Tnz_control.csv", index_col = 0)
    pheno_df = pd.read_csv("data/phenos_to_predict_SL.csv", index_col = 0)


    parser = argparse.ArgumentParser()
    parser.add_argument("-f", action="store_true")
    args = parser.parse_args()

    if (args.f):
        # run linear regression with only fixed l1_ratio
        print("Running elastic nets with fixed l1_ratio ...")

        print("Flowering time mean with all genes. l1_ratio = 0.5")
        meanft_enet_output = run_model_with_parallel_loocv(
            expr_df,
            pheno_df["FT_mean"],
            functools.partial(
                elastic_net_model_fn,
                l1_list = [.06]
            )
        )
        pickle.dump(meanft_enet_output, open("outputs/meanft_enet_05.sav", "wb"))

        
        print("Flowering time standard error size with all genes, l1_ratio = 0.5")
        ftse_enet_output = run_model_with_parallel_loocv(
            expr_df,
            pheno_df["FT_se"],
            functools.partial(
                elastic_net_model_fn,
                l1_list = [.5]
            )
        )
        pickle.dump(ftse_enet_output, open("outputs/ftse_enet_05.sav", "wb"))

       
        print("mean number of leaves with all genes, l1_ratio = 0.5")
        leavesmean_enet_output = run_model_with_parallel_loocv(
            expr_df,
            pheno_df["leaves_mean"],
            functools.partial(
                elastic_net_model_fn,
                l1_list = [.15]
            )
        )
        pickle.dump(leavesmean_enet_output, open("outputs/leavesmean_enet_05.sav", "wb"))

        print("standard error number of leaves with all genes, l1_ratio = 0.5")
        leavesse_enet_output = run_model_with_parallel_loocv(
            expr_df,
            pheno_df["leaves_se"],
            functools.partial(
                elastic_net_model_fn,
                l1_list = [.5]
            )
        )
        pickle.dump(leavesse_enet_output, open("outputs/leavesse_enet_05.sav", "wb"))

       
       
        print("mean biomass with all genes, l1_ratio = 0.5")
        biomassmean_enet_output = run_model_with_parallel_loocv(
            expr_df,
            pheno_df["weights_mean"],
            functools.partial(
                elastic_net_model_fn,
                l1_list = [.005]
            )
        )
        pickle.dump(biomassmean_enet_output, open("outputs/biomassmean_enet_05.sav", "wb"))

        print("standard error of biomass with all genes, l1_ratio = 0.5")
        biomassse_enet_output = run_model_with_parallel_loocv(
            expr_df,
            pheno_df["weights_se"],
            functools.partial(
                elastic_net_model_fn,
                l1_list = [.5]
            )
        )
        pickle.dump(biomassse_enet_output, open("outputs/biomassse_enet_05.sav", "wb"))

       


    else: 
        # run linear regression in all possible combos
        print("Running elastic nets with varying l1_ratio ...")

        print("Mean flowering time with all genes")
        meanft_enet_output = run_model_with_parallel_loocv(
            expr_df,
            pheno_df["FT_mean"],
            elastic_net_model_fn
        )
        pickle.dump(meanft_enet_output, open("outputs/meanft_enet.sav", "wb"))

        print("standard error flowering time with all genes")
        ftse_enet_output = run_model_with_parallel_loocv(
            expr_df,
            pheno_df["FT_se"],
            elastic_net_model_fn
        )
        pickle.dump(ftse_enet_output, open("outputs/ftse_enet.sav", "wb"))

        print("Mean number of leaves with all genes")
        leavesmean_enet_output = run_model_with_parallel_loocv(
            expr_df,
            pheno_df["leaves_mean"],
            elastic_net_model_fn
        )
        pickle.dump(leavesmean_enet_output, open("outputs/leavesmean_enet.sav", "wb"))

        print("standard error of number of leaves with all genes")
        leavesse_enet_output = run_model_with_parallel_loocv(
            expr_df,
            pheno_df["leaves_se"],
            elastic_net_model_fn
        )
        pickle.dump(leavesse_enet_output, open("outputs/leavesse_enet.sav", "wb"))

        print("Mean biomass with all genes")
        biomassmean_enet_output = run_model_with_parallel_loocv(
            expr_df,
            pheno_df["weights_mean"],
            elastic_net_model_fn
        )
        pickle.dump(biomassmean_enet_output, open("outputs/biomassmean_enet.sav", "wb"))

        print("Standard error biomass with all genes")
        biomassse_enet_output = run_model_with_parallel_loocv(
            expr_df,
            pheno_df["weights_se"],
            elastic_net_model_fn
        )
        pickle.dump(biomassse_enet_output, open("outputs/biomassse_enet.sav", "wb"))

