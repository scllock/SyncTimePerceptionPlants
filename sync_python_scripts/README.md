# The genetic basis for synchronised time perception in plant populations - Python scripts

This repository contains the code used to train elastic net models on the gene expression, variant, and physiology data sets. Based from https://github.com/stressedplants/SinglePlantOmics/tree/main

## Folder structure

This folder contains this README and a Conda environment file `flowering_variability_env.yml`, which can be used to reproduce the conda environment used for this project.


-'DEG_Tnz_ws_fixed' contains sub folders for each of the conditions 'Control' and 'Treatment'. Each of these folders then contain: the Python file used to train the models on gene expression data 'train_elastic_nets.py' (Also is bash script and log of script running); the input data files under folder 'data'; and outputs from the training procedure in 'outputs'. (Note that some of the outputs contain the suffix `_smaller_save`, since they have been reduced in size by the analysis jupyter notebook).
- `analysis` contains: a Jupyter notebook to create figures; and, all of the figures and supplemental data produced by the analysis.

-'DEG_Tnz_ws_vary' contains sub folders for each of the conditions 'Control' and 'Treatment'. Each of these folders then contain: the Python file used to train the models on gene expression data 'train_elastic_nets.py' (Also is bash script and log of script running); the input data files under folder 'data'; and outputs from the training procedure in 'outputs'. (Note that some of the outputs contain the suffix `_smaller_save`, since they have been reduced in size by the analysis jupyter notebook).
- `analysis` contains: a Jupyter notebook to create figures; and, all of the figures and supplemental data produced by the analysis.




## Reproducing the Conda environment

Please see [this guide](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#activating-an-environment) for how to use the `.yml` file.

## Training the models 

Data and scripts are included to train the models. They can be found in each 'Control' or 'Treatment' folder. Make sure your working directory is set to the correct folder before running. 


To train models on the gene expression data, *with varying values for `l1_ratio`*, you should run the command `python train_elastic_nets.py`. To train models with *fixed* `l1_ratio` values, you should run `python train_elastic_nets.py -f`. Note, the fixed values for l1_ratio have been hard coded, so they will need to be edited in the code if you want to change them. The training parallelizes the different LOOCV folds. To make easier the two inital sub folders were used 'DEG_Tnz_Ws_fixed' for fixed l1 ratio and 'DEG_Tnz_Ws_vary' for varied l1 ratio.

For the paper, these commands were run on the [Viking computing cluster](https://vikingdocs.york.ac.uk/). For the most time consuming stage (running models on gene expression with varying `l1_ratio`), we used 16 cores and approximately 200GB of RAM, for a total run time of about 1 hour.

## Analysing the outputs

A Jupyter notebook has been included to create the plots for the paper. After installing the conda environment, you should be able to run the command `jupyter notebook` and open this folder in your web browser. Make sure the environment has access to both the `Control ` and `Treatment' folders. 

*Note: the logistic regression outputs have been reduced in size by running `remove_coefs_from_dict` in `analysis_elastic_nets.ipynb`, which removes the large `coefs_paths_` attribute.*
