# Titer model

## About the titer model  

Antigenic distances between pairs of dengue viruses are experimentally measured via neutralization titers. These titer values are prone to noise, and there is a limited amount of available titer data. If the antigenic heterogeneity observed in the raw data is truly the result of an underlying evolutionary process, we expect that changes in antigenic phenotype correspond to underlying changes in viral genotype.

The titer model maps changes in antigenic phenotype (titer drops) to specific branches in the viral phylogeny as described below. This allows us to directly quantify the extent to which observed phenotypic variation is explained by an underlying genetic evolutionary process. We can also use variations of the model formulation to directly compare competing hypotheses about the nature of dengue antigenic evolution.

The titer model pipeline has three main steps.
First, we build a phylogeny of dengue virus sequences to establish the genetic relationships between viruses.
Next, we infer how much antigenic change has occurred along each branch of the phylogeny by mapping titer changes to individual branches.
This assigns each branch _b_ an antigenic distance _d_<sub>_b_</sub>.
With this in hand, we estimate the antigenic distance between all pairs of viruses by tracing the path between them in the phylogeny, summing branch-specific distances $d_b$ as we go.

To learn these values of _d_<sub>_b_</sub>, we first split our dataset into training (random 90% of measurements) and test data (the remaining 10% of values).
We take the training data and fit _d_<sub>_b_</sub> for each branch in the tree, subject to regularization.
Parsimoniously, we expect that antigenic change is more likely to occur through larger changes on a few branches than through small changes on many branches; correspondingly, our prior expectation of values of _d_<sub>_b_</sub> is exponentially distributed such that most values of _d_<sub>_b_</sub> = 0.
This is analogous to lasso regression to identify a few parameters with positive weights and set other parameters to 0.
Additionally, some viruses have greater binding avidity, and some sera are more potent than others; these 'row' and 'column' effects, respectively, are normally distributed and are taken into account when estimating titers.
The model uses convex optimization to learn the values of _d_<sub>_b_</sub> that minimize the sum of squared errors (SSE) between observed and predicted titers in the training data.

## Model variations

![Titer model variation](https://raw.githubusercontent.com/blab/dengue-antigenic-dynamics/master/figures/png/titer_model_performance.png)

**A.** The 'interserotype model' only allows branches that lie between serotypes to contribute to antigenic evolution.
All other branches are assigned _d_<sub>_b_</sub> = 0.
**B.** The 'full tree model' allows any branch in the phylogeny to contribute to antigenic evolution (_d_<sub>_b_</sub> &ge; 0).
**C,D.** Predictive performance of each model on the test dataset (aggregated from 10-fold cross-validation).

## Running the titer model

This model was originally published in [Neher et al (PNAS, 2016)](http://dx.doi.org/10.1073/pnas.1525578113) and implemented as part of the Nextstrain [augur package](https://github.com/nextstrain/augur). The relevant portions of the repository have been reproduced here under [`implementation-nextstrain-augur`](implementation-nextstrain-augur/). Documentation, including install instructions, for the full augur pipeline can be found [here](https://github.com/nextstrain/augur/tree/6d9f7088d8792196e5021c67b876d9de1d2a13dd).

### 1 - Prepare input files  

`dengue.prepare.py` handles subsampling and other basic dataset config. You can edit the function `make_config` here with any desired changes.  

```
cd dengue-antigenic-dynamics/titer_model/implementation-nextstrain-augur/builds/dengue/  
python dengue.prepare.py
mv ./prepared/dengue_all.json ./prepared/dengue_config.json
```

**NB: For this dataset, step 1 has been run for you; see `./dengue_config.json`**

### 2 - Run the titer pipeline (and estimate clade frequencies)  

`dengue.process.py` handles the actual analysis; parameter settings, etc. can be changed in the `make_config` function here.  

Run `python dengue.process.py --titer_model full_tree` for the "full tree" model.  

Run `python dengue.process.py --titer_model interserotype` for the "interserotype" model.  

**NB: For this dataset, step 2 has been run for you; see `./full-tree-model-output/` and `./interserotype-model-output/`**

### 3 - Check out your results

Results are output in JSON format, found in `./processed/`. You can parse and examine results using the notebooks found [here](../figures/).
