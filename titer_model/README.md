# Titer model

## About the titer model

The titer model seeks to assign discrete units of antigenic change to specific mutations across the dengue phylogeny.
Full model details can be found in [Neher et al.](https://www.pnas.org/content/113/12/E1701) and our [preprint](https://bedford.io/papers/bell-dengue-antigenic-dynamics/).

## Running the titer model

This model was originally published in [Neher et al (PNAS, 2016)](http://dx.doi.org/10.1073/pnas.1525578113) and implemented as part of the Nextstrain [augur package](https://github.com/nextstrain/augur). The relevant portions of the repository have been reproduced here under [`implementation-nextstrain-augur`](implementation-nextstrain-augur/). Documentation, including install instructions, for the full augur pipeline can be found [here](https://github.com/nextstrain/augur/tree/6d9f7088d8792196e5021c67b876d9de1d2a13dd). Note that for future reuse, I would actually recommend running the model via the newly modularized implementation [here](https://github.com/nextstrain/augur/blob/master/augur/titers.py).

### 1 - Prepare input files

`dengue.prepare.py` handles subsampling and other basic dataset config. You can edit the function `make_config` here with any desired changes.

```
cd dengue-antigenic-dynamics/titer_model/implementation-nextstrain-augur/builds/dengue/
python dengue.prepare.py #--options
```

**NB: For this dataset, step 1 has been run for you; see `./implementation-nextstrain-augur/dengue/prepared/titered.json`**

### 2 - Run the titer pipeline (and estimate clade frequencies)

`dengue.process.py` handles the actual analysis; parameter settings, etc. can be changed in the `make_config` function here.  
`python dengue.process.py -j path_to_prepared_json`

**NB: For this dataset, step 2 has been run for you; see `./titered_output/`**

### 3 - Check out your results

Results are output in JSON format, found in `./processed/`. You can parse and examine results using the notebook found [here](../data_wrangling_scripts/).
