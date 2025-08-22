# Four Metrics of the Parametric Bootstrap Test

## Environment Setup

Install dependencies using Conda:

```bash
conda env create -f environment.yml
conda activate pbt_new_metrics
```

## Four Metrics for Model Adequacy Evaluation

### 1. Mean Diversity

***The average number of distinct amino acid types*** observed at each site within a simulated dataset or empirical dataset.
 
The closer the mean diversity of the simulated dataset is to that of the empirical dataset, the better the model is considered to capture the across site heterogeneity of amino acid composition in the real data.

Evaluation results for this metric are written to the file: `mean_diversity_results.txt`

> This method was proposed in [Giacomelli et al., 2025](https://doi.org/10.1093/gbe/evae273).

### 2. Mean Entropy

***The average uncertainty of the amino acid frequency distribution (Shannon entropy)*** at each site within a simulated dataset.

The closer the mean entropy of the simulated dataset is to that of the empirical dataset, the better the model is considered to describe the across site heterogeneity of amino acid composition in the real data.

Evaluation results for this metric are written to the file: `mean_entropy_results.txt`

### 3. CvM Test on Sitewise Diversity

For each simulated or empirical dataset, the diversity at all sites is calculated to generate a sitewise diversity distribution. *A two-sample Cramér–von Mises test* is used to compare the distribution from simulated data with that from the empirical data.

The ***W² statistic*** from the Cramér–von Mises test is used as the core metric to quantify the difference between the two distributions. A smaller W² value (closer to 0) indicates that the simulated distribution is more similar to the empirical distribution.

Evaluation results for this metric are written to the file: `cvm_diversity_results.txt`

### 4. CvM Test on Sitewise Entropy

Same as above, but applied to the comparison of sitewise entropy distributions.

Evaluation results for this metric are written to the file: `cvm_entropy_results.txt`

## Step One: Model Fitting

This part of the phylogenetic analysis is performed using the [IQ-TREE](https://github.com/iqtree/iqtree3) software package.
Below is an example pipeline using the **LG-C60opt-R6-PMSF** model:

1) Run the LG+C60+R6 model on the protein alignment file `alignment.nex`, using a guide tree, and optimize the profile weights of C60:

    ```
    iqtree3 -s alignment.nex -m LG+C60+R6 -te guide.treefile -me 0.99 -safe -mwopt -pre LG_C60opt -prec 10
    ```

2) Use the `extract_profiles_weights.py` function to extract the optimized weights from the `LG_C60opt.iqtree` file generated in the first step. The weights will be saved as a `.txt` file in the same directory.

    ```
    python extract_profiles_weights.py <input_directory>
    ```

   Here, <input_directory> refers to the folder containing the `LG_C60opt.iqtree` file.


3) Define a new model based on the optimized weights and save it as `C60opt.nex`. Then run the first phase of the PMSF method in IQ-TREE to estimate the mixture model parameters and infer the site-specific frequency profile.

    ```
    iqtree3 -s alignment.nex -ft guide.treefile -m C60opt -mdef C60opt.nex -safe -pre LG_C60opt_PMSF_step1 -n 0
    ```

4) Run the second phase of the PMSF method, performing phylogenetic inference using the site-frequency profile obtained from the previous step.

    ```
    iqtree3 -s alignment.nex -fs LG_C60opt_PMSF_step1.sitefreq -m C60opt -mdef C60opt.nex -safe -pre infer/LG_C60opt_PMSF_step2 --wbtl -wsr -bb 1000
    ```

5) Use the `combine_sitefq_tate.py` function to generate a `LG_C60opt_PMSF.nex` file. This file defines a partition model where each site is treated as a separate partition, which is required for sequence simulation.

    ```
    python combine_sitefq_tate.py <input_folder>
    ```

    Note:
    
    To use `combine_sitefq_tate.py`, please make sure the files `XX_step1.sitefreq`, `XX_step2.rate`, and `XX_step2.iqtree` are all in the same folder, and share the same filename prefix.

6) Use the [Alisim](https://doi.org/10.1093/molbev/msac092) function in IQ-TREE to simulate sequence alignments based on the inferred partition model and tree.

   ```
   iqtree3 --alisim LG_C60opt_PMSF/seq --seqtype AA -p LG_C60opt_PMSF.nex -t LG_C60opt_PMSF_step2.contree --length [seq_length] --out-format fasta --num-alignments 100
   ```

For more details, please refer to the `data/tutorial.pdf` file.


## Step Two: Starting Parametric Bootstrap Test

To start the parametric bootstrap test, please make sure your input folder is structured as follows:

```
—— data
   ├── model_1        <- Contains m simulated datasets in FASTA format
   │   ├── simulate_1.fa
   │   ├── simulate_2.fa
   │   ├── ...
   │   └── simulate_m.fa
   ├── model_2
   ├── ...
   └── model_n
```

After completing the simulation of datasets, run the following command to perform model evaluation:

```
python PBT.py <SimRootFolder> <OriginalAlignment>  [options]
```
- `SimRootFolder`: the path to the folder containing simulated datasets for n models.

- `OriginalAlignment`: the path to the original empirical dataset.

- `[options]`: optional flags to control which metrics are computed and the number of CPUs to use.

   Available options:

- `-mdiv` : run mean-diversity pipeline

- `-ment` : run mean-entropy pipeline

- `-cvmdiv` : run CvM-diversity pipeline

- `-cvment` : run CvM-entropy pipeline

- `-T N`: number of CPUs to use per model (default: 1)

The output folder structure after running PBT.py is as follows:

```
—— data
   ├── model_1        <- Contains m simulated datasets in FASTA format
   │   ├── simulate_1.fa
   │   ├── simulate_2.fa
   │   ├── ...
   │   └── simulate_m.fa
   ├── model_2
   ├── ...
   ├── model_n
   ├── cvm_diversity_results.txt        <- The result file of the CvM test on sitewise diversity
   ├── cvm_entropy_results.txt        <- The result file of the CvM test on sitewise entropy
   ├── mean_diversity_results.txt        <- The result file of the mean diversity
   ├── mean_entropy_results.txt        <- The result file of the mean entropy
   ├── original.sitewise_diversity.txt        <- Per-site diversity records from the empirical dataset
   ├── original.sitewise_entropy.txt        <- Per-site entropy records from the empirical dataset
   └── PBT.log        <- A complete log file recording all steps, messages, and model selection results of the analysis.
```