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

<details>
<summary> Tips</summary>

For more details, please refer to the `data/tutorial.pdf` file.
</details>


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
python PBT.py <SimRootFolder> <OriginalAlignment> [NumCPUs]
```
- Here, `SimRootFolder` is the path to the folder containing simulated datasets for n models, and `OriginalAlignment` is the path to the original empirical dataset.

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