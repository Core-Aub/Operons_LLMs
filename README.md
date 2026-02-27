# Operons Prediction With Large Language Models

**"Transformer-Based Operon Prediction Using Textual Representations of Gene Pairs"** by Rida Assaf and Basel Fakhri.
<p align="center">
  <img src="LLM.png" alt="Operon prediction pipeline" width="750">
</p>

# Reference
-To be added when published 

Please cite this paper if you find our work useful.

# Description
This repository contains code and data to reproduce the experiments from the paper **"Transformer-Based Operon Prediction Using Textual Representations of Gene Pairs"**.
## Repository structure

- `Data/`
  - `data_generation_pipeline/`
    - `serialize_gene_pairs.py`: generates serialized datasets by converting adjacent gene pairs into a natural-language description.
  - `features/`
    - `patric_features/`: PATRIC genome annotations (`*.PATRIC.features.tab`) used to list genes in genomic order.
    - `extracted_features/`: precomputed features used in the textual descriptions
      - `GC_ratio.csv`: per-gene GC ratio
      - `codon_bias.csv`: per-gene codon bias vectors (stored as a dictionary per gene)
      - `conservation_scores.tsv`: tab-separated file with conservation scores between gene families
      - `length.csv`: per-gene length (kept for completeness)
      - `string_scores.csv`: Gene-pair STRING interaction scores.
  - `generated_datasets/`: ready-to-train datasets (CSV) with columns `genome_id`, `gene1`, `gene2`, `text`, `label`
      This directory contains three subfolders corresponding to different training settings:

    - `BSUBTILIS/`  
      Generated datasets corresponding to labels from `main_dataset/LOSO_bsubtilis.csv`.  
      Used for the LOSO setting where *Bacillus subtilis* is held out for testing.  
      Contains multiple `.csv` files generated with different feature configurations.
  
    - `ECOLI/`  
      Generated datasets corresponding to labels from `main_dataset/LOSO_ecoli.csv`.  
      Used for the LOSO setting where *Escherichia coli* is held out for testing.  
      Contains multiple `.csv` files generated with different feature configurations.
  
    - `ROUND_ROBIN/`  
      Generated datasets corresponding to labels from `main_dataset/LOSO_roundrobin.csv`.  
      Used for the random split training setting.  
      Contains multiple `.csv` files generated with different feature configurations.
  
    Each subfolder includes several dataset variants built using different combinations of features  
    (e.g., vanilla, +function, +family, +conservation, +STRING).   
  - `main_dataset/`: Base labeled gene pairs.
    - `LOSO_bsubtilis.csv`  
      Used for the Leave-One-Species-Out (LOSO) training setting when testing on *Bacillus subtilis*.  
      - Training positive examples are gathered from the Operon DataBase (ODB) [1].  
      - Testing data for *Bacillus subtilis* are obtained from RegulonDB [2] and DBTBS [3].  
    - `LOSO_ecoli.csv`  
      Used for the Leave-One-Species-Out (LOSO) training setting when testing on *Escherichia coli*.  
      - Training positive examples are gathered from the Operon DataBase (ODB) [1].  
      - Testing data for *Escherichia coli* are obtained from RegulonDB [2] and DBTBS [3].  
    - `roundrobin.csv`  
      Used for the random split training setting.  
      - All positive examples are gathered from the Operon DataBase (ODB) [1].

---

- `Training/`
  - `RandomSplit_pipeline/Operons_RandomSplit.ipynb`: fine-tuning and evaluation with a random stratified train/test split.
  - `LOSO_pipeline/Operons_LOSO.ipynb`: 
    Leave-One-Species-Out (LOSO) evaluation by holding out one or more `genome_id`s.  
    This notebook:
    - Trains models while excluding the selected species from the training set.  
    - Evaluates generalization performance on the held-out species.  
    - Enables systematic cross-species validation experiments.  
    - Includes a dedicated section for testing model resilience across varying feature sets  
      (e.g., full features, -function, -family, -conservation),  
      allowing controlled analysis of how different biological feature combinations affect performance and robustness.

- `Classifier/`
  - `roberta_operons_classifier.ipynb`:  
    It includes a variety of models trained under different LOSO settings and across multiple feature configurations (e.g., vanilla, +function, +family, +conservation, +STRING), enabling comparative evaluation and robustness analysis.:            notebook for loading the final trained model for inference, as well as for further improvement and fine-tuning on new or existing datasets.


## Quickstart

### 1) Install dependencies

The notebooks were written for Google Colab, but they also run locally.

```bash
pip install -U transformers datasets evaluate torch scikit-learn pandas numpy matplotlib
```
## 2) Script Usage (`serialize_gene_pairs.py`)

To run the script:

```bash
python3 serialize_gene_pairs.py <genome_dir> <gc.csv> <labels.csv> <function> <family> <conservation.tsv?> <conservation+string?> <conservation+string+codon?> <output.csv>
```

### Arguments

- `<genome_dir>`  
  Path to the directory containing genome files.

- `<gc.csv>`  
  GC content file.

- `<labels.csv>`  
  Labels file.

- `<function>`  
  Insert `1` to include functional features, `0` otherwise.

- `<family>`  
  Insert `1` to include family features, `0` otherwise.

- `<conservation.tsv?>`  
  Optional. Provide a conservation file to include conservation features.

- `<conservation+string?>`  
  Optional. Use this if you want to add STRING scores along with conservation.

- `<conservation+string+codon?>`  
  Optional. Use this if you want to add codon features in addition to conservation and STRING scores.

- `<output.csv>`  
  Desired output file path.


### References

[1] Operon DataBase (ODB):  
Pertea, M., Ayanbule, K., Smedinghoff, M., & Salzberg, S. L. (2009).  
OperonDB: a comprehensive database of predicted operons in microbial genomes.  
*Nucleic Acids Research*, 37(Database issue), D479–D482.

[2] RegulonDB:  
Santos-Zavaleta, A., et al. (2019).  
RegulonDB v10.5: tackling challenges to unify classic and high throughput knowledge of gene regulation in *Escherichia coli* K-12.  
*Nucleic Acids Research*, 47(D1), D212–D220.

[3] DBTBS (Database of Transcriptional Regulation in *Bacillus subtilis*):  
Sierro, N., et al. (2008).  
DBTBS: a database of transcriptional regulation in *Bacillus subtilis*.  
*Nucleic Acids Research*, 36(Database issue), D93–D96.

