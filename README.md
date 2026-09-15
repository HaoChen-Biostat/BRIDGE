**BRIDGE:** The code implements the proposed method (BRIDGE) in "A Statistical Framework for Integrative Imaging Genomics with Biclustering and Ensemble Penalized Regression in Alzheimer’s Disease" by Hao Chen, Yong He, Lei Hou,  Lei Liu and Chencheng Ma.

Imaging genomics provides a powerful paradigm for decoding the complex interplay between molecular variation and brain functional phenotypes in neurodegenerative disorders. In this study, the molecular modality is specifically gene expression data rather than time-invariant genetic variants. However, statistically integrating these high-dimensional modalities remains challenging, as existing methods often rely on incomplete biological priors or overlook the latent modular structure of the data. To address this, we propose BRIDGE (Biclustering and ensemble penalized Regression for Integrative Data in Genomics and nEuroimaging), a unified statistical framework that synergistically integrates biclustering with ensemble penalized regression. Our approach proceeds in three interdependent stages: (1) accurate estimation of individual brain functional networks via sparse precision matrices; (2) data-driven discovery of regulatory modules—defined as subsets of gene expression features co-varying with subsets of brain connections—using a sparse biclustering algorithm; and (3) stable disease classification and biomarker identification via dimensionality reduction integrated with an ensemble penalized logistic regression model. Simulation studies demonstrate that BRIDGE outperforms ablation variants of the proposed framework and external competing benchmark methods in both feature selection performance and classification accuracy. In an application to the Alzheimer’s Disease Neuroimaging Initiative (ADNI) cohort, the framework achieved robust out-of-sample classification performance and uncovered biologically interpretable disease-related modules. These results highlight BRIDGE as a rigorous tool for integrative imaging genomics, capable of identifying pathogenic factors without reliance on a priori knowledge.
Keywords:Imaging genomics; Gene expression; High-Dimensional Feature Selection; Sparse Biclustering; Ensemble Penalized Regression; Functional Connectivity Network.

## Descriptions

The `Simulation` folder contains files for implementing the **proposed method (BRIDGE)** :

-`SW function.R`: functions for generating small-world networks

-`bridge_functions.R`: core functions for imaging data generation, multivariate sparse regression, module identification, module-wise PCA, and BRIDGE classification

-`Evaluation_functions`.R: functions for evaluating classification, feature selection, and module recovery, and for generating result plots

-`example.R`: example usage
