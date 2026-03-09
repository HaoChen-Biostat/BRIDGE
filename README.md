# BRIDGE
BRIDGE: The code implements the proposed method (BRIDGE) in "A Statistical Framework for Integrative Imaging Genomics with Biclustering and Ensemble Penalized Regression in Alzheimer’s Disease" by Hao Chen, Yong He, Lei Hou,  Lei Liu and Chencheng Ma.

Imaging genomics provides a powerful paradigm for decoding the complex interplay between genetics and brain functional phenotypes in neurodegenerative disorders. However, statistically integrating these high-dimensional modalities remains challenging, 
as existing methods often rely on incomplete biological priors or overlook the latent modular structure of the data. To address this, we propose BRIDGE (Biclustering and Regression for Integrative Data in Genomics and nEuroimaging), a unified statistical 
framework that synergies biclustering with ensemble penalized regression. Our approach proceeds in three interdependent stages: (1) accurate estimation of individual brain functional networks via sparse precision matrices; (2) data-driven discovery of 
regulatory modules—defined as subsets of genes co-varying with subsets of brain connections—using a sparse biclustering algorithm; and (3) stable disease classification and biomarker identification via dimensionality reduction integrated with an ensemble 
penalized logistic regression model. Simulation studies demonstrate that BRIDGE significantly outperforms alternative approaches in both feature selection stability and classification accuracy. In an application to the Alzheimer’s Disease Neuroimaging 
Initiative (ADNI) cohort, the framework achieved robust out-of-sample classification performance and uncovered biologically interpretable disease-related modules. These results highlight BRIDGE as a rigorous tool for integrative imaging genomics, 
capable of identifying pathogenic factors without reliance on a priori knowledge.
Keywords: Imaging genomics; Multi-omics integration; Alzheimer’s disease; Biclustering; Penalized regression; Brain functional connectivity.

## Descriptions

The `simulation` folder contains files for reproducing the simulation studies:
The `BRIDGE` folder: Files implementing the proposed method (BRIDGE) 
- `SW function.R`: functions to generate small-world network
- `functions.R`: functions required for multivariate sparse regression (theta estimation) and module identification
- `example.R`: example usage
