# metabolomics
Analysis pipelines for metabolomics data

## Template 1.0
An analysis pipeline with a detailed instruction manual designed to be used by people with little-to-no experience with R or programming in general. It includes example .csv files that the pipeline references. This allows the user to set up experiment information more easily. It was designed for use in Kyu Rhee's lab at Weill Cornell Medicine. 

The pipeline uses xcms for untargeted peak picking and feature grouping. It has optimized parameters for the Rhee lab's Agilent QTOF and TOF instruments, which run an aqueous normal phase liquid chromatography method using a diamond hydride column. It also includes some QC to confirm appropriate peak picking, using metabolites easily observed from Mycobacterium tuberculosis lysates. The pipeline then filters features by m/z corresponding to expected \[M+H\]+ or \[M-H\]- from a set of  unique formulae from KEGG Mycobacterium tuberculosis metabolic pathways. Finally it creates volcano plots of features, and prints chromatograms of significantly changing features. For more details, please refer to the Instructions for XCMS Pipeline V1.0. 
