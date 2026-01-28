# rabies_seqmod
Repo for the output of rigid body docking of target RABV G protein and inferred host protein proteins performed on Clustpro. Target CDSs of study rabv were shared. 
*Source: region X/Borneo

*metadata column fields:
| **Seq-RUN-ID**|	**Sample-ID**	|  **sample_type**|  **species** |  **locality** |  **year**  |
|---------------|---------------|-----------------|--------------|---------------|------------|

# rabies_serostats
##Overview
This R script performs descriptive and inferential analyses of longitudinal antibody titre data collected from repeated sampling of individual animals over time. The workflow is designed to summarize temporal trends in antibody responses and to model non-linear dynamics using generalized additive models (GAMs).

##Key Features
*Import and preprocessing of longitudinal antibody titre datasets
*Descriptive statistics of antibody titres by sampling time, individual, and group (e.g., vaccination status)
*Visualization of temporal antibody trajectories using line plots and smooth curves
*Fitting of generalized additive models (GAMs) to capture non-linear changes in antibody titres over time
* Model diagnostics and goodness-of-fit assessment

##Analytical Approach
*Antibody titres are summarized using appropriate central tendency and dispersion metrics (e.g., mean, median, interquartile range).
*Longitudinal patterns are explored using smoothing techniques to account for repeated measurements.
*GAMs are fitted with time as a smooth term, with optional inclusion of covariates such as age, vaccination history, or sampling interval.
