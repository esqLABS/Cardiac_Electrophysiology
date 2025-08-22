README

The current scripts allow simulating the effect of benchmark CiPA drugs on a cardiac electrophysiology model. The scripts were translated from MATLAB and adapted. The R implementation was validated by reproducing the original results for the benchmark CiPA drugs. 

Original publication: 
Llopis-Lorente J, Baroudi S, Koloskoff K, Mora MT, Basset M, Romero L, Benito S, Dayan F, Saiz J, Trenor B. Combining pharmacokinetic and electrophysiological models for early prediction of drug-induced arrhythmogenicity. Comput Methods Programs Biomed. 2023 Dec;242:107860. doi: 10.1016/j.cmpb.2023.107860. Epub 2023 Oct 11. PMID: 37844488.

The virtual population was generated as part of that initial publication and is simply used as input in the current R workflow. 
The input Free plasma concentration (FPC) for each CiPA drugs is the one from the original publication in order to reproduce the results and validate the R implementation. Those FPCs were computed for each drug and each scenario with non-compartmental pop-PK models in the original paper. 

R workflow and scripts:
- Master file: modelRunnerDrugs.R (Adapt the name of drugs to be simulated within the script): uses parallel computing (individuals run in parallel) 
- Results are saved in 'Results'
- For postprocessing, biomarker calculation and plotting, run 'SimResults_Analysis.R' (needs to be adapted for each input result file) 

Example simulation outputs are available in the folder 'Results' with a limited population size.