sub-directory that has all MATLAB scripts to process rat data

IndCell_ExtractData.m — creates 15 mat files for each rat from spike data (8 rats for no drug, 4 rats for bicuculline, 3 rats for muscimol), named 
Rat[j]_IndCell_EB_[NoDrug/Bic/Mus].mat (has raw data of spike counts after minimal processing)

OptDecodeTwin.m — parse & saves results, creates IndCellThrsh_OptWin_EB_[1/2/3]Spon.mat  (has decoding accuracies by rat, can ignore p-vals)

normFRfiles.m — creates IndCell_normFR_EB.mat, uses rat[j]_IndCell_EB_[NoDrug/Bic/Mus].mat. 
IndCell_normFR_EB.mat — has NormStruct variable, a 3x2 structure, row1=ND, row2=Bic, row3=Mus, col1=orth, col2=retro. Each entry contains a matrix of net firing rates: #cells X 10 (trials), combined across rats.  E.g., NormStruct{1,1} is 913x10 matrix of net firing rate, all 913 MCs (no drug) by 10 trials with ortho. 

make_DistMeans.m — creates file DistMeans.m, used for distance between |orth-retro| (trial-avg), IndCell_normFR_EB.mat
make_DecodAccMat.m — quick script to create DecodAcc.mat, combining all decoding accuracies across rats from IndCellThrsh_OptWin_EB_2sSpon.mat

calc_DistDecode.m — plots decoding acc as fcn of distance |muO-muR| and displays correlation (Pearson’s & Spearman’s Rank)

plot_Hists_Decod.m — loads prev saved decode accuracies (DecodAcc.mat & DecAcc_mnGreater.mat), shows histograms of avg MC , and conditioned on muO<muR & vice-versa.  Has p-values for significance tests too!
plot_Selectivity.m — shows proportion of retro/orth select cells based on trial-avg firing (horizontal bar graphs), calls Rat[j]*.mat (Fig 1Aiii)

get_plot_PSTH.m — plots various pop PSTH, in Fig 1Aii and Fig 2B, as well as p-values (time series) for when drug diff significant

TrialAvg_SigTests.m — does comparison of population averages of trial-to-trial spike rate variance, using IndCell_normFR_EB.mat. Shows results in 2 Tables for population mean of trial-to-trial var.

Population Decoding (Fig 5): 
get_popDecod_PCA.m — implements PCA+LDA, creates PopDecod_pcaLDA.mat (if save_flag=1) used to plot histogram via plot_Popul_PCA.m
get_popDecod_SVM.m  — implements SVM, crease PopDecod_Bayes_SVM.mat (if save_flag=1) used to plot histogram via plot_Popul_SVM.m 


get_find_Optim_TspTev.m — Code for Figure S1 — p-values with different windows, using mat files: IndCellThrsh_OptWin_EB_[1/2/3]Spon.mat
