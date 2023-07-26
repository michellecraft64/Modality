sub-directory that has all MATLAB scripts to process rat data

IndCell_ExtractData.m — creates 15 mat files for each rat from spike data (8 rats for no drug, 4 rats for bicuculline, 3 rats for muscimol), named 
Rat[j]_IndCell_EB_[NoDrug/Bic/Mus].mat (has raw data of spike counts after minimal processing)

OptDecodeTwin_RawFR.m — parse & saves results, creates IndCellRAW_OptWin_EB_0sSpon.mat (has decoding accuracies by rat)

rawFRfiles.m — creates IndCell_rawFR_EB.mat, uses uses rat[j]_IndCell_EB_[NoDrug/Bic/Mus].mat. 
IndCell_rawFR_EB.mat — has NormStruct variable, a 3x2 structure, row1=ND, row2=Bic, row3=Mus, col1=orth, col2=retro. Each entry contains a matrix of raw firing rates: #cells X 10 (trials), combined across rats.  E.g., NormStruct{1,1} is 913x10 matrix of net firing rate, all 913 MCs (no drug) by 10 trials with ortho. 

make_DistMeans_raw.m — creates file DistMeans_raw.mat, used for distance between |orth-retro| (trial-avg), IndCell_normFR_EB.mat

make_DecodAccMat_raw.m — quick script to create DecodAcc_raw.mat, combining all decoding accuracies across rats from  IndCellRAW_OptWin_EB_0sSpon.mat ; comment out parts to get IndCellRAW3shift_OptWin_EB_0sSpon.mat

calc_DistDecode_raw.m — plots decoding acc as fcn of distance |muO-muR| and displays correlation (Pearson’s & Spearman’s Rank), also with |muO/stdO-muR/stdR|

plot_Hists_Decod_raw.m — loads prev saved decode accuracies (DecodAcc.mat & DecAcc_mnGreater.mat), shows histograms of avg MC , and conditioned on muO<muR & vice-versa.  Has p-values for significance tests, as well as Effect Size calc.

plot_Selectivity.m — shows proportion of retro/orth select cells based on trial-avg firing (horizontal bar graphs), calls Rat[j]*.mat (Fig 1Aiii)

get_plot_PSTH.m — plots various pop PSTH, in Fig 1Aii, Fig 2, Fig 3B, as well as p-values (time series) for when drug diff significant. 
Contains code to plot the heatmap for trial-avg firing rate and std dev across trials (Fig 2), tables for correlation between trial-var and firing rate (Fig A3, Table A2). 

get_plot_PSTH_extndTime.m — extended time so returns to baseline (Fig A1)

TrialAvg_SigTests_raw.m — same as above but on raw firing rates (IndCell_rawFR_EB.mat), , as well as Effect Size calc.

Population Decoding (Fig 6): 
get_popDecod_PCA.m — implements PCA+LDA, creates PopDecod_pcaLDA.mat (if save_flag=1) used to plot histogram via plot_Popul_PCA.m
get_popDecod_SVM.m  — implements SVM, crease PopDecod_Bayes_SVM.mat (if save_flag=1) used to plot histogram via plot_Popul_SVM.m 

Fig A2:
avgDecod_varyWind.m — gets and plots avg +/- std of decoding accuracy as time windows vary.

Fig A2:
get_find_Optim_Tev_raw.m — Code for Figure A2 — p-values with different windows, using mat files: IndCellRAW_OptWin_EB_0sSpon.mat 
and IndCellRAW3shift_OptWin_EB_0sSpon.mat

offResponse_anal.m — script to produce table of percentage of “Off Response” cells that fire more after odor is removed. Displays Table A1.
