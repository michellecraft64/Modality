Contains MATLAB code to  ``Modality of odor perception is transmitted to cortical brain regions from the olfactory bulb“  by Michelle F. Craft et al. for data analysis and models.

Sub-directories:


#~/data_analysis/

plot_Hists_Decod.m — loads prev saved decode accuracies (DecodAcc.mat & DecAcc_mnGreater.mat), shows histograms of avg MC , and conditioned on muO<muR & vice-versa
plot_Selectivity.m — shows proportion of retro/orth select cells based on trial-avg firing (horizontal bar graphs), calls Rat[j]*.mat

Population Decoding: 
get_popDecod_PCA.m — implements PCA+LDA, creates PopDecod_pcaLDA.mat (if save_flag=1) used to plot histogram via plot_Popul_PCA.m
get_popDecod_SVM.m  — implements SVM, crease PopDecod_Bayes_SVM.mat (if save_flag=1) used to plot histogram via plot_Popul_SVM.m 


#~/firing_rate/model/

FRmodel_plots.m — driver to run model and show plots, calls WCnet_EIcells.m
WCnet_EIcells.m — implementing firing rate model E,I pair with Runge-Kutta (4th order). 
getDcdAcc_LONG.m — script to vary means in ortho/retro on large 775 grid of points, using same rho to fit mean decodAcc, creates dJustDeAcc.mat
plot_Hists_Decod.m —  plots results from getDcdAcc_LONG.m, showing decodAcc across drug preps, as fcn of |muO-muR|
