sub-directory that contains code to implement, simulate, and plot firing rate model

FRmodel_plots.m — driver to run model and show plots, calls WCnet_EIcells.m . Also displays table of relative diff in firing rate (comparing synDepr), and table of simulated trial var in model. 

WC_EIcells.m — implementing firing rate model E,I pair with Runge-Kutta (4th order). 

getDcdAcc_LONG.m — script to vary means in ortho/retro on large 775 grid of points, using same rho to fit mean decodAcc, creates dJustDeAcc.mat

plotJustDeAcc.m —  plots results from getDcdAcc_LONG.m, showing decodAcc across drug preps, as fcn of |muO/sigO-muR/sigR|.  
