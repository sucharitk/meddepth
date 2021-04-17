%% add paths
load_medriv_project
addpath(genpath('/Users/skatyal/OneDrive/Projects/NeuroMatlabToolboxes/HRVTool/'))

%% 1 - initialize experiment parameters

exp_medepth = medriv_experiment;

%%% add questionnaire data
exp_medepth = medriv_addquestdata(exp_medepth);

%% 2 - extract frequency bands, heart rate, respiration, etc. and save data in R file for stats

chanset = 0;
remove_artif = false;

windowsize = 60; % if positive gives the size of the windows to break the epochs in, if negative remove that much from the initial
noverlap = 20; % no overlap between windows
freq_range = [4 7; 7 13; 15 30; 30 50];
filename = 'icacomprem_med_medriv_physio';

medriv_med_eegfreq(exp_medepth, freq_range, filename, chanset, ...
    remove_artif, windowsize, noverlap);

%% 3 - topographic plots of interactions and main effects from R stats

plot_type = 01; % 0-interaction plot for two groups combined, 1-interaction between alpha/theta and qnum, 2-five individual qnums separately
grayscale = false;
gpname = {'all', 'ltm', 'ctl'};

pthresh = 0.025; % control condition only CH for LTM (Supplementary fig 3)
do_fdr = false;
pltgp = 1:2;


nfreqs = 2;
for group = pltgp % 0-both, 1-long-term meditators, 2-controls
    medriv_med_topos(plot_type, gpname{group+1}, pthresh, nfreqs, grayscale,...
        do_fdr)
end


