%% add paths
load_medriv_project
addpath(genpath('/Users/skatyal/OneDrive/Projects/NeuroMatlabToolboxes/HRVTool/'))
addpath(genpath('/Users/skatyal/OneDrive/Projects/NeuroMatlabToolboxes/xjview/'))
rmpath('~/google_drive/Projects/eeglab2019_1/plugins/Fieldtrip-lite20210601/external/signal/')

%% 1 - initialize experiment parameters

exp_meddepth = meddepth_experiment;

%%% add questionnaire data
exp_meddepth = medriv_addquestdata(exp_meddepth);

%% 2 - extract frequency bands, heart rate, respiration, etc. and save data in R file for stats

chanset = 0;
remove_artif = false;

windowsize = 60; % if positive gives the size of the windows to break the epochs in, if negative remove that much from the initial
noverlap = 20; % no overlap between windows
% windowsize = 20; % if positive gives the size of the windows to break the epochs in, if negative remove that much from the initial
% noverlap = 20; % no overlap between windows
freq_range = [4 7; 7 13; 15 30; 30 50];
filename = 'icacomprem_med_medriv_physio';

medriv_med_eegfreq(exp_meddepth, freq_range, filename, chanset, ...
    remove_artif, windowsize, noverlap);

%% 5.1 - Figure 3A - Interaction of alpha/theta for two groups combined

plot_type = 0; % 0-interaction plot for two groups combined, 1-interaction between alpha/theta and qnum, 2-five individual qnums separately
grayscale = false;
gpname = {'all', 'ltm', 'ctl'};

pthresh = 2.5E-8/26; % .001 corresponds to boneferroni .025 in 26 channels
do_fdr = false;
pltgp = 0;

nfreqs = 2;
for group = pltgp % 0-both, 1-long-term meditators, 2-controls
    medriv_med_topos(plot_type, gpname{group+1}, pthresh, nfreqs, grayscale,...
        do_fdr)
end

%% 5.2 - Figure 3B and 3C - Interaction fo alpha/theta for two groups separate

plot_type = 01; %  1-interaction between alpha/theta and qnum

pthresh = .025/26; % .001 corresponds to boneferroni .025 in 26 channels
do_fdr = false;
pltgp = 1:2;

nfreqs = 2;
for group = pltgp % 0-both, 1-long-term meditators, 2-controls
    medriv_med_topos(plot_type, gpname{group+1}, pthresh, nfreqs, grayscale,...
        do_fdr)
end

%% 5.3 - Supp Figure 1 - Interaction alpha/theta, separate groups, last 6 minutes of M1 and M2

plot_type = 02; %  1-interaction between alpha/theta and qnum

pthresh = .025/26; % .001 corresponds to boneferroni .025 in 26 channels
do_fdr = false;
pltgp = 1:2;

nfreqs = 2;
for group = pltgp % 0-both, 1-long-term meditators, 2-controls
    medriv_med_topos(plot_type, gpname{group+1}, pthresh, nfreqs, grayscale,...
        do_fdr)
end

%% 5.4 - Supp Fig 2 - Interaction alpha/theta, separate groups, no chanting

plot_type = 03; %  1-interaction between alpha/theta and qnum
grayscale = false;
gpname = {'all', 'ltm', 'ctl'};

pthresh = .025/26; % .001 corresponds to boneferroni .025 in 26 channels
do_fdr = false;
pltgp = 1:2;

nfreqs = 2;
for group = pltgp % 0-both, 1-long-term meditators, 2-controls
    medriv_med_topos(plot_type, gpname{group+1}, pthresh, nfreqs, grayscale,...
        do_fdr)
end

%% 5.5 - Supp Fig 3 - Interaction alpha/theta, separate groups, chanting only

plot_type = 04; %  1-interaction between alpha/theta and qnum
grayscale = false;
gpname = {'all', 'ltm', 'ctl'};

pthresh = .025; % .001 corresponds to boneferroni .025 in 26 channels
do_fdr = false;
pltgp = 1;

nfreqs = 2;
for group = pltgp % 0-both, 1-long-term meditators, 2-controls
    medriv_med_topos(plot_type, gpname{group+1}, pthresh, nfreqs, grayscale,...
        do_fdr)
end


%% 5.6 - Supp Fig 7AB - theta difference BL CH

plot_type = 05; %  1-interaction between alpha/theta and qnum
grayscale = false;
gpname = {'all', 'ltm', 'ctl'};

pthresh = .05/26; % .001 corresponds to b oneferroni .025 in 26 channels
do_fdr = false;
pltgp = 1:2;

nfreqs = 2;
for group = pltgp % 0-both, 1-long-term meditators, 2-controls
    medriv_med_topos(plot_type, gpname{group+1}, pthresh, nfreqs, grayscale,...
        do_fdr)
end

%% 5.6 - Supp Fig 7C - gamma difference BL CH in lTM

plot_type = 06; %  1-interaction between alpha/theta and qnum
grayscale = false;
gpname = {'all', 'ltm', 'ctl'};

pthresh = .05/26; % .001 corresponds to boneferroni .025 in 26 channels
do_fdr = false;
pltgp = 1;

nfreqs = 2;
for group = pltgp % 0-both, 1-long-term meditators, 2-controls
    medriv_med_topos(plot_type, gpname{group+1}, pthresh, nfreqs, grayscale,...
        do_fdr)
end

