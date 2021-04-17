function exp_medriv = medriv_experiment

% exp_medriv.session_dir = 'C:\g\Projects\Experiments\Meditation_Rivalry\Data\';
exp_medriv.session_dir = '~/OneDrive/Projects/Experiments/Meditation_Rivalry/Data';

exp_medriv.eeg_datadir = 'Data/eeglab';
exp_medriv.eeg_epochdir = 'Epochs';

subj_codes = [4140 4143 4144 4145 4146 4147 4148 4149 4150 4151 4152 ...
    4153 4154 4155 4156 4158 4159 4160 4161 4162 4163 ...
    4164 4165 4167 4168 4169 4170 4174 4173 4181 4186 ...
    4190 4191 4192 4193];

years_practice = [43 52 NaN 44 46 43 40 17 16 38 27,...
    0 0 35 15 38 11 0 29 38 20, ...
    0 0 0 0 0 28 10 8 0 0, ...
    0 0 0 0];

income = [NaN NaN NaN 30000 120000 NaN 6000 NaN NaN 5000 NaN,...
    85000 NaN 50000 35000 6800 3600 NaN NaN 2000 NaN, ...
    90000 NaN NaN 96000 NaN 6500 30000 49100 63672 48000, ...
    40000 NaN 45000 30000];

nsubj = numel(subj_codes);

subj_age = [67 63 67 68 66 69 64 37 68 58 54,...
    51 60 65 33 54 43 43 62.5 62 42 ...
    73 29 57 60 32 62.5 30 34 62 56 ...
    65 29 52 68];

subj_sex = [0 zeros(1, 7) ones(1, 3), ...
    1 0 0 1 0 1 1 0 0 1, ...
    0 1 1 0 0 1 1 0 1 0, ...
    0 1 0 0]; % 0 males, 1 females

subj_edu = [20 12 20 16 20 12 18 20 4150 12 18 ... % check 4145
    18 4154 18 16 18 18 18 14 14 18 ...
    18 14 18 16 14 16 16 21 18 20 ...
    14 14 14 16];

subj_valid = true(1, nsubj);

subj_valid_riv = ones(1, nsubj); % for most subjects where rivalry is valid 
% subj_valid_riv(subj_codes==4143 | subj_codes==4147 | subj_codes==4164 | ...
% %    subj_codes==4155) = 0; % 4143-didn't see rivalry, 4147-pressed the wrong buttons, 4155-didn't have enough rivalry runs as he was in a hurry, 4164-highly disproportional eye dominance
subj_valid_riv(subj_codes==4143 | subj_codes==4147) = 0; % 4143-didn't see rivalry, 4147-pressed the wrong buttons, 4155-didn't have enough rivalry runs as he was in a hurry
subj_valid_riv(subj_codes==4164) = 0; % 4164-highly disproportional eye dominance
subj_valid_riv(subj_codes==4148) = 1; % for one participant whose dominant is not valid but mixed are
subj_valid_riv(subj_codes==4191) = 0; % highly disproportional eye dominance duration

subj_valid_sret = 2*ones(1, nsubj);
subj_valid_sret(subj_codes==4140 | subj_codes==4174) = 0; % for participants who ran the task before we had sret (also their baseline is eyes-closed not listening to story)
subj_valid_sret(subj_codes==4143 | subj_codes==4144 | subj_codes==4145) = 1; % 50 trial subjects

riv_runs = cell(1, nsubj);
riv_runs{subj_codes==4140} = {1:5, 6:9, 10:13};
riv_runs{subj_codes==4148} = {2:5, 6:10, 11:15}; % subject saw a lot of mixed, and for the very first run didn't report any alternations so get rid of that one
riv_runs{subj_codes==4153} = {[1:4 6], 7:11, 12:16};
riv_runs{subj_codes==4155} = {1:4, 5:8, 9:10};
riv_runs{subj_codes==4174} = {1:4, 5:8, 9:12};

bad_riv_eeg_runs = cell(1, nsubj);
bad_riv_eeg_runs{subj_codes==4160} = {[1 5], [], []};
bad_riv_eeg_runs{subj_codes==4161} = {[3 5], [], [2]};
bad_riv_eeg_runs{subj_codes==4169} = {[3], [], [4]};

% bad_resp_runs = cell(1, nsubj);
% % bad_resp_runs{subj_codes==4144} = {1;2;3;4;5;6;7;8;9;10;11;12;13;14;15};
% bad_resp_runs{subj_codes==4144} = {[1 0 12]; [4 0 10]; [6 0 20]; [7 0 28];...
%     [9 15 28]; [11 12 30]; 12; 13; 14; 15};
% bad_resp_runs{subj_codes==4145} = {[1 0 25]; [2 0 12.6]; 3; [4]; [5 0 45]; ...
%     [6]; [7]; 8; [9 ]; [10 ]; 11; 12; 13; 14; 15};
% bad_resp_runs{subj_codes==4146} = {[3 40 50]; [4 0 15]; 6; 7; 8; 9; 10; ...
%     [13 0 10 62.3 70]};
% bad_resp_runs{subj_codes==4148} = {[4 0 10]; [5 0 20]; [10 0 20]};
% bad_resp_runs{subj_codes==4149} = {[1 0 7.5]; [2 0 10]; 3; 4; [5 45 -1];...
%     [6 0 30]; 7; 8; 9 ;10; 11; 12; 13; [14 0 5]; 15};
% bad_resp_runs{subj_codes==4150} = {[5 0 9.5]; [11 20 25]; [12 60 -1]; ...
%     [13 25 -1]; 14; 15};
% bad_resp_runs{subj_codes==4151} = {[1 0 10 60 70]; [2 0 10]; [3 0 15 50 60];...
%     [4 0 20]; 5; 6; [7 0 20]; 8; [9 0 30 50 57.5]; 10; [11 0 40]; [12 60 -1]; ...
%     [13 25 -1]; 14; 15};
% bad_resp_runs{subj_codes==4152} = {[11 0 18.5]; [12 16.3 20.5]; [13 0 30]; [15 0 30]};
% bad_resp_runs{subj_codes==4153} = {[1 0 10]; [2 0 10]; [3 0 5];...
%     [4 0 10]; [5 0 10]; [6 0 20]; [7 0 15]; [8 0 15]; [9 0 15]; [10 10 20];...
%     [11 0 20]; [12 0 30]; [13 0 30]; [15 0 30]};
% bad_resp_runs{subj_codes==4154} = {[1 0 20]; [2 0 10]; [3 0 10];...
%     [4 0 10]; [5 0 10]; [6 0 30]; [11 0 20]};
% bad_resp_runs{subj_codes==4155} = {[1 0 10]; [2 0 10]; [3 0 10];...
%     [4 0 10]; [5 0 25]; [6 0 20]; [7 0 20]; [8 0 20]; [9 0 20]; [10]};
% bad_resp_runs{subj_codes==4156} = {[2 0 5]; [3 0 10];...
%     [5 0 10]; [6 0 15]; [7 0 17.5]; [9 0 18]; [10 0 20];...
%     [11 0 20]};
% bad_resp_runs{subj_codes==4158} = {[2 0 5]; [3 0 6];...
%     [4 0 6.5]; [5 0 6]; [6 0 25.5]; [7]; [8 0 45]; [10 0 16];...
%     [11 0 22]; [12]; [13]; [14 15.5 59.1]; [15 0 46]};
% bad_resp_runs{subj_codes==4159} = {[2 0 12.5]; [3 0 17.7];...
%     [4 0 15]; [5 0 18.5]; [6 0 20]; [8 0 35]; [9 0 30]; [10 0 25];...
%     [11 0 30]; [12 0 31]; [13 0 40]};
% bad_resp_runs{subj_codes==4160} = {[2 0 10]};
% bad_resp_runs{subj_codes==4161} = {[2 0 10]; ...
%     [4 0 10]; [6 ]; [7]; [10 ];...
%     [11 0 20]; [12]; [15 ]};
% bad_resp_runs{subj_codes==4162} = {[5 50 -1]; [6];7;8;9;10;11;12;13;14;15};
% bad_resp_runs{subj_codes==4163} = {[4 0 5.7]; [6 20 35]; [15 20 29]};
% bad_resp_runs{subj_codes==4165} =  {[1 47 64]; [2 0 30]; [3 0 30];...
%     [11 0 30]};
% bad_resp_runs{subj_codes==4167} =  {[2 0 16]; ...
%     [4 0 7.5]; [5 0 8.8]; [6 0 20]; [7]; [8 0 20]; [9 0 20]; [10 0 20];...
%     [11 0 25]; [12 0 24]; [13 0 23]; [14 0 22]; [15 ]};
% bad_resp_runs{subj_codes==4168} =  {[1]; [2 0 10]; [3];...
%     [4 ]; [5 ]; [6 ]; [7]; [8 0 30]; 9; [10 ];...
%     [11 ]; [12]; [13]; [14 ]; [15 21 25 ]};
% bad_resp_runs{subj_codes==4170} =  {[1 0 10]; [2 0 15]; [3 ];...
%     [4 ]; [5 ]; [6 ]; [7 ]; [8 ]; 9; [10 ];...
%     [11 ]; [12]; [13]; [14 ]; [15 ]};
% bad_resp_runs{subj_codes==4173} =  {[1 0 6.4]; [2 0 10]; [6 0 16.8];...
%     [7]; [8 57 -1]; 9; [11 0 28]; [12]; [13 0 20]; [14 0 30]; [15 ]};
% 
% bad_resp_runs{subj_codes==4192} = {[4 0 11];[6 0 32];[7 7.3 31.7];...
%     [8 10 22]; [9 0 26.6]; [11 0 50.8]; [12 14.8 25.3]; [14 13.9 24.8];...
%     [15 10.8 29.4]};
% bad_resp_runs{subj_codes==4193} = {1;2;3;4;5;6;7;8;9;10;11;12;13;14;15};

paf_qual = zeros(1, nsubj); % 1= not good alpha everywhere; 
paf_qual(subj_codes==4144) = 1; % 2= not good alpha everywhere, not good 1/f;
paf_qual(subj_codes==4154) = 2; %  3= not good occipital alpha, but good frontal alpha; 
paf_qual(subj_codes==4156) = 3; %4= not good 1/f % 5= freq spectrum highly variable across runs
paf_qual(subj_codes==4159) = 3;
paf_qual(subj_codes==4160) = 4;
paf_qual(subj_codes==4161) = 5;
paf_qual(subj_codes==4162) = 3;
paf_qual(subj_codes==4168) = 3;
paf_qual(subj_codes==4186) = 3;

paf_qual_med = zeros(1, nsubj); % 1= not good alpha everywhere
paf_qual_med(subj_codes==4154) = 1;
paf_qual_med(subj_codes==4168) = 1;
paf_qual_med(subj_codes==4174) = 1;

exp_medriv.nsubj = nsubj;

cd(exp_medriv.session_dir)
datadirs = dir;
datadirs = datadirs([datadirs.isdir]);
datadirs = datadirs(3:end);

alldirs = {datadirs.name};

exp_medriv.groups = {'AM', 'HC', 'LT'};
exp_medriv.med_trigs = {'S108', 'S109'; 'S110', 'S111'; 'S112', 'S113'; 'S114', 'S115'};
exp_medriv.riv_trigs = {'S101', 'S102'};

exp_medriv.chans.locs{1} = {'O1', 'Oz', 'O2'};
exp_medriv.chans.locs{2} = {'CP5', 'P3', 'P7'};
exp_medriv.chans.locs{3} = {'CP6', 'P4', 'P8'};
exp_medriv.chans.locs{4} = {'FC1', 'FC5', 'F3'};
exp_medriv.chans.locs{5} = {'FC2', 'FC6', 'F4'};
exp_medriv.chans.locs{6} = {'Fz', 'FC1', 'FC2', 'Cz'};
exp_medriv.chans.group_name = {'occip', 'lpar', 'rpar', 'lfront', 'rfront', ...
    'mfront'};

for ii = 1:nsubj
    
    findirs = strfind(alldirs, [num2str(subj_codes(ii)) '_']);
    dirname = alldirs{~isemptycell(findirs)};

    exp_medriv.data(ii).subj_code = subj_codes(ii);
    exp_medriv.data(ii).subj_age = subj_age(ii);
    exp_medriv.data(ii).subj_edu = subj_edu(ii);
    exp_medriv.data(ii).subj_sex = subj_sex(ii);
    exp_medriv.data(ii).years_practice = years_practice(ii);
    exp_medriv.data(ii).subj_valid = subj_valid(ii);
    exp_medriv.data(ii).subj_valid_riv = subj_valid_riv(ii);
    exp_medriv.data(ii).subj_valid_sret = subj_valid_sret(ii);
    exp_medriv.data(ii).riv_runs = riv_runs{ii};
    exp_medriv.data(ii).bad_riv_eeg_runs = bad_riv_eeg_runs{ii};
    %     exp_medriv.data(ii).bad_resp_runs = bad_resp_runs{ii};

    exp_medriv.data(ii).paf_quality = paf_qual(ii);
    exp_medriv.data(ii).paf_quality_med = paf_qual_med(ii);

    exp_medriv.data(ii).dir_name = dirname;
    sdate = dirname(6:11);
    sdate = sdate([3 4 1 2 5 6]);
    exp_medriv.data(ii).session_date = sdate;
    exp_medriv.data(ii).group = find(strcmp(exp_medriv.groups, dirname(end-1:end)));

end
end