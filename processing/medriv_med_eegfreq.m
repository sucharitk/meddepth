function  medriv_med_eegfreq(exp_medriv, freq_range, filename, chanset, ...
    remov_artif, windowsize, noverlap)

group_n = zeros(1, 3);

cd(exp_medriv.session_dir)

nsubj = exp_medriv.nsubj;
do_plots = false;

alpharange = [7 13];
onebyfrange = [.5 5.5];
alpharange_aroundpeak = .5;
sub1byflag = true;

nfreq = size(freq_range, 1);
powmeth = 1; %1- fft, 2-pwelch, 3-fft^2, 4, relative fft power
powmethname = {'fft', 'pwelch', 'fft2', 'relfftpow'};

blks_to_an = [1 2 3 4]; % 4 blocks of measurement
nblocks = numel(blks_to_an);

fs = 250;


freqpow_subj = cell(nsubj, nblocks);
timefrend_subj = cell(nsubj, nblocks);
timefrbeg_subj = cell(nsubj, nblocks);
paf_subj = cell(nsubj, nblocks);
pafsnr_subj = paf_subj;
rrate_subj = cell(nsubj, nblocks);
iqrr_subj = rrate_subj;
hrate_subj = rrate_subj;
hrv_subj = rrate_subj;

start_subj = 1;

for ns = start_subj:nsubj
    
    subj_data = exp_medriv.data(ns);
    
    if subj_data.subj_valid && any(~isnan(subj_data.questresp(:)))
        if subj_data.group==1 || subj_data.group==2
            
            cd(fullfile(exp_medriv.session_dir, subj_data.dir_name))
            group_n(subj_data.group) = group_n(subj_data.group)+1;
            
            fprintf('\nanalyzing subject %g\n', subj_data.subj_code)
            epload = load(['Epochs/' filename]);
            epochs = epload.epochs;
            artif = epload.artif;
            chanlocs = epload.chanlocs;
            
            if ~exist('chanlab32', 'var')
                
                if chanset
                    chanlab32 = chans_to_an;
                else
                    chanlab32 = {chanlocs.labels};
                    chanlab32 = chanlab32(1:32);
                end
                %                 nchanan = 32;
            end
            
            chanlabels = {chanlocs.labels};
            if chanset
                chanloc_inds = get_channels_from_labels(chanlabels, chans_to_an);
                nbchan = sum(chanloc_inds);

                chanmiss = get_channels_from_labels(chanlab32, chans_to_an);
                chanmiss = ~chanmiss;
            else
                % average over all available channels
                nbchan = numel([chanlocs.X]);
                chanloc_inds = ~isemptycell({chanlocs.X});
                
                rrate_ch = strcmp(chanlabels,'RSPRate');
                hrate_ch = strcmp(chanlabels, 'HRBVP');
                bvp_ch = strcmp(chanlabels, 'BVP');
                resp_ch = strcmp(chanlabels, 'RSP');
                if sum(rrate_ch) || sum(hrate_ch) || sum(bvp_ch) || sum(resp_ch)
                    phys = true;
                else
                    phys = false; % subject does not have physio data
                    fprintf('no phys data on %d\n', subj_data.subj_code)
                end
                
                chanlabels = chanlabels(chanloc_inds);
                chanmiss = get_channels_from_labels(chanlab32, chanlabels);

            end
            
            fprintf('block number: ')
            parfor bn = 1:nblocks
                fprintf('%g ',bn)
                nb = blks_to_an(bn);

                if ~isempty(epochs{nb})
                    if remov_artif
                        valid_epoch = epochs{nb}(chanloc_inds, artif{nb});
                    else
                        valid_epoch = epochs{nb}(chanloc_inds, :);
                    end
                    
                    if windowsize>0
                        [epout, timefrend, timefrbeg] = ...
                            break_epochs_withoverlap(valid_epoch, fs, ...
                            windowsize, noverlap);
                    else
                        if bn==1
                            epout = {valid_epoch(:, -windowsize/10*fs:end)};
                        else
                            epout = {valid_epoch(:, -windowsize*fs:end)};
                        end
                    end
                    
                    if phys
                        %                         rrate = epochs{nb}(rrate_ch, :);
                        resp = epochs{nb}(resp_ch, :);
                        bvp = epochs{nb}(bvp_ch, :);
                        
                        zresp = ~resp;
                        zbvp = ~bvp;
                        bvp(zbvp) = NaN;
                        resp(zresp) = NaN;
                        
                        
                        if windowsize>0
                            respout = break_epochs_withoverlap(resp, fs, ...
                                windowsize, noverlap);
                            bvpout = break_epochs_withoverlap(bvp, fs, ...
                                windowsize, noverlap);
                        else
                           % to work on, if needed for some analysis
                        end

                    end
                    
                    nepsize = numel(epout);
                    freqpow = NaN(nepsize, nfreq, nbchan);
                    paf = NaN(nepsize,nbchan); pafsnr = paf; %slope = paf;
                    iqrr2 = NaN(1, nepsize); resp2 = iqrr2; hrv2 = iqrr2; hr2=iqrr2;
                    for neps = 1:nepsize
                        
                        [fpow, ft_ep, freqs] = ...
                            freqrangepow(epout{neps}, fs, freq_range, powmeth);
                        
                        freqpow(neps, :, chanmiss) = fpow;
                        
                        [paf(neps, chanmiss),~,~,pafsnr(neps,chanmiss)] =...
                            peak_alpha_frequency(...
                            ft_ep, freqs, alpharange, ...
                            onebyfrange, sub1byflag, alpharange_aroundpeak);
                        
                        
                        if phys
                            % if physiology data is available for the
                            % subject
                            
                            lepoch = numel(bvpout{neps});
                            if sum(isnan(bvpout{neps}))<0.25*lepoch
                                % for one subject technical problems with nexus-II 
                                [~, hf,hr2(neps)]= hrv_calc(bvpout{neps}', fs);
                                hrv2(neps) = hf;
                                
                            else
                                hrv2(neps) = NaN;
                                hr2(neps) = NaN;
                            end
                            
                            if any(isnan(respout{neps}))
                                resp2(neps) = NaN;
                                iqrr2(neps) = NaN;
                            else
                                [resp2(neps),iqrr2(neps)] = ...
                                    filtresp(double(respout{neps}), fs,...
                                    do_plots);
                            end
                            
                        end
                    end
                    freqpow_subj{ns, bn} = freqpow;
                    paf_subj{ns, bn} = paf;
                    pafsnr_subj{ns, bn} = pafsnr;
                    %                     slope_subj{ns, bn} = slope;
                    rrate_subj{ns, bn} = resp2;
                    iqrr_subj{ns,bn} = iqrr2;
                    hrate_subj{ns, bn} = hr2;
                    hrv_subj{ns, bn} = hrv2;
                    timefrend_subj{ns, bn} = timefrend;
                    timefrbeg_subj{ns, bn} = timefrbeg(:,2);
                else
                    fprintf('err')
                end
                
            end
            fprintf('\n')
        end
    end
    
end
cd(exp_medriv.session_dir)

%%% define variables to save
all_theta = [];
all_alpha = [];
all_beta = [];
all_gamma = [];
all_block = [];
all_group = [];
all_subj = [];
all_chans = [];
all_yearspract = [];
all_hourspract = [];
all_quest1 = [];
all_quest2 = [];
all_quest3 = [];
all_quest4 = [];
all_quest5 = [];
all_rrate = [];
all_iqrr = [];
all_hrate = [];
all_hrv = [];
all_btime = [];
all_scode = [];
all_timefrend1 = [];
all_timefrend2 = [];
all_timefrbeg = [];
all_paf = [];
all_pafsnr = [];
all_medi = [];
all_age = [];

age = [exp_medriv.data.subj_age];
group = [exp_medriv.data.group];
scode = [exp_medriv.data.subj_code];
subj_data = [exp_medriv.data];
years = [subj_data.years_practice];
hours = [subj_data.hours_practice];
% pafqual = [subj_data.paf_quality_med];

questresp = [subj_data.questresp];
questresp = reshape(questresp, [4 14 nsubj]);
questclust = exp_medriv.questclust;
nclust = max(questclust);
clustq = NaN(nclust, 4, nsubj);
for cq = 1:nclust
    clustq(cq, :, :) = nanmean(questresp(:, questclust==cq, :), 2);
end

medi = [subj_data.medi];
medi = reshape(medi, [4 nsubj]);


nbchan = 32;

for ns = 1:nsubj
    for nb = 1:nblocks
        if ~isempty(freqpow_subj{ns, nb})
            
            neps = size(freqpow_subj{ns, nb},1);
            
            fs = (squeeze(freqpow_subj{ns, nb}(:, 1, :)));
            all_theta = [all_theta; fs(:)];
            fs = (squeeze(freqpow_subj{ns, nb}(:, 2, :)));
            all_alpha = [all_alpha; fs(:)];
            fs = (squeeze(freqpow_subj{ns, nb}(:, 3, :)));
            all_beta = [all_beta; fs(:)];
            fs = (squeeze(freqpow_subj{ns, nb}(:, 4, :)));
            all_gamma = [all_gamma; fs(:)];
            
            ps = paf_subj{ns,nb};
            all_paf = [all_paf; ps(:)];
            
            ps = pafsnr_subj{ns,nb};
            all_pafsnr = [all_pafsnr; ps(:)];
            %
            
            ps = rrate_subj{ns,nb};
            ps = squeeze(repmat(ps, [1 1 nbchan]));
            all_rrate = [all_rrate; ps(:)];
            
            ps = iqrr_subj{ns,nb};
            ps = squeeze(repmat(ps, [1 1 nbchan]));
            all_iqrr = [all_iqrr; ps(:)];

            ps = hrate_subj{ns,nb};
            ps = squeeze(repmat(ps, [1 1 nbchan]));
            all_hrate = [all_hrate ;ps(:)];

            ps = hrv_subj{ns,nb};
            ps = squeeze(repmat(ps, [1 1 nbchan]));
            all_hrv = [all_hrv ;ps(:)];

            ps = repmat(nb, [neps nbchan]);
            all_block = [all_block ;ps(:)];
            
            ps = repmat(group(ns), [neps nbchan]);
            all_group = [all_group; ps(:)];
            
            ps = repmat(age(ns), [neps nbchan]);
            all_age = [all_age; ps(:)];
            
            ps = repmat(scode(ns), [neps nbchan]);
            all_scode = [all_scode; ps(:)];

            ps = repmat(ns, [neps nbchan]);
            all_subj = [all_subj; ps(:)];
            
            ps = squeeze(repmat(1:nbchan, [1 1 1 neps]))';
            all_chans = [all_chans ;ps(:)];

            
            ps = timefrend_subj{ns,nb};
            ps = permute(repmat(ps, [1 1 nbchan]), [1 3 2]);
            ps2 = ps(:,:,1);
            all_timefrend1 = [all_timefrend1; ps2(:)];
            ps2 = ps(:,:,2);
            all_timefrend2 = [all_timefrend2; ps2(:)];
            
            ps = timefrbeg_subj{ns,nb};
            ps = permute(repmat(ps, [1 1 nbchan]), [1 3 2]);
            all_timefrbeg = [all_timefrbeg; ps(:)];

            ps = repmat(years(ns), [neps nbchan]);
            all_yearspract = [all_yearspract; ps(:)];
            ps = repmat(hours(ns), [neps nbchan]);
            all_hourspract = [all_hourspract; ps(:)];
            
            ps = repmat(clustq(1,nb,ns), [neps nbchan]);
            all_quest1 = [all_quest1 ;ps(:)];
            ps = repmat(clustq(2,nb,ns), [ neps nbchan]);
            all_quest2 = [all_quest2 ;ps(:)];
            ps = repmat(clustq(3,nb,ns), [ neps nbchan]);
            all_quest3 = [all_quest3 ;ps(:)];
            ps = repmat(clustq(4,nb,ns), [ neps nbchan]);
            all_quest4 = [all_quest4 ;ps(:)];
            ps = repmat(clustq(5,nb,ns), [ neps nbchan]);
            all_quest5 = [all_quest5 ;ps(:)];
            
            ps = repmat(medi(nb,ns), [neps nbchan]);
            all_medi = [all_medi; ps(:)];
            
            ps = squeeze(repmat(1:neps, [1 1 nbchan]));
            all_btime = [all_btime; ps(:)];
            
        end
    end
end


dr = double([all_theta, all_alpha, all_beta, all_gamma,...
    all_rrate, all_iqrr, all_hrate, all_hrv, all_block,...
    all_group, all_subj, all_chans, all_medi, all_paf, all_pafsnr, ...
    all_quest1, all_quest2, all_quest3, all_quest4, all_quest5, all_hourspract...
    all_yearspract, all_btime, all_scode, all_timefrend1, all_timefrend2,...
    all_timefrbeg, all_age]);

oufilename = [filename '_freqs_' powmethname{powmeth} '_' num2str(windowsize) '.csv'];
csvwrite(oufilename, dr)

header_text = 'theta,alpha,beta,gamma,rrate,iqrr,hrate,hrv,block,group,';
header_text = [header_text 'subj,chans,medi,paf,pafsnr,q1,q2,q3,q4,q5,hourspract,yearspract,time,'];
header_text = [header_text 'subjcode,timefrend1,timefrend2,timefrbeg,age'];

add_header_csv(oufilename, header_text)

fprintf('saved %s\n', oufilename)
end

