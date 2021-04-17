function exp_medriv = medriv_addquestdata(exp_medriv)
cd(exp_medriv.session_dir)

% questionnaire data

load('medeq')
subjids = medeq{:,1};
szs = size(medeq);
nsi = szs(1);
subj_code = [exp_medriv.data.subj_code];
scode = cell(1, numel(subj_code));
exp_medriv.questclust = [3, 1, 5, 2, 3, 4, 1, 4, 4, 4, 5, 4, 5, 4]; %
questdir = exp_medriv.questclust-2; 
questdir(questdir>=0) = 1;
for ns = 1:numel(subj_code)
    scode{ns} = num2str(subj_code(ns));
    exp_medriv.data(ns).questresp(1:4, 1:szs(2)-1) = NaN;
    exp_medriv.data(ns).medi(1:4) = NaN;
end
for ns = 1:nsi
    subjnum = ~isemptycell(strfind(scode, subjids{ns}(3:end)));
    if any(subjnum)
        qresp = NaN(1, szs(2)-1);
        epnum = str2double(subjids{ns}(1));
        %         if epnum~=2
        %             if epnum>1, epnum = epnum-1; end
        for qn = 2:szs(2)
            qresp(qn-1) = medeq{ns, qn};
        end
        exp_medriv.data(subjnum).questresp(epnum, :) = qresp;
        medi = qresp*questdir';
        exp_medriv.data(subjnum).medi(epnum) = medi;
        %         end
    end
end

% cluster according to the medeq paper
% exp_medriv.questclust = 4*ones(1, 14);
% exp_medriv.questclust([1 5]) = 1; % cluster according to the question of present moment awareness
% exp_medriv.questclust([3 11 13]) = 2; % cluster according to the question of present moment awareness
% exp_medriv.questclust([2]) = 3; % cluster according to the question of present moment awareness

fprintf('questionnaire data saved\n')
end