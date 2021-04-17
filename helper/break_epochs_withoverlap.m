function [epout, timefrend] = break_epochs_withoverlap(epoch, fs, windowsize, noverlap)

% break an EEG epoch into smaller windows starting from the end

windowsize = windowsize*fs; % size of the each epoch in seconds
noverlap = noverlap*fs; % how much to not overlap between epochs

szep = size(epoch);
numdims = numel(szep);
pdims = 1:numdims; pdims(1) = numdims; pdims(end) = 1; % swap the first and last dimension
epoch = permute(epoch, pdims);
szep = szep(end);

epsizes = NaN(600,2);
for ii=1:size(epsizes, 1)
    epsizes(ii, :) = [szep-(ii-1)*noverlap szep-(ii-1)*noverlap-windowsize];
end
if epsizes(end,2)>0, disp('not capturing everything'); end
epsizes = epsizes(epsizes(:, 2)>0, :);
nepsize = size(epsizes,1);
epout = cell(1, nepsize);

for neps = 1:nepsize
     epo = epoch(epsizes(neps,2):epsizes(neps,1), :, :, :, :, :);
     epout{neps} = permute(epo, pdims); % reswap first and last dimensions
end

timefrend = (szep-epsizes)/fs;
end