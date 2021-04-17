function [chan_inds, chanind_pos] = ...
    get_channels_from_labels(chanlabs_from, chanlabs_select)
%
% from a larger chanlocs structure get a subset of chanloc indices based on
% the chanlabels
%

nc = numel(chanlabs_from);
chan_inds = false(1, nc);
nc2 = numel(chanlabs_select);
chanind_pos = [];

for cn = 1:nc2
    fcn = strcmp(chanlabs_from, chanlabs_select{cn});
    %     if any(fcn)
    chan_inds = chan_inds | fcn;
    if nargout>1
        chanind_pos = [chanind_pos find(fcn)];
    end
    %     end
end

end