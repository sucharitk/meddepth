function [freqpow, ft_ep, freqs] = freqrangepow(epoch, fs, freq_range, method)
%
% function freqpow = freqrangepow(epoch, fs, freq_range)
%
% calculate power for different freq bands for an epoch
%

% method = 1; %1- fft, 2-pwelch, 3-fft^2, 4, relative fft power
nfreqs = size(freq_range, 1);
if ~exist('freq_range', 'var') || isempty(freq_range)
    freq_range = [.5 3; 4 7; 7.5 13; 15 25; 30 49];
end

nbchans = size(epoch, 1);
freqpow = NaN(nfreqs, nbchans);

switch method
    case {1, 3, 4}
        NFFT = 2^16;
        freqs = linspace(0, fs/2, NFFT/2+1);
        ft_ep = AbsFFT(epoch, fs, NFFT);
        
        if method>=3
            ft_ep = ft_ep.^2;
        end
        if method ==4
            ft_ep = ft_ep./repmat(sum(ft_ep,2), 1, NFFT/2+1);
        end
        for nf = 1:nfreqs
            ff_range = freqs>=freq_range(nf, 1) & freqs<=freq_range(nf, 2);
            freqpow(nf, :) = mean(ft_ep(:, ff_range), 2);
        end
    case 2
        [ft_ep,freqs] = pwelch(epoch',200,100,500,fs);
        ft_ep = ft_ep';
        for nf = 1:nfreqs
            ff_range = freqs>=freq_range(nf, 1) & freqs<=freq_range(nf, 2);
            freqpow(nf, :) = mean(ft_ep(:, ff_range), 2);
        end
end

