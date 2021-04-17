function [rr1,iqrr] = filtresp(D, Fs,do_plots)

if ~exist('do_plots', 'var'),do_plots=false; end

min_peak_dist = 2; % sec
Fn = Fs/2;                                                  % Nyquist Frequency
Ts = 1/Fs;                                                  % Sampling Time (sec)
L = numel(D);
% tdelay = 0;
% transtimes=transtimes-tdelay;
t = linspace(0, L, L)*Ts;                                   % Time Vector (sec)

% FTD = fft(D-mean(D))/L;                                     % Fourier Transform
% Fv = linspace(0, 1, fix(L/2)+1)*Fn;                         % Frequency Vector
% Iv = 1:numel(Fv);                                           % Index Vector

Wp = [0.040 0.35]/Fn;                                        % Passband Frequency (Normalised)
Ws = [0.035 0.40]/Fn;                                        % Stopband Frequency (Normalised)
Rp =   1;                                                   % Passband Ripple (dB)
Rs =  50;                                                   % Stopband Ripple (dB)
[n,Ws]  = cheb2ord(Wp,Ws,Rp,Rs);                            % Filter Order
[z,p,k] = cheby2(n,Rs,Ws);                                  % Filter Design, Sepcify Bandpass
[sos,g] = zp2sos(z,p,k);                                    % Convert To Second-Order-Section For Stability
% figure(3)
% freqz(sos, 2^16, Fs)                                        % Filter Bode Plot

meand = mean(D);
stdd = std(D);
D2 = D;
dlims = [meand-stdd*2.5 meand+stdd*2.5];
D2(D<dlims(1)) = dlims(1);
D2(D>dlims(2)) = dlims(2);


dfilt = filtfilt(sos, g, D2);                           % Filter Signal

[pks,locs] = findpeaks(dfilt, 'MinPeakDist',min_peak_dist*Fs);
locs = locs(pks>0);
pks = pks(pks>0);
% ntt = numel(transtimes); ttrans = NaN(1, ntt);
% for nt = 1:ntt
%     [~,ttrans(nt)] = min(abs(t-transtimes(nt)));
% end

tdif = diff([0 t(locs)]);                                   % Time Difference Between Peaks (sec)
drate = 60./tdif;                                            % Frequency (Respirations/Minute)
rr2 = (numel(locs)/t(end))*60;
rr1 = median(drate);

% absftd = abs(FTD(Iv))*2;
% [~,mxind]=max(absftd);

% rr2 = Fv(mxind)*60;
% phase = angle(detrend(hilbert(dfilt)));
% phas_trans = phase(ttrans);
iqrr = iqr(drate);
if do_plots
    figure(2)
    subplot(211)
    plot(t, dfilt, t(locs), pks, '^r')
    grid
    xlabel('Time')
    ylabel('Amplitude')
    title(sprintf('rr1 = %.4g,      rr2 = %0.4g,     iqr = %.4g',rr1, rr2, iqrr))
    
    subplot(212)
    plot(t, D, t, D2)


end
if iqrr>10, rr1 = NaN; end
end