% For the CENOGRID dataset we perform a windowed analysis with the spike 
% powerspectrum, in order to obtain a wavelet-like spectrum.

% Here: TRAJECTORY BASED ON CONCATENATION OF O18 AND C13 (FILTERED TIME SERIES)


clear, clc

data = load("../data/detrended.txt");

t = data(:,1);
O18_raw = flipud(data(:,2));
C13_raw = flipud(data(:,3));
t = flipud(t);

O18 = load("../data/O18_filtered_time_reverse.txt");
C13 =load("../data/C13_filtered_time_reverse.txt");


%% Set parameters

epsilon = 0.05; % recurrence threshold
windowsize2 = 1000; % windowsize for trajectory
windowsize1 = 100; % windowsize for the spectrum
ws = 1; % windowstep
N = length(t); % length of the time series
taus = 0:100; % delay range passed to PECUZAL
lambda = 0.003;

M = length(1:ws:N-windowsize2);

spectrum_raw = zeros(windowsize1,M); 
spectrum_filtered = zeros(windowsize1,M); 

%% Compute spike powerspectra in windowed analysis

cnt = 1;
tic
for i = 1:ws:M
    display(i)
    % extract time series
    x = O18_raw(i:i+windowsize2);
    y = C13_raw(i:i+windowsize2);
    Yx = horzcat(x,y);
    x = O18(i:i+windowsize2);
    y = C13(i:i+windowsize2);
    Yy = horzcat(x,y);
    
    % compute RP and tau-rr
    RPx = rp(Yx,epsilon,'var');
    RPy = rp(Yy,epsilon,'var');
    tauRRx = tau_recurrence_rate(RPx);
    tauRRy = tau_recurrence_rate(RPy);
    % compute spike powerspectra
    spectrum_raw(:,cnt) = get_spike_spectrum(tauRRx(1:windowsize1), lambda);
    spectrum_filtered(:,cnt) = get_spike_spectrum(tauRRy(1:windowsize1), lambda);
    
    cnt = cnt + 1;
    
end
toc
    
save("mixed_spectrum_raw.mat", "spectrum_raw")
save("mixed_spectrum_filtered.mat", "spectrum_filtered")

%% Load data

clear, clc

spectrum_raw = load("mixed_spectrum_raw.mat");
spectrum_raw = spectrum_raw.spectrum_raw;
spectrum_filtered = load("mixed_spectrum_filtered.mat");
spectrum_filtered = spectrum_filtered.spectrum_filtered;


%% Plot the wavelet-like spectrograms

tt = 5000:5000:500000;

% Transform into logspace

len = 90;

[x1,t1] = logimage(spectrum_raw(1:len,:),tt(1:len));
[x2,~] = logimage(spectrum_filtered(1:len,:),tt(1:len));

fs = 22;

% new t-axis (y-axis)
t2 = 10.^(t1);

figure('Units','normalized','Position',[.01 .01 .99 .99])
subplot(211)
imagesc(t(1:12421),t2,x1)
colormap(parula)
caxis([0 0.013])
title('\delta^{18}O & \delta^{13}C composite spike wavelet (raw)')
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
set(gca,'Yscale','log')
set(gca,'FontSize',fs)
set(gca,'XDir','reverse')
set(gca,'YDir','normal')
grid on

subplot(212)
imagesc(t(1:12421),t2,x2)
colormap(parula)
caxis([0 0.013])
title('\delta^{18}O & \delta^{13}C composite spike wavelet (filtered)')
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
set(gca,'Yscale','log')
set(gca,'FontSize',fs)
set(gca,'XDir','reverse')
set(gca,'YDir','normal')
grid on

