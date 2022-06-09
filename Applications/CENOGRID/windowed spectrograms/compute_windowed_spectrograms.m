% For the CENOGRID dataset we perform a windowed analysis with the spike 
% powerspectrum, in order to obtain a wavelet-like spectrum.


% Here: TRAJECTORY BASED ON PECUZAL, MCDTS & NO EMBEDDING (RAW TIME SERIES)


clear, clc

data = load("../data/detrended.txt");

t = data(:,1);
O18 = data(:,2);
C13 = data(:,3);

t = flipud(t);
O18 = flipud(O18);
C13 = flipud(C13);

%%%%%%%%% CAUTION %%%%%%%%%%%%% O18 and C13 are swapped!!!!! 
%% ---> actually O18 is C13 in this case and vice versa %%%%%%%%%%%%%%

%% Set parameters

epsilon = 0.05; % recurrence threshold
windowsize2 = 1000; % windowsize for trajectory
windowsize1 = 100; % windowsize for the spectrum
ws = 1; % windowstep
N = length(t); % length of the time series
taus = 0:100; % delay range passed to PECUZAL
lambda = 0.003;

M = length(1:ws:N-windowsize2);

spectrum_O18_pec = zeros(windowsize1,M); 
spectrum_C13_pec = zeros(windowsize1,M);
spectrum_O18_mcdts = zeros(windowsize1,M); 
spectrum_C13_mcdts = zeros(windowsize1,M);
spectrum_O18_ne = zeros(windowsize1,M); 
spectrum_C13_ne = zeros(windowsize1,M);

%% Compute spike powerspectra in windowed analysis

cnt = 1;
tic
for i = 1:ws:M
    display(i)
    % extract time series
    x = O18(i:i+windowsize2);
    y = C13(i:i+windowsize2);
    % 1) Pecuzal embedding gained from whole dataset analysis
    Yx = genembed(x,[0,13,7,10,56,27,3,5,77,42,20,17,15], ones(1,13));
    Yy = genembed(y,[0,13,7,10,56,27,3,5,77,42,20,17,15], ones(1,13));
    % compute RP and tau-rr
    RPx = rp(Yx,epsilon,'var');
    RPy = rp(Yy,epsilon,'var');
    tauRRx = tau_recurrence_rate(RPx);
    tauRRy = tau_recurrence_rate(RPy);
    % compute spike powerspectra
    spectrum_O18_pec(:,cnt) = get_spike_spectrum(tauRRx(1:windowsize1), lambda);
    spectrum_C13_pec(:,cnt) = get_spike_spectrum(tauRRy(1:windowsize1), lambda);
    
    % 2) MCDTS embedding gained from whole dataset for prediction
    Yx = embed(x,3,1);
    Yy = embed(y,3,1);
    % compute RP and tau-rr
    RPx = rp(Yx,epsilon,'var');
    RPy = rp(Yy,epsilon,'var');
    tauRRx = tau_recurrence_rate(RPx);
    tauRRy = tau_recurrence_rate(RPy);
    % compute spike powerspectra
    spectrum_O18_mcdts(:,cnt) = get_spike_spectrum(tauRRx(1:windowsize1), lambda);
    spectrum_C13_mcdts(:,cnt) = get_spike_spectrum(tauRRy(1:windowsize1), lambda);
    
    % 2) No embedding
    Yx = x;
    Yy = y;
    % compute RP and tau-rr
    RPx = rp(Yx,epsilon,'fan');
    RPy = rp(Yy,epsilon,'fan');
    tauRRx = tau_recurrence_rate(RPx);
    tauRRy = tau_recurrence_rate(RPy);
    % compute spike powerspectra
    spectrum_O18_ne(:,cnt) = get_spike_spectrum(tauRRx(1:windowsize1), lambda);
    spectrum_C13_ne(:,cnt) = get_spike_spectrum(tauRRy(1:windowsize1), lambda);
   
    cnt = cnt + 1;
    
end
toc
    
% save("spectrum_O18_pec.mat", "spectrum_O18_pec")
% save("spectrum_O18_mcdts.mat", "spectrum_O18_mcdts")
% save("spectrum_O18_ne.mat", "spectrum_O18_ne")
% save("spectrum_C13_pec.mat", "spectrum_C13_pec")
% save("spectrum_C13_mcdts.mat", "spectrum_C13_mcdts")
% save("spectrum_C13_ne.mat", "spectrum_C13_ne")

%% Load data

spectrum_O18_pec = load("spectrum_O18_pec.mat");
spectrum_O18_pec = spectrum_O18_pec.spectrum_O18_pec;
spectrum_O18_mcdts_mse = load("spectrum_O18_mcdts.mat");
spectrum_O18_mcdts_mse = spectrum_O18_mcdts_mse.spectrum_O18_mcdts;
spectrum_O18_ne = load("spectrum_O18_ne.mat");
spectrum_O18_ne = spectrum_O18_ne.spectrum_O18_ne;
spectrum_C13_pec = load("spectrum_C13_pec.mat");
spectrum_C13_pec = spectrum_C13_pec.spectrum_C13_pec;
spectrum_C13_mcdts_mse = load("spectrum_C13_mcdts.mat");
spectrum_C13_mcdts_mse = spectrum_C13_mcdts_mse.spectrum_C13_mcdts;
spectrum_C13_ne = load("spectrum_C13_ne.mat");
spectrum_C13_ne = spectrum_C13_ne.spectrum_C13_ne;


tt = 5000:5000:500000;

%% Plot the wavelet-like spectrograms

% Transform into logspace

len = 90;

[x1,t1] = logimage(fliplr(spectrum_O18_pec(1:len,:)),tt(1:len));
[x2,~] = logimage(fliplr(spectrum_O18_mcdts_mse(1:len,:)),tt(1:len));
[x3,~] = logimage(fliplr(spectrum_O18_ne(1:len,:)),tt(1:len));

[y1,~] = logimage(fliplr(spectrum_C13_pec(1:len,:)),flipud(tt(1:len)));
[y2,~] = logimage(fliplr(spectrum_C13_mcdts_mse(1:len,:)),flipud(tt(1:len)));
[y3,~] = logimage(fliplr(spectrum_C13_ne(1:len,:)),flipud(tt(1:len)));


fs = 22;
lw2 = 1;

len = 90;

% new t-axis (y-axis)
t2 = 10.^(t1);

cax_max = 0.013;

figure('Units','normalized','Position',[.01 .01 .99 .99])
subplot(211)
imagesc(t(1:12421),t2,x1)
xline(13.9,'r--','LineWidth',lw2)
xline(17,'r--','LineWidth',lw2)
xline(34,'r--','LineWidth',lw2)
xline(47,'r--','LineWidth',lw2)
xline(56,'r--','LineWidth',lw2)
xline(61,'r--','LineWidth',lw2)
yline(21000,'k--','LineWidth',lw2)
yline(41000,'k--','LineWidth',lw2)
yline(100000,'k--','LineWidth',lw2)
yline(405000,'k--','LineWidth',lw2)
colormap(parula)
caxis([0 cax_max])
title('\delta^{18}O (PECUZAL)')
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
set(gca,'FontSize',fs)
set(gca,'YScale','log')
set(gca,'XDir','reverse')
set(gca,'YDir','normal')
grid on
subplot(212)
imagesc(t(1:12421),t2,y1)
xline(13.9,'r--','LineWidth',lw2)
xline(17,'r--','LineWidth',lw2)
xline(34,'r--','LineWidth',lw2)
xline(47,'r--','LineWidth',lw2)
xline(56,'r--','LineWidth',lw2)
xline(61,'r--','LineWidth',lw2)
yline(21000,'k--','LineWidth',lw2)
yline(41000,'k--','LineWidth',lw2)
yline(100000,'k--','LineWidth',lw2)
yline(405000,'k--','LineWidth',lw2)
colormap(parula)
caxis([0 cax_max])
title('\delta^{13}C (PECUZAL)')
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
set(gca,'FontSize',fs)
set(gca,'YScale','log')
set(gca,'XDir','reverse')
set(gca,'YDir','normal')
grid on

figure('Units','normalized','Position',[.01 .01 .99 .99])
subplot(211)
imagesc(t(1:12421),t2,x2)
xline(13.9,'r--','LineWidth',lw2)
xline(17,'r--','LineWidth',lw2)
xline(34,'r--','LineWidth',lw2)
xline(47,'r--','LineWidth',lw2)
xline(56,'r--','LineWidth',lw2)
xline(61,'r--','LineWidth',lw2)
yline(21000,'k--','LineWidth',lw2)
yline(41000,'k--','LineWidth',lw2)
yline(100000,'k--','LineWidth',lw2)
yline(405000,'k--','LineWidth',lw2)
colormap(parula)
caxis([0 cax_max])
title('\delta^{18}O (MCDTS MSE)')
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
set(gca,'FontSize',fs)
set(gca,'YScale','log')
set(gca,'XDir','reverse')
set(gca,'YDir','normal')
grid on
subplot(212)
imagesc(t(1:12421),t2,y2)
xline(13.9,'r--','LineWidth',lw2)
xline(17,'r--','LineWidth',lw2)
xline(34,'r--','LineWidth',lw2)
xline(47,'r--','LineWidth',lw2)
xline(56,'r--','LineWidth',lw2)
xline(61,'r--','LineWidth',lw2)
yline(21000,'k--','LineWidth',lw2)
yline(41000,'k--','LineWidth',lw2)
yline(100000,'k--','LineWidth',lw2)
yline(405000,'k--','LineWidth',lw2)
colormap(parula)
caxis([0 cax_max])
title('\delta^{13}C (MCDTS MSE)')
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
set(gca,'FontSize',fs)
set(gca,'YScale','log')
set(gca,'XDir','reverse')
set(gca,'YDir','normal')
grid on


figure('Units','normalized','Position',[.01 .01 .99 .99])
subplot(211)
imagesc(t(1:12421),t2,x3)
xline(13.9,'r--','LineWidth',lw2)
xline(17,'r--','LineWidth',lw2)
xline(34,'r--','LineWidth',lw2)
xline(47,'r--','LineWidth',lw2)
xline(56,'r--','LineWidth',lw2)
xline(61,'r--','LineWidth',lw2)
yline(21000,'k--','LineWidth',lw2)
yline(41000,'k--','LineWidth',lw2)
yline(100000,'k--','LineWidth',lw2)
yline(405000,'k--','LineWidth',lw2)
colormap(parula)
caxis([0 cax_max])
title('\delta^{18}O (no embedding)')
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
set(gca,'FontSize',fs)
set(gca,'YScale','log')
set(gca,'XDir','reverse')
set(gca,'YDir','normal')
grid on
subplot(212)
imagesc(t(1:12421),t2,y3)
xline(13.9,'r--','LineWidth',lw2)
xline(17,'r--','LineWidth',lw2)
xline(34,'r--','LineWidth',lw2)
xline(47,'r--','LineWidth',lw2)
xline(56,'r--','LineWidth',lw2)
xline(61,'r--','LineWidth',lw2)
yline(21000,'k--','LineWidth',lw2)
yline(41000,'k--','LineWidth',lw2)
yline(100000,'k--','LineWidth',lw2)
yline(405000,'k--','LineWidth',lw2)
colormap(parula)
caxis([0 cax_max])
title('\delta^{13}C (no embedding)')
set(gca,'FontSize',fs)
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
set(gca,'YScale','log')
set(gca,'XDir','reverse')
set(gca,'YDir','normal')
grid on
 

%% Plot FFT and spike spectra wavelets together

dt = abs(mean(diff(t))) * 100000; % Sampling period
Fs = 1/dt; % Sampling frequency
                
spec_raw = load("periodograms_raw.mat");
f = spec_raw.f;
spec_O18_raw = spec_raw.pxx_raw;
spec_C13_raw = spec_raw.pxy_raw;

spec_filtered = load("periodograms_filtered.mat");
spec_O18_filtered = spec_filtered.pxx_filtered;
spec_C13_filtered = spec_filtered.pxy_filtered;

%% FFT- and spike spectra for raw data
fs = 26;
lw1 = 1;
lw2 = 2.5;

factor = 0.01;

xlines = [13.9 17 34 47 56 61];
ylines1 = [1/21000 1/41000 1/100000 1/405000];
ylines2 = [21000 41000 100000 405000];
ylimms1 = [1/500000 6/100000];
ylimms2 = [10000 450000];

figure('Units','normalized','Position',[.01 .01 .99 .99])
subplot(211)
imagesc(t(1:12421),f(2:end),spec_C13_raw(2:end,:))
for i = 1:length(xlines)
    xline(xlines(i),'r--','LineWidth',lw2)
end
for i = 1:length(ylines1)
    yline(ylines1(i),'k--','LineWidth',lw1)
end
colormap(parula)
caxis([0 factor*(max(spec_C13_raw(:)))])
title('\delta^{13}C (evolutionary FFT spectrogram)')
ylabel('Frequency [1/yrs] ')
ylim(ylimms1)
set(gca,'FontSize',fs,'XDir','reverse')
grid on

subplot(212)
imagesc(t(1:12421),t2,y1)
for i = 1:length(xlines)
    xline(xlines(i),'r--','LineWidth',lw2)
end
for i = 1:length(ylines2)
    yline(ylines2(i),'k--','LineWidth',lw1)
end
colormap(parula)
caxis([0 cax_max])
title('\delta^{13}C (evolutionary spike spectrogram based on PECUZAL embedding)')
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
ylim(ylimms2)
set(gca,'FontSize',fs,'YScale','log','XDir','reverse','YDir','normal')
grid on


figure('Units','normalized','Position',[.01 .01 .99 .99])
subplot(211)
imagesc(t(1:12421),f(2:end),spec_O18_raw(2:end,:))
for i = 1:length(xlines)
    xline(xlines(i),'r--','LineWidth',lw2)
end
for i = 1:length(ylines1)
    yline(ylines1(i),'k--','LineWidth',lw1)
end
colormap(parula)
caxis([0 factor*(max(spec_O18_raw(:)))])
title('\delta^{18}O (evolutionary FFT spectrogram)')
ylabel('Frequency [1/yrs] ')
ylim(ylimms1)
set(gca,'FontSize',fs,'XDir','reverse')
grid on

subplot(212)
imagesc(t(1:12421),t2,x1)
for i = 1:length(xlines)
    xline(xlines(i),'r--','LineWidth',lw2)
end
for i = 1:length(ylines2)
    yline(ylines2(i),'k--','LineWidth',lw1)
end
colormap(parula)
caxis([0 cax_max])
title('\delta^{18}O (evolutionary spike spectrogram based on PECUZAL embedding)')
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
ylim(ylimms2)
set(gca,'FontSize',fs,'YScale','log','XDir','reverse','YDir','normal')
grid on


figure('Units','normalized','Position',[.01 .01 .99 .99])
subplot(211)
imagesc(t(1:12421),f(2:end),spec_C13_raw(2:end,:))
for i = 1:length(xlines)
    xline(xlines(i),'r--','LineWidth',lw2)
end
for i = 1:length(ylines1)
    yline(ylines1(i),'k--','LineWidth',lw1)
end
colormap(parula)
caxis([0 factor*(max(spec_C13_raw(:)))])
title('\delta^{13}C (evolutionary FFT spectrogram)')
ylabel('Frequency [1/yrs] ')
ylim(ylimms1)
set(gca,'FontSize',fs,'XDir','reverse')
grid on

subplot(212)
imagesc(t(1:12421),t2,y2)
for i = 1:length(xlines)
    xline(xlines(i),'r--','LineWidth',lw2)
end
for i = 1:length(ylines2)
    yline(ylines2(i),'k--','LineWidth',lw1)
end
colormap(parula)
caxis([0 cax_max])
title('\delta^{13}C (evolutionary spike spectrogram based on MCDTS embedding)')
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
ylim(ylimms2)
set(gca,'FontSize',fs,'YScale','log','XDir','reverse','YDir','normal')
grid on


figure('Units','normalized','Position',[.01 .01 .99 .99])
subplot(211)
imagesc(t(1:12421),f(2:end),spec_O18_raw(2:end,:))
for i = 1:length(xlines)
    xline(xlines(i),'r--','LineWidth',lw2)
end
for i = 1:length(ylines1)
    yline(ylines1(i),'k--','LineWidth',lw1)
end
colormap(parula)
caxis([0 factor*(max(spec_O18_raw(:)))])
title('\delta^{18}O (evolutionary FFT spectrogram)')
ylabel('Frequency [1/yrs] ')
ylim(ylimms1)
set(gca,'FontSize',fs,'XDir','reverse')
grid on

subplot(212)
imagesc(t(1:12421),t2,x2)
for i = 1:length(xlines)
    xline(xlines(i),'r--','LineWidth',lw2)
end
for i = 1:length(ylines2)
    yline(ylines2(i),'k--','LineWidth',lw1)
end
colormap(parula)
caxis([0 cax_max])
title('\delta^{18}O (evolutionary spike spectrogram based on MCDTS embedding)')
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
ylim(ylimms2)
set(gca,'FontSize',fs,'YScale','log','XDir','reverse','YDir','normal')
grid on


figure('Units','normalized','Position',[.01 .01 .99 .99])
subplot(211)
imagesc(t(1:12421),f(2:end),spec_C13_raw(2:end,:))
for i = 1:length(xlines)
    xline(xlines(i),'r--','LineWidth',lw2)
end
for i = 1:length(ylines1)
    yline(ylines1(i),'k--','LineWidth',lw1)
end
colormap(parula)
caxis([0 factor*(max(spec_C13_raw(:)))])
title('\delta^{13}C (evolutionary FFT spectrogram)')
ylabel('Frequency [1/yrs] ')
ylim(ylimms1)
set(gca,'FontSize',fs,'XDir','reverse')
grid on

subplot(212)
imagesc(t(1:12421),t2,y3)
for i = 1:length(xlines)
    xline(xlines(i),'r--','LineWidth',lw2)
end
for i = 1:length(ylines2)
    yline(ylines2(i),'k--','LineWidth',lw1)
end
colormap(parula)
caxis([0 cax_max])
title('\delta^{13}C (evolutionary spike spectrogram based on non embedded time series)')
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
ylim(ylimms2)
set(gca,'FontSize',fs,'YScale','log','XDir','reverse','YDir','normal')
grid on



figure('Units','normalized','Position',[.01 .01 .99 .99])
subplot(211)
imagesc(t(1:12421),f(2:end),spec_O18_raw(2:end,:))
for i = 1:length(xlines)
    xline(xlines(i),'r--','LineWidth',lw2)
end
for i = 1:length(ylines1)
    yline(ylines1(i),'k--','LineWidth',lw1)
end
colormap(parula)
caxis([0 factor*(max(spec_O18_raw(:)))])
title('\delta^{18}O (evolutionary FFT spectrogram)')
ylabel('Frequency [1/yrs] ')
ylim(ylimms1)
set(gca,'FontSize',fs,'XDir','reverse')
grid on

subplot(212)
imagesc(t(1:12421),t2,x3)
for i = 1:length(xlines)
    xline(xlines(i),'r--','LineWidth',lw2)
end
for i = 1:length(ylines2)
    yline(ylines2(i),'k--','LineWidth',lw1)
end
colormap(parula)
caxis([0 cax_max])
title('\delta^{18}O (evolutionary spike spectrogram based on non embedded time series)')
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
ylim(ylimms2)
set(gca,'FontSize',fs,'YScale','log','XDir','reverse','YDir','normal')
grid on
