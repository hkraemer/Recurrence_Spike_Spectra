% For the CENOGRID dataset we perform a windowed analysis with the spike 
% powerspectrum, in order to obtain a wavelet-like spectrum.


% Here: TRAJECTORY BASED ON KENNEL, MCDTS & NO EMBEDDING (FILTERED TIME SERIES)

clear, clc

data = load("../data/detrended.txt");

t = data(:,1);
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

spectrum_O18_ken = zeros(windowsize1,M); 
spectrum_C13_ken = zeros(windowsize1,M);
spectrum_O18_mcdts_mse = zeros(windowsize1,M); 
spectrum_C13_mcdts_mse = zeros(windowsize1,M);
spectrum_O18_mcdts_mse_KL = zeros(windowsize1,M); 
spectrum_C13_mcdts_mse_KL = zeros(windowsize1,M);
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
    X = horzcat(x,y);
    % 1) Kennel embedding gained from whole dataset analysis
    Yx = genembed(x,[0,51,102,153,204,255],[1,1,1,1,1,1]);
    Yy = genembed(y,[0,51,102,153,204,255],[1,1,1,1,1,1]);
    % compute RP and tau-rr
    RPx = rp(Yx,epsilon,'var');
    RPy = rp(Yy,epsilon,'var');
    tauRRx = tau_recurrence_rate(RPx);
    tauRRy = tau_recurrence_rate(RPy);
    % compute spike powerspectra
    spectrum_O18_ken(:,cnt) = get_spike_spectrum(tauRRx(1:windowsize1), lambda);
    spectrum_C13_ken(:,cnt) = get_spike_spectrum(tauRRy(1:windowsize1), lambda);
    
    % 2) MCDTS embedding gained from whole dataset for prediction
    Yx = genembed(x,[0,1,2],[1,1,1]);
    Yy = genembed(y,[0,1,2],[1,1,1]);
    % compute RP and tau-rr
    RPx = rp(Yx,epsilon,'var');
    RPy = rp(Yy,epsilon,'var');
    tauRRx = tau_recurrence_rate(RPx);
    tauRRy = tau_recurrence_rate(RPy);
    % compute spike powerspectra
    spectrum_O18_mcdts_mse(:,cnt) = get_spike_spectrum(tauRRx(1:windowsize1), lambda);
    spectrum_C13_mcdts_mse(:,cnt) = get_spike_spectrum(tauRRy(1:windowsize1), lambda);
    
    % 3) MCDTS embedding gained from whole dataset for prediction
    Yx = genembed(X,[0,8,83,135,67,131],[1,2,1,1,1,1]); %MCDTS-C-MSE-KL (m)
    Yy = genembed(y,[0,67,92,116,15],[1,1,1,1,1]); %MCDTS-C-MSE-KL
    % compute RP and tau-rr
    RPx = rp(Yx,epsilon,'var');
    RPy = rp(Yy,epsilon,'var');
    tauRRx = tau_recurrence_rate(RPx);
    tauRRy = tau_recurrence_rate(RPy);
    % compute spike powerspectra
    spectrum_O18_mcdts_mse_KL(:,cnt) = get_spike_spectrum(tauRRx(1:windowsize1), lambda);
    spectrum_C13_mcdts_mse_KL(:,cnt) = get_spike_spectrum(tauRRy(1:windowsize1), lambda);
    
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
    
save("filtered_spectrum_O18_ken.mat", "spectrum_O18_ken")
save("filtered_spectrum_O18_mcdts_mse.mat", "spectrum_O18_mcdts_mse")
save("filtered_spectrum_O18_mcdts_mse_KL.mat", "spectrum_O18_mcdts_mse_KL")
save("filtered_spectrum_O18_ne.mat", "spectrum_O18_ne")
save("filtered_spectrum_C13_ken.mat", "spectrum_C13_ken")
save("filtered_spectrum_C13_mcdts_mse.mat", "spectrum_C13_mcdts_mse")
save("filtered_spectrum_C13_mcdts_mse_KL.mat", "spectrum_C13_mcdts_mse_KL")
save("filtered_spectrum_C13_ne.mat", "spectrum_C13_ne")

%% Load data

spectrum_O18_ken = load("filtered_spectrum_O18_ken.mat");
spectrum_O18_ken = spectrum_O18_ken.spectrum_O18_ken;
spectrum_O18_mcdts_mse = load("filtered_spectrum_O18_mcdts_mse.mat");
spectrum_O18_mcdts_mse = spectrum_O18_mcdts_mse.spectrum_O18_mcdts_mse;
spectrum_O18_mcdts_mse_KL = load("filtered_spectrum_O18_mcdts_mse_KL.mat");
spectrum_O18_mcdts_mse_KL = spectrum_O18_mcdts_mse_KL.spectrum_O18_mcdts_mse_KL;
spectrum_O18_ne = load("filtered_spectrum_O18_ne.mat");
spectrum_O18_ne = spectrum_O18_ne.spectrum_O18_ne;
spectrum_C13_ken = load("filtered_spectrum_C13_ken.mat");
spectrum_C13_ken = spectrum_C13_ken.spectrum_C13_ken;
spectrum_C13_mcdts_mse = load("filtered_spectrum_C13_mcdts_mse.mat");
spectrum_C13_mcdts_mse = spectrum_C13_mcdts_mse.spectrum_C13_mcdts_mse;
spectrum_C13_mcdts_mse_KL = load("filtered_spectrum_C13_mcdts_mse_KL.mat");
spectrum_C13_mcdts_mse_KL = spectrum_C13_mcdts_mse_KL.spectrum_C13_mcdts_mse_KL;
spectrum_C13_ne = load("filtered_spectrum_C13_ne.mat");
spectrum_C13_ne = spectrum_C13_ne.spectrum_C13_ne;

tt = 5000:5000:500000;

%% Plot the wavelet-like spectrograms

% Transform into logspace

len = 90;

[x1,t1] = logimage(spectrum_O18_ken(1:len,:),tt(1:len));
[x2,~] = logimage(spectrum_O18_mcdts_mse(1:len,:),tt(1:len));
[x3,~] = logimage(spectrum_O18_mcdts_mse_KL(1:len,:),tt(1:len));
[x4,~] = logimage(spectrum_O18_ne(1:len,:),flipud(tt(1:len)));

[y1,~] = logimage(spectrum_C13_ken(1:len,:),flipud(tt(1:len)));
[y2,~] = logimage(spectrum_C13_mcdts_mse(1:len,:),flipud(tt(1:len)));
[y3,~] = logimage(spectrum_C13_mcdts_mse_KL(1:len,:),flipud(tt(1:len)));
[y4,~] = logimage(spectrum_C13_ne(1:len,:),flipud(tt(1:len)));


fs = 22;
lw2 = 2;

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
yline(21000,'r.-','LineWidth',lw2)
yline(41000,'r.-','LineWidth',lw2)
yline(100000,'r.-','LineWidth',lw2)
yline(405000,'r.-','LineWidth',lw2)
colormap(parula)
caxis([0 cax_max])
title('\delta^{18}O (Kennel)')
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
yline(21000,'r.-','LineWidth',lw2)
yline(41000,'r.-','LineWidth',lw2)
yline(100000,'r.-','LineWidth',lw2)
yline(405000,'r.-','LineWidth',lw2)
colormap(parula)
caxis([0 cax_max])
title('\delta^{13}C (Kennel)')
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
yline(21000,'r.-','LineWidth',lw2)
yline(41000,'r.-','LineWidth',lw2)
yline(100000,'r.-','LineWidth',lw2)
yline(405000,'r.-','LineWidth',lw2)
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
yline(21000,'r.-','LineWidth',lw2)
yline(41000,'r.-','LineWidth',lw2)
yline(100000,'r.-','LineWidth',lw2)
yline(405000,'r.-','LineWidth',lw2)
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
yline(21000,'r.-','LineWidth',lw2)
yline(41000,'r.-','LineWidth',lw2)
yline(100000,'r.-','LineWidth',lw2)
yline(405000,'r.-','LineWidth',lw2)
colormap(parula)
caxis([0 cax_max])
title('\delta^{18}O (MCDTS MSE KL)')
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
yline(21000,'r.-','LineWidth',lw2)
yline(41000,'r.-','LineWidth',lw2)
yline(100000,'r.-','LineWidth',lw2)
yline(405000,'r.-','LineWidth',lw2)
colormap(parula)
caxis([0 cax_max])
title('\delta^{13}C (MCDTS MSE KL)')
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
set(gca,'FontSize',fs)
set(gca,'YScale','log')
set(gca,'XDir','reverse')
set(gca,'YDir','normal')
grid on

figure('Units','normalized','Position',[.01 .01 .99 .99])
subplot(211)
imagesc(t(1:12421),t2,x4)
xline(13.9,'r--','LineWidth',lw2)
xline(17,'r--','LineWidth',lw2)
xline(34,'r--','LineWidth',lw2)
xline(47,'r--','LineWidth',lw2)
xline(56,'r--','LineWidth',lw2)
xline(61,'r--','LineWidth',lw2)
yline(21000,'r.-','LineWidth',lw2)
yline(41000,'r.-','LineWidth',lw2)
yline(100000,'r.-','LineWidth',lw2)
yline(405000,'r.-','LineWidth',lw2)
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
imagesc(t(1:12421),t2,y4)
xline(13.9,'r--','LineWidth',lw2)
xline(17,'r--','LineWidth',lw2)
xline(34,'r--','LineWidth',lw2)
xline(47,'r--','LineWidth',lw2)
xline(56,'r--','LineWidth',lw2)
xline(61,'r--','LineWidth',lw2)
yline(21000,'r.-','LineWidth',lw2)
yline(41000,'r.-','LineWidth',lw2)
yline(100000,'r.-','LineWidth',lw2)
yline(405000,'r.-','LineWidth',lw2)
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

spec_raw = load("periodograms_raw.mat");
% periods = spec_raw.periods;
% periods = fliplr(periods);
% periods = periods(1:503);
spec_O18_raw = spec_raw.pxx_raw;
spec_C13_raw = spec_raw.pxy_raw;
spec_O18_raw = spec_O18_raw(1:503,:);
spec_C13_raw = spec_C13_raw(1:503,:);

spec_filtered = load("periodograms_filtered.mat");
periods = spec_filtered.periods;
periods = fliplr(periods);
periods = periods(1:503);
spec_O18_filtered = spec_filtered.pxx_filtered;
spec_C13_filtered = spec_filtered.pxy_filtered;
spec_O18_filtered = spec_O18_filtered(1:503,:);
spec_C13_filtered = spec_C13_filtered(1:503,:);

%% FFT- and spike spectra for raw data
fs = 22;
lw2 = 2;

figure('Units','normalized','Position',[.01 .01 .99 .99])
subplot(211)
imagesc(t(1:12421),periods,flipud(spec_C13_raw))
xline(13.9,'r--','LineWidth',lw2)
xline(17,'r--','LineWidth',lw2)
xline(34,'r--','LineWidth',lw2)
xline(47,'r--','LineWidth',lw2)
xline(56,'r--','LineWidth',lw2)
xline(61,'r--','LineWidth',lw2)
yline(21000,'r.-','LineWidth',lw2)
yline(41000,'r.-','LineWidth',lw2)
yline(100000,'r.-','LineWidth',lw2)
yline(405000,'r.-','LineWidth',lw2)
colormap(parula)
caxis([0 0.005*(max(spec_C13_raw(:)))])
title('\delta^{13}C (normal FFT FILTERED)')
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
ylim([10000 450000])
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
yline(21000,'r.-','LineWidth',lw2)
yline(41000,'r.-','LineWidth',lw2)
yline(100000,'r.-','LineWidth',lw2)
yline(405000,'r.-','LineWidth',lw2)
colormap(parula)
caxis([0 cax_max])
title('\delta^{13}C (Kennel FILTERED)')
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
ylim([10000 450000])
set(gca,'FontSize',fs)
set(gca,'YScale','log')
set(gca,'XDir','reverse')
set(gca,'YDir','normal')
grid on


figure('Units','normalized','Position',[.01 .01 .99 .99])
subplot(211)
imagesc(t(1:12421),periods,flipud(spec_O18_raw))
xline(13.9,'r--','LineWidth',lw2)
xline(17,'r--','LineWidth',lw2)
xline(34,'r--','LineWidth',lw2)
xline(47,'r--','LineWidth',lw2)
xline(56,'r--','LineWidth',lw2)
xline(61,'r--','LineWidth',lw2)
yline(21000,'r.-','LineWidth',lw2)
yline(41000,'r.-','LineWidth',lw2)
yline(100000,'r.-','LineWidth',lw2)
yline(405000,'r.-','LineWidth',lw2)
colormap(parula)
caxis([0 0.01*(max(spec_O18_raw(:)))])
title('\delta^{18}O (normal FFT FILTERED)')
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
ylim([10000 450000])
set(gca,'FontSize',fs)
set(gca,'YScale','log')
set(gca,'XDir','reverse')
set(gca,'YDir','normal')
grid on
subplot(212)
imagesc(t(1:12421),t2,x1)
xline(13.9,'r--','LineWidth',lw2)
xline(17,'r--','LineWidth',lw2)
xline(34,'r--','LineWidth',lw2)
xline(47,'r--','LineWidth',lw2)
xline(56,'r--','LineWidth',lw2)
xline(61,'r--','LineWidth',lw2)
yline(21000,'r.-','LineWidth',lw2)
yline(41000,'r.-','LineWidth',lw2)
yline(100000,'r.-','LineWidth',lw2)
yline(405000,'r.-','LineWidth',lw2)
colormap(parula)
caxis([0 cax_max])
title('\delta^{18}O (Kennel FILTERED)')
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
ylim([10000 450000])
set(gca,'FontSize',fs)
set(gca,'YScale','log')
set(gca,'XDir','reverse')
set(gca,'YDir','normal')
grid on


figure('Units','normalized','Position',[.01 .01 .99 .99])
subplot(211)
imagesc(t(1:12421),periods,flipud(spec_C13_raw))
xline(13.9,'r--','LineWidth',lw2)
xline(17,'r--','LineWidth',lw2)
xline(34,'r--','LineWidth',lw2)
xline(47,'r--','LineWidth',lw2)
xline(56,'r--','LineWidth',lw2)
xline(61,'r--','LineWidth',lw2)
yline(21000,'r.-','LineWidth',lw2)
yline(41000,'r.-','LineWidth',lw2)
yline(100000,'r.-','LineWidth',lw2)
yline(405000,'r.-','LineWidth',lw2)
colormap(parula)
caxis([0 0.005*(max(spec_C13_raw(:)))])
title('\delta^{13}C (normal FFT FILTERED)')
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
ylim([10000 450000])
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
yline(21000,'r.-','LineWidth',lw2)
yline(41000,'r.-','LineWidth',lw2)
yline(100000,'r.-','LineWidth',lw2)
yline(405000,'r.-','LineWidth',lw2)
colormap(parula)
caxis([0 cax_max])
title('\delta^{13}C (MCDTS MSE FILTERED)')
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
ylim([10000 450000])
set(gca,'FontSize',fs)
set(gca,'YScale','log')
set(gca,'XDir','reverse')
set(gca,'YDir','normal')
grid on


figure('Units','normalized','Position',[.01 .01 .99 .99])
subplot(211)
imagesc(t(1:12421),periods,flipud(spec_O18_raw))
xline(13.9,'r--','LineWidth',lw2)
xline(17,'r--','LineWidth',lw2)
xline(34,'r--','LineWidth',lw2)
xline(47,'r--','LineWidth',lw2)
xline(56,'r--','LineWidth',lw2)
xline(61,'r--','LineWidth',lw2)
yline(21000,'r.-','LineWidth',lw2)
yline(41000,'r.-','LineWidth',lw2)
yline(100000,'r.-','LineWidth',lw2)
yline(405000,'r.-','LineWidth',lw2)
colormap(parula)
caxis([0 0.01*(max(spec_O18_raw(:)))])
title('\delta^{18}O (normal FFT FILTERED)')
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
ylim([10000 450000])
set(gca,'FontSize',fs)
set(gca,'YScale','log')
set(gca,'XDir','reverse')
set(gca,'YDir','normal')
grid on
subplot(212)
imagesc(t(1:12421),t2,x2)
xline(13.9,'r--','LineWidth',lw2)
xline(17,'r--','LineWidth',lw2)
xline(34,'r--','LineWidth',lw2)
xline(47,'r--','LineWidth',lw2)
xline(56,'r--','LineWidth',lw2)
xline(61,'r--','LineWidth',lw2)
yline(21000,'r.-','LineWidth',lw2)
yline(41000,'r.-','LineWidth',lw2)
yline(100000,'r.-','LineWidth',lw2)
yline(405000,'r.-','LineWidth',lw2)
colormap(parula)
caxis([0 cax_max])
title('\delta^{18}O (MCDTS MSE FILTERED)')
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
ylim([10000 450000])
set(gca,'FontSize',fs)
set(gca,'YScale','log')
set(gca,'XDir','reverse')
set(gca,'YDir','normal')
grid on


figure('Units','normalized','Position',[.01 .01 .99 .99])
subplot(211)
imagesc(t(1:12421),periods,flipud(spec_C13_raw))
xline(13.9,'r--','LineWidth',lw2)
xline(17,'r--','LineWidth',lw2)
xline(34,'r--','LineWidth',lw2)
xline(47,'r--','LineWidth',lw2)
xline(56,'r--','LineWidth',lw2)
xline(61,'r--','LineWidth',lw2)
yline(21000,'r.-','LineWidth',lw2)
yline(41000,'r.-','LineWidth',lw2)
yline(100000,'r.-','LineWidth',lw2)
yline(405000,'r.-','LineWidth',lw2)
colormap(parula)
caxis([0 0.005*(max(spec_C13_raw(:)))])
title('\delta^{13}C (normal FFT FILTERED)')
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
ylim([10000 450000])
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
yline(21000,'r.-','LineWidth',lw2)
yline(41000,'r.-','LineWidth',lw2)
yline(100000,'r.-','LineWidth',lw2)
yline(405000,'r.-','LineWidth',lw2)
colormap(parula)
caxis([0 cax_max])
title('\delta^{13}C (MCDTS MSE KL FILTERED)')
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
ylim([10000 450000])
set(gca,'FontSize',fs)
set(gca,'YScale','log')
set(gca,'XDir','reverse')
set(gca,'YDir','normal')
grid on


figure('Units','normalized','Position',[.01 .01 .99 .99])
subplot(211)
imagesc(t(1:12421),periods,flipud(spec_O18_raw))
xline(13.9,'r--','LineWidth',lw2)
xline(17,'r--','LineWidth',lw2)
xline(34,'r--','LineWidth',lw2)
xline(47,'r--','LineWidth',lw2)
xline(56,'r--','LineWidth',lw2)
xline(61,'r--','LineWidth',lw2)
yline(21000,'r.-','LineWidth',lw2)
yline(41000,'r.-','LineWidth',lw2)
yline(100000,'r.-','LineWidth',lw2)
yline(405000,'r.-','LineWidth',lw2)
colormap(parula)
caxis([0 0.01*(max(spec_O18_raw(:)))])
title('\delta^{18}O (normal FFT FILTERED)')
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
ylim([10000 450000])
set(gca,'FontSize',fs)
set(gca,'YScale','log')
set(gca,'XDir','reverse')
set(gca,'YDir','normal')
grid on
subplot(212)
imagesc(t(1:12421),t2,x3)
xline(13.9,'r--','LineWidth',lw2)
xline(17,'r--','LineWidth',lw2)
xline(34,'r--','LineWidth',lw2)
xline(47,'r--','LineWidth',lw2)
xline(56,'r--','LineWidth',lw2)
xline(61,'r--','LineWidth',lw2)
yline(21000,'r.-','LineWidth',lw2)
yline(41000,'r.-','LineWidth',lw2)
yline(100000,'r.-','LineWidth',lw2)
yline(405000,'r.-','LineWidth',lw2)
colormap(parula)
caxis([0 cax_max])
title('\delta^{18}O (MCDTS MSE KL FILTERED)')
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
ylim([10000 450000])
set(gca,'FontSize',fs)
set(gca,'YScale','log')
set(gca,'XDir','reverse')
set(gca,'YDir','normal')
grid on



figure('Units','normalized','Position',[.01 .01 .99 .99])
subplot(211)
imagesc(t(1:12421),periods,flipud(spec_C13_raw))
xline(13.9,'r--','LineWidth',lw2)
xline(17,'r--','LineWidth',lw2)
xline(34,'r--','LineWidth',lw2)
xline(47,'r--','LineWidth',lw2)
xline(56,'r--','LineWidth',lw2)
xline(61,'r--','LineWidth',lw2)
yline(21000,'r.-','LineWidth',lw2)
yline(41000,'r.-','LineWidth',lw2)
yline(100000,'r.-','LineWidth',lw2)
yline(405000,'r.-','LineWidth',lw2)
colormap(parula)
caxis([0 0.005*(max(spec_C13_raw(:)))])
title('\delta^{13}C (normal FFT FILTERED)')
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
ylim([10000 450000])
set(gca,'FontSize',fs)
set(gca,'YScale','log')
set(gca,'XDir','reverse')
set(gca,'YDir','normal')
grid on
subplot(212)
imagesc(t(1:12421),t2,y4)
xline(13.9,'r--','LineWidth',lw2)
xline(17,'r--','LineWidth',lw2)
xline(34,'r--','LineWidth',lw2)
xline(47,'r--','LineWidth',lw2)
xline(56,'r--','LineWidth',lw2)
xline(61,'r--','LineWidth',lw2)
yline(21000,'r.-','LineWidth',lw2)
yline(41000,'r.-','LineWidth',lw2)
yline(100000,'r.-','LineWidth',lw2)
yline(405000,'r.-','LineWidth',lw2)
colormap(parula)
caxis([0 cax_max])
title('\delta^{13}C (no embedding FILTERED)')
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
ylim([10000 450000])
set(gca,'FontSize',fs)
set(gca,'YScale','log')
set(gca,'XDir','reverse')
set(gca,'YDir','normal')
grid on


figure('Units','normalized','Position',[.01 .01 .99 .99])
subplot(211)
imagesc(t(1:12421),periods,flipud(spec_O18_raw))
xline(13.9,'r--','LineWidth',lw2)
xline(17,'r--','LineWidth',lw2)
xline(34,'r--','LineWidth',lw2)
xline(47,'r--','LineWidth',lw2)
xline(56,'r--','LineWidth',lw2)
xline(61,'r--','LineWidth',lw2)
yline(21000,'r.-','LineWidth',lw2)
yline(41000,'r.-','LineWidth',lw2)
yline(100000,'r.-','LineWidth',lw2)
yline(405000,'r.-','LineWidth',lw2)
colormap(parula)
caxis([0 0.01*(max(spec_O18_raw(:)))])
title('\delta^{18}O (normal FFT FILTERED)')
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
ylim([10000 450000])
set(gca,'FontSize',fs)
set(gca,'YScale','log')
set(gca,'XDir','reverse')
set(gca,'YDir','normal')
grid on
subplot(212)
imagesc(t(1:12421),t2,x4)
xline(13.9,'r--','LineWidth',lw2)
xline(17,'r--','LineWidth',lw2)
xline(34,'r--','LineWidth',lw2)
xline(47,'r--','LineWidth',lw2)
xline(56,'r--','LineWidth',lw2)
xline(61,'r--','LineWidth',lw2)
yline(21000,'r.-','LineWidth',lw2)
yline(41000,'r.-','LineWidth',lw2)
yline(100000,'r.-','LineWidth',lw2)
yline(405000,'r.-','LineWidth',lw2)
colormap(parula)
caxis([0 cax_max])
title('\delta^{18}O (no embedding FILTERED)')
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
ylim([10000 450000])
set(gca,'FontSize',fs)
set(gca,'YScale','log')
set(gca,'XDir','reverse')
set(gca,'YDir','normal')
grid on