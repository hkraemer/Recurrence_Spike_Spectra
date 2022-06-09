% For the CENOGRID dataset we perform a windowed analysis with the spike 
% powerspectrum, in order to obtain a wavelet-like spectrum.


% Here: TRAJECTORY BASED ON PECUZAL, MCDTS & NO EMBEDDING (RAW TIME SERIES)


clear, clc

data = load("../data/detrended.txt");
ecc = load("../data/Laskar/ecc_Laskar.txt");

t = data(:,1);
C13 = data(:,2);
O18 = data(:,3);

t = flipud(t);
O18 = flipud(O18);
C13 = flipud(C13);

%% Set parameters
% 
% epsilon = 0.05; % recurrence threshold
% windowsize2 = 1000; % windowsize for trajectory
% windowsize1 = 100; % windowsize for the spectrum
% ws = 1; % windowstep
% N = length(t); % length of the time series
% taus = 0:100; % delay range passed to PECUZAL
% lambda = 0.003;
% 
% M = length(1:ws:N-windowsize2);
% 
% spectrum_O18_pec = zeros(windowsize1,M); 
% spectrum_C13_pec = zeros(windowsize1,M);
% spectrum_O18_mcdts = zeros(windowsize1,M); 
% spectrum_C13_mcdts = zeros(windowsize1,M);
% spectrum_O18_ne = zeros(windowsize1,M); 
% spectrum_C13_ne = zeros(windowsize1,M);
% 
% %% Compute spike powerspectra in windowed analysis
% 
% cnt = 1;
% tic
% for i = 1:ws:M
%     display(i)
%     % extract time series
%     x = O18(i:i+windowsize2);
%     y = C13(i:i+windowsize2);
%     % 1) Pecuzal embedding gained from whole dataset analysis
%     Yx = genembed(x,[0,13,7,10,56,27,3,5,77,42,20,17,15], ones(1,13));
%     Yy = genembed(y,[0,13,7,10,56,27,3,5,77,42,20,17,15], ones(1,13));
%     % compute RP and tau-rr
%     RPx = rp(Yx,epsilon,'var');
%     RPy = rp(Yy,epsilon,'var');
%     tauRRx = tau_recurrence_rate(RPx);
%     tauRRy = tau_recurrence_rate(RPy);
%     % compute spike powerspectra
%     spectrum_O18_pec(:,cnt) = get_spike_spectrum(tauRRx(1:windowsize1), lambda);
%     spectrum_C13_pec(:,cnt) = get_spike_spectrum(tauRRy(1:windowsize1), lambda);
%     
%     % 2) MCDTS embedding gained from whole dataset for prediction
%     Yx = embed(x,3,1);
%     Yy = embed(y,3,1);
%     % compute RP and tau-rr
%     RPx = rp(Yx,epsilon,'var');
%     RPy = rp(Yy,epsilon,'var');
%     tauRRx = tau_recurrence_rate(RPx);
%     tauRRy = tau_recurrence_rate(RPy);
%     % compute spike powerspectra
%     spectrum_O18_mcdts(:,cnt) = get_spike_spectrum(tauRRx(1:windowsize1), lambda);
%     spectrum_C13_mcdts(:,cnt) = get_spike_spectrum(tauRRy(1:windowsize1), lambda);
%     
%     % 2) No embedding
%     Yx = x;
%     Yy = y;
%     % compute RP and tau-rr
%     RPx = rp(Yx,epsilon,'fan');
%     RPy = rp(Yy,epsilon,'fan');
%     tauRRx = tau_recurrence_rate(RPx);
%     tauRRy = tau_recurrence_rate(RPy);
%     % compute spike powerspectra
%     spectrum_O18_ne(:,cnt) = get_spike_spectrum(tauRRx(1:windowsize1), lambda);
%     spectrum_C13_ne(:,cnt) = get_spike_spectrum(tauRRy(1:windowsize1), lambda);
%    
%     cnt = cnt + 1;
%     
% end
% toc
    
% save("spectrum_O18_pec.mat", "spectrum_O18_pec")
% save("spectrum_O18_mcdts.mat", "spectrum_O18_mcdts")
% save("spectrum_O18_ne.mat", "spectrum_O18_ne")
% save("spectrum_C13_pec.mat", "spectrum_C13_pec")
% save("spectrum_C13_mcdts.mat", "spectrum_C13_mcdts")
% save("spectrum_C13_ne.mat", "spectrum_C13_ne")

%% Load data: ISS of isotopes and of eccentricity

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

spectrum_ecc_pec = load("spectrum_ecc_pec2.mat");
spectrum_ecc_pec = spectrum_ecc_pec.spectrum_ecc_pec2;

%% Load data: Fourierspectra of isotopes

dt = abs(mean(diff(t))) * 100000; % Sampling period
Fs = 1/dt; % Sampling frequency
                
spec_raw = load("periodograms_raw.mat");
f = spec_raw.f;
spec_C13_raw = spec_raw.pxx_raw;
spec_O18_raw = spec_raw.pxy_raw;
spec_ecc = spec_raw.p_ecc;

spec_filtered = load("periodograms_filtered.mat");
spec_C13_filtered = spec_filtered.pxx_filtered;
spec_O18_filtered = spec_filtered.pxy_filtered;


%% Plot the wavelet-like spectrograms for both isotope-time series

method = 2; % 1 = PECUZAL, 2 = MCDTS, 3 = no embedding

m_string = ["(PECUZAL)", "(MCDTS)", "(no embedding)"];

fs = 22;
lw2 = 1;

len = 90;
cax_max = 0.013;

if method == 1
    data1 = spectrum_C13_pec;
    data2 = spectrum_O18_pec;
elseif method == 2
    data1 = spectrum_C13_mcdts_mse;
    data2 = spectrum_O18_mcdts_mse; 
elseif method == 3
    data1 = spectrum_C13_ne;
    data2 = spectrum_O18_ne;
end

tt = t(1:size(data1,2));
dt = abs(mean(diff(t))).* 1e6;
ttt = (1:len).*dt;

figure('Units','normalized','Position',[.01 .01 .99 .99])

subplot(4,1,[1,2])
surf(tt,ttt,fliplr(data1(1:len,:)))
view(2)
shading flat
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
title(strcat('\delta^{18}O','{ }',m_string(method)))
xticks([])
ylabel('Period [yrs] ')
ylim([ttt(1) ttt(end)])
xlim([tt(end) tt(1)])
set(gca,'FontSize',fs)
set(gca,'YScale','log')
set(gca,'XDir','reverse')
grid on

subplot(4,1,[3,4])
surf(tt,ttt,fliplr(data2(1:len,:)))
view(2)
shading flat
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
title(strcat('\delta^{13}C','{ }',m_string(method)))
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
ylim([ttt(1) ttt(end)])
xlim([tt(end) tt(1)])
set(gca,'FontSize',fs)
set(gca,'YScale','log')
set(gca,'XDir','reverse')
grid on


%% Plot FFT and spike spectra wavelets together


%% FFT- and ISS for raw data

method = 1; % 1 = PECUZAL, 2 = MCDTS, 3 = no embedding

m_string = ["(PECUZAL)", "(MCDTS)", "(no embedding)"];

len = 90;
cax_max = 0.013;

if method == 1
    data1 = spectrum_C13_pec;
    data2 = spectrum_O18_pec;
elseif method == 2
    data1 = spectrum_C13_mcdts_mse;
    data2 = spectrum_O18_mcdts_mse; 
elseif method == 3
    data1 = spectrum_C13_ne;
    data2 = spectrum_O18_ne;
end

tt = t(1:size(data1,2));
dt = abs(mean(diff(t))).* 1e6;
ttt = (1:len).*dt;

fs = 26;
lw1 = 1;
lw2 = 2.5;

factor = 0.01;

xlines = [13.9 17 34 47 56 61];
ylines1 = [1/21000 1/41000 1/100000 1/405000];
ylines2 = [21000 41000 100000 405000];
ylimms1 = [1/500000 6/100000];
ylimms2 = [10000 450000];

ttf = t(1:size(spec_C13_raw,2));
tttf = 1./f(2:end);

figure('Units','normalized','Position',[.01 .01 .99 .99])

subplot(4,1,[1,2])
surf(ttf,tttf,spec_C13_raw(2:end,:))
view(2)
shading flat
for i = 1:length(xlines)
    xline(xlines(i),'r--','LineWidth',lw2)
end
for i = 1:length(ylines1)
    yline(ylines2(i),'k--','LineWidth',lw1)
end
colormap(parula)
caxis([0 factor*(max(spec_C13_raw(:)))])
title('\delta^{13}C evolutionary FFT spectrogram')
ylabel('Period [yrs] ')
ylim([10000 ttt(end)])
xlim([ttf(end) ttf(1)])
xticklabels([])
set(gca,'FontSize',fs,'LineWidth',2,'YScale','log','XDir','reverse')
grid on

subplot(4,1,[3,4])
surf(tt,ttt,fliplr(data2(1:len,:)))
view(2)
shading flat
for i = 1:length(xlines)
    xline(xlines(i),'r--','LineWidth',lw2)
end
for i = 1:length(ylines2)
    yline(ylines2(i),'k--','LineWidth',lw1)
end
colormap(parula)
caxis([0 cax_max])
title(strcat('\delta^{13}C evolutionary ISS of \tau-RR','{ }',m_string(method)))
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
ylim([10000 ttt(end)])
xlim([tt(end) tt(1)])
set(gca,'FontSize',fs,'LineWidth',2,'YScale','log','XDir','reverse')
grid on


figure('Units','normalized','Position',[.01 .01 .99 .99])

subplot(4,1,[1,2])
surf(ttf,tttf,spec_O18_raw(2:end,:))
view(2)
shading flat
for i = 1:length(xlines)
    xline(xlines(i),'r--','LineWidth',lw2)
end
for i = 1:length(ylines1)
    yline(ylines2(i),'k--','LineWidth',lw1)
end
colormap(parula)
caxis([0 factor*(max(spec_O18_raw(:)))])
title('\delta^{18}O evolutionary FFT spectrogram')
ylabel('Period [yrs] ')
ylim([10000 ttt(end)])
xlim([ttf(end) ttf(1)])
xticklabels([])
set(gca,'FontSize',fs,'LineWidth',2,'YScale','log','XDir','reverse')
grid on

subplot(4,1,[3,4])
surf(tt,ttt,fliplr(data1(1:len,:)))
view(2)
shading flat
for i = 1:length(xlines)
    xline(xlines(i),'r--','LineWidth',lw2)
end
for i = 1:length(ylines2)
    yline(ylines2(i),'k--','LineWidth',lw1)
end
colormap(parula)
caxis([0 cax_max])
title(strcat('\delta^{18}O evolutionary ISS of \tau-RR','{ }',m_string(method)))
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
ylim([10000 ttt(end)])
xlim([tt(end) tt(1)])
set(gca,'FontSize',fs,'LineWidth',2,'YScale','log','XDir','reverse')
grid on

%% Eccentricity & ISS

figure('Units','normalized','Position',[.01 .01 .99 .99])

subplot(4,1,[1,2])
surf(tt,ttt,spectrum_ecc_pec(1:len,:))
view(2)
shading flat
for i = 1:length(xlines)
    xline(xlines(i),'r--','LineWidth',lw2)
end
for i = 1:length(ylines1)
    yline(ylines2(i),'k--','LineWidth',lw1)
end
colormap(parula)
caxis([0 0.05])
title(strcat('Evolutionary inter spike spectrum of Eccentricity','{ }',m_string(method)))
ylabel('Period [yrs] ')
ylim([10000 ttt(end)])
xlim([tt(end) tt(1)])
xticks([])
set(gca,'FontSize',fs,'LineWidth',2,'YScale','log','XDir','reverse')
grid on

subplot(4,1,[3,4])
surf(tt,ttt,fliplr(data2(1:len,:)))
view(2)
shading flat
for i = 1:length(xlines)
    xline(xlines(i),'r--','LineWidth',lw2)
end
for i = 1:length(ylines2)
    yline(ylines2(i),'k--','LineWidth',lw1)
end
colormap(parula)
caxis([0 cax_max])
title(strcat('\delta^{13}C evolutionary ISS of \tau-RR','{ }',m_string(method)))
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
ylim([10000 ttt(end)])
xlim([tt(end) tt(1)])
set(gca,'FontSize',fs,'LineWidth',2,'YScale','log','XDir','reverse')
grid on


figure('Units','normalized','Position',[.01 .01 .99 .99])

subplot(4,1,[1,2])
surf(tt,ttt,spectrum_ecc_pec(1:len,:))
view(2)
shading flat
for i = 1:length(xlines)
    xline(xlines(i),'r--','LineWidth',lw2)
end
for i = 1:length(ylines1)
    yline(ylines2(i),'k--','LineWidth',lw1)
end
colormap(parula)
caxis([0 0.05])
title(strcat('Evolutionary inter spike spectrum of Eccentricity','{ }',m_string(method)))
ylabel('Frequency [1/yrs] ')
ylim([10000 ttt(end)])
xlim([tt(end) tt(1)])
xticks([])
set(gca,'FontSize',fs,'LineWidth',2,'YScale','log','XDir','reverse')
grid on

subplot(4,1,[3,4])
surf(tt,ttt,fliplr(data1(1:len,:)))
view(2)
shading flat
for i = 1:length(xlines)
    xline(xlines(i),'r--','LineWidth',lw2)
end
for i = 1:length(ylines2)
    yline(ylines2(i),'k--','LineWidth',lw1)
end
colormap(parula)
caxis([0 cax_max])
title(strcat('\delta^{18}O evolutionary ISS of \tau-RR','{ }',m_string(method)))
xlabel('time [Mio yrs BP]')
ylabel('Period [yrs] ')
ylim([10000 ttt(end)])
xlim([tt(end) tt(1)])
set(gca,'FontSize',fs,'LineWidth',2,'YScale','log','XDir','reverse')
grid on