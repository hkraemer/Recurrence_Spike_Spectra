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


%% Load data: ISS of isotopes and of eccentricity

spectrum_O18_mcdts_mse = load("./computed spectra/spectrum_O18_mcdts.mat");
spectrum_O18_mcdts_mse = spectrum_O18_mcdts_mse.spectrum_O18_mcdts;
spectrum_O18_cao = load("./computed spectra/spectrum_O18_cao.mat");
spectrum_O18_cao = spectrum_O18_cao.spectrum_O18_cao;
spectrum_O18_ne = load("./computed spectra/spectrum_O18_ne.mat");
spectrum_O18_ne = spectrum_O18_ne.spectrum_O18_ne;
spectrum_C13_mcdts_mse = load("./computed spectra/spectrum_C13_mcdts.mat");
spectrum_C13_mcdts_mse = spectrum_C13_mcdts_mse.spectrum_C13_mcdts;
spectrum_C13_cao = load("./computed spectra/spectrum_C13_cao.mat");
spectrum_C13_cao = spectrum_C13_cao.spectrum_C13_mcdts;
spectrum_C13_ne = load("./computed spectra/spectrum_C13_ne.mat");
spectrum_C13_ne = spectrum_C13_ne.spectrum_C13_ne;

spectrum_ecc_mcdts = load("./computed spectra/spectrum_ecc_mcdts.mat");
spectrum_ecc_mcdts = spectrum_ecc_mcdts.spectrum_ecc_mcdts;
spectrum_ecc_ne = load("./computed spectra/spectrum_ecc_ne.mat");
spectrum_ecc_ne = spectrum_ecc_ne.spectrum_ecc_ne;
spectrum_ecc_fft = load("./computed spectra/spectrum_ecc_fft.mat");
spectrum_ecc_fft = spectrum_ecc_fft.spectrum_ecc_fft;


%% Load data: Fourierspectra of isotopes and eccentricity
% computed in script `compute_windowed_fft_periodograms.m`

dt = abs(mean(diff(t))) * 100000; % Sampling period
% dt = abs(mean(diff(t))); % Sampling period
Fs = 1/dt; % Sampling frequency
                
spec_filtered = load("./computed spectra/periodograms_short_time_fourier.mat");
spec_C13_filtered = spec_filtered.pxx_filtered;
spec_C13_filtered = spec_filtered.pxx_raw;
spec_O18_filtered = spec_filtered.pxy_filtered;
spec_O18_filtered = spec_filtered.pxy_raw;
spec_ecc = spec_filtered.p_ecc;
f = spec_filtered.f;


%% Plot the wavelet-like spectrograms for both isotope-time series

method = 2; % 1 = MCDTS, 2 = CAO, 3 = no embedding

m_string = ["(MCDTS)", "(CAO)", "(no embedding)"];

fs = 22;
lw2 = 1;

len = 90;
cax_max = 0.013;
cax_max = 0.9;

if method == 1
    data1 = spectrum_C13_mcdts_mse;
    data2 = spectrum_O18_mcdts_mse; 
elseif method == 2
    data1 = spectrum_C13_cao;
    data2 = spectrum_O18_cao;
elseif method == 3
    data1 = spectrum_C13_ne;
    data2 = spectrum_O18_ne;
end

tt = t(1:size(data1,2));
dt = abs(mean(diff(t))).* 1e6;
ttt = (1:len).*dt;
tf = (1:len).*Fs.*0.001;
ttf = 1./tf;

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

method = 1; % 1 = MCDTS, 2 = CAO, 3 = no embedding

m_string = ["(MCDTS)", "(CAO)", "(no embedding)"];

len = 90;
cax_max = 0.013;

tt = t(1:size(spectrum_ecc_mcdts,2)); % for ISS
tt2 = t(1:size(spec_ecc,2)); % for short time fourier
dt = abs(mean(diff(t))).* 1e6;
ttt = (1:len).*dt;
tttf = 1./f(2:end);

xlines = [13.9 17 34 47 56 61];
ylines1 = [1/21000 1/41000 1/100000 1/405000];
ylines2 = [21000 41000 100000 405000];
ylimms1 = [1/500000 6/100000];
ylimms2 = [10000 450000];


fs = 26;
lw1 = 1;
lw2 = 2.5;

factor = 0.01;

%%
if method == 1
    data1 = spectrum_C13_mcdts_mse;
    data2 = spectrum_O18_mcdts_mse; 
elseif method == 2
    data1 = spectrum_C13_cao;
    data2 = spectrum_O18_cao;
elseif method == 3
    data1 = spectrum_C13_ne;
    data2 = spectrum_O18_ne;
end


ttf = t(1:size(spec_C13_filtered,2));


figure('Units','normalized','Position',[.01 .01 .99 .99])

subplot(4,1,[1,2])
surf(ttf,tttf,spec_C13_filtered(2:end,:))
view(2)
shading flat
for i = 1:length(xlines)
    xline(xlines(i),'r--','LineWidth',lw2)
end
for i = 1:length(ylines1)
    yline(ylines2(i),'k--','LineWidth',lw1)
end
colormap(parula)
caxis([0 factor*(max(spec_C13_filtered(:)))])
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
surf(ttf,tttf,spec_O18_filtered(2:end,:))
view(2)
shading flat
for i = 1:length(xlines)
    xline(xlines(i),'r--','LineWidth',lw2)
end
for i = 1:length(ylines1)
    yline(ylines2(i),'k--','LineWidth',lw1)
end
colormap(parula)
caxis([0 factor*(max(spec_O18_filtered(:)))])
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
surf(tt,ttt,spectrum_ecc_mcdts(1:len,:))
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

%% Eccentricity

figure('Units','normalized','Position',[.01 .01 .99 .99])

h1 = subplot(6,1,[1,2]);
surf(tt2(1:length(tt)),tttf,spec_ecc(2:end,1:length(tt)))
view(2)
shading flat
% for i = 1:length(xlines)
%     xline(xlines(i),'r--','LineWidth',lw2)
% end
for i = 1:length(ylines1)
    yline(ylines2(i),'k--','LineWidth',lw1)
end
colormap(parula)
caxis([0 0.05])
title('Evolutionary FFT of Eccentricity (short time FT Hamming window)')
ylabel('Period [yrs] ')
ylim([10000 ttt(end)])
xlim([tt(end) tt(1)])
xticks([])
set(gca,'FontSize',fs,'LineWidth',2,'YScale','log','XDir','reverse')
h1_pos = get(h1,'Position');
grid on

h2 = subplot(6,1,[3,4]);
surf(tt,ttt,spectrum_ecc_mcdts(1:len,:))
view(2)
shading flat
% for i = 1:length(xlines)
%     xline(xlines(i),'r--','LineWidth',lw2)
% end
for i = 1:length(ylines1)
    yline(ylines2(i),'k--','LineWidth',lw1)
end
colormap(parula)
caxis([0 0.05])
title({" ";'Evolutionary inter spike spectrum of \tau-RR of embedded Eccentricity'})
ylabel('Period [yrs] ')
ylim([10000 ttt(end)])
xlim([tt(end) tt(1)])
xticks([])
set(gca,'FontSize',fs,'LineWidth',2,'YScale','log','XDir','reverse')
h2_pos=get(h2,'Position'); 
set(h2,'Position',[h1_pos(1) h1_pos(2)-(1.2*h1_pos(4)) h2_pos(3:end)])
grid on

h3 = subplot(6,1,[5,6]);
surf(tt,ttf,spectrum_ecc_fft(1:len,:))
view(2)
shading flat
% for i = 1:length(xlines)
%     xline(xlines(i),'r--','LineWidth',lw2)
% end
for i = 1:length(ylines1)
    yline(ylines2(i),'k--','LineWidth',lw1)
end
colormap(parula)
% caxis([0 0.05])
title(strcat('Evolutionary FFT spectrum of \tau-RR of embedded Eccentricity'))
ylabel('Period [yrs] ')
xlabel('time [Mio yrs BP]')
ylim([10000 ttt(end)])
xlim([tt(end) tt(1)])
set(gca,'FontSize',fs,'LineWidth',2,'YScale','log','XDir','reverse')
h3_pos=get(h2,'Position'); 
set(h3,'Position',[h1_pos(1) h2_pos(2)-(1.25*h1_pos(4)) h2_pos(3:end)])
grid on



