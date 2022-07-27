% For the LASKAR dataset (eccentricity time series from Laskar 2011) we
% plot the results for the evolutionary spectra (computed in
% `compute_windowed_spectrograms_embedded_eccentricity.m` and in script 
% `compute_windowed_fft_periodograms.m`)

clear, clc

ecc = load("../data/Laskar/ecc_Laskar.txt");

dt = 5000; % Sampling period/time
Fs = 1/dt; % Sampling frequency

t = linspace(0,6.7100e+07,length(ecc)); % time vector

%% Load computed spectra

% embedded time series

% inter spike spectra
spectrum_ecc_mcdts = load("./computed spectra/spectrum_ecc_mcdts.mat");
spectrum_ecc_mcdts = spectrum_ecc_mcdts.spectrum_ecc_mcdts;
% fourier spectra
spectrum_ecc_fft = load("./computed spectra/spectrum_ecc_fft.mat");
spectrum_ecc_fft = spectrum_ecc_fft.spectrum_ecc_fft;

% not embedded time series

% inter spike spectra
spectrum_ecc_ne = load("./computed spectra/spectrum_ecc_ne.mat");
spectrum_ecc_ne = spectrum_ecc_ne.spectrum_ecc_ne;
% fourier spectra
spectrum_ecc_fft_ne = load("./computed spectra/spectrum_ecc_fft_ne.mat");
spectrum_ecc_fft_ne = spectrum_ecc_fft_ne.spectrum_ecc_fft_ne;

% FFT based on raw time series (without any RP - tau-RR computations)
spec_filtered = load("./computed spectra/periodograms_short_time_fourier.mat");
spec_ecc = spec_filtered.p_ecc;
f = spec_filtered.f;


%% Plot the wavelet-like spectrograms 

method = 1; % 1 = MCDTS, 2 = no embedding

m_string = ["(MCDTS)", "(no embedding)"];

fs = 22;
lw1 = 1;
lw2 = 2.5;

len = 90;
cax_max = 0.013;

if method == 1
    data1 = spectrum_ecc_mcdts;
    data2 = spectrum_ecc_fft; 
elseif method == 2
    data1 = spectrum_ecc_ne;
    data2 = spectrum_ecc_fft_ne;
end

% prepare time vector for plotting
tt = t(1:length(data1))./1000000; % for plotting in Mio yrs

ttt = t(1:len);
window_tau = 200; % length of the input signal to the FFT, see script `compute_windowed_spectrograms_embedded_eccentricity.m`
ff = Fs*(0:(window_tau/2))/window_tau;
tf = 1./ff(2:end); % periods for windowed tau-RR
ttf = 1./f(2:end); % periods for short FT of raw time series

%% normalization:

ylim1 = 10000;
ylim2 = ttt(end); % corresponds to 445000

% indices for FFT of raw time series time frame corresponding to ylim1 & ylim2
idx1 = find(ttf==ylim1);
idx2 = find(ttf>ylim2,1,'last');

% indices for FFT of tau-RR time frame corresponding to ylim1 & ylim2
idx1_ = find(tf==ylim1);
idx2_ = find(tf>ylim2,1,'last');

M1 = length(idx2:idx1); % length of y-vals for ecc-fft
M2 = length(1:len); % length of y-values for ISS-spectograms
M3 = length(idx2_:idx1_); % length of y-vals for FFT of tauRR
M4 = length(1:length(tt)); % length of x-vals for all spectra

spec1 = zeros(M1, M4); % ecc FFT-spectro
spec2 = zeros(M2, M4); % ISS of tau-RR of ecc
spec3 = zeros(M3, M4); % FFT of tau-RR of ecc


for i = 1:M4
    spec1(:,i) = spec_ecc(idx2:idx1,i) ./ sum(spec_ecc(idx2:idx1,i));
    spec2(:,i) = data1(1:len,i) ./ sum(data1(1:len,i));
    spec3(:,i) = data2(idx2_:idx1_,i) ./ sum(data2(idx2_:idx1_,i));
end

% for i = 1:length(spec_ecc)
%     spec_ecc(:,i) = spec_ecc(:,i) ./ sum(spec_ecc(:,i));
% end
% for i = 1:length(spectrum_ecc_ne)
%     spectrum_ecc_ne(:,i) = spectrum_ecc_ne(:,i) ./ sum(spectrum_ecc_ne(:,i));
%     spectrum_ecc_mcdts(:,i) = spectrum_ecc_mcdts(:,i) ./ sum(spectrum_ecc_mcdts(:,i));
%     spectrum_ecc_fft(:,i) = spectrum_ecc_fft(:,i) ./ sum(spectrum_ecc_fft(:,i));
%     spectrum_ecc_fft_ne(:,i) = spectrum_ecc_fft_ne(:,i) ./ sum(spectrum_ecc_fft_ne(:,i));
% end



%%
ylines2 = [95000 124000 405000];
ylimms2 = [10000 450000];

c_ax_upper = 0.05;

figure('Units','normalized','Position',[.001 .001 .99 .99])

% h1 = subplot(6,1,[1,2]);
tiledlayout(3,1)
nexttile
surf(tt,ttf,spec_ecc(2:end,1:length(tt)))
surf(tt,ttf(idx2:idx1),spec1)
view(2)
shading flat
for i = 1:length(ylines2)
    yline(ylines2(i),'k--','LineWidth',lw1)
end
colormap(parula)
caxis([0 c_ax_upper])
title('Evolutionary FFT of Eccentricity (short time FT Hamming window)')
ylabel('Period [yrs] ')
ylim([10000 ttt(end)])
xlim([tt(1) tt(end)])
xticks([])
set(gca,'FontSize',fs,'LineWidth',2,'YScale','log','XDir','reverse')
% h1_pos = get(h1,'Position');
grid on


% h2 = subplot(6,1,[3,4]);
nexttile
% surf(tt,ttt,data1(1:len,:))
surf(tt,ttt,spec2)
view(2)
shading flat
for i = 1:length(ylines2)
    yline(ylines2(i),'k--','LineWidth',lw1)
end
colormap(parula)
caxis([0 c_ax_upper])
title({" ";'Evolutionary inter spike spectrum of \tau-RR of embedded Eccentricity'})
ylabel('Period [yrs] ')
ylim([10000 ttt(end)])
xlim([tt(1) tt(end)])
xticks([])
set(gca,'FontSize',fs,'LineWidth',2,'YScale','log','XDir','reverse')
% h2_pos=get(h2,'Position'); 
% set(h2,'Position',[h1_pos(1) h1_pos(2)-(1.2*h1_pos(4)) h2_pos(3:end)])
grid on

nexttile
% h3 = subplot(6,1,[5,6]);
% surf(tt,tf,data2)
surf(tt,tf(idx2_:idx1_),spec3)
view(2)
shading flat
for i = 1:length(ylines2)
    yline(ylines2(i),'k--','LineWidth',lw1)
end
colormap(parula)
caxis([0 c_ax_upper])
title(strcat('Evolutionary FFT spectrum of \tau-RR of embedded Eccentricity'))
ylabel('Period [yrs] ')
xlabel('time [Mio yrs BP]')
ylim([10000 ttt(end)])
xlim([tt(1) tt(end)])
set(gca,'FontSize',fs,'LineWidth',2,'YScale','log','XDir','reverse')
% h3_pos=get(h2,'Position'); 
% set(h3,'Position',[h1_pos(1) h2_pos(2)-(1.25*h1_pos(4)) h2_pos(3:end)])
grid on

cb = colorbar;
cb.Layout.Tile = 'south';
