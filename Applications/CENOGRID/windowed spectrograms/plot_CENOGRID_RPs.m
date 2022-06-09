% For the CENOGRID dataset we perform a windowed analysis with the spike 
% powerspectrum, in order to obtain a wavelet-like spectrum.


% Here: Plot RPs (RAW TIME SERIES)


clear, clc

data = load("../data/Loess20.txt");

t = data(:,1);
O18 = data(:,2);
C13 = data(:,3);

%% first interpolate the data
% t2 = 0:0.005:67.1;
% data_interp1 = interp1(t, O18, t2, 'spline');
% data_interp2 = interp1(t, C13, t2, 'spline');
% 
% clear data
% data(:,1) = t2;
% data(:,2) = data_interp1;
% data(:,3) = data_interp2;
% 
% save("../data/Loess_20_int.txt","data","-ascii")

%%
clear, clc

data = load("../data/Loess_20_int.txt");

t = data(:,1);
C13 = data(:,2);
O18 = data(:,3);

t = flipud(t);
O18 = flipud(O18);
C13 = flipud(C13);

%%

[RP1, D1] = rp(C13,0.08,'fan');
[RP2, D2] = rp(O18,0.08,'fan');

%%
RP1 = D1<0.08;
RP2 = D2<0.08;

fs = 22; % fontsize

figure('Units','normalized','Position',[.01 .01 .99 .99])

subplot(121)
imagesc(t,t,RP1), colormap([1 1 1;0 0 0]), axis xy square
set(gca,'YDir','reverse')
set(gca,'XDir','reverse')
set(gca,'LineWidth',2)
set(gca,'FontSize',fs)
xlabel("time [Mio yrs BP]")
xticks(0:5:65)
yticks(0:5:65)
grid on

subplot(122)
imagesc(t,t,RP2), colormap([1 1 1;0 0 0]), axis xy square
set(gca,'YDir','reverse')
set(gca,'XDir','reverse')
set(gca,'YAxisLocation','right')
set(gca,'LineWidth',2)
set(gca,'FontSize',fs)
xlabel("time [Mio yrs BP]")
xticks(0:5:65)
yticks(0:5:65)
grid on

%%

figure('Units','normalized','Position',[.01 .01 .99 .3])

subplot(121)
plot(t,C13,'LineWidth',2)
set(gca,'XDir','reverse')
set(gca,'LineWidth',2)
set(gca,'FontSize',fs)
ylabel("\delta^{13}C [‰]")
title("benthic \delta^{13}C stable isotope time series")
xticks(0:5:65)
xlim([0 67.1])
grid on

subplot(122)
plot(t,O18,'LineWidth',2)
set(gca,'YDir','reverse')
set(gca,'XDir','reverse')
set(gca,'LineWidth',2)
set(gca,'FontSize',fs)
ylabel("\delta^{18}O [‰]")
title("benthic \delta^{18}O stable isotope time series")
xticks(0:5:65)
xlim([0 67.1])
ylim([-3 6])
grid on


