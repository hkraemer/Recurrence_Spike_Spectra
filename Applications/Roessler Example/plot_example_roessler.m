clear, clc
cd("./computed data")

tau_rr1 = load('tau_rr_1.csv');
tau_rr2 = load('tau_rr_2.csv');
tau_rr3 = load('tau_rr_3.csv');
spectrum1(1,:)= load('tauRR_spectrum1_roessler_0.9.csv');
spectrum1(2,:) = load('tauRR_spectrum1_roessler_0.95.csv');
spectrum1(3,:)= load('tauRR_spectrum1_roessler_0.99.csv');
spectrum2(1,:)= load('tauRR_spectrum2_roessler_0.9.csv');
spectrum2(2,:) = load('tauRR_spectrum2_roessler_0.95.csv');
spectrum2(3,:)= load('tauRR_spectrum2_roessler_0.99.csv');
spectrum3(1,:)= load('tauRR_spectrum3_roessler_0.9.csv');
spectrum3(2,:) = load('tauRR_spectrum3_roessler_0.95.csv');
spectrum3(3,:)= load('tauRR_spectrum3_roessler_0.99.csv');
RP_1 = load('RP_1.csv');
RP_2 = load('RP_2.csv');
RP_3 = load('RP_3.csv');
Y1 = load('Y_1.csv');
Y2 = load('Y_2.csv');
Y3 = load('Y_3.csv');

tau_rr1_n = load('tau_rr_1_n.csv');
tau_rr2_n = load('tau_rr_2_n.csv');
tau_rr3_n = load('tau_rr_3_n.csv');
spectrum1_n(1,:)= load('tauRR_spectrum1_n_roessler_0.9.csv');
spectrum1_n(2,:) = load('tauRR_spectrum1_n_roessler_0.95.csv');
spectrum1_n(3,:)= load('tauRR_spectrum1_n_roessler_0.99.csv');
spectrum2_n(1,:)= load('tauRR_spectrum2_n_roessler_0.9.csv');
spectrum2_n(2,:) = load('tauRR_spectrum2_n_roessler_0.95.csv');
spectrum2_n(3,:)= load('tauRR_spectrum2_n_roessler_0.99.csv');
spectrum3_n(1,:)= load('tauRR_spectrum3_n_roessler_0.9.csv');
spectrum3_n(2,:) = load('tauRR_spectrum3_n_roessler_0.95.csv');
spectrum3_n(3,:)= load('tauRR_spectrum3_n_roessler_0.99.csv');
RP_1_n = load('RP_1_n.csv');
RP_2_n = load('RP_2_n.csv');
RP_3_n = load('RP_3_n.csv');
Y1_n = load('Y_1_n.csv');
Y2_n = load('Y_2_n.csv');
Y3_n = load('Y_3_n.csv');

%% Compute FFT-periodograms of the tau-RRs 

freq = 1;
[perio1, fs1] = periodogram(tau_rr1,hamming(length(tau_rr1)),length(tau_rr1),freq);
[perio2, fs2] = periodogram(tau_rr2,hamming(length(tau_rr2)),length(tau_rr2),freq);
[perio3, fs3] = periodogram(tau_rr3,hamming(length(tau_rr3)),length(tau_rr3),freq);

[perio1_n, fs1_n] = periodogram(tau_rr1_n,hamming(length(tau_rr1)),length(tau_rr1),freq);
[perio2_n, fs2_n] = periodogram(tau_rr2_n,hamming(length(tau_rr2)),length(tau_rr2),freq);
[perio3_n, fs3_n] = periodogram(tau_rr3_n,hamming(length(tau_rr3)),length(tau_rr3),freq);


%%
lw = 3;
lw2 = 2;
fs = 14;

figure('Units','normalized','Position',[.01 .01 .6 .3])
subplot(131)
plot3(Y1(:,1),Y1(:,2),Y1(:,3))
grid on
xticklabels([])
yticklabels([])
zticklabels([])
set(gca,'LineWidth',2)
set(gca,'Fontsize',fs)
title('Trajectory of Rössler system (a = 0.36)')

subplot(132)
plot3(Y2(:,1),Y2(:,2),Y2(:,3))
grid on
xticklabels([])
yticklabels([])
zticklabels([])
set(gca,'LineWidth',2)
set(gca,'Fontsize',fs)
title('Trajectory of Rössler system (a = 0.41)')

subplot(133)
plot3(Y3(:,1),Y3(:,2),Y3(:,3))
grid on
xticklabels([])
yticklabels([])
zticklabels([])
set(gca,'LineWidth',2)
set(gca,'Fontsize',fs)
title('Trajectory of Rössler system (a = 0.428)')

%% Pick the threshold val and plot spectros

i = 3; % (1: 0.9, 2: 0.95, 3: 0.99)

figure('Units','normalized','Position',[.01 .01 .6 .99])

subplot(331)
imagesc(RP_1), colormap([0 0 0 ; 1 1 1])
%set(gca, 'YDir', 'normal')
set(gca,'LineWidth',2)
set(gca,'Fontsize',fs)
title('RP')
xlabel('time [a.u.]')
ylabel('time [a.u.]')

subplot(332)
imagesc(RP_2), colormap([0 0 0 ; 1 1 1])
%set(gca, 'YDir', 'normal')
set(gca,'LineWidth',2)
set(gca,'Fontsize',fs)
title('RP')
xlabel('time [a.u.]')
ylabel('time [a.u.]')

subplot(333)
imagesc(RP_3), colormap([0 0 0 ; 1 1 1])
%set(gca, 'YDir', 'normal')
set(gca,'LineWidth',2)
set(gca,'Fontsize',fs)
title('RP')
xlabel('time [a.u.]')
ylabel('time [a.u.]')

subplot(334)
plot((1:length(tau_rr1)),tau_rr1, 'Linewidth', lw)
grid on
title('\tau-RR')
set(gca,'LineWidth',2)
set(gca,'Fontsize',fs)
ylabel('\tau-RR')
xlabel('period [a.u.]')
xlim([0 1000])

subplot(335)
plot((1:length(tau_rr1)),tau_rr2, 'Linewidth', lw)
grid on
title('\tau-RR')
set(gca,'LineWidth',2)
set(gca,'Fontsize',fs)
ylabel('\tau-RR')
xlabel('period [a.u.]')
xlim([0 1000])

subplot(336)
plot((1:length(tau_rr1)),tau_rr3, 'Linewidth', lw)
grid on
title('\tau-RR')
set(gca,'LineWidth',2)
set(gca,'Fontsize',fs)
ylabel('\tau-RR')
xlabel('period [a.u.]')
xlim([0 1000])

subplot(337)
plot((1:length(spectrum1)),spectrum1(i,:), '.-', 'Linewidth', lw)
grid on
title('inter spike spectrum of \tau-RR')
set(gca,'LineWidth',2)
set(gca,'Fontsize',fs)
ylabel('rel. power/period')
xlabel('period [a.u.]')

subplot(338)
plot((1:length(spectrum1)),spectrum2(i,:), '.-', 'Linewidth', lw)
grid on
title('inter spike spectrum of \tau-RR')
set(gca,'LineWidth',2)
set(gca,'Fontsize',fs)
ylabel('rel. power/period')
xlabel('period [a.u.]')

subplot(339)
plot((1:length(spectrum1)),spectrum3(i,:), '.-', 'Linewidth', lw)
grid on
title('inter spike spectrum of \tau-RR')
set(gca,'LineWidth',2)
set(gca,'Fontsize',fs)
ylabel('rel. power/period')
xlabel('period [a.u.]')


figure('Units','normalized','Position',[.01 .01 .6 .99])

% subplot(331)
% imagesc(RP_1), colormap([1 1 1 ; 0 0 0])
% set(gca, 'YDir', 'normal')
% set(gca,'LineWidth',2)
% set(gca,'Fontsize',fs)
% title('RP')
% xlabel('time [a.u.]')
% ylabel('time [a.u.]')
% 
% subplot(332)
% imagesc(RP_2), colormap([1 1 1 ; 0 0 0])
% set(gca, 'YDir', 'normal')
% set(gca,'LineWidth',2)
% set(gca,'Fontsize',fs)
% title('RP')
% xlabel('time [a.u.]')
% ylabel('time [a.u.]')
% 
% subplot(333)
% imagesc(RP_3), colormap([1 1 1 ; 0 0 0])
% set(gca, 'YDir', 'normal')
% set(gca,'LineWidth',2)
% set(gca,'Fontsize',fs)
% title('RP')
% xlabel('time [a.u.]')
% ylabel('time [a.u.]')
% 
% subplot(334)
% plot((1:length(spectrum1)),tau_rr1, 'Linewidth', lw)
% grid on
% title('\tau-RR')
% set(gca,'LineWidth',2)
% set(gca,'Fontsize',fs)
% ylabel('\tau-RR')
% xlabel('period [a.u.]')
% xlim([0 1000])
% 
% subplot(335)
% plot((1:length(spectrum2)),tau_rr2, 'Linewidth', lw)
% grid on
% title('\tau-RR')
% set(gca,'LineWidth',2)
% set(gca,'Fontsize',fs)
% ylabel('\tau-RR')
% xlabel('period [a.u.]')
% xlim([0 1000])
% 
% subplot(336)
% plot((1:length(spectrum3)),tau_rr3, 'Linewidth', lw)
% grid on
% title('\tau-RR')
% set(gca,'LineWidth',2)
% set(gca,'Fontsize',fs)
% ylabel('\tau-RR')
% xlabel('period [a.u.]')
% xlim([0 1000])

subplot(337)
plot(1./fs1,10*log10(perio1),'.-', 'Linewidth', lw2)
grid on
title('Fourier powerspectrum of \tau-RR')
set(gca,'LineWidth',2)
set(gca,'Fontsize',fs)
ylabel('rel. power/period')
xlabel('period [a.u.]')
xlim([0 500])

subplot(338)
plot(1./fs1,10*log10(perio2), '.-', 'Linewidth', lw2)
grid on
title('Fourier powerspectrum of \tau-RR')
set(gca,'LineWidth',2)
set(gca,'Fontsize',fs)
ylabel('rel. power/period')
xlabel('period [a.u.]')
xlim([0 500])

subplot(339)
plot(1./fs1,10*log10(perio1), '.-', 'Linewidth', lw2)
grid on
title('Fourier powerspectrum of \tau-RR')
set(gca,'LineWidth',2)
set(gca,'Fontsize',fs)
ylabel('rel. power/period')
xlabel('period [a.u.]')
xlim([0 500])

%%

%% Noise 1

figure('Units','normalized','Position',[.01 .01 .6 .3])
subplot(131)
plot3(Y1_n(:,1),Y1_n(:,2),Y1_n(:,3))
grid on
xticklabels([])
yticklabels([])
zticklabels([])
set(gca,'LineWidth',2)
set(gca,'Fontsize',fs)
title('Trajectory of Rössler system 5% add. noise (a = 0.36)')

subplot(132)
plot3(Y2_n(:,1),Y2_n(:,2),Y2_n(:,3))
grid on
xticklabels([])
yticklabels([])
zticklabels([])
set(gca,'LineWidth',2)
set(gca,'Fontsize',fs)
title('Trajectory of Rössler system 5% add. noise (a = 0.41)')

subplot(133)
plot3(Y3_n(:,1),Y3_n(:,2),Y3_n(:,3))
grid on
xticklabels([])
yticklabels([])
zticklabels([])
set(gca,'LineWidth',2)
set(gca,'Fontsize',fs)
title('Trajectory of Rössler system 5% add. noise (a = 0.428)')


%% Noise 2

figure('Units','normalized','Position',[.01 .01 .6 .99])

subplot(331)
imagesc(RP_1_n), colormap([0 0 0 ; 1 1 1])
%set(gca,'Ydir','normal')
set(gca,'LineWidth',2)
set(gca,'Fontsize',fs)
title('RP')
xlabel('time [a.u.]')
ylabel('time [a.u.]')

subplot(332)
imagesc(RP_2_n), colormap([0 0 0 ; 1 1 1])
%set(gca,'Ydir','normal')
set(gca,'LineWidth',2)
set(gca,'Fontsize',fs)
title('RP')
xlabel('time [a.u.]')
ylabel('time [a.u.]')

subplot(333)
imagesc(RP_3_n), colormap([0 0 0 ; 1 1 1])
%set(gca, 'YDir', 'normal')
set(gca,'LineWidth',2)
set(gca,'Fontsize',fs)
title('RP')
xlabel('time [a.u.]')
ylabel('time [a.u.]')

subplot(334)
plot((1:length(tau_rr1_n)),tau_rr1_n, 'Linewidth', lw)
grid on
title('\tau-RR')
set(gca,'LineWidth',2)
set(gca,'Fontsize',fs)
ylabel('\tau-RR')
xlabel('period [a.u.]')
xlim([0 1000])

subplot(335)
plot((1:length(tau_rr1_n)),tau_rr2_n, 'Linewidth', lw)
grid on
title('\tau-RR')
set(gca,'LineWidth',2)
set(gca,'Fontsize',fs)
ylabel('\tau-RR')
xlabel('period [a.u.]')
xlim([0 1000])

subplot(336)
plot((1:length(tau_rr1_n)),tau_rr3_n, 'Linewidth', lw)
grid on
title('\tau-RR')
set(gca,'LineWidth',2)
set(gca,'Fontsize',fs)
ylabel('\tau-RR')
xlabel('period [a.u.]')
xlim([0 1000])

subplot(337)
plot((1:length(spectrum1_n)),spectrum1_n(i,:), 'Linewidth', lw)
grid on
title('inter spike spectrum of \tau-RR')
set(gca,'LineWidth',2)
set(gca,'Fontsize',fs)
ylabel('rel. power/period')
xlabel('period [a.u.]')

subplot(338)
plot((1:length(spectrum1_n)),spectrum2_n(i,:), 'Linewidth', lw)
grid on
title('inter spike spectrum of \tau-RR')
set(gca,'LineWidth',2)
set(gca,'Fontsize',fs)
ylabel('rel. power/period')
xlabel('period [a.u.]')

subplot(339)
plot((1:length(spectrum1_n)),spectrum3_n(i,:), 'Linewidth', lw)
grid on
title('inter spike spectrum of \tau-RR')
set(gca,'LineWidth',2)
set(gca,'Fontsize',fs)
ylabel('rel. power/period')
xlabel('period [a.u.]')


figure('Units','normalized','Position',[.01 .01 .6 .99])
% 
% subplot(331)
% imagesc(RP_1), colormap([1 1 1 ; 0 0 0])
% set(gca, 'YDir', 'normal')
% set(gca,'LineWidth',2)
% set(gca,'Fontsize',fs)
% title('RP')
% xlabel('time [a.u.]')
% ylabel('time [a.u.]')
% 
% subplot(332)
% imagesc(RP_2), colormap([1 1 1 ; 0 0 0])
% set(gca, 'YDir', 'normal')
% set(gca,'LineWidth',2)
% set(gca,'Fontsize',fs)
% title('RP')
% xlabel('time [a.u.]')
% ylabel('time [a.u.]')
% 
% subplot(333)
% imagesc(RP_3), colormap([1 1 1 ; 0 0 0])
% set(gca, 'YDir', 'normal')
% set(gca,'LineWidth',2)
% set(gca,'Fontsize',fs)
% title('RP')
% xlabel('time [a.u.]')
% ylabel('time [a.u.]')
% 
% subplot(334)
% plot((1:length(spectrum1)),tau_rr1, 'Linewidth', lw)
% grid on
% title('\tau-RR')
% set(gca,'LineWidth',2)
% set(gca,'Fontsize',fs)
% ylabel('\tau-RR')
% xlabel('period [a.u.]')
% xlim([0 1000])
% 
% subplot(335)
% plot((1:length(spectrum2)),tau_rr2, 'Linewidth', lw)
% grid on
% title('\tau-RR')
% set(gca,'LineWidth',2)
% set(gca,'Fontsize',fs)
% ylabel('\tau-RR')
% xlabel('period [a.u.]')
% xlim([0 1000])
% 
% subplot(336)
% plot((1:length(spectrum3)),tau_rr3, 'Linewidth', lw)
% grid on
% title('\tau-RR')
% set(gca,'LineWidth',2)
% set(gca,'Fontsize',fs)
% ylabel('\tau-RR')
% xlabel('period [a.u.]')
% xlim([0 1000])

subplot(337)
plot(1./fs1_n,10*log10(perio1_n),'.-', 'Linewidth', lw2)
grid on
title('Fourier powerspectrum of \tau-RR')
set(gca,'LineWidth',2)
set(gca,'Fontsize',fs)
ylabel('rel. power/period')
xlabel('period [a.u.]')
xlim([0 500])

subplot(338)
plot(1./fs1_n,10*log10(perio2_n), '.-', 'Linewidth', lw2)
grid on
title('Fourier powerspectrum of \tau-RR')
set(gca,'LineWidth',2)
set(gca,'Fontsize',fs)
ylabel('rel. power/period')
xlabel('period [a.u.]')
xlim([0 500])

subplot(339)
plot(1./fs1_n,10*log10(perio1_n), '.-', 'Linewidth', lw2)
grid on
title('Fourier powerspectrum of \tau-RR')
set(gca,'LineWidth',2)
set(gca,'Fontsize',fs)
ylabel('rel. power/period')
xlabel('period [a.u.]')
xlim([0 500])

