clear, clc

Y = load("Y.csv");
tau_rr = load("tau_rr.csv");
tau_pec = load("tau_pec.csv");
RP = load("RP.csv");
perio_power = load("perio_power.csv");
perio_freq = load("perio_freq.csv");
perio_power2 = load("perio_power2.csv");
perio_freq2 = load("perio_freq2.csv");
data = load("data.csv");

%% Plot periodogram of x-component

periodogram(data(:,1))

%% Plot periodogram of tau-RR

periodogram(tau_rr)

%%
subset = 1:length(tau_rr);

t = 1:length(tau_rr(subset));

lw = 3;
fs = 20;

figure('Units','normalized','Position',[.1 .1 .7 .8])
subplot(3,3,[1,4])
plot(tau_rr(subset),t, 'LineWidth', lw)
set(gca,'XDir','reverse');
set(gca,'LineWidth',2)
set(gca,'FontSize',fs)
xlabel("\tau-RR")
ylabel("time [\Deltat]")
xlim([0 1])
yticklabels([])
grid on

subplot(3,3,[2,3,5,6])
imagesc(t,t,RP), colormap([0 0 0; 1 1 1])
title("Recurrence Plot (8% rec.-rate)")
xticklabels([])
yticklabels([])
set(gca,'LineWidth',2)
set(gca,'FontSize',fs)
grid on

subplot(3,3,[8,9])
plot(t,tau_rr(subset), 'LineWidth', lw)
set(gca,'YDir','reverse');
set(gca,'LineWidth',2)
set(gca,'FontSize',fs)
ylabel("\tau-RR")
ylim([0 1])
xlabel("time [\Deltat]")
xticklabels([])
grid on


figure
subplot(211)
plot(perio_freq, perio_power.^2, 'LineWidth', lw)
xlabel("frequency [1/\Deltat]")
ylabel("Power")
title("Periodogram of x-component")
set(gca,'LineWidth',2)
set(gca,'FontSize',fs)
%set(gca, 'YScale', 'log')
xlim([0 4])
grid on
subplot(212)
plot(perio_freq2, perio_power2.^2, 'LineWidth', lw)
xlabel("frequency [1/\Deltat]")
ylabel("Power")
title("Periodogram of \tau-RR")
xlim([0 4])
ylim([0 0.1])
set(gca,'LineWidth',2)
set(gca,'FontSize',fs)
%set(gca, 'YScale', 'log')
grid on

figure
subplot(311)
plot(data(:,1), 'LineWidth', lw)
xlim([0 6000])
xlabel("time [\Deltat]")
title("x-component of Lorenz system")
set(gca,'LineWidth',2)
set(gca,'FontSize',fs)
grid on

subplot(3,1,[2,3])
plot3(Y(:,1), Y(:,2), Y(:,3))
xlabel("x(t)")
ylabel(strcat("x(t+",num2str(tau_pec(2)),")"))
zlabel(strcat("x(t+",num2str(tau_pec(3)),")"))
title("Reconstruction from x-component of Lorenz system")
set(gca,'LineWidth',2)
set(gca,'FontSize',fs)
grid on

%% Roessler limit cycle example
clear, clc
Y = load("data_roessler.csv");

RP = rp(Y,0.15, 'var');
RP_d = rp_diagonal(RP);

%%
% save("RP_limit_roessler.csv","RP","-ascii")
% save("RP_skeleton_limit_roessler.csv","RP_d","-ascii")

%% compute tau-RRs

% compute tau recurrence rate
N = size(RP,1);
tauRR = zeros(1,N);
tauRR_d = zeros(1,N);
cnt = 1;
for k = 0:N-1
    d = diag(RP,k);
    d2 = diag(RP_d,k);
    tauRR(cnt) = length(find(d))/length(d);
    tauRR_d(cnt) = length(find(d2))/length(d2);
    clear d, clear d2
    cnt = cnt + 1;
end

%%
% save("tauRR_limit_roessler.csv","tauRR","-ascii")
% save("tauRR_skeleton_limit_roessler.csv","tauRR_d","-ascii")

%% Plot 

fs = 20;
lw = 3;

perio_freq = load("perio_roessler_freq.csv");
perio_power = load("perio_roessler_power.csv");
perio_freq_d = load("perio_d_roessler_freq.csv");
perio_power_d = load("perio_d_roessler_power.csv");
tau_RR_spectrum = load("tauRR_spectrum_roessler.csv");
tau_RR_spectrum_d = load("tauRR_spectrum_d_roessler.csv");


subset = 1:500;

figure('Units','normalized','Position',[.1 .1 .9 .5])
subplot(131)
plot3(Y(:,1),Y(:,2),Y(:,3))
title("Rössler limit cycle attractor")
xticklabels([])
yticklabels([])
zticklabels([])
set(gca, 'LineWidth',2)
set(gca, 'FontSize',fs)
grid on

subplot(132)
imagesc(subset, subset, RP(subset,subset)), colormap([1 1 1; 0 0 0]), axis xy square
title("RP (subset)")
xticklabels([])
yticklabels([])
set(gca, 'LineWidth',2)
set(gca, 'FontSize',fs)


subplot(133)
imagesc(subset, subset, RP_d(subset,subset)), colormap([1 1 1; 0 0 0]), axis xy square
title("skeletonized RP (subset)")
xticklabels([])
yticklabels([])
set(gca, 'LineWidth',2)
set(gca, 'FontSize',fs)

figure('Units','normalized','Position',[.1 .1 .7 .7])
subplot(321)
plot(tauRR(subset))
title("\tau-RR of RP")
set(gca, 'LineWidth',2)
set(gca, 'FontSize',fs)
xlabel("time [\Deltat]")
grid on

subplot(322)
plot(tauRR_d(subset))
title("\tau-RR of skeletonized RP")
set(gca, 'LineWidth',2)
set(gca, 'FontSize',fs)
xlabel("time [\Deltat]")
grid on

subplot(323)
plot(1./(perio_freq*0.2), perio_power, 'LineWidth', lw)
title("Periodogram of \tau-RR")
set(gca, 'LineWidth',2)
set(gca, 'FontSize',fs)
xlim([0 100])
xlabel("time [\Deltat]")
grid on

subplot(324)
plot(1./(perio_freq_d*0.2), perio_power_d, 'LineWidth', lw)
title("Periodogram of skeletonized \tau-RR")
set(gca, 'LineWidth',2)
set(gca, 'FontSize',fs)
xlim([0 100])
xlabel("time [\Deltat]")
grid on

subplot(325)
plot(tau_RR_spectrum, 'LineWidth', lw)
title("Spectrum of \tau-RR")
set(gca, 'LineWidth',2)
set(gca, 'FontSize',fs)
xlim([0 100])
xlabel("time [\Deltat]")
grid on

subplot(326)
plot(tau_RR_spectrum_d, 'LineWidth', lw)
title("Spectrum of skeletonized \tau-RR")
set(gca, 'LineWidth',2)
set(gca, 'FontSize',fs)
xlim([0 100])
xlabel("time [\Deltat]")
grid on

