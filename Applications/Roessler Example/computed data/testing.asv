clear, clc

a = load("tauRR_spectrum1_roessler_0.9_1.csv");
b = load("tauRR_spectrum1_roessler_0.9_50.csv");

subplot(211)
plot(a)
subplot(212)
plot(b)


tau_1 = load("tau_rr_1_all.csv");

%%
subplot(211)
plot(tau_1(1,:))
subplot(212)
plot(tau_1(51,:))

s1 = inter_spike_spectrum(tau_1(1,:), 'type', "normal", 'threshold', 0.9);
s2 = inter_spike_spectrum(tau_1(51,:), 'type', "normal", 'threshold', 0.9);

%%
subplot(211)
plot(s1)
subplot(212)
plot(s2)