% We plot the distances of all noisy spectra with respect to the noise free
% case for all three considered systems (limit-2, limit-3, chaos). We
% considered 50 different noise levels ranging from 1%-50% additive noise.

% The noisy time series along with their recurrence plots and their
% tau-recurrence rate have been computed in 
% `compute_roessler_trajectories_and_tau_rr_for_different_noise_levels.m`

% The corresponding inter spike spectra, using three different
% regularization thresholds of 0.9, 0.95 and 0.99 have been computed in 
% `compute_roessler_spike_spectra_distances.jl`. We used a normal lasso
% regression.

% The distances of all noisy spectra compared to the noise-free spectrum
% have been computed in `compute_roessler_spike_spectra_distances.jl`

clear, clc

num_levels = 50; % number of noise levels (excluding 0)
sigmas = linspace(0.01,0.5,num_levels).*100; % noise levels

fs = 25; % fontsize for plotting
lw = 2.5; % linewidth for plotting
sz = 40; % scatter size
% the different considered systems
systems = ["Rössler system (a = 0.36)", "Rössler system (a = 0.41)", "Rössler system (a = 0.428)"];

%% Cosine distance

for i = 1:3
    dist1 = load(strcat("./computed data/dist_cosine_rho_0.9_system_", num2str(i), ".csv"));
    dist2 = load(strcat("./computed data/dist_cosine_rho_0.95_system_", num2str(i), ".csv"));
    dist3 = load(strcat("./computed data/dist_cosine_rho_0.99_system_", num2str(i), ".csv"));
    
    figure('Units', 'normalized', 'Position', [.2 .2 .6 .6])
    p1 = plot(sigmas, dist1, 'LineWidth', lw, 'LineStyle', '-'); hold on
    scatter(sigmas, dist1, sz, p1.Color, 'filled'), hold on
    p2 = plot(sigmas, dist2, 'LineWidth', lw, 'LineStyle', '-'); hold on
    scatter(sigmas, dist2, sz, p2.Color, 'filled'), hold on
    p3 = plot(sigmas, dist3, 'LineWidth', lw, 'LineStyle', '-'); hold on
    scatter(sigmas, dist3, sz, p3.Color, 'filled')
    xlabel("percentage of additive noise")
    ylabel("Cosine Distance")
    title(systems(i))
    legend([p1(1), p2(1), p3(1)], "reg. threshold = 0.90", "reg. threshold = 0.95", "reg. threshold = 0.99")
    legend('Location','northwest')
    grid on
    set(gca, 'LineWidth', 2, 'FontSize', fs)
    ylim([0 0.6])
end

%% KL Divergence

for i = 1:3
    dist1 = load(strcat("./computed data/dist_kl_div_rho_0.9_system_", num2str(i), ".csv"));
    dist2 = load(strcat("./computed data/dist_kl_div_rho_0.95_system_", num2str(i), ".csv"));
    dist3 = load(strcat("./computed data/dist_kl_div_rho_0.99_system_", num2str(i), ".csv"));
    
    figure('Units', 'normalized', 'Position', [.2 .2 .6 .6])
    p1 = plot(sigmas, dist1, 'LineWidth', lw, 'LineStyle', '-'); hold on
    scatter(sigmas, dist1, sz, p1.Color, 'filled'), hold on
    p2 = plot(sigmas, dist2, 'LineWidth', lw, 'LineStyle', '-'); hold on
    scatter(sigmas, dist2, sz, p2.Color, 'filled'), hold on
    p3 = plot(sigmas, dist3, 'LineWidth', lw, 'LineStyle', '-'); hold on
    scatter(sigmas, dist3, sz, p3.Color, 'filled')
    xlabel("percentage of additive noise")
    ylabel("KL divergence")
    title(systems(i))
    legend([p1(1), p2(1), p3(1)], "reg. threshold = 0.90", "reg. threshold = 0.95", "reg. threshold = 0.99")
    legend('Location','northwest')
    grid on
    set(gca, 'LineWidth', 2, 'FontSize', fs)
    ylim([0 0.6])
end

%% Hellinger distance

for i = 1:3
    dist1 = load(strcat("./computed data/dist_hellinger_rho_0.9_system_", num2str(i), ".csv"));
    dist2 = load(strcat("./computed data/dist_hellinger_rho_0.95_system_", num2str(i), ".csv"));
    dist3 = load(strcat("./computed data/dist_hellinger_rho_0.99_system_", num2str(i), ".csv"));
    
    figure('Units', 'normalized', 'Position', [.2 .2 .6 .6])
    p1 = plot(sigmas, dist1, 'LineWidth', lw, 'LineStyle', '-'); hold on
    scatter(sigmas, dist1, sz, p1.Color, 'filled'), hold on
    p2 = plot(sigmas, dist2, 'LineWidth', lw, 'LineStyle', '-'); hold on
    scatter(sigmas, dist2, sz, p2.Color, 'filled'), hold on
    p3 = plot(sigmas, dist3, 'LineWidth', lw, 'LineStyle', '-'); hold on
    scatter(sigmas, dist3, sz, p3.Color, 'filled')
    xlabel("percentage of additive noise")
    ylabel("Hellinger distance")
    title(systems(i))
    legend([p1(1), p2(1), p3(1)], "reg. threshold = 0.90", "reg. threshold = 0.95", "reg. threshold = 0.99")
    legend('Location','northwest')
    grid on
    set(gca, 'LineWidth', 2, 'FontSize', fs)
    ylim([0 0.6])
end

%% Wasserstein

for i = 1:3
    dist1 = load(strcat("./computed data/dist_wasserstein_rho_0.9_system_", num2str(i), ".csv"));
    dist2 = load(strcat("./computed data/dist_wasserstein_rho_0.95_system_", num2str(i), ".csv"));
    dist3 = load(strcat("./computed data/dist_wasserstein_rho_0.99_system_", num2str(i), ".csv"));
    
    figure('Units', 'normalized', 'Position', [.2 .2 .6 .6])
    p1 = plot(sigmas, dist1, 'LineWidth', lw, 'LineStyle', '-'); hold on
    scatter(sigmas, dist1, sz, p1.Color, 'filled'), hold on
    p2 = plot(sigmas, dist2, 'LineWidth', lw, 'LineStyle', '-'); hold on
    scatter(sigmas, dist2, sz, p2.Color, 'filled'), hold on
    p3 = plot(sigmas, dist3, 'LineWidth', lw, 'LineStyle', '-'); hold on
    scatter(sigmas, dist3, sz, p3.Color, 'filled')
    xlabel("percentage of additive noise")
    ylabel("Wasserstein distance")
    title(systems(i))
    legend([p1(1), p2(1), p3(1)], "reg. threshold = 0.90", "reg. threshold = 0.95", "reg. threshold = 0.99")
    legend('Location','northwest')
    grid on
    set(gca, 'LineWidth', 2, 'FontSize', fs)

end