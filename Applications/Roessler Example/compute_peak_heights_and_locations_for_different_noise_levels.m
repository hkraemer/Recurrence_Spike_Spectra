% For the obtained tau-recurrence spectra for the different noise levels of
% all three Roessler setups (computed in 
% `compute_roessler_spike_spectra_for_different_noise_levels.jl`) 
% we find the peaks and their positions.

clear, clc

num_levels = 51; % number of noise levels
thresholds = [0.9; 0.95; 0.99]; % regularization thresholds


% peak prominence parameter for findpeaks-function
peak_prominence = 0.01;

% preallocation
peak_heights = cell(3, 3, num_levels);
peak_locs = cell(3, 3, num_levels);

% for each spectrum compute the peak heights and their locations
for system = 1:3
    for i = 1:length(thresholds)
        rho = thresholds(i);
        for k = 1:num_levels
            spectrum = load(strcat("./computed data/tauRR_spectrum",num2str(system),"_roessler_",num2str(rho),"_",num2str(k),".csv"));
            [peak_heights{system, i, k}, peak_locs{system, i, k}] = findpeaks(spectrum, 'MinPeakProminence', peak_prominence);
        end
    end
end

%% Plot the output

rho = 2; % corresponds to thresholds(rho)


factor = 300; % scatter size
c = lines(100); % matlabs colormap
% c = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4660 0.6740 0.1880];

sigmas = linspace(0,0.5,num_levels).*100; % noise levels
fs = 20; % fontsize for plotting

% the different considered systems
systems = ["Rössler system (a = 0.36)", "Rössler system (a = 0.41)", "Rössler system (a = 0.428)"];


figure('Units', 'normalized', 'Position', [.001 .2 .999 .6])

for i = 1:3
    heights = squeeze(peak_heights(i, rho, :));
    locs = squeeze(peak_locs(i, rho, :));

    subplot(1,3,i)
    for j = 1:length(heights)
        loc_values = locs{j};
        height_values = heights{j};
    
        for k = 1:length(loc_values)
            scatter(loc_values(k), j, height_values(k)*factor, c(k,:), 'filled'), hold on
        end
    end
    yticklabels([0 10 20 30 40 50])
    ylabel("percentage of additive noise")
    xlabel("peak position")
    title(systems(i))
    set(gca, 'LineWidth', 2, 'FontSize', fs)
    grid on
    ylim([0 50])
    xlim([0 500])
end

