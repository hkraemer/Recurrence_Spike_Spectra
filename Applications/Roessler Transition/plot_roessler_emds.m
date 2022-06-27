% After computing the EMDs of the spectra in
% `compute_emd_roessler_transitions.jl`, we here plot those pairwise
% distances as a recurrence plot.

clear, clc

cd('/Users/hkraemer/Documents/Git/Recurrence_Spike_Spectra/Applications/Roessler Transition')

emd1 = load("./results/emd1.csv");
emd2 = load("./results/emd2.csv");
emd3 = load("./results/emd3.csv");
emd4 = load("./results/emd4.csv");

emd1t = load("./results/emd1t.csv");
emd2t = load("./results/emd2t.csv");
emd3t = load("./results/emd3t.csv");
emd4t = load("./results/emd4t.csv");

FFT1 = load("./results/emd_FFT1.csv");
FFT2 = load("./results/emd_FFT2.csv");
FFT3 = load("./results/emd_FFT3.csv");

%% Plot Lyapunov spectrum and mark the transitions

as = 0.36:0.0001:0.43;
lyap = load("./results/Lyaps_Roessler_0_36_to_0_43.csv");
lyaps = load("./results/All_lyaps_Roessler_0_36_to_0_43.csv");

figure
plot(as,lyap)
grid on

figure
plot(as,lyaps(1,:)), hold on
plot(as,lyaps(2,:)), hold on
plot(as,lyaps(3,:)), hold on
legend("\lambda_1","\lambda_2","\lambda_3")
grid on
xlabel("control parameter a")

%%
transitions = [0.3747 0.3833 0.3856 0.3858 0.3859 0.3868 0.3874 0.3889 0.3907 0.391 0.3924 0.3967 0.4 0.4004 0.4022 0.4069 0.4092 0.4109 0.4118 0.412 0.4183 0.4184 0.4217 0.4229 0.4249 0.4251 0.43]; 


%% Compute RPs an perform a running window for computing RR and RT

% params for windowed analysis
w = 40;
ws = 1;
% recurrence selection method
method = "var";
% recurrence threshold
e = 0.1;

[RR_emd1, RT_emd1, RP_emd1] = compute_rp_and_rqa(emd1, w, ws, method, e, as);
[RR_emd2, RT_emd2, RP_emd2] = compute_rp_and_rqa(emd2, w, ws, method, e, as);
[RR_emd3, RT_emd3, RP_emd3] = compute_rp_and_rqa(emd3, w, ws, method, e, as);
[RR_emd4, RT_emd4, RP_emd4] = compute_rp_and_rqa(emd4, w, ws, method, e, as);

[RR_emd1t, RT_emd1t, RP_emd1t] = compute_rp_and_rqa(emd1t, w, ws, method, e, as);
[RR_emd2t, RT_emd2t, RP_emd2t] = compute_rp_and_rqa(emd2t, w, ws, method, e, as);
[RR_emd3t, RT_emd3t, RP_emd3t] = compute_rp_and_rqa(emd3t, w, ws, method, e, as);
[RR_emd4t, RT_emd4t, RP_emd4t] = compute_rp_and_rqa(emd4t, w, ws, method, e, as);

[RR_fft1, RT_fft1, RP_fft1] = compute_rp_and_rqa(FFT1, w, ws, method, e, as);
[RR_fft2, RT_fft2, RP_fft2] = compute_rp_and_rqa(FFT2, w, ws, method, e, as);
[RR_fft3, RT_fft3, RP_fft3] = compute_rp_and_rqa(FFT3, w, ws, method, e, as);


%% Plot the final results

lw1 = 0.5;

figure('Units','normalized','Position',[.01 .01 .99 .99])
subplot(3,4,[1,2,5,6])
imagesc(as,as,RP_emd1), colormap([1 1 1; 0 0 0]), axis xy square
% for j = 1:length(transitions)
%     xline(transitions(j),'r--','linewidth',lw1)
% end
title("EMD 1")
grid on

subplot(3,4,[3,4,7,8])
imagesc(as,as,RP_emd2), colormap([1 1 1; 0 0 0]), axis xy square
% for j = 1:length(transitions)
%     xline(transitions(j),'r--','linewidth',lw1)
% end
title("EMD 2")
grid on

subplot(3,4,[9,10])
plot(as,RR_emd1), hold on
% for j = 1:length(transitions)
%     xline(transitions(j),'r--','linewidth',lw1)
% end
xlim([as(1) as(end)])
grid on
subplot(3,4,[11,12])
plot(as,RR_emd2), hold on
% for j = 1:length(transitions)
%     xline(transitions(j),'r--','linewidth',lw1)
% end
xlim([as(1) as(end)])
grid on


figure('Units','normalized','Position',[.01 .01 .99 .99])
subplot(3,4,[1,2,5,6])
imagesc(as,as,RP_emd3), colormap([1 1 1; 0 0 0]), axis xy square
% for j = 1:length(transitions)
%     xline(transitions(j),'r--','linewidth',lw1)
% end
title("EMD 3")
grid on

subplot(3,4,[3,4,7,8])
imagesc(as,as,RP_emd4), colormap([1 1 1; 0 0 0]), axis xy square
% for j = 1:length(transitions)
%     xline(transitions(j),'r--','linewidth',lw1)
% end
title("EMD 4")
grid on

subplot(3,4,[9,10])
plot(as,RR_emd3), hold on
% for j = 1:length(transitions)
%     xline(transitions(j),'r--','linewidth',lw1)
% end
xlim([as(1) as(end)])
grid on
subplot(3,4,[11,12])
plot(as,RR_emd4), hold on
% for j = 1:length(transitions)
%     xline(transitions(j),'r--','linewidth',lw1)
% end
xlim([as(1) as(end)])
grid on

figure('Units','normalized','Position',[.01 .01 .99 .99])
subplot(3,4,[1,2,5,6])
imagesc(as,as,RP_fft1), colormap([1 1 1; 0 0 0]), axis xy square
% for j = 1:length(transitions)
%     xline(transitions(j),'r--','linewidth',lw1)
% end
title("FFT 1")
grid on

subplot(3,4,[3,4,7,8])
imagesc(as,as,RP_fft3), colormap([1 1 1; 0 0 0]), axis xy square
% for j = 1:length(transitions)
%     xline(transitions(j),'r--','linewidth',lw1)
% end
title("FFT 3")
grid on

subplot(3,4,[9,10])
plot(as,RR_fft1), hold on
% for j = 1:length(transitions)
%     xline(transitions(j),'r--','linewidth',lw1)
% end
xlim([as(1) as(end)])
grid on
subplot(3,4,[11,12])
plot(as,RR_fft3), hold on
% for j = 1:length(transitions)
%     xline(transitions(j),'r--','linewidth',lw1)
% end
xlim([as(1) as(end)])
grid on

%% Helper function

function [recrate,rectimes,y] = compute_rp_and_rqa(P, w, ws, method, e, as)

len = length(P);
% compute recurrence matrix
if strcmp(method,"var")
    epsilon = quantile(P(:),e);
    y = double(P < epsilon); 
else
    q = quantile(P,e); % distance that corresponds to the fraction e of rec. points per column
    thresholds = repmat(q,len,1); % q has to be applied for each row in d
    % apply individual thresholds
    % apply threshold(s)
    y=double(P<thresholds);
end

numRec = (w+1)*(w+1);
% apply a running window and measure RR
rr_win = zeros(1,ceil((len-w)/ws)); % preallocation
rt_win = zeros(1,ceil((len-w)/ws));
for k = 1:ws:len-w
    
    R = y(k:(k+w),k:(k+w));
    
    rr_win(k) = sum(R(:))/numRec;
    
    N = size(R);
    % recurrence times ("white" vertical lines)
    rt_hist = zeros(1,N(1)); % allocate vector
    for i = 1:N(1)
       cnt = 0;

       % boolean variable to avoid counting white lines at the edges of RP
       first_flag = false;

       for j = 1:N(2)
          if ~R(j,i) % are we on a white line?
             if first_flag % line does not cross the RP's edges
                cnt = cnt + 1; % count number of points along the vertical line
             end
          else % we meet a recurrence point
             first_flag = true; % we are for sure within the RP
             if cnt
                 rt_hist(cnt) = rt_hist(cnt) + 1; % store line length
             end
             cnt = 0;
          end
       end
    end
    rt_win(k) = sum(rt_hist .* (1:N(1))) / sum(rt_hist);
end
recrate = NaN*ones(size(as));
recrate(w/2+1:end-w/2) = rr_win;
rectimes = NaN*ones(size(as));
rectimes(w/2+1:end-w/2) = rt_win;
end