function Z = acorr(varargin)
% Input:
% Minimum input-arguments : 1
% Maximum input-arguments : 3
% Z = acorr(Y,tau_max,Show)
% 
%
% This function computes the autocorrelationvalues of timeseries Y versus the 
% time delayed version of itself. If not specified in the second input-argument 
% the maximum  timelag value 'tau_max' is set to 10.
% The optional output contains the timelag-values (first column) and
% the corresponding autocorrelation-values (second column)
% If input 'Show' is set to 1 (Default is 0) the autocorrelation as a
% function of the time lag is plotted along with a horizontal line marking
% the 1/exp-decay.
%
% K.H.Kraemer, 2019

%% Assign input

y = varargin{1};

try
    tau_max = varargin{2};
catch 
    tau_max = 10;
end

try
    show = varargin{3};
catch 
    show = 0;
end

N = length(y);


%% Check input
narginchk(1,3)
nargoutchk(0,1)

% Check if tau-values are positive integers
if rem(tau_max,1)~=0
   error('tau-values must be positive integers.')
end

if tau_max<1
    error('tau-values must be positive integers.')
end
   
% Check if input-vector is a column- or line-vector
if size(y,1)<size(y,2)
    y=y';
end

if show~=0 && show~=1
    warning('input show needs to be 1 (display figure) or 0 (no figure displayed). Now set to 0.')
    show = 0;
end

tau_min = 0;

%% Calculate Autocorrelation

mu=nanmean(y);
v=nanvar(y);

auto=zeros(1,length(tau_min:tau_max));
cnt = 1;
for i=tau_min:tau_max
    
    % lag the time series 

    N2 = N-i;
    y_ = y(1:N2,:);
    y_lag = y(1+i:N2+i,:);
        
    auto(cnt) = ((y_lag-mu)'*(y_-mu)/N2) / v;
    cnt = cnt + 1; 
end

Z(:,1)= tau_min:tau_max;
Z(:,2)= auto(:);

%% Plotting Autocorrelationfunction against lag tau
if show
    figure
    plot(Z(:,1),Z(:,2),'-.*','LineWidth',2), hold on    
    thres1 = Z(1,2)*exp(-1);
    yline(thres1);
    grid on
    title('Autocorrelationfunction')
    xlabel('timelag \tau')
    ylabel('Autocorrelation')
    set(gca,'LineWidth',2)
    set(gca,'FontSize',12)
end



end