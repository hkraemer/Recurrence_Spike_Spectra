function [t2,x2] = Rossler(varargin)
% This function numerically integrates the Roessler-System with input
% paramters a, b and c starting with input inital condition 
% x0 = [x(t0); y(t0); z(t0)]. If you wish to receive a Phase Space Portrait
% of the integrated System set the optional input parameter Show_PS = 1.
%
% [t,x] = Rossler(a,b,c,x0,N,dt,transients,Show_PS)
%
% t is the time vector containing just the times, where the samping took
% place. The corresponding values to these t2-times are stored in x.
%
% Here the sampling rate is dt and the output is a vector-series of 
% length N (input parameter) and `transients` samples removed.


% bind input
a = varargin{1};
b = varargin{2};
c = varargin{3};
x0= varargin{4};
N = varargin{5};
dt = varargin{6};
transients = varargin{7};
try
    Show_PS = varargin{8};
catch
    Show_PS = 0;
end


% Define your sampling time
time_sample = dt;

time_interval = (N+transients)* time_sample;     % k+1 is the number of columns in the solution matrix M

% This results in a corresponding time vector
t2 = 0:time_sample:time_interval;

 
f = @(t,x) [-x(2)-x(3); x(1)+a*x(2); b+x(3)*(x(1)-c)];  % Rossler system
SOL = ode45(f,[t2(1) t2(end)],x0);
x2 = deval(SOL,t2);

% Store sampled time series    
x2 = x2(:,transients+1:end-1);
t2 = t2(1:end-transients+1);

% Plot Phase Space
if Show_PS == 1
    % Plot sampled Phase Space
    figure
    plot3(x2(1,:),x2(2,:),x2(3,:))
    xlabel('x_1')
    ylabel('x_2')
    zlabel('x_3')
    title('numerically integrated Roessler Attractor')
    grid on

end

end
