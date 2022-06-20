function out = matcont_model_4DOF
out{1} = @init;
out{2} = @fun_eval;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t, y, m, mu, kx, ky, kz1, kz2, klin1, klin2, knl1, knl2, cx, cy, cz1, cz2, clin1, clin2)

% call the "one and only" ODE definition in order to stay consistent
[dydt] = ode_4DOF(t, y, m, mu, kx, ky, kz1, kz2, klin1, klin2, knl1, knl2, cx, cy, cz1, cz2, clin1, clin2);

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
y0=[1, 0, 0, 0, 0, 0, 0, 0];
handles = feval(matcont_model_4DOF);
% options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
options = odeset('RelTol', 1e-6);
tspan = [0 100];

% --------------------------------------------------------------------------
function jac = jacobian(t, y, m, mu, kx, ky, kz1, kz2, klin1, klin2, knl1, knl2, cx, cy, cz1, cz2, clin1, clin2)
% --------------------------------------------------------------------------
function jacp = jacobianp(t, y, m, mu, kx, ky, kz1, kz2, klin1, klin2, knl1, knl2, cx, cy, cz1, cz2, clin1, clin2)
% --------------------------------------------------------------------------
function hess = hessians(t, y, m, mu, kx, ky, kz1, kz2, klin1, klin2, knl1, knl2, cx, cy, cz1, cz2, clin1, clin2)
% --------------------------------------------------------
function hessp = hessiansp(t, y, m, mu, kx, ky, kz1, kz2, klin1, klin2, knl1, knl2, cx, cy, cz1, cz2, clin1, clin2)
%---------------------------------------------------------------------------
function tens3  = der3(t, y, m, mu, kx, ky, kz1, kz2, klin1, klin2, knl1, knl2, cx, cy, cz1, cz2, clin1, clin2)
%---------------------------------------------------------------------------
function tens4  = der4(t, y, m, mu, kx, ky, kz1, kz2, klin1, klin2, knl1, knl2, cx, cy, cz1, cz2, clin1, clin2)
%---------------------------------------------------------------------------
function tens5  = der5(t, y, m, mu, kx, ky, kz1, kz2, klin1, klin2, knl1, knl2, cx, cy, cz1, cz2, clin1, clin2)
% -------------------------------------------------------------------------

