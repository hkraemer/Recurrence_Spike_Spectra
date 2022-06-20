function dydt = ode_4DOF(t, y, m, mu, kx, ky, kz1, kz2, klin1, klin2, knl1, knl2, cx, cy, cz1, cz2, clin1, clin2)

% definition of the extended Hoffmann-Gaul minimal model according to

% Stender, Merten, Merten Tiedemann, and Norbert Hoffmann. "Characterization 
% of complex states for friction?excited systems." PAMM 17.1 (2017): 45-46.

% kz2 = kz2-(t*1/600)*100;

% transformation matrices for extension masses 1 and 2
T1 = [1/2 1/2 -1/sqrt(2) 0;
    1/2 1/2 -1/sqrt(2) 0;
    -1/sqrt(2) -1/sqrt(2) 1 0;
    0 0 0 0];

T2 = [1/2 1/2 0 -1/sqrt(2);
    1/2 1/2 0 -1/sqrt(2);
    0 0 0 0;
    -1/sqrt(2) -1/sqrt(2) 0 1];

M = diag([m, m, m, m]); 
C = diag([cx, cy, cz1, cz2])+T1.*clin1+T2.*clin2;
K = diag([kx, ky, kz1, kz2])+[0 -mu*ky 0 0; zeros(3,4)]+T1.*klin1+T2.*klin2;

% displacement in the joint
u1 = (y(3)-(y(1)+y(2))/sqrt(2))^3;
u2 = (y(4)-(y(1)+y(2))/sqrt(2))^3;

% cubic force terms
fnl = [-1/sqrt(2).*ones(2,2); eye(2,2)]*[knl1.*u1; knl2.*u2];

% definition of the first order system of equations
dydt = [zeros(4), eye(4); 
    -K/M, -C/M]*y-[zeros(4,1); M\fnl];

end

