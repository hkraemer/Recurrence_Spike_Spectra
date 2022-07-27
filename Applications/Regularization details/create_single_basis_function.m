function cs = create_single_basis_function(N, period)
% CREATE_SINGLE_BASIS_FUNCTION generates a matrix containing all (shifted)
% basis functions of a spike train (Dirac-Comb) of length N.
%    
%    basis = generate_basis_functions(N, period) generates
%    (period-by-N)-matrix containing spike-trains (Dirac-Combs) with
%    interspike-interval `period` in all possible shifted variants.
%
%    Further reading:
%    H. K. Kraemer et al., … 2021
%
%    See also ....

% Copyright (c) 2021
% K. Hauke Kraemer
% Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% This program is free software and runs under MIT licence.
%% check input
narginchk(2,2)
nargoutchk(1,1)
assert(isscalar(N),"Input N must be a postive integer number")
assert(mod(N,1) == 0,"Input N must be a postive integer number")
assert(N > 0,"Input N must be a postive integer number")
assert(isscalar(period),"Input period must be a postive integer number")
assert(mod(period,1) == 0,"Input period must be a postive integer number")
assert(period > 0,"Input period must be a postive integer number")
assert(N>=period, "N must be larger or equal to the given priod")
%%
cs = zeros(period,N);
for i = 1:period:N
    cs(1,i) = 1;
end

css = horzcat(cs(1,:),zeros(1,period-1));
for j = 1:period-1
    inp_array = circshift(css,j);
    cs(j+1,:) = inp_array(1:N);
end

end