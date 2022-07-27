function basis = generate_basis_functions(N)
% GENERATE_BASIS_FUNCTION generates a matrix containing all (shifted)
% basis functions of a spike train (Dirac-Comb) of length N.
%    
%    basis = GENERATE_BASIS_FUNCTION(N) generates spike-trains (Dirac-Combs) 
%    of all possible periods (=interspike interval), including all possible
%    time shifts of those spike trains. This is a (M-by-N)-matrix, with M
%    being `sum(1:ceil(N/2))+1`.
%
%    Further reading:
%    H. K. Kraemer et al., … 2022
%
%    See also ....

% Copyright (c) 2022
% K. Hauke Kraemer
% Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% This program is free software and runs under MIT licence.
%%
narginchk(1,1)
nargoutchk(1,1)
assert(isscalar(N),"Input N must be a postive integer number")
assert(mod(N,1) == 0,"Input N must be a postive integer number")
assert(N > 0,"Input N must be a postive integer number")

%%
num_of_basis_functions = sum(1:ceil(N/2))+1;
basis = zeros(num_of_basis_functions, N);
cnt = 1;
for i = 1:ceil(N/2)
   basis(cnt:cnt+i-1,:) = create_single_basis_function(N, i);
   cnt = cnt + i;
end
temp = create_single_basis_function(N, ceil(N/2)+1);
basis(end,:) = temp(1,:);
end