function reg = regenerate_signal(Theta, coefs)
% regenerate a decomposed signal from basis functions and its coefficients

assert(size(Theta,2) == length(coefs),"Number of basis functions must match number of coefficients")
reg =  Theta*coefs;

end