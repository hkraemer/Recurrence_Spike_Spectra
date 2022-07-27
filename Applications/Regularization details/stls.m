function coeffs = stls(s, Theta, lambda)
% Sequential Thresholded Least Squares sparse regression method

max_iterations = 10;
tol = 1e-5;

coeffs = pinv(Theta)*s;     % initial guess least squares
biginds_old = coeffs;

k = 1;
while k < max_iterations
    smallinds(:) = (abs(coeffs)<lambda);  % find small coefficients 
    coeffs(smallinds) = 0;  % threshold these coeffs
    biginds(:) = ~smallinds;
    % check convergence
    if max(abs(biginds - biginds_old')) < tol
        break
    end
    % regress onto remaining terms
    coeffs(biginds) = pinv(Theta(:,biginds))*s;
    biginds_old(:) = biginds;
    k = k + 1;
end
end