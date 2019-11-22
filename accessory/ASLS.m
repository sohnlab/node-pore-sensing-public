function z = ASLS(y, lambda, p, maxIter)
% applies asymmetric least squares to estimate baseline, z.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this code is taken from Eilers, Boelens 'Baseline Correction with 
% Asymmetric Least Squares Smoothing' 2005
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% .001 < p < .1
% 10e2 < lambda < 10^9

m = length(y);
D = diff(speye(m), 2);
w = ones(m,1);

for it = 1:maxIter
    W = spdiags(w, 0, m, m);
    C = chol(W + lambda * (D.' * D));
    z = C \ (C' \ (w .* y));
       
    w = p * (y < z) + (1 - p) * (y > z);
end

end
