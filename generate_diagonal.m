function [Matrix] = generate_diagonal(vector)
% Input vector with (L+1)*1
% Output Matrix with diag(vector)
L = length(vector)-1;
vectorlong(1,1) = vector(1);
for k=2:L+1
    vectoradd = vector(k)*ones(2*k-1,1);
    vectorlong = [vectorlong;vectoradd];
end
Matrix = diag(vectorlong);