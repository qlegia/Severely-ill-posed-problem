function Phi = generate_Phi(X,Y,eta)
N = length(X);
M = length(Y);
for j=1:M
    test = Y(j,:);
    for i=1:N
        distance(i) = norm(X(i,:)-test,2);
        Phi(j,i) = max(1-distance(i)/eta,0)^2;%*(4*distance(i)/eta+1);
    end
end

Phi = 1/eta^2*Phi;
