function Y = spherical_harmonics(X,L)
% Generate the basis of spherical harmonics Y_l^m of S^2
% with -l<=m<=l
zx = X(:,1);
zy = X(:,2);
zz = X(:,3);
[theta1,phi1,r1] = cart2sph(zx,zy,zz);
theta = pi/2-phi1;
phi = theta1;

    Y = [];
    Ybelow =[];
    Yupper = [];
for l=1:L+1
    legendre_lm = legendre(l-1,zz);
    Kl0 = sqrt((2*(l-1)+1)/(4*pi));
    Yzero(1,:) = Kl0*legendre_lm(1,:);
    for m=1:l-1
        Klm= Kl0*sqrt(factorial(l-1-m)/factorial(l-1+m)); 
        Yupper(m,:) =sqrt(2)*Klm*sin(m*phi)'.*legendre_lm(m+1,:);
        Ybelow(m,:) = sqrt(2)*Klm*cos(m*phi)'.*legendre_lm(m+1,:);
    end
   
    Y = [Y;Ybelow(end:-1:1,:);Yzero;Yupper];
     clear Yupper;clear Ybelow
end

