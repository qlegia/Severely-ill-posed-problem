function phi=cmp_sbf2(t,eta)
% the local supported function on the sphere, based on Wendland function
% of smoothness of order 2 in R^3
%      (1-r)^4_+ (4r+1) where r = sqrt(2-2*t)
% this kernel generates H^{5/2}(S^2)
% ----------------------------------------------------------------
r = real(sqrt(2.-2.*t));
mask = (abs(r)/eta < 1);
phi = mask .*(1-r/eta).^2.*(4.*r/eta+1);
phi = phi/eta^2;    
  
     
