function phi0=Wend_sbf0(t,eta)
% the spherical basis function based on Wendland function
% restricted to the sphere
% phi(t) = (1-r)^2_+  with r = sqrt(2-2t)
% this kernel generates H^{3/2}(S^2)
r = real(sqrt(2.0-2.0*t));
mask = (abs(r)/eta < 1.0);
phi0 = mask.*(1.0-r/eta).^2/eta^2;
