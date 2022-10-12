function phi=Wend_sbf4(t)
% this kernel generates H^{7/2} (S^2)
r = real(sqrt(2.0-2*t));
mask = (abs(r)<1.0);
phi = mask .* (1.0-r).^6.*(35.0*r.^2+18.0*r+3.0);
