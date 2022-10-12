function phi=Leg_Wend(m,x,el,eta)
% compute the product of a P_el(t) phi(t)/ norm^2(P_el)
% with norm^2(P_el) = 2/(2el+1) 
% where P_el is the Legendre function and
% phi_m(t) is the Wendland based spherical basis function
if (m == 0)
  y =Wend_sbf0(x,eta);
elseif (m == 2)
  y =Wend_sbf2(x,eta);
elseif (m == 4)
  y =Wend_sbf4(x);
elseif (m == 6)
  y =Wend_sbf6(x);
end
ps  =legendre(el,x);
pl  = (el+0.5)*ps(1,:); % normalize
phi =y.*pl;  % phi = SBF(x)*P_n(x) 
