function al = calculate_fourier(L,eta)
format long;
Lmax = L; 
al=zeros(Lmax+1,1);
fnct=zeros(Lmax+1,1);
tic
for el=0:Lmax
  [al(el+1),ct] = quadgk(@(x)Leg_Wend(0,x,el,eta),-1,1,'RelTol',1e-8,'AbsTol',1e-14);
  al(el+1) = al(el+1)*4*pi/(2*el+1); % which generates \widehat\phi_\ell
  el;
  al(el+1);
  ct;
  fnct(el+1) = ct;
end