function X = halton_pts(N)

b1 = 3;
b2 = 5;

z = 2*vdcorput(N,b2)-1;
varphi = 2*pi*vdcorput(N,b1)-1;

for k=1:N
  point = [sqrt(1- z(k)^2)*cos(varphi(k)); 
           sqrt(1- z(k)^2)*sin(varphi(k));
           z(k)];
  X(:,k) = point;
end
	
