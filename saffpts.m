function [L,hL,qL]=saffpts(N, n)
% generates m points on sphere according to Ed Saff's algorithm
% [L, hL, qL] = saffpts(N,Z)
%    N = number of points
%    n = number of zones
% return L which is a (3:n) matrix containing cartesian coordinates 
% of the points, hL is the mesh norm, qL is the separate norm
% ----------------------------------------------------------------
% Author : Q.T Le Gia
% Last modified : May 10, 2002.
if (nargin < 2)
    n = floor(sqrt(pi*N)/2)+1;
end
beta = 4/sqrt(N);
step_theta = (pi-beta)/(n-2);
theta1 = [0:step_theta:(n-3)*step_theta] + beta/2;
theta2 = [step_theta:step_theta:(n-2)*step_theta] + beta/2;
mbar = ones(1,n);
mbar(2:n-1) = N*(cos(theta1)-cos(theta2))/2;
alpha = 0;
for i=1:n
    if (mbar(i)-floor(mbar(i)+alpha)<0.5)
        m(i) = floor(mbar(i)+alpha);
    else
        m(i) = floor(mbar(i)+alpha) + 1;
    end
    alpha = alpha + mbar(i) - m(i);
end    
L(:,1) = [0,0,1]';
cur_index = 2;
q = ones(1,n);
for i=2:n-1
    % generate points on one partition
    t = (i-3/2)*step_theta + beta/2;
    r = sin(t); h = cos(t);
    alpha = mod((i-2),2)/2;
    s = (2*pi*[0:m(i)-1]+alpha)/m(i);
    L(:,cur_index:cur_index+m(i)-1) = [r*cos(s); r*sin(s); h*ones(1,m(i))];
    % compute the separate norm by taking the geodesic distance with previous level
    tmp = L(:,cur_index)'* L(:,cur_index-m(i-1):cur_index-1);
    q(i) = min(acos([tmp L(:,cur_index)'*L(:,cur_index+1)])); 
    % advance one more partition
    cur_index = cur_index+m(i);
end    
L(:,cur_index)=[0 0 -1]';
tmp = L(:,cur_index)'*L(:,cur_index-m(n-1):cur_index-1);
q(n) = min(acos(tmp));
qL = min(q)/2.0;
% -----------------------------------------------------------------------------
% compute the mesh norm
% the following part of code is incorrect, I think 
h = zeros(1,n);
for i=2:n-1
    if (i==2)
        t = step_theta/4 + beta/2;
    else    
        t = (i-3/2-1/2)*step_theta + beta/2;
    end    
    r = sin(t); hz = cos(t);
    alpha = mod((i-2),2)/2;
    % ?????
    %% s = (2*pi*[0:m(i)-1]+alpha)/m(i);
    midpt = [r*cos((alpha+pi)/m(i)); r*sin((alpha+pi)/m(i)); hz];
    htmp = midpt'*L(:,cur_index-m(i-1):cur_index-1);
    h(i) = min(acos([htmp midpt'*L(:,cur_index)]));
end
hL = max(h);
