clear all;
tic;
close all;
L=[15;18;21;24]+10;
Lend = 50;
nv = 2.51;
rho=1.2;  % rho = R in the paper
%n = [32, 125, 500, 2000]; %% ??? [32, 250, 400, 2000]
%n = [32, 250, 400, 2000];
%h = [0.4371;0.2288;0.1327;0.0613];
n = [32, 128, 512, 2048];
h = [0.4404; 0.2212; 0.1108; 0.0553];
%L=34;
%n=2000;
%h=0.0613;
level = length(n);
x = zeros(n(end),level);y=x;z=x;
f_observe = zeros(n(end),level);

for el=1:Lend+1
    a(el,1) = rho^(-(el-1));                 % symbols of the operator
    c(el,1) = a(el,1)*((el-1)+0.5)^(-nv);    % c_ell*g_{ell,k} = Fourier coeffs of RHS 
    d(el,1) = ((el-1)+0.5)^(-nv);            % =hat{u}_{ell,k} / g_{ell,k}
    
end 
g = load('g.txt');      % g_{ell,k} for ell=0,...,Lend=50
%g = load('g1.txt');    % g_{ell,k} for ell=0,...,Lend=100
D = generate_diagonal(d);
C = generate_diagonal(c);
%[Z,lat,lon] = saffsph_deg3(10000);  Z = Z';
[Z,hz,qz] = saffpts(10000);  Z = Z';
NN = length(Z);
zx = Z(:,1);
zy = Z(:,2);
zz = Z(:,3);
Y = spherical_harmonics(Z,Lend);
f_exact = Y'*D*g; % which is denoted by u(x) in section 6

%figure;
%scatter3(zx,zy,zz,15,f_exact,'filled'); 
%[K,fh] = plotfunc(f_exact,Z','Exact function');
%set(gca,'CLim',[min(f_exact),max(f_exact)]); colormap(jet);

f_appro = zeros(10000,1);
% saving the results, for plotting later
f_appro_s = zeros(10000,level);
err_s = zeros(10000,level);

delta = 0.01; % noise level
noise = load('noise.txt');

% no-noise
noise = zeros(NN,level);

for k=1:level
    %[X,lat,lon] = saffsph_deg3(n(k)); X=X';%h(k,1) = msh_norm(X');
    [X,hx,qx] = saffpts(n(k)); X = X';
    x(1:n(k),k) = X(:,1);
    y(1:n(k),k) = X(:,2);
    z(1:n(k),k) = X(:,3);
    Y = spherical_harmonics(X,Lend);
    YY= spherical_harmonics(X,L(k));  % choosing the set Y = set X ??
    f_observe(1:n(k),k) = Y'*C*g/rho+noise(1:n(k),k);
    %eta(k,1) = 0.5*h(k)/(a(L(k)+1,1))
    eta(k,1) = 1/5*h(k)/(a(L(k)+1,1))^(1/1.5);
    %eta=0.2;
    residual = f_observe(1:n(k),k);
    %epsilon(k,1) =0.1*a(L(k)+1,1)+delta; % -> Result_1.mat
    %epsilon(k,1) =0.01*a(L(k)+1,1)+delta; % -> Result_2.mat
    %epsilon(k,1) = 0*a(L(k)+1,1)+delta; % -> Result_3.mat
    %epsilon(k,1) = 0.7*delta; % -> Result_4.mat
    %epsilon(k,1) = 0.1*delta; % -> Result_5.mat
    epsilon(k,1) = 0.1*a(L(k)+1,1); % -> Result_6.mat
    lambda(k,1) =0.01^2*a(L(k)+1,1)^2*0.5^3;
    al = calculate_fourier(L(k),eta(k,1));
    F = generate_diagonal(al);
    A = generate_diagonal(a(1:L(k)+1,1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X=[x(1:n(k),k),y(1:n(k),k),z(1:n(k),k)];
    Phi = generate_Phi(X,X,eta(k,1));
    NN = YY'*A*F*YY/rho;
    for m=1:k-1
        X_previous = [x(1:n(m),m),y(1:n(m),m),z(1:n(m),m)];
        Y_previous = spherical_harmonics(X_previous,L(k));
        al_previous = calculate_fourier(L(k),eta(m,1));
        F_previous = generate_diagonal(al_previous);
        residual = residual-YY'*A*F_previous*Y_previous*alpha(1:n(m),m)/rho;
    end
    toc;
    CC = [NN', -eye(n(k)), zeros(n(k),n(k)); -NN', zeros(n(k),n(k)), -eye(n(k)) ];
    d = [residual+epsilon(k)*ones(n(k),1); -residual+epsilon(k)*ones(n(k),1)];
    T = [lambda(k)*Phi, zeros(n(k),2*n(k)); zeros(2*n(k),n(k)), eye(2*n(k),2*n(k))];
    Aeq=[]; beq=[]; lb = [-Inf*ones(n(k),1);zeros(2*n(k),1)]; ub = Inf*ones(3*n(k),1); %x0=zeros(3*n,1);
    options = optimset('MaxIter',50000,'LargeScale','off','TolFun',1e-14,'Algorithm','interior-point-convex');
    [result,output,exitflag] = quadprog(2*T,[],CC,d,Aeq,beq,lb,ub,[],options);
    toc;
    %[result,output,exitflag] = linprog(c,C,d,Aeq,beq,lb,ub,[]);
    alpha(1 : n(k), k) = result(1 : n(k), 1);
    Phi1 = generate_Phi(X,Z,eta(k,1));
    f_appro = f_appro + Phi1 * alpha(1 :n(k), k);
    err  = f_appro - f_exact;
    f_appro_s(:,k) = f_appro;
    err_s(:,k) = err; 
    %err_norm (k, :)= [norm(err,inf)/norm(f_exact,inf),sqrt(sum(err.^2)*4*pi/NN)/sqrt(sum(f_exact.^2)*4*pi/NN) ]
    err_norm (k, :)= [norm(err,inf)/norm(f_exact,inf),norm(err,2)/norm(f_exact,2) ]
    % (k, :)= sqrt(sum(err.^2)*4*pi/NN)/sqrt(sum(f_exact.^2)*4*pi/NN)
    %figure;
    %scatter3(zx,zy,zz,15,f_exact,'filled');
    
    %figure;
    %scatter3(zx,zy,zz,15,f_appro,'filled');
    %[K,fh] = plotfunc(f_appro,Z','Reconstruction');
    %colormap(jet); colorbar; set(gca,'CLim',[min(f_exact),max(f_exact)]);  
    %drawnow;
    
    %figure; 
    %scatter3(zx,zy,zz,15,err,'filled'); 
    %[K,fh] = plotfunc(err,Z','Error');
    %colormap(jet); colorbar;set(gca,'CLim',[-0.2,0.2]); 

    %drawnow;
end
%figure;
%scatter3(zx,zy,zz,15,f_exact,'filled');
err4 = max(max(abs(f_exact-f_appro)))
toc; 
%save Result_1 f_exact f_appro f_appro_s err_s f_observe epsilon lambda eta err_norm err4 a c d n h nv L noise;
%save Result_2 f_exact f_appro f_appro_s err_s f_observe epsilon lambda eta err_norm err4 a c d n h nv L noise;
%save Result_3 f_exact f_appro f_appro_s err_s f_observe epsilon lambda eta err_norm err4 a c d n h nv L noise;
%save Result_4 f_exact f_appro f_appro_s err_s f_observe epsilon lambda eta err_norm err4 a c d n h nv L noise;
%save Result_5 f_exact f_appro f_appro_s err_s f_observe epsilon lambda eta err_norm err4 a c d n h nv L noise;
save Result_6 f_exact f_appro f_appro_s err_s f_observe epsilon lambda eta err_norm err4 a c d n h nv L noise;
