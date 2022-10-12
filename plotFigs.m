clear all;
close all;
% load results

[Z,hz,qz] = saffpts(10000);  Z = Z';
NN = length(Z);

%load Result_1;
%load Result_2;
%load Result_3;
%load Result_4;
%load Result_5;
load Result_6;

level = size(f_appro_s,2);

% plot exact function
[K,fh] = plotfunc(f_exact,Z','Exact function',1);
set(gca,'CLim',[min(f_exact),max(f_exact)]); colormap(jet);

f_appro = zeros(NN,1);
err = zeros(NN,1);

for  k=1:level
    f_appro = f_appro_s(:,k);
    err = err_s(:,k);
    figure;
    [K,fh] = plotfunc(f_appro,Z','Reconstruction',2*(k-1)+2);
    colormap(jet); colorbar; set(gca,'CLim',[min(f_exact),max(f_exact)]);  
    drawnow;
    
    figure; 
    [K,fh] = plotfunc(err,Z','Error',2*(k-1)+3);
    colormap(jet); colorbar;set(gca,'CLim',[-0.2,0.2]); 
end

