function [xf,Pf] = fuse_measurements(xa,P,za,Ra, ind)
% Obtain maximum aposteriori estimate using Bayesian Fusion
% Inputs:
% x: prior mean
% P: prior covariance
% z: measurements
% R: measurement variances
% ind: measurement indices in grid map
% ---
% Outputs:
% xf: map with fused measurements
% Pf: covariance with fused measurements
%均值
[m,n] = size(xa);
x = reshape(xa,[],1);
%子图均值
[r,s] = size(za);
z = reshape(za,[],1);
R = diag(reshape(Ra,[],1));
%
H = zeros(r*s,m*n);
for i=1:length(ind)
    H(i,ind(i)) = 1;
end
%x为观测矩阵，z为观测值，H为状态观测矩阵，是一个单位矩阵，z-H*x为测量余量
%根据观测值预测
[x,Pf] = KF_update_cholesky(x,P,z-H*x,R,H);

xf = reshape(x,m,n);
% Debugging
disp(['Matrix trace before update: ', num2str(trace(P))])
disp(['Matrix trace after update: ', num2str(trace(Pf))])

end