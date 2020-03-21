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
% ---
% Teresa Vidal Calleja 12/2016
%
%��ֵ
[m,n] = size(xa);
x = reshape(xa,[],1);
%��ͼ��ֵ
[r,s] = size(za);
z = reshape(za,[],1);
R = diag(reshape(Ra,[],1));
%
H = zeros(r*s,m*n);
j=1;
for i=1:length(ind)
    H(j,ind(i)) = 1;
    j=j+1;
end
%xΪ�۲����zΪ�۲�ֵ��HΪ״̬�۲������һ����λ����z-H*xΪ��������
%���ݹ۲�ֵԤ��
[x,Pf] = KF_update_cholesky(x,P,z-H*x,R,H);

xf = reshape(x,m,n);

end