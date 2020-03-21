function [x,P]= KF_update_cholesky(x,P,v,R,H)
% Calculate the KF (or EKF) update given the prior state [x,P], the innovation v, the 
% observe uncertainty R, and the (linearised) observation model H. The result is calculated 
% using Cholesky factorisation, which is more numerically stable than a naive implementation.
% 求卡尔曼参数PHt*Si
PHt = P*H';
s = H*PHt + R;
S = (s+s')*0.5; % ensure S is symmetric 

Si=inv(S);
%更新x，卡尔曼滤波估计值
x = x + PHt*Si*v; % update 
P = P - PHt*Si*PHt';
end