function [x,P]= KF_update_cholesky(x,P,v,R,H)
% Calculate the KF (or EKF) update given the prior state [x,P], the innovation v, the 
% observe uncertainty R, and the (linearised) observation model H. The result is calculated 
% using Cholesky factorisation, which is more numerically stable than a naive implementation.
% 求卡尔曼参数PHt*Si
PHt = P*H';
s = H*PHt + R;
S = (s+s')*0.5; % ensure S is symmetric 
[Sc,p]  = chol(S);  % note: S = Sc'*Sc,矩阵A必须是对称正定的，分解出的矩阵R是上三角的
%引入cholesky，简化计算
if ~p
    Sci = inv(Sc);  % note: inv(S) = Sci*Sci'
    Wc = PHt * Sci;
    W  = Wc * Sci';

    x = x + W*v; % update 
    P = P - Wc*Wc';
else
    Si=inv(S);
    x = x + PHt*Si*v; % update 
    P = P - PHt*Si*H*P;
    % x = x + PHt/S*v; % update 
    % P = P - PHt/S*PHt';
end
end