function [varargout]= KF_update_cholesky(varargin)
% Calculate the KF (or EKF) update given the prior state [x,P], the innovation v, the 
% observe uncertainty R, and the (linearised) observation model H. The result is calculated 
% using Cholesky factorisation, which is more numerically stable than a naive implementation.
% �󿨶�������PHt*Si
% Covariance and mean
if (nargin > 3)
    x = varargin{1};
    P = varargin{2};
    v = varargin{3};
    R = varargin{4};
    H = varargin{5};
    cov_only = 0;
% Covariance only
else
    P = varargin{1};
    R = varargin{2};
    H = varargin{3};
    cov_only = 1;
end

PHt = P*H';
s = H*PHt + R;
S = (s+s')*0.5; % ensure S is symmetric 
[Sc,p]  = chol(S);  % note: S = Sc'*Sc,����A�����ǶԳ������ģ��ֽ���ľ���R�������ǵ�
%����cholesky���򻯼���
if ~p
    Sci = inv(Sc);  % note: inv(S) = Sci*Sci'
    Wc = PHt * Sci;
    W  = Wc * Sci';
    if (~cov_only)
        x = x + W*v; % update
    end
    P = P - Wc*Wc';
else
    Si=inv(S);
    if (~cov_only)
        x = x + PHt*Si*v; % update 
    end
    P = P - PHt*Si*H*P;
    % x = x + PHt/S*v; % update 
    % P = P - PHt/S*PHt';
end

if (~cov_only)
    varargout{1} = x;
    varargout{2} = P;
else
    varargout{1} = P;
end

end