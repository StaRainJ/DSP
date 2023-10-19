function [g_gric] = compute_gric(sqr_res,sigma,lambda1, lambda2, model,heta)
%COMPUTE_GRIC
if(nargin<6)
    heta = 1;
end
n = numel(sqr_res);
if(nargin ==1)
    if(n<10)
        lambda1 = 2;
        lambda2 = 4;
    else
        lambda1 = 1;
        lambda2 = 5;
    end
    sigma = std(sqr_res);
end



switch lower(model)
    case {'general','fundamental'}
        
        k = 7; % number of parameters
        d = 3; % dimension of the manifold
        r = 4; % dimension of the data (pairs of image points)
    case 'affine'
        
        k = 4; % number of parameters
        d = 3; % dimension of the manifold
        r = 4; % dimension of the data (pairs of image points)
        
    case 'homography'
        
        k = 8; % number of parameters
        d = 2; % dimension of the manifold
        r = 4; % dimension of the data (pairs of image points)
        
    case 'affinity'
        
        k = 6; % number of parameters
        d = 2; % dimension of the manifold
        r = 4; % dimension of the data (pairs of image points)
    case 'circle'
        k = 3; % number of parameters
        d = 1; % dimension of the manifold
        r = 2;
    case 'line'
        k = 2; % number of parameters
        d = 1;
        r = 2;
    case 'plane'
        k = 4; % number of parameters
        d = 2;
        r = 3;
    case 'cylinder'
        k = 7; % number of parameters
        d = 2;
        r = 3;
    case 'parabola'
        k = 3; % number of parameters
        d = 1;
        r = 2;
    otherwise
        disp('Unknown model.')
end

g_gric = sum(min(sqr_res./sigma^2, 2*(r-d))) + lambda1*n*d + lambda2*k*heta;
%g_gric = sum(sqr_res./sigma^2) + lambda1*n*d + lambda2*k;
%g_gric = sum(min(sqr_res./sigma^2, 2*(r-d))) + lambda2*d;

end



