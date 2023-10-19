function opts = Parm_Init1(N,basis_num)
% Parameters setting of  our DSP Solver for basis X_{basis_num} :
%     [X^*, E^*]= argmin { 0.5||M*X-E||_2^2 + lambda*||X||_1 + mu*||E||_1 
%     s.t.   X'*X = I 

switch basis_num
    case 1
        opts.lambda = 0.005*log(4*N).*[1,1,0.5,1,1,0.5,0.5,0.5,0.1]'; 
% opts.lambda = 0.05*log(4*N).*[1,1,.1,1,1,.1,.1,.1,.01]'; 
        opts.mu     = 0.06; % 
        opts.tol = 1e-7;
        opts.max_iter = 2000;
    case 2
        opts.lambda = 0.005.*log(4*N)*[1,1,0.5,1,1,0.5,0.5,0.5,0.1]'; 
%         opts.lambda = 0.05*log(4*N).*[1,1,.1,1,1,.1,.1,.1,.01]'; 
        opts.mu  = 0.03;  
        opts.tol = 1e-6;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
        opts.max_iter = 1000;  
end

opts.ifplotLoss = 0; % judge if showing loss curve
opts.ifplotMtX = 0;  % judge if showing M*X 
