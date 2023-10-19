function [X, iter,cost,E]= DSP_Optim1(M, Xinit,Einit,opts)
% The Dual Sparsity Solver for Problem:
%     [X^*, E^*]= argmin { 0.5||M*X-E||_2^2 + lambda*||X||_1 + mu*||E||_1 
%     s.t.   X'*X = 1 
%
% where:
%     M : Embedded data matrix with shape [n_samples, n_dim] 
%     X : optimization  basis  
% 
%     ||.||_1 : L1 norm for any vector A is defined by
% 
%         ||A||_1 = \\sum_i |A_i|
% 
% We solve the unknown model reasoning and fitting problem with our DSP method, which is proposed in
% our NeurIPS 2023 paper (Spotlight):
%
%     Jiang Xingyu, Ma Jiayi. Robust Model Reasoning and Fitting via Dual Sparsity Pursuit. NeurIPS 2023 
% 
% Please refer to the paper for details, and kindly cite our work if you find it is useful.
% 
% PARAMETERS
% ----------
% 
%     See Parm.Init1.m for details
% 
% OUTPUT
% ------
% 
%     X :  Single basis of length {n_dim}
% 
%     E :  Error vector of length {n_samples}
%
%     cost :  final objective value
% 
%     iter :  number of used iterations 
% 
%% initionalization
    if ~isfield(opts, 'max_iter')
        opts.max_iter = 500;
    end
    if ~isfield(opts, 'regul')
        opts.regul = 'l1';
    end     
    if ~isfield(opts, 'pos')
        opts.pos = false;
    end
    
    if ~isfield(opts, 'tol')
        opts.tol = 1e-6;
    end
    
    if ~isfield(opts, 'verbose')
        opts.verbose = false;
    end
    
    if ~isfield(opts, 'ifplotLoss')
        opts.ifplotLoss = 0;
    end
   
lambda = opts.lambda; % for basis term
mu = opts.mu; % for error term
N = size(M,1);
%% Construc cost function
    %% cost f
    function cost = calc_f(X,E)
        G = M*X-E;
        cost = 1/2 *normF2(G);
    end 
    %% cost function 
    function cost = calc_F(X,E,lambda,mu)
            cost = calc_f(X,E) + norm1(lambda.*X)+mu*norm1(E);
    end 
%% gradient: f'_x = M'*(M*X-E)
    DtM = M'*M; % dxd 
    Mt = M';
    function res = grad(X,E) 
        DtE = Mt*E;
        res = DtM*X -DtE;
    end

    function res1 = U_E(X)
        res1 = M*X;
    end
%% L-Lipschitz constant 
    EE = eig(DtM);
    L = EE(end);
%% sparse optimizing    
    [X,iter,cost,E] = DSP_Solving1(@grad, @proj_l1,@U_E, Xinit,Einit, L, opts, @calc_F);    
    
end