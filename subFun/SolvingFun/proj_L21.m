function X = proj_L21(U, opts)
%% solving  X = argmin_X 0.5*||X - U||_F^2 + mu||X||_2,1 
% function X = proj_l1(U, opts)
% Description:
% Inputs: U: double dense matrix d x k 
%         opts: struct 
%             opts.pos: positive constraint (default = false)
% Outputs: X: a full matrix in d x k

    if ~isfield(opts, 'pos')
        opts.pos = false;
    end 
%%   
    mu = opts.mu;
    
    if numel(mu) > 1 && size(mu, 2) == 1 % column vector 
        mu = repmat(mu, 1, size(U, 2));
    end 
    
 %%
 U2 = sum(U.^2,2);
 UU = max(0, U2 - mu); 
 X = UU.*U./(U2);

end