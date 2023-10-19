function X = proj_l12(U, opts)
% -----------------------------------------------
	if nargin == 0 
		d = 1000;
		n = 1000;
		U = normc(rand(n, d))';        
		lambda = 1;
        opts.lambda = lambda;
	end 
 	X = proj_l2(U', opts)';
    if nargin == 0 
        imagesc(X);
        X = [];
    end 
end
	