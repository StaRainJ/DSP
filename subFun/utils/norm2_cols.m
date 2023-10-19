function res = norm2_cols(X)
% function res = norm2_cols(X)
% return norm2 of each column of X 
% ******************************************************************************


    res = sqrt(sum(X.^2, 1));
end
