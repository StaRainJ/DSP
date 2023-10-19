function res = normF2(A)
    % compatible with multi-channel
	%% ================== File info ==========================
	%% ================== end File info ==========================
    B = A.^2;
    B = B(:);
    res = sum(B);
	% res = norm(A, 'fro')^2;
end 