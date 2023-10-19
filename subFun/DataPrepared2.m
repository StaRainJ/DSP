function [D] =  DataPrepared2(X, Y)

N = size(X,2);

[xs1, ~] = normalise2dpts(X);
[xs2, ~] = normalise2dpts(Y);
Xt = xs1; Yt= xs2;
%%
D1 = zeros(9,N);
for nn = 1:N
    D1(:,nn) = kron(Yt(:,nn),Xt(:,nn)); % 
end
nml = sqrt(sum(D1.^2)); % D/norm2£¨D£©
D1 = D1*diag(1./nml);
%%
D = D1' ;