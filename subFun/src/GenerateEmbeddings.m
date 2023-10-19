function [D1, D2, D3, D4] = GenerateEmbeddings(Xn, Yn, Xt, Yt, N)
D2 = []; D3 = [];D4 = [];
%% general  fundamental embedding F
D1 = zeros(9,N);
for nn = 1:N
    D1(:,nn) = kron(Yt(:,nn),Xt(:,nn)); % 基础矩阵的表达形式
end
nml = sqrt(sum(D1.^2)); % D/norm2（D）
D1 = D1*diag(1./nml);

%% homography  H 
% D2 = [];
% ooo  = zeros(1,3);
% for k=1:N
%   p1 = Xt(:,k);
%   p2 = Yt(:,k);
%   D2 = [ D2;
%     p1'*p2(3) ooo -p1'*p2(1)
%     ooo p1'*p2(3) -p1'*p2(2)
%    ];
% end
% D2 = D2';
% nml = sqrt(sum(D2.^2));
% D2 = D2*diag(1./nml);
% 
 %% affine  H_A
% D3 = [];
% ooo  = zeros(1,3);
% for k=1:N
%   p1 = Xt(:,k);
%   p2 = Yt(:,k); 
%   D3 = [ D3;
%     p1' ooo -p2(1)
%     ooo p1' -p2(2)
%    ];
% end
% D3 = D3';
% nml = sqrt(sum(D3.^2));
% D3 = D3*diag(1./nml);
% 

%% for orthographic camera E_o or affine camera  F_A
% D4 = [Xn,Yn,ones(N,1)]';
% nml = sqrt(sum(D4.^2));
% D4 = D4*diag(1./nml);