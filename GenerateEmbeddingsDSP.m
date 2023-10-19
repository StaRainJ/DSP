function  [D,X,Y] = GenerateEmbeddingsDSP(X,Y,GenerateType)

N = size(X,1);

% [xs1, ~] = normalise2dpts(X);
% [xs2, ~] = normalise2dpts(Y);
% Xt = xs1; Yt= xs2;
[X, Y, Xn, Yn, Xt, Yt] =  DataPreprocessing(X, Y);
switch GenerateType
    case 1
%% general  fundamental embedding F
D1 = zeros(9,N);
for nn = 1:N
    D1(:,nn) = kron(Yt(:,nn),Xt(:,nn)); 
end
nml = sqrt(sum(D1.^2)); % D/norm2£¨D£©
D1 = D1*diag(1./nml);
D = D1';
%% homography  H 
    case 2
D2 = [];
ooo  = zeros(1,3);
for k=1:N
  p1 = Xt(:,k);
  p2 = Yt(:,k);
  D2 = [ D2;
    p1'*p2(3) ooo -p1'*p2(1)
    ooo p1'*p2(3) -p1'*p2(2)
   ];
end
D2 = D2';
nml = sqrt(sum(D2.^2));
D2 = D2*diag(1./nml);
D = D2';
 %% affine  H_A
    case 3
D3 = [];
ooo  = zeros(1,3);
for k=1:N
  p1 = Xt(:,k);
  p2 = Yt(:,k); 
  D3 = [ D3;
    p1' ooo -p2(1)
    ooo p1' -p2(2)
   ];
end
D3 = D3';
nml = sqrt(sum(D3.^2));
D3 = D3*diag(1./nml);

D = D3';
%% for orthographic camera E_o or affine camera  F_A
    case 4
D4 = [Xn,Yn,ones(N,1)]';
nml = sqrt(sum(D4.^2));
D4 = D4*diag(1./nml);
D = D4';
    otherwise
        disp('Error: please set valid embedding type, using 1 to 4')
end