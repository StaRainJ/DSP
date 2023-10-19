function [C0] = DataInitialization(X,Y)

% % Parameters Setting
gamma = 2;          % tradeoff in distance calculation : d=1+gama*exp(-min(d1,d2)); small gamma for ratation and scaling
Bound = [3,30];      % constraint the range of MinPts
pct   = 0.025;        % the  rate to choose MinPts: 0.02-0.05
Mu    = [0.15,0.08];  % the  rate to choose epsilon:

N      = size(X,1);
NoiseN = round(max(min(0.05.*N,5),10));
NoiseX = [X;rand(NoiseN,2).*max(X)];
NoiseY = [Y;rand(NoiseN,2).*max(Y)];
X = NoiseX; Y = NoiseY;

X1 = [X,Y,Y-X];
D  = myDist(X1,gamma);% 20  h=2;


X1=mapminmax(X1',0,1)'; 
% D = pdist2(X1,X1);
D = pdist2(X1(:,1:2),X1(:,1:2))+pdist2(X1(:,3:4),X1(:,3:4))+pdist2(X1(:,5:6),X1(:,5:6));

% ******* iteration1 *****
C0 = 1:size(X,1);% CorrectIndex;%
mu = Mu(1);
[D1,epsilon,MinPts,CorePoints] = GetParamater(D,C0,pct,mu,Bound,N);
%%
% [IDX1,~] = DBSCAN(D1,MinPts,epsilon);  
% IDX1=IDX1(1:N);

C0 = CorePoints;%find( IDX1>0 );  %  CorrectIndex;% 