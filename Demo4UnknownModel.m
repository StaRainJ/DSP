
%% ************* Demo for Unknown model fitting via DSP **************************
% The Dual Sparsity Solver for model reasoning (unknown model type):
%     [X^*, E^*, r^*]= argmin { 0.5||M*X-E||_2^2 + lambda*||X||_1 + mu*||E||_1 
%     s.t.   X'*X = I
%
% ------
% We solve the robust model reasoning and fitting problem with our DSP method, which is proposed in
% our NeurIPS 2023 paper (Spotlight):
%
%     Jiang Xingyu, Ma Jiayi. Robust Model Reasoning and Fitting via Dual Sparsity Pursuit. NeurIPS 2023 
% 
% Please refer to the paper for details, and kindly cite our work if you find it is useful.
% ------

close all;
clear;
clc;

addpath(genpath('data'));
addpath(genpath('subFun'));


ifplotMatches = 0; %  set 1 to see inlier detection


DataName = 'EVD'; % 
Mgt = 2;

% DataName = 'kusvod2GT'; % 
% Mgt = 1;



ModelType = 'unKnown'; %   

DataDir = ['./data/',DataName,'/'];
str = dir([DataDir,'*.mat']);


for iii =1:length(str) 
close all
clear CorrectIndex inliers   X1 X2

%% read data
load([DataDir,str(iii).name]); 

X0 = X(:,1:2); Y0 =Y(:,1:2); inliers = CorrectIndex;

N0 = size(X0,1);

% data prepare for better visibility
% [X0,Y0,inliers,I1,I2] = DataPrepared(X0,Y0,inliers,I1,I2);

tic
%% data Embedding for optimization

[D0,X,Y] = GenerateEmbeddingsDSP(X0,Y0,1);  % X0--Nx2, X--3xN

%% Pre-filtering outliers
IfPreProcess = 2; %%
[C0] = PreProcedure(X0,Y0,IfPreProcess); 
N = length(C0);
D1 = D0(C0,:);

%%  Solving Basis1
disp('******* Solving Basis1 ************')

[v,~]= eig(D1'*D1);
Xinit = v(:,1);
Einit = zeros(N,1);

opts = Parm_Init1(N,1);
[X1,MaxItr1,~,E1] = DSP_Optim1(D1, Xinit,Einit,opts);

% Finding inliers from basis X1; 
DX1    = D0*X1;    

% ====== adaptive threshold setting =========== 
% jj = 1; aa = [];
% opts_proj2.lambda = 0.02;
% for tt = 0.05:-0.002:0.01
%    aa(jj) = sum( abs(DX1) < tt);
% 
%    if aa(jj)./aa(1)<=0.95
%        opts_proj2.lambda = tt+0.002;
%        break;
%    end
%       
%           jj = jj+1;
% end

% ======================================= 

opts_proj2.lambda = 0.02;
E11 = feval(@proj_l1, DX1, opts_proj2);
InlierIndex1 = find(abs(E11) <= 0.00001);


if ifplotMatches >=2
  figure(11),plot_matches(I1, I2, X(1:2,:)', Y(1:2,:)', InlierIndex1, inliers);
end

%% solving Basis2
disp('******* Solving Basis2 ************')

D2 = D0(InlierIndex1,:);
% D2 = D0;
N2 = size(D2,1);
% [v,~]= eig(D2'*D2);
Yinit = v(:,2); % ones(6,1);

Yinit = OrthProj(X1,Yinit);
Yinit = SphereProj(Yinit);
Einit = zeros(N2,1);

opts2 = Parm_Init1(N2,2);
opts2.mu = opts_proj2.lambda;

[X2,MaxItr2,~,E2] = DSP_Optim2(D2, Yinit,Einit,opts2,X1);
InlierIndex2 = InlierIndex1((abs(E2) <= 0.0001));

if length(InlierIndex2) < 20
    InlierIndex2 = InlierIndex1;
end

if ifplotMatches >= 2
   figure(22),plot_matches(I1, I2, X(1:2,:)', Y(1:2,:)', InlierIndex2, inliers); %% basis 2 results
end

%% ----- compute L(M,x,e) = 0.5||M'x -e||^2 + lambda ||x||_1 + gamma ||e||_1
tt = 0.04;
Lambda1 =  opts.lambda;
SE1 = norm1(truncLoss(E11,tt));

DX2    = D0*X2;    
Lambda2 =  opts2.lambda;
E22 = feval(@proj_l1, DX2, opts_proj2);
SE2 = norm1( truncLoss(E22,tt) )+ norm1( truncLoss( E11(logical(E11)&~logical(E22)),tt ));


L1 = 0.5*normF2(DX1-E11) + norm1(Lambda1.*X1) + 0.02 * SE1;
L2 = 0.5*normF2(DX2-E22) + norm1(Lambda1.*X2) + 0.02 * SE2;


ModelRatio = L2./L1;


disp('******* Succesfully solving unKnown Model Fitting via DSP ! ************')

%% Model reasoning 
if strcmp(ModelType,'unKnown')
    tao = 1.2;
    if ModelRatio < tao
      MT = 'Homo';
      Mest = 2;
    else 
      MT = 'Fund';
      Mest = 1;
    end
end

disp(['Identified Model Type: ' MT]) 
disp(['Model Rato: ' num2str([ModelRatio])])

%% PostProcessing

 thF= 1; thH = 2; 

if strcmp(MT,'Fund')
       [VFCIndexF, F_est, d] = Post_F2(X, Y, 100, thF, InlierIndex1); 
        e1 = mean(SampsonDistanceF(X(1:3,inliers),Y(1:3,inliers), F_est));
        ee = e1;
        FinalIndex =  VFCIndexF;
%         FinalIndex = find(d<2);
elseif strcmp(MT,'Homo')  
        [VFCIndexH, H_est, d] = Post_H2(X, Y, 100, thH, InlierIndex2);
        e2 = mean(SampsonDistanceH(X(1:3,inliers),Y(1:3,inliers),H_est));
        ee = e2;
        FinalIndex =  VFCIndexH;
%         FinalIndex = find(d<3);
end

Time =toc;

%% show results

if ifplotMatches >= 1
    % plot original data
    figure(100), plot_matches(I1, I2, X(1:2,:)', Y(1:2,:)', 1:N0, inliers);  
    % plot final results
    figure(111), plot_matches(I1, I2, X(1:2,:)', Y(1:2,:)', FinalIndex , inliers);  
end

GE(iii) = ee;
RT(iii) = Time;
MoTy{iii} = MT;
M_judge(iii) = (Mest == Mgt);
MR(iii,:) = [ModelRatio];
end
clc
disp(['******* Final Overall Results of ',DataName,' ************'])
disp(['The Median GE is ',num2str(median(GE))]);
disp(['The Mean Run Time is ',num2str(mean(RT)),' s']);
disp(['The Failure Ratio is ',num2str(sum(GE>5|~M_judge)./length(GE))]);

