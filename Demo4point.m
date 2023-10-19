%==================================================================
% The Dual Sparsity Solver for model reasoning and fitting from 2-D point data:
%     [X^*, E^*]= argmin { 0.5||M*x-e||_2^2 + lambda*||x||_1 + mu*||e||_1 
%     s.t.   x'*x = 1
% Note : Test potential models are 'Line'  'Parabola' 'Ellipse'...
%
% ------
% We solve the robust model reasoning and fitting problem with our DSP method, which is proposed in
% our NeurIPS 2023 paper (Spotlight):
%
%     Jiang Xingyu, Ma Jiayi. Robust Model Reasoning and Fitting via Dual Sparsity Pursuit. NeurIPS 2023 
% 
% Please refer to the paper for details, and kindly cite our work if you find it is useful.
%==================================================================

close all;
clear;

addpath(genpath('data'));
addpath(genpath('subFun'));

Num_in = 300;
Num_out = 300;
L = Num_in + Num_out;
ModelType = 'Ellipse'; % 'Line'  'Parabola' 'Ellipse'
%% para setting
opts.lambda = log(2*L).*[0.01,0.1,0.1,1,1]'; 
opts.tol = 1e-7;
opts.mu     = 0.03;
opts.max_iter = 2000;

disp('******************* Generating Data *************************')
%% Line
if strcmp(ModelType,'Line')
xa = 0.2; ya= 0.1;
xb = 0.8; yb =0.9;
sigma = 0.01;
[x, y,gp] = generateLine(xa,ya,xb,yb,sigma,Num_in, Num_out);

xx = [x,y];
xx = xx';

L = size(x,1);
D = [ones(L,1),x, y, x.^2, y.^2];
nml = sqrt(sum(D'.^2)); % D/norm2£¨D£©
D = (D'*diag(1./nml))'; 
Din = D(gp>0,1:3);
[u0,s0,v0]= svd(Din);
W0 = v0(:,end);
W0 = [W0;0;0];
W0 = W0 ./ sqrt(W0'*W0);
%%   Parabola
elseif strcmp(ModelType,'Parabola')
[x, y,gp] = generateParabola( Num_in, Num_out);

L = size(x,1);
D = [ones(L,1),x, y, x.^2, y.^2];
nml = sqrt(sum(D'.^2)); %  D/norm2£¨D£©
D = (D'*diag(1./nml))'; 
Din = D(gp>0,1:4);
[u0,s0,v0]= svd(Din);
W0 = v0(:,end);
W0 = [W0;0];
W0 = W0 ./ sqrt(W0'*W0);
figure,plot(D*W0);

%% Ellipse 
elseif strcmp(ModelType,'Ellipse')
num = 1;
radius = 1;
epsilon = 0.01;
[x,y,gp] = generate_Circles( num, Num_in, Num_out, radius, epsilon);
%
L = size(x,1);
D = [ones(L,1),x, y, x.^2, y.^2];
nml = sqrt(sum(D'.^2)); % D/norm2£¨D£©
D = (D'*diag(1./nml))'; 
Din = D(gp>0,:);
[u0,s0,v0]= svd(Din);
W0 = v0(:,end);
W0 = W0 ./ sqrt(W0'*W0);
end
%% Data Show

figure(11),plot(x,y, 'b.','MarkerSize',5);

[v,~] = eig(D'*D);

Xinit = v(:,1); %
% Xinit = rand(5,1);
Xinit = Xinit ./ sqrt(Xinit'*Xinit);
Einit = zeros(L,1);


%% RUN DSP
disp('******************* Runing DSP  *************************')
tic 

[W,MaxItr,cost,E]= DSP_Optim1(D, Xinit,Einit,opts);
if size(cost,1) >1
    cost = cost';
end

time =toc;
W1 = W./sqrt(W'*W);
%% Plot results
CorrectIndex = zeros(L,1);
% CorrectIndex(abs(E) < 0.001.*max(abs(E))) = 1;
CorrectIndex(abs(E) < 0.01) = 1;
CorrectIndex = logical(CorrectIndex);
figure(22), plot(x(CorrectIndex), y(CorrectIndex),'b*','MarkerSize',5),hold on;
plot(x(~CorrectIndex), y(~CorrectIndex),'ro'),hold off;
legend('Inliers','outliers');
axis equal

%%
disp('************* Solution for problem: [1, x, y, x.^2, y.^2]*W = 0  *******')
disp('The optimal solution of W: ')
disp(W0');
disp(['Estimated W by DSP for ',ModelType,' model: '])
disp(W1')
disp('RunTime (s):'); disp(time)

disp('Cost:'); disp(cost)

disp('*****************************************************************')