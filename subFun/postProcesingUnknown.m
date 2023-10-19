function [FinalIndex,ModelType,par,ee,e1,e2] = postProcesingUnknown(opts,X, Y,th, C0,D,X1,X2,InlierIndex1, InlierIndex2, inliers)
LL1 = length(InlierIndex1);
LL2 = length(InlierIndex2);

if LL2 < 10
    InlierIndex2 = InlierIndex1;
end


%% PostProcessing
[VFCIndexH, H_est, d] = Post_H2(X, Y, 300, th, InlierIndex2);
 e1 = mean(SampsonDistanceH(X(1:3,inliers),Y(1:3,inliers),H_est));
 dH = d(C0(InlierIndex1)); 
%  dH = d(C0(InlierIndex1(InlierIndex2))); 
% MH= mean(dH);
 dH(dH>5) = 5;
 ddH = mean(dH);

[VFCIndexF, F_est, d] = Post_F2(X, Y, 300, th, InlierIndex1); 
e2 = mean(SampsonDistanceF(X(1:3,inliers),Y(1:3,inliers), F_est));
dF = d(C0(InlierIndex1));
% dF(dF>th) = th;
% MF= mean(dF);
dF(dF>5) = 5;
ddF = mean(dF);

 %% Model selection  Using RSBC
[ccc,Esum,CC,Sparsity1,Sparsity2 ] = ModelSelection(D,X1,X2,InlierIndex1,InlierIndex2);
par = ddH./ccc./ddF;
if ddH < 3*ccc*ddF    % 2  4
    Model = 'H';
else
    Model = 'F';
end

%%
if strcmp(Model,'H')
    FinalIndex = VFCIndexH;
ee = e1;
ModelType = 1;
else 
 FinalIndex = VFCIndexF;
ee = e2;
ModelType = 2;
end