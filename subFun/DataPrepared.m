function [X0,Y0,CorrectIndex,I1,I2]=DataPrepared(X0,Y0,CorrectIndex,I1,I2)
%% for better visible
N0 = size(X0,1);
idx1 = CorrectIndex;
idx0 = setdiff(1:N0,CorrectIndex);

% X0 = X0(idx1,:);
% Y0 = Y0(idx1,:); % for testing on only inliers

X00 = X0; Y00 = Y0;

X0 = X00([idx1,idx0],:);
Y0 = Y00([idx1,idx0],:);
inliers = 1:length(CorrectIndex);
CorrectIndex = inliers;

%% uniform image size
if size(I1,3)==1
    I1=repmat(I1,[1 1 3]);
end
if size(I2,3)==1
    I2=repmat(I2,[1 1 3]);
end