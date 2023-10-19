function [X, Y, Xn, Yn, Xt, Yt] =  DataPreprocessing(X, Y)

 X0 = X;
 Y0 = Y;
 N = size(X,1);
 X = [X0';ones(1,N)];
 Y = [Y0';ones(1,N)];
 
 [Xn,m1,s1] = DataNorm(X0);
 [Yn,m2,s2] = DataNorm(Y0);
 Xt = [Xn';ones(1,N)];
 Yt = [Yn';ones(1,N)];