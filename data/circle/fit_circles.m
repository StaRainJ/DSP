function [ M ] = fit_circles(X,C)
%FIT_CIRCLE Summary of this function goes here
%
%   fits a circle  in x,y plane in a more accurate
%   (less prone to ill condition )
%  procedure than circfit2 but using more memory
%  x,y are column vector where (x(i),y(i)) is a measured point
%
%  result is center point (yc,xc) and radius R
%  an optional output is the vector of coeficient a
% describing the circle's equation
%
%   x^2+y^2+a(1)*x+a(2)*y+a(3)=0
%
label=unique(C);
N=length(label); %numero di cluster;
M=[];



for i=1:N
    
    L  = label(i);
    P = X(:,C==L);
    x=P(1,:)'; y=P(2,:)';
    a=[x y ones(size(x))]\[-(x.^2+y.^2)];
    xc = -.5*a(1);
    yc = -.5*a(2);
    R  =  sqrt((a(1)^2+a(2)^2)/4-a(3));
    
    M=[M,[xc;yc;R]];
    
end

