function [ H ] = hpCircles( X,S )
%GENERATEHYPOTESISLINES Summary of this function goes here
%  


m=size(S,1);
H=zeros(3,m);
for i=1:m
    x = X(1:2,S(i,:));
    M = zeros(3,3);
    b = zeros(3,1);
    for j=1:3
        M(j,:) = [x(1,j), x(2,j), 1 ];
        b(j,:)= - x(1,j)^2 - x(2,j)^2;
    end
    y= M\b;
    
    H(1,i) = -y(1)/2;                      % center x
    H(2,i) = -y(2)/2;                      % center y
    H(3,i) = sqrt(H(1,i)^2+H(2,i)^2-y(3)); % radius
     
    
end


