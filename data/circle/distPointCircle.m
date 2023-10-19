function [ d ] = distPointCircle( x,C )
%DISTPOINTLINE calcola la distanza euclidea tra un punto x e un piano H
%   Detailed explanation goes here

d=abs(norm(x-C(1:2))-C(3));

end

