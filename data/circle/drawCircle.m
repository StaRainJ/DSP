function [ output_args ] = drawCircle( x,y,r )
%DRAWCIRCLE draw a circle of given center (x,y) and radius r
%   Detailed explanation goes here
t=[0:0.1:2*pi];
X=x+r*cos(t);
Y=y+r*sin(t);
plot(X,Y,'b.');
axis equal

end

