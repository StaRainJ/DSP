function y = OrthProj(x,y)
% x ,y, z are all dx1 vectors 
y = y - x*x'*y;