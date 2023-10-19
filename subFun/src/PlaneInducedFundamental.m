function F = PlaneInducedFundamental(H, x, y)

xt = H*x;
xt = xt*diag(1./xt(3,:));

d = [xt(1:2,1) y(1:2,1) xt(1:2,2) y(1:2,2)];
e = linlinintersect(d');
et = [0 -1 e(2);1 0 -e(1);-e(2) e(1) 0];
F = et*H;
F = F/F(3,3);

end
