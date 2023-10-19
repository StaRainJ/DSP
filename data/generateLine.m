function [x, y,groups ]=generateLine(xa,ya,xb,yb,sigma,Ni,No)

LINE = [xa,ya,xb,yb,sigma,Ni,];
 
NOISE = [1, 1, No];

groups = [];
x= []; y = [];
%%
for i = 1
    x1 = LINE(i,1); y1 = LINE(i,2); x2 = LINE(i,3); y2 = LINE(i,4); sigma = LINE(i,5); NN = LINE(i,6);
    Vx = x2 -x1;  Vy = y2 -y1;
    t = rand (NN,1);
    xnew = x1 + t.*Vx + sigma.*rand(NN,1);
    ynew = y1 + t.*Vy + sigma.*rand(NN,1);
    x = [x;xnew];
    y = [y;ynew];
    groups = i.*ones(NN,1);
end



%%  generate outliers
range_x = NOISE(1);
range_y = NOISE(2);
amount = NOISE(3);

NoiseX = rand(amount,1)*range_x;
NoiseY = rand(amount,1)*range_y;
x= [x;NoiseX]; y = [y;NoiseY];
groups = [groups;zeros(amount,1)];


%%
% figure, plot(x,y,'b.');