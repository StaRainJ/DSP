function [x ,y, groups] = generate_Circles( num, Ni,No, radius, epsilon )
%GENERATE_CIRCLES generate circles: number of points is uniform in density
%   Detailed explanation goes here

X = zeros(2,num*Ni);
C = [1*rand(1,num);1*rand(1,num)]; %centers
% C = [0;0];
T = rand(num,Ni)*2*pi;
cont = 1;
for i = 1:num
    for j = 1:Ni
    X(:,cont) = C(:,i) + [radius*cos(T(i,j)); radius*sin(T(i,j))];
    cont = cont+1;
    end
  groups = i.*ones(Ni,1);
end
X=X+randn(2,num*Ni).*epsilon;  

x0 = X(1,:);
y0 = X(2,:);
% x = X(1,:);
% y = X(2,:);
x = mapminmax(x0,-1,1);
y = mapminmax(y0,-1,1);
%%  generate outliers
range_x = 2.2;
range_y = 2.2;
amount = No;

NoiseX = rand(amount,1)*range_x-1.0;
NoiseY = rand(amount,1)*range_y-1.0;
x= [x';NoiseX]; y = [y';NoiseY];
groups = [groups;zeros(amount,1)];



end
