function [ X, G ] = generate_Circles_uniform_density( num_circle,num_points, epsilon )
%GENERATE_CIRCLES generate circles: number of points is uniform in density
%   Detailed explanation goes here

%%
C = [10*rand(1,num_circle);10*rand(1,num_circle)];
radius=3*rand(1,num_circle)+.5;%centers
points=ceil(num_points*radius);
X=zeros(2,sum(points));
G=nan(size(X,2),1);
%%
cont = 1;
for i = 1:num_circle
    for j = 1:points(i)
        if(mod(cont,2))
            T = rand(1,1)*2*pi+pi;
        else
            T = rand(1,1)*pi;
        end
    X(:,cont) = C(:,i) + [radius(i)*cos(T); radius(i)*sin(T)];
    G(cont)=i;
    cont = cont+1;
    end
end
%%
X=X+randn(2,size(X,2)).*epsilon;  

end


