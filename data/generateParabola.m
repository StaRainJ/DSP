function [x, y,groups] = generateParabola(Ni,No)

LINE = [200, 500, 0.004, 200, 2, Ni];
gp = 0;
x =[];y = []; groups = [];
%%
for i = 1 %:maxDisplayStructure
    x1 = LINE(i,1); y1 = LINE(i,2); a = LINE(i,3); d = LINE(i,4); sigma =LINE(i,5);  NN = LINE(i,end);
    t = rand (NN,1);
    xnew = rand(NN,1).*600+50;
    ynew =  a.*( xnew - x1).*(xnew - y1)+d + sigma.*rand(NN,1);
    x = [x;xnew];
    y = [y;ynew];
     gp = gp+1;
    groups = [groups; gp.*ones(NN,1)];
end

% figure, plot(x,y,'b.');
%  X = [x,y]' ;
 
x = mapminmax(x',-0.9,0.9)';
y = mapminmax(y',-0.8,0.8)';
 

%%  generate outliers
NOISE = [2, 2, No];

range_x = NOISE(1);
range_y = NOISE(2);
amount = NOISE(3);

NoiseX = rand(amount,1)*range_x-1;
NoiseY = rand(amount,1)*range_y-1;

x= [x;NoiseX]; y = [y;NoiseY];
groups = [groups;zeros(amount,1)];
