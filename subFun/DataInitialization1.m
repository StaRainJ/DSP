function [C0] = DataInitialization1(Xt,Yt,idx00)
[m,n] = size(Xt);
if m > 2
    Xt = Xt';
    Yt = Yt';
end
idx = 1:m;
if m >=50    
K0=min(max( round(sqrt(m./10)),4) ,15); 
kdtreeX = vl_kdtreebuild(Xt(:,idx00));
kdtreeY = vl_kdtreebuild(Yt(:,idx00));
[neighborX, ~] =vl_kdtreequery(kdtreeX, Xt(:,idx00), Xt, 'NumNeighbors', K0+1) ;
[neighborY, ~] =vl_kdtreequery(kdtreeY, Yt(:,idx00), Yt, 'NumNeighbors', K0+1) ;

bb=[neighborX;neighborY];
as=diff(sort(bb));
dd=(double(as==zeros(size(as,1),size(as,2))));
count_both_point=sum(double(dd))-1;
pp_inlier=count_both_point./K0;
tau = std(pp_inlier);
idx=find(pp_inlier>min(max(tau-0.1,0.1),0.5));
% idx=find(pp_inlier>0.3);
idx11=length(idx);

if idx11>=50
 K1= min(max(round(sqrt(idx11./10)),5) ,15);  
kdtreeX = vl_kdtreebuild(Xt(:,idx));
kdtreeY = vl_kdtreebuild(Yt(:,idx));
[neighborX, ~] =vl_kdtreequery(kdtreeX, Xt(:,idx), Xt, 'NumNeighbors', K1+1) ;
[neighborY, ~] =vl_kdtreequery(kdtreeY, Yt(:,idx), Yt, 'NumNeighbors', K1+1) ;
neighborX=idx(neighborX);neighborY=idx(neighborY);
bb=[neighborX;neighborY];
as=diff(sort(bb));
dd=(double(as==zeros(size(as,1),size(as,2))));
count_both_point=sum(double(dd))-1;
pp_inlier=count_both_point./K1;
% tau = std(pp_inlier);
idx=find(pp_inlier>min(max(tau-0.05,0.2),0.8));
% idx=find(pp_inlier>0.5);
if length(idx) <20
  idx  =  find(pp_inlier>0.3);
end

end

end
C0 =idx;