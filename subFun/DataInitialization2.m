function [C0] = DataInitialization2(Xt,Yt,idx)
if size(Xt,1) > 2
    Xt = Xt';
    Yt = Yt';
end
C0 = idx;
idx11 = length(idx);
if idx11 >100
K1= min(max(sqrt(idx11./3),5) ,15);  
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
% tau = std(1-pp_inlier);
C0 = find(pp_inlier>0.6);
end
% if idx22<idx11
%     idx=idx2;
% end   
end

