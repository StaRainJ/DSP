function bestInlier = DegeneracyUpdate(pts1h, pts2h, th, bestInlier, H)

d = SampsonDistanceH(pts1h,pts2h,H);
ind = find(d>2.0*th);
ind0 = find(d<=th);

if length(ind) >= 2 && length(ind0) >= 6
    bestS = length(bestInlier);
    for i = 1:100
        sp = randsample(ind,2);
        % F = PlaneInducedFundamental(H, pts1h(:,sp), pts2h(:,sp));
        % sp0 = randsample(ind0,6);
        sp0 = ind0;
        tmp = [sp0,sp];
        F = norm8Point(pts1h(:,tmp),pts2h(:,tmp));
        d = SampsonDistanceF(pts1h,pts2h,F);
        curS = length(find(d<=th));
        if curS > bestS
            bestS = curS;
            bestInlier = find(d<=th);
        end
    end
end