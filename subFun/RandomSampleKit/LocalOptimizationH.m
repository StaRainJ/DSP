function [bestInliers,bestH] = LocalOptimizationH(curInliers,pts1h,pts2h,t)
H = norm4Point(pts1h(:, curInliers), pts2h(:, curInliers));
bestInliers = curInliers;
bestNInliers = length(bestInliers);
bestH = H;
irlsSteps = 10;
s = 4;

th_multiplier = 4; th_step_size = (th_multiplier*t - t)./irlsSteps;

for loirls = 0:irlsSteps
    d = SampsonDistanceH(pts1h, pts2h, H);
    loind = find(d<=(th_multiplier*t - th_step_size*loirls));
    if length(loind) >= 4
        loind2 = randsample(loind, min(s*7, length(loind)));
        w = 1./(1+3*d(loind2)/t);
        H = weightedNorm4Point(pts1h(:, loind2), pts2h(:, loind2), w);
        d = SampsonDistanceH(pts1h, pts2h, H);
        loind = find(d<=t);
        if length(loind) > bestNInliers
            bestInliers = loind;
            bestNInliers = length(bestInliers);
            bestH = H;
        end
    end
end
end