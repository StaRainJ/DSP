function [bestInliers,bestf] = LocalOptimizationF(curInliers, pts1h, pts2h, t)
f = norm8Point(pts1h(:, curInliers), pts2h(:, curInliers));
bestInliers = curInliers;
bestNInliers = length(bestInliers);
bestf = f;
irlsSteps = 10;
s = 8;

th_multiplier = 4*sqrt(2); th_step_size = (th_multiplier*t - t)./irlsSteps;

for loirls = 0:irlsSteps
    d = SampsonDistanceF(pts1h, pts2h, f);
    loind = find(d<=(th_multiplier*t - th_step_size*loirls));
    if length(loind)>=8
        loind2 = randsample(loind, min(s*7, length(loind)));
        w = 1./(1+3*d(loind2)/t);
        f = weightedNorm8Point(pts1h(:, loind2), pts2h(:, loind2), w);
        d = SampsonDistanceF(pts1h, pts2h, f);
        loind = find(d<=t);
        if length(loind) > bestNInliers
            bestInliers = loind;
            bestNInliers = length(bestInliers);
            bestf = f;
        end
    end
end