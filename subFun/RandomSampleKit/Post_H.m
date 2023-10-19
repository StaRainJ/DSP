function [inliers, H, d] = Post_H(pts1h, pts2h, nTrials, th)

N = size(pts1h,2);
inliers = [];
maxNTrials = nTrials;
curNTrials = 0;
bestScore = -N;

while curNTrials < maxNTrials
    [~, curInliers, ~] = MinimalSample_H(pts1h, pts2h, N, th);
    curScore = length(curInliers);
    if bestScore < curScore    
        inliers = curInliers;
        [inliers,~] = LocalOptimizationH(inliers, pts1h, pts2h, th);
        bestScore = length(inliers);
    end
    curNTrials = curNTrials + 1;
end
if length(inliers)>=4
    H = norm4Point(pts1h(:, inliers), pts2h(:, inliers));
    d = SampsonDistanceH(pts1h, pts2h, H);
else
    H = rand(3,3);
    d = SampsonDistanceH(pts1h, pts2h, H);
end
inliers = find(d<=th);
end