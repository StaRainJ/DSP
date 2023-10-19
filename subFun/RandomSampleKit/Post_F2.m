function [inliers, F, d] = Post_F2(pts1h, pts2h, nTrials, th, idx)

N = size(pts1h,2);
inliers = [];
maxNTrials = nTrials;
curNTrials = 0;
bestScore = -N;

while curNTrials < maxNTrials
    [~, curInliers, ~] = MinimalSample2_F(pts1h, pts2h, idx, th);
    curScore = length(curInliers);
    if bestScore < curScore    
%         [H, degeneracy] = DegeneracyCheck(pts1h, pts2h, indices, th);
%         if degeneracy
%             inliers = DegeneracyUpdate(pts1h, pts2h, th, inliers, H);
%         else
%             inliers = curInliers;
%         end
        inliers = curInliers;
        [inliers, ~] = LocalOptimizationF(inliers, pts1h, pts2h, th);
        bestScore = length(inliers);
        
%         maxNTrials = min(updateNumTrials(1/N, log(1-0.99), length(inliers), maxNTrials), nTrials);
%         maxNTrials = max(maxNTrials, 100);
    end
    curNTrials = curNTrials + 1;
end

bestScore = [];


if length(inliers)>=8
    F = norm8Point(pts1h(:, inliers), pts2h(:, inliers)); % diag(pts2h'*F*pts1h)  = 0;
    d = SampsonDistanceF(pts1h, pts2h, F);
else
    F = rand(3,3);
    d = SampsonDistanceF(pts1h, pts2h, F);
end
inliers = find(d<=th);

end