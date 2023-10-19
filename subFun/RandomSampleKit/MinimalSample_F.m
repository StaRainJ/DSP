function [bestF, bestInlier, indices] = MinimalSample_F(pts1h, pts2h, nPts, th)

indices = randperm(nPts, 7);
F = norm7Point(pts1h(:, indices), pts2h(:, indices));
bestF = [];
bestS = -1;
bestInlier = [];
for i = 1:size(F,3)
    d = SampsonDistanceF(pts1h, pts2h, squeeze(F(:,:,i)));
    inlier = find(d<=th);
    if length(inlier) > bestS
        bestS = length(inlier);
        bestF = squeeze(F(:,:,i));
        bestInlier = inlier;
    end
end

end