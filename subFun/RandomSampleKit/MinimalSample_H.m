function [H, Inlier, indices] = MinimalSample_H(pts1h, pts2h, nPts, th)

indices = randperm(nPts, 4);
H = norm4Point(pts1h(:, indices), pts2h(:, indices));
d = SampsonDistanceH(pts1h, pts2h, H);
Inlier = find(d<=th);

end