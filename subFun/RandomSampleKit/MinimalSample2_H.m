function [H, Inlier, indices] = MinimalSample2_H(pts1h, pts2h, idx, th)

indices = randsample(idx, 4);
H = norm4Point(pts1h(:, indices), pts2h(:, indices));
d = SampsonDistanceH(pts1h, pts2h, H);
Inlier = find(d<=th);

end