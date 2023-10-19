function f = norm8Point(pts1h, pts2h)
% Normalize the points
num = size(pts1h, 2);

[pts1h, T1] = normalise2dpts(pts1h);
[pts2h, T2] = normalise2dpts(pts2h);
pts1h(isnan(pts1h)) = 1;
pts2h(isnan(pts2h)) = 1;

% Compute the constraint matrix
m = zeros(num, 9);
for idx = 1: num
  m(idx,:) = [...
    pts1h(1,idx)*pts2h(1,idx), pts1h(2,idx)*pts2h(1,idx), pts2h(1,idx), ...
    pts1h(1,idx)*pts2h(2,idx), pts1h(2,idx)*pts2h(2,idx), pts2h(2,idx), ...
                 pts1h(1,idx),              pts1h(2,idx), 1];
end

% Find out the eigen-vector corresponding to the smallest eigen-value.
[~, ~, vm] = svd(m, 0);
f = reshape(vm(:, end), 3, 3)';

% Enforce rank-2 constraint
[u, s, v] = svd(f);
s(end) = 0;
f = u * s * v';

f = f/norm(f(:));

f = T2'*f*T1;

end