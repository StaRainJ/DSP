function [inliers, H, d] = HomographyLORANSAC2(pts1h, pts2h, idx, nTrials, threshold)

inliers = [];

maxNTrials = nTrials;
curNTrials = 0;
bestNInliers = 0;

while curNTrials < maxNTrials
    d = estTFormDistance(pts1h, pts2h, idx);    
    [curInliers, curNInliers] = findInliers(d, threshold);
    if bestNInliers < curNInliers
        bestNInliers = curNInliers;
        inliers = curInliers;
        if bestNInliers >= 4
            [curInliers,~] = LocalOptimization(curInliers, pts1h, pts2h, threshold);
            curNInliers = length(curInliers);
            if bestNInliers < curNInliers
                bestNInliers = curNInliers;
                inliers = curInliers;
            end
        end
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
inliers = find(d<=threshold);
end


function [d, H] = estTFormDistance(pts1h, pts2h, idx)

indices = randperm(length(idx), 4);
indices = idx(indices);
H = norm4Point(pts1h(:, indices), pts2h(:, indices));
d = SampsonDistanceH(pts1h, pts2h, H);
end

function H = norm4Point(xs1,xs2)
[~,c] = size(xs1);

if (size(xs1) ~= size(xs2))
 error ('Input point sets are different sizes!')
end

if (size(xs1,1) == 2)
  xs1 = [xs1 ; ones(1,size(xs1,2))];
  xs2 = [xs2 ; ones(1,size(xs2,2))];
end

[xs1, T1] = normalise2dpts(xs1);
[xs2, T2] = normalise2dpts(xs2);
xs1(isnan(xs1)) = 1;
xs2(isnan(xs2)) = 1;

D = [];
ooo  = zeros(1,3);
for k=1:c
  p1 = xs1(:,k);
  p2 = xs2(:,k);
  D = [ D;
    p1'*p2(3) ooo -p1'*p2(1)
    ooo p1'*p2(3) -p1'*p2(2)
   ];
end

[~,~,v] = svd(D, 0);
h = v(:,9);
H = reshape(h,3,3)';

H = T2^(-1)*H*T1;

end

function [inliers, nInliers] = findInliers(distance, threshold)
inliers = find(distance <= threshold);
nInliers = length(inliers);
end

function [bestInliers,bestH] = LocalOptimization(curInliers,pts1h,pts2h,t)

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
        if length(loind)>bestNInliers
            bestInliers = loind;
            bestNInliers = length(bestInliers);
            bestH = H;
        end
    end
end
end


function H = weightedNorm4Point(xs1,xs2,w)
[~,c] = size(xs1);

if (size(xs1) ~= size(xs2))
 error ('Input point sets are different sizes!')
end

if (size(xs1,1) == 2)
  xs1 = [xs1 ; ones(1,size(xs1,2))];
  xs2 = [xs2 ; ones(1,size(xs2,2))];
end

[xs1, T1] = normalise2dpts(xs1);
[xs2, T2] = normalise2dpts(xs2);
xs1(isnan(xs1)) = 1;
xs2(isnan(xs2)) = 1;

D = [];
ooo  = zeros(1,3);
for k=1:c
  p1 = xs1(:,k);
  p2 = xs2(:,k);
  D = [ D;
    p1'*p2(3) ooo -p1'*p2(1)
    ooo p1'*p2(3) -p1'*p2(2)
   ];
end
w = kron(diag(sqrt(w)),[1 0;0 1]);
[~,~,v] = svd(w*D, 0);
h = v(:,9);
H = reshape(h,3,3)';

H = T2^(-1)*H*T1;
end