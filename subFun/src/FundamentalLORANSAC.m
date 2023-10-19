function [inliers, f, d] = FundamentalLORANSAC(pts1h, pts2h, nTrials, threshold, confidence)

disType = 'sampson';
integerClass = 'int32';
nPts = size(pts1h,2);
outputClass = 'double';


inliers = false(1, nPts);

maxNTrials = nTrials;
curNTrials = zeros(1, integerClass);
bestNInliers = zeros(1, integerClass);
logOneMinusConf = log(ones(1, outputClass) - confidence);
oneOverNPts = ones(1, outputClass) / cast(nPts, outputClass);

while curNTrials < maxNTrials
    d = estTFormDistance(disType, pts1h, pts2h, nPts, outputClass,...
        integerClass);
    
    [curInliers, curNInliers] = findInliers(d, nPts, threshold);
    
    if bestNInliers < curNInliers
        bestNInliers = curNInliers;
        inliers = curInliers;
        if bestNInliers >= 8
            [curInliers,~] = LocalOptimization(curInliers, pts1h, pts2h, threshold, disType);
            curNInliers = length(curInliers);
            if bestNInliers < curNInliers
                bestNInliers = curNInliers;
                inliers = curInliers;
            end
        end
        % Update the number of trials
        maxNTrials = updateNumTrials(oneOverNPts, logOneMinusConf, ...
            outputClass, integerClass, bestNInliers, maxNTrials);
    end
    curNTrials = curNTrials + 1;
end
if length(inliers)>8
    f = norm8Point(pts1h(:, inliers), pts2h(:, inliers), outputClass,...
        integerClass);
    d = computeDistance(disType, pts1h, pts2h, f);
else
    f = rand(3,3);
    d = computeDistance(disType, pts1h, pts2h, f);
end
inliers = find(d<=threshold);
end

function maxNTrials = updateNumTrials(oneOverNPts, logOneMinusConf, ...
  outputClass, integerClass, curNInliers, maxNTrials)

ratioOfInliers = cast(curNInliers, outputClass) * oneOverNPts;
if ratioOfInliers > cast(1, outputClass) - eps(outputClass)
  newNum = zeros(1, integerClass);
else
  ratio8 = ratioOfInliers^8;
  if ratio8 > eps(ones(1, outputClass))
    logOneMinusRatio8 = log(ones(1, outputClass) - ratio8);
    newNum = cast(ceil(logOneMinusConf / logOneMinusRatio8), integerClass);
  else
    newNum = intmax(integerClass);
  end
end

if maxNTrials > newNum
  maxNTrials = newNum;
end
end

function [d, f] = estTFormDistance(disType, pts1h, pts2h, nPts,...
  outputClass, integerClass)

indices = cast(randperm(nPts, 8), integerClass);
f = norm8Point(pts1h(:, indices), pts2h(:, indices), outputClass,...
  integerClass);
d = computeDistance(disType, pts1h, pts2h, f);
end

function f = norm8Point(pts1h, pts2h, outputClass, integerClass)
% Normalize the points
num = cast(size(pts1h, 2), integerClass);
% [pts1h, t1] = vision.internal.normalizePoints(pts1h, 2, outputClass);
% [pts2h, t2] = vision.internal.normalizePoints(pts2h, 2, outputClass);

[pts1h, T1] = normalise2dpts(pts1h);
[pts2h, T2] = normalise2dpts(pts2h);
pts1h(isnan(pts1h)) = 1;
pts2h(isnan(pts2h)) = 1;
% Compute the constraint matrix
m = coder.nullcopy(zeros(num, 9, outputClass));
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

function [inliers, nInliers] = findInliers(distance, nPts, threshold)
inliers = find(distance <= threshold);
nInliers = cast(length(inliers), 'like', nPts);
end

function d = computeDistance(disType, pts1h, pts2h, f)
pfp = (pts2h' * f)';
pfp = pfp .* pts1h;
d = sum(pfp, 1) .^ 2;

if strcmp(disType, 'sampson')
  epl1 = f * pts1h;
  epl2 = f' * pts2h;
  d = d ./ (epl1(1,:).^2 + epl1(2,:).^2 + epl2(1,:).^2 + epl2(2,:).^2);
end
d = sqrt(d);
end

function [bestInliers,bestf] = LocalOptimization(curInliers, pts1h, pts2h, t, disType)

integerClass = 'int32';
outputClass = 'double';

f = norm8Point(pts1h(:, curInliers), pts2h(:, curInliers), outputClass,...
    integerClass);

bestInliers = curInliers;
bestNInliers = length(bestInliers);
bestf = f;
irlsSteps = 10;
s = 8;

th_multiplier = 4*sqrt(2); th_step_size = (th_multiplier*t - t)./irlsSteps;

for loirls = 0:irlsSteps
    d = computeDistance(disType, pts1h, pts2h, f);
    loind = find(d<=(th_multiplier*t - th_step_size*loirls));
    if length(loind)>8
    loind2 = randsample(loind, min(s*7, length(loind)));
    w = 1./(1+3*d(loind2)/t);
    f = weightedNorm8Point(pts1h(:, loind2), pts2h(:, loind2), w, outputClass, integerClass);
    d = computeDistance(disType, pts1h, pts2h, f);
    loind = find(d<=t);
    if length(loind)>bestNInliers
        bestInliers = loind;
        bestNInliers = length(bestInliers);
        bestf = f;
    end
    end
end


end

function f = weightedNorm8Point(pts1h, pts2h, w, outputClass, integerClass)

num = cast(size(pts1h, 2), integerClass);

[pts1h, T1] = normalise2dpts(pts1h);
[pts2h, T2] = normalise2dpts(pts2h);
pts1h(isnan(pts1h)) = 1;
pts2h(isnan(pts2h)) = 1;

% Compute the constraint matrix
m = coder.nullcopy(zeros(num, 9, outputClass));
for idx = 1: num
  m(idx,:) = [...
    pts1h(1,idx)*pts2h(1,idx), pts1h(2,idx)*pts2h(1,idx), pts2h(1,idx), ...
    pts1h(1,idx)*pts2h(2,idx), pts1h(2,idx)*pts2h(2,idx), pts2h(2,idx), ...
                 pts1h(1,idx),              pts1h(2,idx), 1];
end

m = diag(sqrt(w))*m;
% Find out the eigen-vector corresponding to the smallest eigen-value.
[~, ~, vm] = svd(m, 0);
f = reshape(vm(:, end), 3, 3)';

% Enforce rank-2 constraint
[u, s, v] = svd(f);
s(end) = 0;
f = u * s * v';

% Normalize the fundamental matrix.
f = f / norm(f(:));

f = T2'*f*T1;

end