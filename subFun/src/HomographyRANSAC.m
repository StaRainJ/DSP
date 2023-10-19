function [inliers, H, d] = HomographyRANSAC(pts1h, pts2h, nTrials, threshold, confidence)

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
    d = estTFormDistance(pts1h, pts2h, nPts, integerClass);
    
    [curInliers, curNInliers] = findInliers(d, nPts, threshold);
    
    if bestNInliers < curNInliers
        bestNInliers = curNInliers;
        inliers = curInliers;
        
        % Update the number of trials
        maxNTrials = updateNumTrials(oneOverNPts, logOneMinusConf, ...
            outputClass, integerClass, curNInliers, maxNTrials);
    end
    curNTrials = curNTrials + 1;
end
H = norm4Point(pts1h(:, inliers), pts2h(:, inliers));
d = computeDistance(pts1h, pts2h, H);
inliers = find(d<=threshold);
end

function maxNTrials = updateNumTrials(oneOverNPts, logOneMinusConf, ...
  outputClass, integerClass, curNInliers, maxNTrials)

ratioOfInliers = cast(curNInliers, outputClass) * oneOverNPts;
if ratioOfInliers > cast(1, outputClass) - eps(outputClass)
  newNum = zeros(1, integerClass);
else
  ratio4 = ratioOfInliers^4;
  if ratio4 > eps(ones(1, outputClass))
    logOneMinusRatio4 = log(ones(1, outputClass) - ratio4);
    newNum = cast(ceil(logOneMinusConf / logOneMinusRatio4), integerClass);
  else
    newNum = intmax(integerClass);
  end
end

if maxNTrials > newNum
  maxNTrials = newNum;
end
end

function [d, H] = estTFormDistance(pts1h, pts2h, nPts, integerClass)

indices = cast(randperm(nPts, 4), integerClass);
H = norm4Point(pts1h(:, indices), pts2h(:, indices));
d = computeDistance(pts1h, pts2h, H);
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



function [inliers, nInliers] = findInliers(distance, nPts, threshold)
inliers = find(distance <= threshold);
nInliers = cast(length(inliers), 'like', nPts);
end

function d = computeDistance(pts1h, pts2h, H)
% V = H*pts1h;
% V = V*(diag(1./V(3,:)));
% d = sqrt(sum((pts2h - V).^2)');
d = SampsonDistanceH(pts1h,pts2h,H);
end







