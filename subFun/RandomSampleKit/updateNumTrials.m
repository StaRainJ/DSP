function maxNTrials = updateNumTrials(oneOverNPts, logOneMinusConf, curNInliers, maxNTrials)

ratioOfInliers = curNInliers * oneOverNPts;
if ratioOfInliers > 1 - eps('double')
  newNum = 1;
else
  ratio8 = ratioOfInliers^8;
  if ratio8 > eps('double')
    logOneMinusRatio8 = log(1 - ratio8);
    newNum = ceil(logOneMinusConf / logOneMinusRatio8);
  else
    newNum = intmax('int32');
  end
end

if maxNTrials > newNum
  maxNTrials = newNum;
end
end