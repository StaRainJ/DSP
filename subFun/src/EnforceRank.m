function fo = EnforceRank(f)

F = [f(1:3)';f(4:6)';f(7:9)'];
[u, s, v] = svd(F);
s(end) = 0;
FPrime = u * s * v';
fo = [FPrime(1,:) FPrime(2,:) FPrime(3,:)]';

