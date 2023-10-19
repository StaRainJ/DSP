function x = SphereProj(x)

Z = 1./sqrt(x'*x);
x = x*Z;