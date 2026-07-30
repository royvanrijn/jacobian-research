needsPackage "Cremona";

R = QQ[x0,x1,x2,x3,x4];
S = QQ[y0,y1,y2,y3,y4];

-- Compactification of grad(x1*x2+x3*x4+x2^3/3).
graphMap = rationalMap map(
    R,
    S,
    {
        x0^2,
        x0*x2,
        x0*x1+x2^2,
        x0*x4,
        x0*x3
    }
);

-- The full polar map has a different first coordinate.
polarMap = rationalMap map(
    R,
    S,
    {
        x1*x2+x3*x4,
        x0*x2,
        x0*x1+x2^2,
        x0*x4,
        x0*x3
    }
);

graphDegrees = projectiveDegrees(graphMap, Certify => true, Verbose => false);
polarDegrees = projectiveDegrees(polarMap, Certify => true, Verbose => false);

print("affine-gradient graph compactification degrees = " | toString graphDegrees);
print("full-polar degrees = " | toString polarDegrees);
assert(graphDegrees == {1,2,2,2,1});
assert(polarDegrees == {1,2,4,4,2});

print("PASS: the top degrees differ (1 versus 2)");

-- Cubic-gradient version of the same triangular constant-Hessian family.
graphMap3 = rationalMap map(
    R,
    S,
    {
        x0^3,
        x0^2*x2,
        x0^2*x1+x2^3,
        x0^2*x4,
        x0^2*x3
    }
);
polarMap3 = rationalMap map(
    R,
    S,
    {
        2*x0*(x1*x2+x3*x4),
        x0^2*x2,
        x0^2*x1+x2^3,
        x0^2*x4,
        x0^2*x3
    }
);

graphDegrees3 = projectiveDegrees(graphMap3, Certify => true, Verbose => false);
polarDegrees3 = projectiveDegrees(polarMap3, Certify => true, Verbose => false);

print("cubic affine-gradient graph degrees = " | toString graphDegrees3);
print("cubic full-polar degrees = " | toString polarDegrees3);
assert(graphDegrees3 == {1,3,3,3,1});
assert(polarDegrees3 == {1,3,6,6,3});

print("PASS: the cubic top degrees differ (1 versus 3)");
