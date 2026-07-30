needsPackage "Cremona";

-- Exact family calibrations for actual affine-map compactifications.
-- These are graph maps [x0^m:F_1^h:...:F_n^h], not full polar maps.

R2 = QQ[c0,cx,cy];
S2 = QQ[a0,a1,a2];

planeTriangular2 = rationalMap map(
    R2,
    S2,
    {c0^2, c0*cx+cy^2, c0*cy}
);
planeTriangularDegrees2 = projectiveDegrees(
    planeTriangular2,
    Certify => true,
    Verbose => false
);
print(
    "plane triangular quadratic degrees = "
    | toString planeTriangularDegrees2
);
assert(planeTriangularDegrees2 == {1,2,1});

planeTriangular3 = rationalMap map(
    R2,
    S2,
    {c0^3, c0^2*cx+cy^3, c0^2*cy}
);
planeTriangularDegrees3 = projectiveDegrees(
    planeTriangular3,
    Certify => true,
    Verbose => false
);
print(
    "plane triangular cubic degrees = "
    | toString planeTriangularDegrees3
);
assert(planeTriangularDegrees3 == {1,3,1});

R4 = QQ[x0,x,y,t,u];
S4 = QQ[b0,b1,b2,b3,b4];

-- Cotangent lift of F=(x+y^2,y), with potential t*(x+y^2)+u*y.
cotangent2 = rationalMap map(
    R4,
    S4,
    {
        x0^2,
        x0*t,
        2*t*y+x0*u,
        x0*x+y^2,
        x0*y
    }
);
cotangentDegrees2 = projectiveDegrees(
    cotangent2,
    Certify => true,
    Verbose => false
);
print(
    "quadratic plane cotangent degrees = "
    | toString cotangentDegrees2
);
assert(cotangentDegrees2 == {1,2,3,2,1});

-- Cotangent lift of F=(x+y^3,y), with potential t*(x+y^3)+u*y.
cotangent3 = rationalMap map(
    R4,
    S4,
    {
        x0^3,
        x0^2*t,
        3*t*y^2+x0^2*u,
        x0^2*x+y^3,
        x0^2*y
    }
);
cotangentDegrees3 = projectiveDegrees(
    cotangent3,
    Certify => true,
    Verbose => false
);
print(
    "cubic plane cotangent degrees = "
    | toString cotangentDegrees3
);
assert(cotangentDegrees3 == {1,3,5,3,1});

R5 = QQ[z0,z1,z2,z3,z4,z5];
S5 = QQ[d0,d1,d2,d3,d4,d5];

-- Quadratic stabilization of the triangular four-variable potential.
stableTriangular2 = rationalMap map(
    R5,
    S5,
    {
        z0^2,
        z0*z2,
        z0*z1+z2^2,
        z0*z4,
        z0*z3,
        z0*z5
    }
);
stableTriangularDegrees2 = projectiveDegrees(
    stableTriangular2,
    Certify => true,
    Verbose => false
);
print(
    "quadratically stabilized triangular quadratic degrees = "
    | toString stableTriangularDegrees2
);
assert(stableTriangularDegrees2 == {1,2,2,2,2,1});

stableTriangular3 = rationalMap map(
    R5,
    S5,
    {
        z0^3,
        z0^2*z2,
        z0^2*z1+z2^3,
        z0^2*z4,
        z0^2*z3,
        z0^2*z5
    }
);
stableTriangularDegrees3 = projectiveDegrees(
    stableTriangular3,
    Certify => true,
    Verbose => false
);
print(
    "quadratically stabilized triangular cubic degrees = "
    | toString stableTriangularDegrees3
);
assert(stableTriangularDegrees3 == {1,3,3,3,3,1});

print("PASS: cotangent lift preserves the top degree but changes the interior list");
print("PASS: quadratic stabilization preserves top degree but appends an interior term");
