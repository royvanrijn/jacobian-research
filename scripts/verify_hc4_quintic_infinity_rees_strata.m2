needsPackage "Cremona";
needsPackage "ReesAlgebra";

-- Generic smooth essential Hessian-rank strata for a quintic top part.
-- Fermat representatives suffice for the exact complete-intersection and
-- linear-type calculation.  The Python companion constructs the universal
-- coefficient forms and verifies their Euler/Hessian/curl identities.

checkRankOne = () -> (
    R1 := QQ[x0,x1,x2,x3,x4];
    T1 := QQ[y0,y1,y2,y3,y4];
    I1 := ideal(x0^4,x1^4);
    assert(codim I1 == 2);
    assert(degree I1 == 16);
    assert(numgens I1 == codim I1);
    assert(isLinearType I1);
    K1 := reesIdeal I1;
    assert(numgens K1 == 1);
    phi1 := rationalMap map(R1,T1,{x0^4,x1^4,0_R1,0_R1,0_R1});
    d1 := projectiveDegrees(phi1, Certify => true, Verbose => false);
    print("rank-one pure-top projective degrees = " | toString d1);
    assert(d1 == {1,4,0,0,0});
);

checkRankTwo = () -> (
    R2 := QQ[z0,z1,z2,z3,z4];
    T2 := QQ[v0,v1,v2,v3,v4];
    I2 := ideal(z0^4,z1^4,z2^4);
    assert(codim I2 == 3);
    assert(degree I2 == 64);
    assert(numgens I2 == codim I2);
    assert(isLinearType I2);
    K2 := reesIdeal I2;
    assert(numgens K2 == 3);
    phi2 := rationalMap map(R2,T2,{z0^4,z1^4,z2^4,0_R2,0_R2});
    d2 := projectiveDegrees(phi2, Certify => true, Verbose => false);
    print("rank-two pure-top projective degrees = " | toString d2);
    assert(d2 == {1,4,16,0,0});
);

checkRankThree = () -> (
    R3 := QQ[w0,w1,w2,w3,w4];
    T3 := QQ[q0,q1,q2,q3,q4];
    I3 := ideal(w0^4,w1^4,w2^4,w3^4);
    assert(codim I3 == 4);
    assert(degree I3 == 256);
    assert(numgens I3 == codim I3);
    assert(isLinearType I3);
    K3 := reesIdeal I3;
    assert(numgens K3 == 6);
    phi3 := rationalMap map(R3,T3,{w0^4,w1^4,w2^4,w3^4,0_R3});
    d3 := projectiveDegrees(phi3, Certify => true, Verbose => false);
    print("rank-three pure-top projective degrees = " | toString d3);
    assert(d3 == {1,4,16,64,0});
);

checkRankOne();
checkRankTwo();
checkRankThree();

print("PASS: generic smooth quintic top ideals are complete intersections");
print("PASS: their Rees ideals are linear type with only Koszul equations");
