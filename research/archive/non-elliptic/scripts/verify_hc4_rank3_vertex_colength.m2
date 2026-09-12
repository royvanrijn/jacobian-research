-- Independent calibration for the HC4 smooth rank-three vertex-colength
-- obstruction.  The universal flatness/socle argument is in the canonical
-- note; this file verifies its complete-intersection and filtered-length
-- numerics on exact Fermat and deformed representatives over QQ.

S = QQ[u1,u2,u3];
J = ideal(u1^4,u2^4,u3^4);
B = S/J;

assert(dim B == 0);
assert(degree B == 64);
degreesToCheck = toList(0..9);
assert(apply(degreesToCheck, index -> hilbertFunction(index,B)) ==
    {1,3,6,10,12,12,10,6,3,1});
assert(hilbertFunction(10,B) == 0);

socleMonomial = u1^3*u2^3*u3^3;
assert(sub(socleMonomial,B) != 0_B);
assert(all({u1,u2,u3}, variable ->
    sub(variable*socleMonomial,B) == 0_B));

R = QQ[e,v1,v2,v3];
I0 = ideal(e^4,v1^4,v2^4,v3^4);
missing0 = e*v1^3;
I0full = I0 + ideal(missing0);
A0 = R/I0;
Q0 = R/I0full;
assert(dim A0 == 0);
assert(degree A0 == 256);

drop0 = degree(A0) - degree(Q0);
assert(drop0 == 48);
assert(drop0 >= 6);

-- An integrable epsilon-deformation from
-- h5=(v1^5+v2^5+v3^5)/5 and h4=t*v1^3+v1*v2^3 on the chart t=1.
use R;
G1 = v1^4 + e*(3*v1^2+v2^3);
G2 = v2^4 + 3*e*v1*v2^2;
G3 = v3^4;
assert(diff(v2,G1) == diff(v1,G2));
assert(diff(v3,G1) == diff(v1,G3));
assert(diff(v3,G2) == diff(v2,G3));
I1 = ideal(e^4,G1,G2,G3);
missing1 = e*v1^3;
I1full = I1 + ideal(missing1);
A1 = R/I1;
Q1 = R/I1full;
assert(dim A1 == 0);
assert(degree A1 == 256);

drop1 = degree(A1) - degree(Q1);
assert(drop1 >= 6);

print("Fermat B Hilbert function = " |
    toString apply(degreesToCheck, index -> hilbertFunction(index,B)));
print("pure associated-graded missing-component length = " | toString drop0);
print("deformed missing-component length = " | toString drop1);
print("PASS: exact rank-three vertex-colength calibrations over QQ");
