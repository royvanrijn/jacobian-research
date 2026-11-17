-- Independent complete-intersection and filtered-length calibrations for
-- the all-dimensional smooth-essential gradient normal-slice theorem.

-- (m,r)=(3,2): B has length 9 and A has length 27.
R32 = QQ[e32,x32,y32];
I32 = ideal(e32^3,x32^3,y32^3);
I32missing = I32 + ideal(e32*x32);
A32 = R32/I32;
Q32 = R32/I32missing;
assert(degree A32 == 27);
assert(degree Q32 == 15);
assert(degree A32-degree Q32 == 12);
-- dim(B*x)=6 and (m-q)*dim(B*x)=2*6=12.

-- (m,r)=(4,3): the HC4 isolated-vertex active algebra.
R43 = QQ[e43,x43,y43,z43];
I43 = ideal(e43^4,x43^4,y43^4,z43^4);
I43missing = I43 + ideal(e43*x43^3);
A43 = R43/I43;
Q43 = R43/I43missing;
assert(degree A43 == 256);
assert(degree A43-degree Q43 == 48);
-- dim(B*x^3)=16 and (m-q)*dim(B*x^3)=3*16=48.

-- (m,r)=(4,2): a unit penultimate term kills epsilon and leaves length 16.
R42 = QQ[e42,x42,y42];
I42 = ideal(e42^4,x42^4,y42^4);
I42unit = I42 + ideal(e42);
A42 = R42/I42;
Q42unit = R42/I42unit;
assert(degree A42 == 64);
assert(degree Q42unit == 16);

print("(m,r)=(3,2) active/drop lengths = 27,12");
print("(m,r)=(4,3) active/drop lengths = 256,48");
print("(m,r)=(4,2) unit-penultimate Segre multiplicity = 16");
print("PASS: all-dimensional normal-slice calibrations over QQ");
