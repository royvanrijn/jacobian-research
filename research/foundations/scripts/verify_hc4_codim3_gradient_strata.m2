-- Independent exact replay for the HC4 codimension-three gradient strata.
-- The canonical note supplies the field-theoretic squarefree divisibility
-- and local Segre arguments.

-- Rank-two quartic-face synchronization after h4|K=t^4.
R = QQ[a_0..a_9,t,w];
F = t^4
    + a_0
    + a_1*w + a_2*t
    + a_3*w^2 + a_4*t*w + a_5*t^2
    + a_6*w^3 + a_7*t*w^2 + a_8*t^2*w + a_9*t^3;
E = ideal(
    24*a_7,
    72*a_6,
    24*a_3+12*a_7*a_9-4*a_8^2,
    36*a_6*a_9-4*a_7*a_8,
    12*a_3*a_9-4*a_4*a_8+4*a_5*a_7,
    12*a_6*a_8-4*a_7^2,
    4*a_3*a_8-4*a_4*a_7+12*a_5*a_6,
    4*a_3*a_5-a_4^2
);
GE = gb E;
assert(a_3^2 % GE == 0);
assert(a_4^3 % GE == 0);
assert(a_6 % GE == 0);
assert(a_7 % GE == 0);
assert(a_8^3 % GE == 0);

-- Generic transverse length on the kernel line.
S = QQ[e,u,v];
Iactive = ideal(e^4,u^4,v^4);
Ifull = Iactive + ideal(e);
Aactive = S/Iactive;
Afull = S/Ifull;
assert(degree Aactive == 64);
assert(degree Afull == 16);
assert(apply(toList(0..6), index ->
    hilbertFunction(index,Afull)) ==
    {1,2,3,4,3,2,1});

-- Ordinary rank-three singularity calibration.
T = QQ[x,y,z];
h5node = x^3*(y^2+z^2)+y^5+z^5;
Cnode = substitute(diff(vars T, transpose diff(vars T,h5node)),
    {x=>1_T,y=>0_T,z=>0_T});
assert(rank Cnode == 2);
assert(Cnode * transpose matrix{{1_T,0_T,0_T}} == 0);

print("rank-two active/full transverse lengths = " |
    toString(degree Aactive) | "," | toString(degree Afull));
print("PASS: rank-two radical synchronization powers");
print("PASS: generic transverse multiplicity sigma_3=16");
print("PASS: nodal rank-three Hessian has radial kernel and rank two");
