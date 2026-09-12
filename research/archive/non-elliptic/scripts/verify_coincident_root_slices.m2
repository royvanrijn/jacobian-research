needsPackage "CoincidentRootLoci";

-- Degree five, partition (3,2).  The package uses t_j for the coefficient of
-- x_0^(n-j)x_1^j, so dehomogenizing x_0=1 identifies t_j with [W^j]M.
X5 = coincidentRootLocus {3,2};
S5 = ring X5;
R5 = QQ[a0,a1,h3,h4,MonomialOrder=>Eliminate 2];
h5 = (-1-h3-2*h4)/3;
h2 = -(h3+h4+h5);
toSlice5 = map(R5,S5,{a0,a1,-h2,-h3,-h4,-h5});
J5 = toSlice5(ideal X5);
Hsecond5 = 2*h2+6*h3+12*h4+20*h5;
K5 = trim eliminate(saturate(J5,ideal(h5*(Hsecond5+2))),{a0,a1});
assert(dim K5 == 3);
assert(codim K5 == 1);

print "PASS: CRL(3,2) projects to a hypersurface of seed dimension one";

-- Degree six, partition (2,2,2).
X6222 = coincidentRootLocus {2,2,2};
S6222 = ring X6222;
R6222 = QQ[a06222,a16222,h36222,h46222,h56222,MonomialOrder=>Eliminate 2];
h66222 = (-1-h36222-2*h46222-3*h56222)/4;
h26222 = -(h36222+h46222+h56222+h66222);
toSlice6222 = map(R6222,S6222,{a06222,a16222,-h26222,-h36222,-h46222,-h56222,-h66222});
J6222 = toSlice6222(ideal X6222);
Hsecond6222 = 2*h26222+6*h36222+12*h46222+20*h56222+30*h66222;
K6222 = trim eliminate(saturate(J6222,ideal(h66222*(Hsecond6222+2))),{a06222,a16222});
assert(dim K6222 == 4);
assert(codim K6222 == 1);
print "PASS: CRL(2,2,2) in degree six has seed dimension two";

-- Degree six, partition (3,3).
X633 = coincidentRootLocus {3,3};
S633 = ring X633;
R633 = QQ[a0633,a1633,h3633,h4633,h5633,MonomialOrder=>Eliminate 2];
h6633 = (-1-h3633-2*h4633-3*h5633)/4;
h2633 = -(h3633+h4633+h5633+h6633);
toSlice633 = map(R633,S633,{a0633,a1633,-h2633,-h3633,-h4633,-h5633,-h6633});
J633 = toSlice633(ideal X633);
Hsecond633 = 2*h2633+6*h3633+12*h4633+20*h5633+30*h6633;
K633 = trim eliminate(saturate(J633,ideal(h6633*(Hsecond633+2))),{a0633,a1633});
assert(dim K633 == 3);
assert(codim K633 == 2);
print "PASS: CRL(3,3) in degree six has seed dimension one";

-- Degree seven, partition (3,2,2).
X7322 = coincidentRootLocus {3,2,2};
S7322 = ring X7322;
R7322 = QQ[a07322,a17322,h37322,h47322,h57322,h67322,MonomialOrder=>Eliminate 2];
h77322 = (-1-h37322-2*h47322-3*h57322-4*h67322)/5;
h27322 = -(h37322+h47322+h57322+h67322+h77322);
toSlice7322 = map(R7322,S7322,{a07322,a17322,-h27322,-h37322,-h47322,-h57322,-h67322,-h77322});
J7322 = toSlice7322(ideal X7322);
Hsecond7322 = 2*h27322+6*h37322+12*h47322+20*h57322+30*h67322+42*h77322;
K7322 = trim eliminate(saturate(J7322,ideal(h77322*(Hsecond7322+2))),{a07322,a17322});
assert(dim K7322 == 4);
assert(codim K7322 == 2);
print "PASS: CRL(3,2,2) in degree seven has seed dimension two";

-- Degree eight is audited on the package parameterization rather than by
-- asking the package to implicitize the entire ambient CRL first.  The latter
-- consumes several gigabytes before the normalized two-coefficient
-- elimination.  Dehomogenizing each marked root gives the standard root
-- normalization; maximal rank on Phi=0 proves that its projected admissible
-- slice has the expected dimension.
X82222 = coincidentRootLocus {2,2,2,2};
f82222 = map X82222;
forms82222 = flatten entries matrix f82222;
A82222 = QQ[r182222,r282222,r382222,r482222];
deh82222 = map(A82222,ring forms82222_0,{1,r182222,1,r282222,1,r382222,1,r482222});
c82222 = apply(forms82222,g -> deh82222(g));
phi82222 = sum(2..8,j -> c82222_j);
d82222 = sum(2..8,j -> j*c82222_j);
second82222 = sum(2..8,j -> j*(j-1)*c82222_j);
jac82222 = jacobian matrix{prepend(phi82222,c82222_{2..8})};
good82222 = c82222_8*d82222*(second82222-2*d82222);
slice82222 = trim saturate(saturate(ideal(phi82222),ideal(good82222)),minors(4,jac82222));
assert(dim X82222 == 4);
assert(codim slice82222 == 1);
print "PASS: parameterized CRL(2,2,2,2) has an admissible projected slice of seed dimension three";

-- Degree eight, partition (3,3,2).
X8332 = coincidentRootLocus {3,3,2};
f8332 = map X8332;
forms8332 = flatten entries matrix f8332;
A8332 = QQ[r18332,r28332,r38332];
deh8332 = map(A8332,ring forms8332_0,{1,r18332,1,r28332,1,r38332});
c8332 = apply(forms8332,g -> deh8332(g));
phi8332 = sum(2..8,j -> c8332_j);
d8332 = sum(2..8,j -> j*c8332_j);
second8332 = sum(2..8,j -> j*(j-1)*c8332_j);
jac8332 = jacobian matrix{prepend(phi8332,c8332_{2..8})};
good8332 = c8332_8*d8332*(second8332-2*d8332);
slice8332 = trim saturate(saturate(ideal(phi8332),ideal(good8332)),minors(3,jac8332));
assert(dim X8332 == 3);
assert(codim slice8332 == 1);
print "PASS: parameterized CRL(3,3,2) has an admissible projected slice of seed dimension two";
