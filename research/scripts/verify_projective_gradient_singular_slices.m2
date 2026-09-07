-- Exact singular-top DVR-profile calibrations over QQ.

-- On the chart y=1, h5=x^3*y^2 has transverse top Jacobian algebra
-- QQ[x]/(x^2), of length two at the repeated root x=0.
R = QQ[e,x];

-- h4=0: the active quotient is free of rank two over QQ[[e]].
Iflat = ideal(e^4,3*x^2,2*x^3);

-- h4=x*y^3: the active quotient has profile
-- R/(e^2) direct-sum R/(e), hence truncated length three.
Imixed = ideal(e^4,3*x^2+e,2*x^3+3*e*x);

-- h4=y^4: both special-fiber generators have order-one e-torsion.
IorderOne = ideal(e^4,3*x^2,2*x^3+4*e);

Qflat = R/Iflat;
assert(degree Qflat == 8);
Qmixed = R/Imixed;
assert(degree Qmixed == 3);
QorderOne = R/IorderOne;
assert(degree QorderOne == 2);

print("binary repeated-root active lengths = 8,3,2");
print("PASS: singular-top flat and torsion profiles over QQ");
