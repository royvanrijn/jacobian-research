\\ Exact PARI/GP certificate for the fixed-map quintic Hasse resolvent.

b = 'b;
y = 'y;
f = 3*(57395628*b^3 - 26749197*b^2 + 4181544*b - 219512)*(49 - 324*b);

W = ellfromeqn(y^2 - f);
if (W != [0, -7996592727, 0, 21315161806412546304, -18938742282153164303832264192], error("unexpected Jacobian model"));

E = ellinit(W);
if (E.j != -91992891/534488452952, error("unexpected j-invariant"));

G = ellglobalred(E);
if (G[1] != 3202146, error("unexpected conductor"));

T = elltors(E);
if (T[1] != 1, error("nontrivial torsion"));

R = ellrank(E, 4);
if (R[1] != 1 || R[2] != 1, error("rank not proved equal to one"));
if (#R[4] < 1, error("missing non-torsion generator"));

\\ Three rational quartic points.  The last two give new Hasse fibers.
pts = [[37/243,19],[3139/19764,137805/3721],[31687/204363,23067311/707281]];
for (i = 1, #pts, if (pts[i][2]^2 != subst(f,b,pts[i][1]), error("point not on quartic")));

x = 'x;
cases = [[19764*x^2-15372*x+3139,-533628*x^3+592920*x^2-216099*x+25429,[2,3,5,61,9187,7765337]],[204363*x^2-158949*x+31687,-1839267*x^3+2043630*x^2-751770*x+89959,[2,3,19,29,389,751,3121,458701]]];

for (i = 1, #cases, q=cases[i][1]; c=cases[i][2]; ps=cases[i][3]; if (!polisirreducible(q) || !polisirreducible(c), error("reducible factor")); if (!issquare(poldisc(c)/poldisc(q)), error("discriminant fields differ")); if (Vec(factor(abs(poldisc(q*c)))[,1]) != ps, error("incomplete ramified-prime list")); for (j=1,#ps,p=ps[j]; if (#polrootspadic(q,p,20) + #polrootspadic(c,p,20) == 0, error("missing local root"))));

print("PASS: resolvent Jacobian has conductor 3202146, trivial torsion, rank one");
print("PASS: two additional rational points give everywhere-local 2+3 fibers");
