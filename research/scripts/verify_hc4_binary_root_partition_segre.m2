-- Active-unit calibrations for singular essential binary quintic tops.

R = QQ[e,x];

-- On y=1 take h5=x^a*y^(5-a) and h4=x*y^3+y^4.
-- The local active generators are
--   a*x^(a-1)+e,
--   (5-a)*x^a+e*(3*x+4).
I2 = ideal(e^4,2*x+e,3*x^2+e*(3*x+4));
I3 = ideal(e^4,3*x^2+e,2*x^3+e*(3*x+4));
I4 = ideal(e^4,4*x^3+e,x^4+e*(3*x+4));

Q2 = R/I2;
assert(degree Q2 == 1);
Q3 = R/I3;
assert(degree Q3 == 2);
Q4 = R/I4;
assert(degree Q4 == 3);

print("active-unit repeated-root lengths for e=2,3,4 are 1,2,3");
print("PASS: HC4 binary root-partition Segre calibrations over QQ");
