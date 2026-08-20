# MW3 component-label enumeration

Exact reduced height Gram:

    H = (1/66) * [[79,17,-1],[17,106,19],[-1,19,259]]

The denominator 66 = 11*3*2 matches the component groups of the semistable
fiber configuration:

    I11 + I3 + I2 + I2 + 6 I1.

For each MW generator enumerate its component labels

    (k11,k3,k2a,k2b)

and use Shioda's formulas:

    <P,P> = 4 + 2(P.O) - sum contr_v(P)

    <P,Q> = 2 + P.O + Q.O - P.Q - sum contr_v(P,Q)

where for I_n

    contr(P) = k(n-k)/n

and

    contr(P,Q) = min(i,j)*(n-max(i,j))/n

for nonidentity components i,j.

Run:

    sage elkies-k3/scripts/enumerate_mw3_component_labels.sage

The best-ranked output minimizes intersections with O and between the three
sections. These configurations should translate directly into local section
valuation conditions at t=infinity,0,1,lambda.
