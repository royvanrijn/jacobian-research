needsPackage "IntegralClosure";

checkPair = (m,r,expected) -> (
    R0 := QQ[P,Q,Z,C,U];
    J := sum(0..r, i ->
        (-1)^(r-i)*binomial(r,i)*P^i*
        sum(0..(m*(r-i)), j ->
            (-1)^j*binomial(m*(r-i),j)*Q^(m*(r-i)-j)*
            U^(r-i+j+1)/(r-i+j+1)));
    A := R0/(C*J-Z*P^(r+1));
    B := integralClosure A;
    pP := minimalPrimes ideal(P_B);
    critical := P_B-U_B*(Q_B-U_B)^m;
    pDelta := minimalPrimes ideal(critical);
    assert(#pP >= 2);
    assert(#pDelta == 2);
    -- The last P-prime is the degree-mr K-cluster.  The first critical
    -- prime is the affine U=P=0 component; the last is E_Delta.
    cluster := pP#(#pP-1);
    delta := pDelta#1;
    contact := cluster+delta;
    assert(dim contact == 2);
    assert(degree contact == expected);
    print concatenate(
        "PASS normalization contact (m,r)=(",
        toString m,
        ",",
        toString r,
        ") has length ",
        toString expected)
    );

checkPair(2,1,4);
checkPair(1,2,2);
checkPair(3,1,9);
checkPair(2,2,8);

print "PASS: independent normalization contacts equal m^2 r";
