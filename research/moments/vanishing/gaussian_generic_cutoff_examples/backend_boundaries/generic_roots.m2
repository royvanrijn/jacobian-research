-- Root extraction over Frac(B/P), with exact normal forms modulo the PRIME P.
-- Coefficient variables are constants for the relative Weyl algebra.
gIntegerRoots = f -> (
    if f==0 then error "zero indicial polynomial";
    lf := listForm f;
    if #lf==0 then return {};
    nums := apply(lf,u -> sub(numerator u#1,gCoefficientQuotient));
    dens := apply(lf,u -> sub(denominator u#1,gCoefficientQuotient));
    cleared := apply(#lf, i -> lift(nums#i * product(select(toList(0..#lf-1),j -> j!=i),j->dens#j),gCoefficientBase));
    mons := unique flatten apply(cleared,p -> apply(listForm p,u -> u#0));
    U := QQ[z];
    polynomials := apply(mons,e -> sum(#lf,i -> coefficient(gCoefficientBase_e,cleared#i)*z^(sum lf#i#0)));
    common := fold(polynomials,(x,y)->gcd(x,y));
    if common==0 then error "zero rational content in generic root test";
    factors := toList factor common;
    linear := select(factors,u -> first degree u#0==1);
    roots := apply(linear,u -> (
        c := listForm u#0;
        a := coefficient(z,u#0);
        b := sub(u#0,{z=>0});
        -lift(b,QQ)/a
        ));
    apply(select(roots,r -> denominator r==1),r -> numerator r)
    );
getIntRoots RingElement := f -> gIntegerRoots f;
