needsPackage "BernsteinSato";
needsPackage "JSON";
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

gEncodeCoefficient = c -> (
    a := lift(sub(numerator c,gCoefficientQuotient),gCoefficientBase);
    b := lift(sub(denominator c,gCoefficientQuotient),gCoefficientBase);
    new HashTable from {
        "numerator"=>apply(listForm a,u->{u#0,toString u#1}),
        "denominator"=>apply(listForm b,u->{u#0,toString u#1})}
    );
gEncodeOperator = f -> apply(listForm f,u->{u#0,gEncodeCoefficient u#1});
gWeight = (f,w) -> max apply(listForm f,u->sum(#w,i->u#0#i*w#i));
gExponents = (n,d) -> (
    if d<0 then return {};
    if n==0 then return {{}};
    flatten apply(toList(0..d),k->apply(gExponents(n-1,d-k),e->prepend(k,e)))
    );

-- Exact Oaku integration presentation, before pruning, over the GENERIC
-- coefficient field. Public bFunction supplies the polynomial; integer roots
-- are extracted above without the library's QQ-only factor formatter.
gGenericSource = (I,n) -> (
    W := ring I;m := n+1;
    wt := join(toList(n:1),{0});ww := join(-wt,wt);
    FI := Fourier I;
    G := gbw(gens FI,ww);
    gs := first entries G;
    degrees := apply(gs,g->gWeight(g,ww));
    b := bFunction(FI,wt);
    roots := gIntegerRoots b;
    if #roots==0 or max roots<0 then error "Gaussian integration has no nonnegative root";
    cutoff := max roots;
    theta := sum(n,i->W_i*W_(m+i));
    beuler := sum(listForm b,u->sub(u#1,W)*theta^(sum u#0));
    initialG := matrix{apply(#gs,i->part(degrees#i,ww,gs#i))};
    liftInitial := matrix{{beuler}} // initialG;
    if initialG*liftInitial != matrix{{beuler}} then error "invalid initial Euler lift";
    lifted := sum(#gs,i->(
        q := liftInitial_(i,0);
        selected := sum(select(listForm q,u->sum(#ww,k->u#0#k*ww#k)==-degrees#i),u->u#1*W_(u#0));
        selected*gs#i));
    originalL := FourierInverse lifted;
    recurrenceWitness := matrix{{originalL}} // gens I;
    if gens I*recurrenceWitness != matrix{{originalL}} then error "invalid original Euler identity";
    basis := gExponents(n,cutoff);
    projected := {};
    unprojected := {};
    scan(#gs,i->scan(gExponents(n,cutoff-degrees#i),e->(
        q := product(n,k->W_(m+k)^(e#k))*gs#i;
        proj := sum(select(listForm q,u->all(toList(0..n-1),k->u#0#k==0)),u->u#1*W_(u#0));
        if proj != 0 then (
            projected=append(projected,FourierInverse proj);
            unprojected=append(unprojected,FourierInverse q);
            );
        )));
    T := (coefficientRing W)[t,dt,WeylAlgebra=>{t=>dt}];
    relations := matrix(T,apply(basis,e->apply(projected,p->sum(
        select(listForm p,u->take(u#0,n)==e),u->(
            if any(take(drop(u#0,m),n),k->k!=0) then error "spatial derivative survived projection";
            u#1*T_0^(u#0#n)*T_1^(last u#0)
            )))));
    M := cokernel relations;
    r := #basis;cols := #projected;
    annihilators := {};witnesses := {};
    scan(r,i->(
        unit := map(T^r,T^1,(j,k)->if j==i then 1_T else 0_T);
        ker := gens kernel map(M,T^1,unit);
        candidates := select(first entries ker,a->a!=0);
        if #candidates==0 then error "coordinate annihilator missing";
        A := first candidates;
        Q := (unit*matrix{{A}})//relations;
        if relations*Q != unit*matrix{{A}} then error "invalid coordinate lift";
        annihilators=append(annihilators,gEncodeOperator A);
        witnesses=append(witnesses,apply(cols,j->gEncodeOperator Q_(j,0)));
        ));
    new HashTable from {
        "rank"=>r,"relation_columns"=>cols,"integration_cutoff"=>cutoff,
        "integration_b_terms"=>apply(listForm b,u->{sum u#0,gEncodeCoefficient u#1}),
        "integration_integer_roots"=>roots,
        "amplitude_exponents"=>basis,
        "relations"=>apply(r,i->apply(cols,j->gEncodeOperator relations_(i,j))),
        "unprojected_relation_operators"=>apply(unprojected,gEncodeOperator),
        "original_recurrence_operator"=>gEncodeOperator originalL,
        "original_recurrence_left_witnesses"=>apply(numgens I,i->gEncodeOperator recurrenceWitness_(i,0)),
        "coordinate_annihilators"=>annihilators,
        "coordinate_annihilator_left_witnesses"=>witnesses}
    );

