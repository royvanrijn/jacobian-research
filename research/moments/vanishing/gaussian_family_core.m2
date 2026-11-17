needsPackage "BernsteinSato";
needsPackage "JSON";

gQuotient = M -> cokernel presentation prune M;
gIsZero = M -> numgens prune M == 0;
gSubset = (I,J) -> all(first entries gens I, f -> f % (gb J) == 0);
gSerializeIdeal = I -> apply(first entries gens gb I, f -> apply(listForm f, u -> {u#0,toString u#1}));

-- Keep the associated graded MODULE, not only its annihilator: conormal
-- multiplicity is the connection rank, and annihilators lose multiplicity.
gCharacteristicData = M -> (
    W := ring M;
    createDpairs W;
    w := apply(numgens W, i -> if member(i,W.dpairInds#1) then 1 else 0);
    G := inw(presentation M,w);
    C := ring G;
    A := (coefficientRing W)(monoid [W.dpairVars#0]);
    toChar := map(C,A,matrix{apply(W.dpairInds#0, i -> C_i)});
    GM := cokernel G;
    ps := if gIsZero GM then {} else minimalPrimes ann GM;
    bs := apply(ps, P -> preimage(toChar,P));
    (G,GM,C,A,ps,bs)
    );

gGenericRank = (M,S) -> (
    if gIsZero M then return (0,new HashTable from {"numerator"=>0,"denominator"=>1});
    d := gCharacteristicData M;
    G := d#0;GM := d#1;C := d#2;A := d#3;ps := d#4;bs := d#5;
    toA := map(A,ring S,vars A);
    SA := toA S;
    chosen := positions(bs, Z -> Z == SA);
    if #chosen == 0 then return (0,new HashTable from {"numerator"=>0,"denominator"=>1});
    if #chosen != 1 then error "conormal component is not unique";
    k := first chosen;P := ps#k;
    f := 1_C;
    scan(#ps, j -> if j != k then (
        cs := select(first entries gens ps#j, q -> q % (gb P) != 0);
        if #cs == 0 then error "failed to isolate conormal component";
        f = f * first cs;
        ));
    N := if f == 1_C then image G else saturate(image G,ideal f);
    numerator := degree cokernel gens N;
    denominator := degree(C^1/P);
    if numerator % denominator != 0 then error "nonintegral conormal multiplicity";
    (numerator // denominator,new HashTable from {
        "numerator"=>numerator,"denominator"=>denominator,
        "conormal_ideal"=>gSerializeIdeal P})
    );

gSourceBad = (M,S) -> (
    B := ring S;
    d := gCharacteristicData M;
    A := d#3;bs := d#5;
    toA := map(A,B,matrix{toList apply(1..numgens B, i -> A_i)});
    SA := toA S;
    bad := {};
    scan(bs, Z -> if Z != SA then (
        J := preimage(toA,Z);
        if J == S then (
            critical := Z + ideal diff(A_0,gens Z);
            J = preimage(toA,critical);
            );
        if J == S then error "no proper noncharacteristic exceptional set";
        if not gSubset(S,J) then error "source exception outside stratum";
        bad = append(bad,J);
        ));
    bad
    );

gConnectionBad = (M,S) -> (
    if gIsZero M then return {};
    B := ring S;
    d := gCharacteristicData M;
    A := d#3;bs := d#5;
    toA := map(A,B,vars A);
    SA := toA S;
    bad := {};
    scan(bs, Z -> if Z != SA then (
        J := preimage(toA,Z);
        if not gSubset(S,J) then error "connection exception outside stratum";
        bad = append(bad,J);
        ));
    bad
    );

-- Input: the EXACT jointly holonomic Gaussian source on (t,parameters).
-- Reuse the ORIGINAL source at every node, not iterated local cohomology.
gFamilySafeIdeal = (M,B,baseIdeal) -> (
    M = gQuotient M;
    W := ring M;
    createDpairs W;
    rr := numgens B;
    if #W.dpairVars#0 != rr+1 then error "source must retain t and every parameter";
    toW := map(W,B,matrix{toList apply(1..rr, i -> W_i)});
    pending := minimalPrimes baseIdeal;
    seen := {};
    safe := {};
    nodes := {};
    while #pending > 0 do (
        S := first pending;pending = drop(pending,1);
        if any(seen,T -> S == T) then continue;
        seen = append(seen,S);
        codim := rr-dim(B/S);
        MS := if codim == 0 then M else gQuotient((localCohom({codim},toW S,M,Strategy=>OaTa))#codim);
        DM := gQuotient Ddual MS;
        NI := DintegrationAll(DM,join({1},toList(rr:0)));
        NC := NI#HomologyModules;
        N0 := gQuotient NC#0;N1 := gQuotient NC#1;
        rankData := gGenericRank(N1,S);
        hrank := rankData#0;
        if hrank != 0 and hrank != 1 then error "Gaussian formal Hom rank exceeds one";
        bads := gSourceBad(MS,S) | gConnectionBad(N0,S) | gConnectionBad(N1,S);
        sing := if codim == 0 then ideal(1_B) else S+minors(codim,jacobian S);
        bads = append(bads,sing);
        bad := intersect bads;
        if bad == S then error "exceptional locus failed to lower dimension";
        children := if hrank == 1 or bad == ideal(1_B) then {} else minimalPrimes bad;
        node := new HashTable from {
            "stratum_closure_ideal"=>gSerializeIdeal S,
            "dimension"=>dim(B/S),"codimension"=>codim,
            "generic_polynomial_Hom_dimension"=>hrank,
            "conormal_multiplicity"=>rankData#1,
            "exceptional_ideal"=>gSerializeIdeal bad,
            "children"=>apply(children,gSerializeIdeal),
            "ambient_restricted_source_presentation"=>toString presentation MS,
            "dual_time_direct_image_H0"=>toString presentation N0,
            "dual_time_direct_image_Hminus1"=>toString presentation N1,
            "whole_closure_safe"=>(hrank==1)};
        nodes = append(nodes,node);
        stdio << "NODE_JSON|" << toJSON node << endl;
        if hrank == 1 then safe = append(safe,S)
        else pending = pending | children;
        );
    J := if #safe == 0 then ideal(1_B) else intersect safe;
    A := B/baseIdeal;
    JA := sub(J,A);
    normalization := if JA == ideal(0_A) then ideal(1_A) else saturate(ideal(0_A),JA);
    torsionPower := JA;
    torsionExponent := 1;
    while normalization*torsionPower != ideal(0_A) do (
        torsionPower = torsionPower*JA;
        torsionExponent = torsionExponent+1;
        );
    if (normalization : JA) != normalization then error "normalization ideal is not colon-stable";
    normalizationUp := baseIdeal + lift(normalization,B);
    new HashTable from {"safe_radical_ideal"=>gSerializeIdeal J,"nodes"=>nodes,
        "safe_closed_components"=>apply(safe,gSerializeIdeal),"status"=>"DECIDED",
        "parameter_ideal"=>gSerializeIdeal baseIdeal,
        "polynomial_Hom_normalization_ideal_lift"=>gSerializeIdeal normalizationUp,
        "normalization_torsion_exponent"=>torsionExponent,
        "normalization_colon_stable"=>true,
        "coefficient_ring_phase_safe"=>(normalization == ideal(1_A))}
    );
