\\ Exact modular-form and finite-matrix data for the level-97 threefold.
check(condition, description) = {
    if(!condition,
        print("Certificate mismatch: ", description);
        quit(1)
    )
};

projective_order(p, trace, determinant) = {
    my(matrix = Mod(1, p) * [0, -determinant; 1, trace]);
    my(power = matid(2));
    for(n = 1, p + 1,
        power *= matrix;
        if(power[1, 2] == 0 && power[2, 1] == 0 && power[1, 1] == power[2, 2],
            return(n)
        )
    );
    error("projective order exceeds bound")
};

{
    level = 97;
    space = mfinit([level, 2], 0);
    cubic_forms = select(form -> poldegree(mfparams(form)[4]) == 3, mfeigenbasis(space));
    check(#cubic_forms == 1, "unique cubic orbit");
    form = cubic_forms[1];
    coefficient_polynomial = y^3 - y^2 - 2*y + 1;
    check(mfparams(form)[4] == coefficient_polynomial, "coefficient field");
    check(poldisc(coefficient_polynomial) == 49 && nfdisc(coefficient_polynomial) == 49,
        "maximal order in the chosen basis");

    coefficients = mfcoefs(form, 1709);
    check(coefficients[2] == 1, "normalized form");
    check(lift(coefficients[3]) == -y - 1, "T_2 generates the full coefficient order");
    hecke_polynomial = charpoly(mfheckemat(space, 2));
    check(mfdim(space) == 7, "dimension of the full newspace");
    check(hecke_polynomial == (x^3 + 4*x^2 + 3*x - 1) * (x^4 - 3*x^3 - x^2 + 6*x - 1),
        "full T_2 characteristic polynomial");
    check(poldisc(hecke_polynomial) == 2^12 * 7^2 * 2777, "T_2 discriminant");
    check(subst(hecke_polynomial, x, 3) == 568, "Eisenstein T_2 eigenvalue");

    examples = [
        [113, -3*y^2 - y + 8, -113, -76433539],
        [1583, 12*y^2 - 23*y - 26, 20579, -168239088703],
        [1709, -4*y^2 + 10*y + 13, -1709, -290525486531]
    ];
    for(i = 1, #examples,
        p = examples[i][1];
        check(isprime(p) && p != level, "good rational prime");
        ap = coefficients[p + 1];
        check(lift(ap) == examples[i][2], "relative Frobenius trace");
        check(norm(ap) == examples[i][3] && valuation(norm(ap), p) == 1,
            "valuation-one trace norm");
        check(norm(ap^2 - 4*p) == examples[i][4], "relative discriminant norm");
        roots = lift(polrootsmod(coefficient_polynomial, p));
        check(#roots == 3, "split coefficient field");
        traces = vector(3, j, lift(Mod(1, p) * subst(lift(ap), y, roots[j])));
        check(vecsum(vector(3, j, traces[j] == 0)) == 1, "exactly one failed coefficient pair");
        check(poldisc(hecke_polynomial) % p != 0 && 568 % p != 0, "no same-level T_2 congruence");
        print([p, roots, traces, norm(ap), norm(ap^2 - 4*p)])
    );
    check(examples[1][4] % examples[2][1] != 0, "distinct CM fields by ramification at 1583");

    p = 113;
    roots = [104, 89, 34];
    quotients = [2, 7, 11];
    expected_traces = [[8, 23, 78], [22, 77, 7], [7, 22, 77]];
    expected_characters = [[1, 1, 1], [1, 1, -1], [-1, -1, 1]];
    expected_orders = [[56, 56, 56], [8, 28, 57], [57, 57, 28]];
    check(vector(3, j, lift(Mod(1, p) * subst(coefficient_polynomial, y, roots[j]))) == [0, 0, 0],
        "three labelled embeddings for 113");
    for(i = 1, #quotients,
        q = quotients[i];
        aq = lift(coefficients[q + 1]);
        traces = vector(3, j, lift(Mod(1, p) * subst(aq, y, roots[j])));
        characters = vector(3, j, kronecker(traces[j]^2 - 4*q, p));
        orders = vector(3, j, projective_order(p, traces[j], q));
        check(traces == expected_traces[i] && characters == expected_characters[i]
            && orders == expected_orders[i], "finite matrix witness");
        print([q, traces, characters, orders])
    );
    projective_traces = vector(3, j, lift(Mod(expected_traces[1][j], p)^2 / 2));
    check(projective_traces == [32, 95, 104], "T_2 separates projective factors");

    remaining_images = [
        [1583, [540, 1257, 1370], [1042, 325, 212], 791, [1496, 1365, 310],
            [[3, 212, 792], [7, 1041, 176], [5, 716, 528]]],
        [1709, [127, 620, 963], [1581, 1088, 745], 1708, [1356, 558, 1509],
            [[11, 1580, 285], [5, 835, 855], [3, 1581, 1710]]]
    ];
    for(i = 1, #remaining_images,
        data = remaining_images[i];
        p = data[1];
        roots = data[2];
        check(vector(3, j, lift(Mod(1, p) * subst(coefficient_polynomial, y, roots[j]))) == [0, 0, 0],
            "three labelled embeddings for the remaining prime");
        a2 = lift(coefficients[3]);
        traces = vector(3, j, lift(Mod(1, p) * subst(a2, y, roots[j])));
        check(traces == data[3], "remaining split trace residues");
        check(vector(3, j, kronecker(traces[j]^2 - 8, p)) == [1, 1, 1], "remaining split witnesses");
        check(vector(3, j, projective_order(p, traces[j], 2)) == vector(3, j, data[4]),
            "remaining split projective orders");
        projective_traces = vector(3, j, lift(Mod(traces[j], p)^2 / 2));
        check(projective_traces == data[5] && #Set(projective_traces) == 3,
            "remaining T_2 witnesses separate the factors");
        print([p, 2, traces, data[4], projective_traces]);
        for(j = 1, 3,
            witness = data[6][j];
            q = witness[1];
            aq = lift(coefficients[q + 1]);
            trace_residue = lift(Mod(1, p) * subst(aq, y, roots[j]));
            check(trace_residue == witness[2] && kronecker(trace_residue^2 - 4*q, p) == -1,
                "remaining nonsplit witness");
            check(projective_order(p, trace_residue, q) == witness[3], "remaining nonsplit projective order");
            print([p, roots[j], q, trace_residue, witness[3]])
        )
    );
    print("PARI/GP ", version(), " verified all level-97 certificates")
}
