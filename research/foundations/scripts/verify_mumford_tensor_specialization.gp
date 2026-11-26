\p 30

check_condition(value, message) = {
  if(!value, error(message));
};

verify_specialization() = {
  /* C_{7,2} of Bouchet, Hanselman, Pieper and Schiavone, (1.1). */
  my(curve = 2*(5/16*x^10 - 567/64*x^9 - 189/2*x^8 - 861/2*x^7
                - 756*x^6 - 378*x^5 + 672*x^4 + 864*x^3 - 448*x));
  my(primes = [13, 29]);
  my(expected = [
    x^8+x^7+15*x^6-41*x^5+69*x^4-533*x^3+2535*x^2+2197*x+28561,
    x^8-x^7-6*x^6+44*x^5-1080*x^4+1276*x^3-5046*x^2-24389*x+707281]);
  my(discriminants = [139546236594961, 234444145241761]);
  for(i = 1, #primes,
    my(p = primes[i], reduction = Mod(1, p)*curve);
    check_condition(poldegree(reduction) == 10, "degree at good prime");
    check_condition(poldegree(gcd(reduction, deriv(reduction))) == 0,
                    "squarefree reduction");

    my(P = hyperellcharpoly(reduction));
    check_condition(P == expected[i] && polisirreducible(P),
                    "irreducible Frobenius polynomial");
    check_condition(polcoeff(P, 4)%p != 0, "ordinary reduction");
    check_condition(nfdisc(P) == discriminants[i], "field discriminant");

    my(degree = poldegree(nfsplitting(P)), pi = Mod(x, P));
    check_condition(degree == 24, "splitting field degree");
    /* Here phi(n) divides degree, and phi(n)^2 >= n/2 bounds the search. */
    for(n = 2, 2*degree^2,
      if(degree%eulerphi(n) == 0,
        check_condition(polisirreducible(charpoly(pi^n, y)),
                        Str("Frobenius power loses degree at ", [p, n]));
      );
    );
    print("Verified absolutely simple ordinary reduction at ", p,
          " with field discriminant ", discriminants[i], ".");
  );
  check_condition(discriminants[1] != discriminants[2], "distinct fields");
};

iferr(verify_specialization(), E, print(E); quit(1));
quit(0);
