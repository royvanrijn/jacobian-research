\\ Exact replay for the Kummer orbit polynomials and primitive formal points.

verify_equal(actual, expected, label) = {
  if(actual != expected, error(label, ": got ", actual, ", expected ", expected))
};

eisenstein_different(f, p, translation, normalized) = {
  my(g = subst(f, x, x + translation), degree = poldegree(f));
  verify_equal(pollead(g), 1, "monic polynomial");
  verify_equal(valuation(polcoeff(g, 0), p), 1, "constant valuation");
  for(i = 1, degree - 1,
    if(polcoeff(g, i) % p != 0, error("nonintegral Eisenstein assertion at ", i))
  );
  verify_equal(valuation(poldisc(f), p), degree * normalized, "different");
  degree * normalized
};

verify_ordinary(p, n, s) = {
  my(t = n - s, D = p^(n-1)*(p+1), q0 = 1+p, q = q0^(p^s));
  my(f, previous, numerator, different = 0, degree = 1);
  my(A = n - 1/(p-1) + 2*p^(1-n)/(p^2-1));
  if(t > 0,
    different += eisenstein_different(x^(p^t)-q0, p, 1, t),
    q0 = 1
  );
  degree += p^t;
  previous = x^(p^t)-q0;
  for(j = 1, s,
    numerator = x^(p^(t+j))-q0^(p^j);
    f = numerator / previous;
    different += eisenstein_different(f, p, 1, t+j-1/(p-1));
    degree += poldegree(f);
    previous = numerator
  );
  for(k = 1, n-1,
    my(m = n-k);
    f = (x^(p^m)-q^(p^k))/(x^(p^(m-1))-q^(p^(k-1)));
    different += eisenstein_different(f, p, 1, m-1/(p-1));
    degree += poldegree(f)
  );
  verify_equal(degree, D, "total degree");
  verify_equal(different / D, A - 2*(p^n-p^t)/((p-1)*D), "line exponent");
};

verify_formal(p, n) = {
  my(predecessor = x, q = p^2, primitive, degree, relative);
  for(j = 1, n,
    primitive = p + predecessor^(q-1);
    degree = poldegree(primitive);
    verify_equal(degree, (q-1)*q^(j-1), "primitive degree");
    eisenstein_different(primitive, p, 0, j-1/(q-1));
    relative = sum(a = 2, p^j-1, if(a % p, p^(2*valuation(a-1,p)), 0));
    verify_equal(relative, p^(j-1)*(p^j-2), "unit-indexed different sum");
    verify_equal(j-1/(q-1)-relative/degree,
                 j-1/(p-1)+2*p^(1-j)/(q-1), "supersingular line exponent");
    predecessor *= primitive
  )
};

verify_all() = {
  foreach([3, 5, 7], p, for(n = 1, 3, for(s = 0, n, verify_ordinary(p, n, s))));
  verify_formal(3, 2);
  verify_formal(5, 1);
  print("Exact Kummer orbit and primitive formal-point differents verified");
};

iferr(verify_all(), failure, print(failure); quit(1));
