\\ Exact PARI/GP certificate for the integral nonproper boundary.
\\
\\ For each generic component, the five recorded valuations are
\\   (v(d), v(t), v(x), v(y), v(z)).
\\ A source point extends over the prime exactly when its source coordinates
\\ are integral and t*q=1 remains a unit identity.

P = x^5 - 3*x^4 + 8*x^3 - 4*x^2 - 8*x + 48;

profile_data(f, p) =
{
  my(nf = nfinit(f), a = Mod(x, f), d, t, xx, yy, zz, elts, dec);
  d = -subst(deriv(P), x, a) / 8;
  t = 1 / d;
  xx = a / (2*d);
  yy = -1 + 3*a - 3*a^2/2 + 5*a^3/8;
  zz = d^2 * (d + yy^2*(1 + 3*t));
  elts = [d,t,xx,yy,zz];
  dec = idealprimedec(nf, p);
  vecsort(vector(#dec, i,
    concat(
      [dec[i].e, dec[i].f],
      vector(5, j, nfeltval(nf, elts[j], dec[i]))
    )
  ));
};

cubic = x^3 - 2*x^2 + 8;
quadratic = x^2 - x + 6;

nf_cubic = nfinit(cubic);
nf_quadratic = nfinit(quadratic);
if (nf_cubic.disc != -23, error("cubic field discriminant must be -23"));
if (nf_quadratic.disc != -23, error("quadratic field discriminant must be -23"));

\\ Rows are (e,f,v(d),v(t),v(x),v(y),v(z)).
if (profile_data(cubic,2) != [[1,3,0,0,0,0,0]], error("unexpected cubic profile at 2"));
if (profile_data(cubic,7) != [[1,1,0,0,0,0,0],[1,2,1,-1,-1,0,1]], error("unexpected cubic profile at 7"));
if (profile_data(cubic,23) != [[1,1,0,0,0,0,0],[2,1,1,-1,-1,0,1]], error("unexpected cubic profile at 23"));

if (profile_data(quadratic,2) != [[1,1,-3,3,2,-3,-12],[1,1,0,0,0,3,0]], error("unexpected quadratic profile at 2"));
if (profile_data(quadratic,7) != [[1,2,1,-1,-1,0,1]], error("unexpected quadratic profile at 7"));
if (profile_data(quadratic,23) != [[2,1,1,-1,-1,0,1]], error("unexpected quadratic profile at 23"));

\\ Identify the retained degree-one branches by their root residues.
alpha = Mod(x,cubic);
beta = Mod(x,quadratic);
cubic_7 = idealprimedec(nf_cubic,7);
cubic_23 = idealprimedec(nf_cubic,23);
quadratic_2 = idealprimedec(nf_quadratic,2);
if (nfmodpr(nf_cubic,alpha,cubic_7[1]) != 1, error("retained cubic branch at 7 must have alpha=1"));
if (nfmodpr(nf_cubic,alpha,cubic_23[2]) != 7, error("retained cubic branch at 23 must have alpha=7"));
if (nfmodpr(nf_quadratic,beta,quadratic_2[2]) != 0, error("retained quadratic branch at 2 must have beta=0"));

\\ Norm(d) shows that no odd prime outside 7 and 23 can support a pole.
a_cubic = Mod(x,cubic);
a_quadratic = Mod(x,quadratic);
d_cubic = -subst(deriv(P),x,a_cubic)/8;
d_quadratic = -subst(deriv(P),x,a_quadratic)/8;
if (nfeltnorm(nf_cubic,d_cubic) != -1127, error("unexpected cubic norm of d"));
if (nfeltnorm(nf_quadratic,d_quadratic) != 1127/8, error("unexpected quadratic norm of d"));
if (1127 != 7^2*23, error("unexpected factorization of 1127"));

print("PASS: both generic component fields have discriminant -23");
print("PASS: the complete boundary valuation profiles at 2, 7, and 23 agree");
print("PASS: no prime outside 2, 7, and 23 supports a reconstruction pole");
print("PASS: retained degree-one branches occur at s=0 mod 2, s=1 mod 7, and s=7 mod 23");
