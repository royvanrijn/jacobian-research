"""Exact bounded controls; unbounded identities are proved in the TeX."""
from pathlib import Path
from math import factorial,comb
import json,hashlib
HERE=Path(__file__).resolve().parent
tests=[]
def check(name,ok):
 if not ok:raise AssertionError(name)
 tests.append(name)
def bc(n,k):return comb(n,k)if 0<=k<=n else 0
def units(p,n,mod=None):
 out=[1];u=1
 for k in range(n):
  for i in range(1,p):
   u*=p*k+i
   if mod:u%=mod
  out.append(u)
 return out
def W(p,n,U):
 return [U[k]*bc(n+(p-1)*k-1,p*k)for k in range(n)]
def H(p,n,U):
 return [p]+[U[k]*(p*bc(n+(p-1)*k-1,p*k)
                   +bc(n+(p-1)*k-1,p*k-1))for k in range(1,n+1)]
for p in [2,3,5,7,11]:
 U=units(p,22)
 for n in range(1,21):
  h=H(p,n,U);w=W(p,n,U);upper=W(p,n+1,U)
  for k in range(n+1):
   direct=p if k==0 else comb(n,k)*factorial(n+(p-1)*k-1)//factorial(n-1)//p**(k-1)
   check("actual_H_coefficient",h[k]==direct)
   lower=w[k]if k<n else 0
   check("primitive_W_coefficient",(n-k)*h[k]==p*n*lower)
   check("exact_contiguous_pair",h[k]==upper[k]+(p-1)*lower)
  check("W_constant",w[0]==1)
for p in [2,3,5,7,11]:
 for K in range(1,5):
  modulus=p**K;L=p**(K-1);U=units(p,2*L+4,modulus)
  multiplier = 1 if p == 2 and K >= 3 else -1
  for k in range(L+5):
   check("unit_factorial_period_multiplier", U[k+L] == multiplier*U[k]%modulus)
  if p != 2:
   check("reject_missing_antiperiod_sign", U[L] != U[0])
# The odd-prime restriction matters.
check("two_adic_antiperiod_failure",units(2,4,8)[4]!=-1%8)
# Verify rational series after clearing the denominator as a bivariate identity.
for p in [2,3,5,7]:
 for K in range(1,4):
  modulus=p**K;L=p**(K-1);mmax=18;U=units(p,mmax+2,modulus)
  multiplier = 1 if p == 2 and K >= 3 else -1
  G=[W(p,m+1,U)for m in range(mmax+1)]
  for m in range(mmax+1):
   lhs=[0]*(m+1)
   for i in range(m+1):
    coeff=(-1)**i*bc(p*L,i)
    for k,val in enumerate(G[m-i]):lhs[k]+=coeff*val
   if m>=L:
    for k,val in enumerate(G[m-L]):lhs[k+L]-=multiplier*val
   rhs=[0]*(m+1)
   for j in range(min(L-1,m)+1):
    rhs[j]=U[j]*(-1)**(m-j)*bc(p*(L-j)-1,m-j)
   for a,b in zip(lhs,rhs):
    check("allprecision_rational_GF",(a-b)%modulus==0)
for p in [5,7,11,13]:
 mod=p**3;w=factorial(p-1)%mod
 S1=sum(pow(i,-1,mod)for i in range(1,p))
 S2=sum(pow(i,-1,mod)**2 for i in range(1,p))
 check("harmonic_one_cubic",S1%(p*p)==0)
 check("harmonic_two_residue",S2%p==0)
 for ell in range(2*p+4):
  block=1
  for i in range(1,p):block=block*(p*ell+i)%mod
  check("uniform_block_cubic",block==w)
  check("reject_modsquare_only_control",(block-(w+p*p))%mod!=0)
 U=units(p,2*p*p+4,mod)
 for k,u in enumerate(U):
  check("unit_factorial_cubic",u==pow(w,k,mod))
 # Cleared cubic GF coefficients for the ENTIRE primitive pair.
 for n in range(1,19):
  actualW=W(p,n,U)
  # K/primitivetrace coefficients replace U[k] by w^k.
  for k,val in enumerate(actualW):
   check("primitive_pair_cubic_coefficient",
         (val-pow(w,k,mod)*bc(n+(p-1)*k-1,p*k))%mod==0)
  actualH=H(p,n,U)
  for k,val in enumerate(actualH):
   trace=p if k==0 else pow(w,k,mod)*(p*bc(n+(p-1)*k-1,p*k)+bc(n+(p-1)*k-1,p*k-1))
   check("physical_trace_cubic_coefficient",(val-trace)%mod==0)
check("three_adic_cubic_counterexample",units(3,2)[2]==40 and (40-2**2)%27!=0)

def centered_coefficients(coefficients, center, degree, modulus):
    """Reduce the Taylor coefficients without constructing a translated polynomial."""
    return [
        sum(
            coefficients[k] * comb(k, j) * pow(center, k-j, modulus)
            for k in range(j, len(coefficients))
        ) % modulus
        for j in range(degree+1)
    ]


def trace_in_jet_quotient(p, precision, degree, center, indices):
    """Use the rank-p reciprocal trace algebra over the truncated Taylor ring."""
    modulus = p**precision
    width = degree+1
    zero = (0,)*width

    def constant(value):
        return (value % modulus,) + zero[1:]

    def add(left, right):
        return tuple((a+b) % modulus for a, b in zip(left, right))

    def scaled(vector, scalar):
        return tuple(scalar*a % modulus for a in vector)

    def multiply(left, right):
        return tuple(
            sum(left[j]*right[i-j] for j in range(i+1)) % modulus
            for i in range(width)
        )

    reciprocal = [constant(comb(p, k)*(-1)**(p-k)) for k in range(p)]
    weighted_center = constant(-factorial(p-1)*center)
    if degree:
        weighted_center = (weighted_center[0], -factorial(p-1) % modulus) + zero[2:]
    reciprocal[p-1] = add(reciprocal[p-1], weighted_center)

    def product(left, right):
        coefficients = [zero]*(2*p-1)
        for i, a in enumerate(left):
            for j, b in enumerate(right):
                coefficients[i+j] = add(coefficients[i+j], multiply(a, b))
        for j in range(2*p-2, p-1, -1):
            for k, coefficient in enumerate(reciprocal):
                index = j-p+k
                coefficients[index] = add(
                    coefficients[index], scaled(multiply(coefficients[j], coefficient), -1)
                )
        return coefficients[:p]

    identity = [constant(1)] + [zero]*(p-1)
    generator = [zero, constant(1)] + [zero]*(p-2)

    def power(exponent):
        result, base = identity, generator
        while exponent:
            if exponent & 1:
                result = product(result, base)
            exponent //= 2
            if exponent:
                base = product(base, base)
        return result

    traces = [constant(p)]
    for j in range(1, p):
        trace = scaled(reciprocal[p-j], -j)
        for i in range(1, j):
            trace = add(trace, scaled(multiply(reciprocal[p-i], traces[j-i]), -1))
        traces.append(trace)

    values = []
    for index in indices:
        trace = zero
        for coefficient, basis_trace in zip(power(index), traces):
            trace = add(trace, multiply(coefficient, basis_trace))
        values.append(trace)
    return values, power


expected = [98,147,245,98,196,336,175,252,252,238,245,315,49,294,250]
actual = centered_coefficients(H(7,145,units(7,145,343)), 4, 14, 343)
check("exact_145_local_coefficients", actual == expected)
check("exact_145_lower_Newton_edges",
      actual[0] % 49 == 0 and actual[0] != 0
      and all(a % 49 == 0 for a in actual[1:5])
      and actual[5] % 7 == 0 and actual[5] % 49 != 0
      and all(a % 7 == 0 for a in actual[6:14]) and actual[14] % 7 != 0)
period = (7**3-1)*7**2
for precision, degree in [(3,0), (2,5), (1,14)]:
    values, power = trace_in_jet_quotient(7, precision, degree, 4, [145,145+period])
    check("reciprocal_quotient_period_identity", power(period) == power(0))
    target = tuple(a % (7**precision) for a in expected[:degree+1])
    check("trace_equals_direct_physical_local_data", values[0] == target)
    check("period_translates_exact_local_data", values[1] == target)


pins = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in [Path(__file__), HERE/"PHYSICAL_PRIME_MOMENT_PAIR_RECURRENCES.tex"]}
receipt = {
    "status": "PASS",
    "controls": len(tests),
    "source_sha256": pins,
    "scope": "Exact coefficient and finite quotient controls. The written proof supplies the prime, index and progression quantifiers.",
    "local_progression": {"prime":7, "center":4, "offset":145, "period":period, "negative_slope_widths":[5,9]},
    "independent_whole_proof_replay": False,
    "formal_verification": False,
    "external_review": False,
}
(HERE/"physical_prime_moment_pair_recurrences_verification.json").write_text(json.dumps(receipt, indent=2)+"\n")
print(json.dumps({"status":"PASS", "controls":len(tests)}))
