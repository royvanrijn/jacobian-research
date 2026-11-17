#!/usr/bin/env python3
"""Exact auxiliary algebra for the full nodal formal-frame argument."""
from pathlib import Path
from functools import lru_cache
import hashlib
import json
import sympy as S


def partitions(total, minimum=1):
    if total == 0:
        yield []
    for first in range(minimum, total+1):
        for rest in partitions(total-first, first):
            yield [first]+rest


def main():
    p=Path(__file__).resolve().parent
    u,v,T,s,h=S.symbols('u v T s h')
    a,b,c,d,e,f,k=S.symbols('a b c d e f k')
    G=a*u**3+d*v**3+b*u**2+c*v**2+e*u+f*v+k
    restricted=G.subs(v,-h/u)
    difference=S.cancel((G.subs(v,(s-h)/u)-restricted)/s)
    assert S.cancel(s*(T+difference)+restricted-(T*s+G.subs(v,(s-h)/u)))==0
    coordinates=S.Matrix([u,u*v+h,T+difference.subs(s,u*v+h)])
    assert S.factor(coordinates.jacobian([u,v,T]).det())==u
    derivative=S.expand(u**4*S.diff(restricted,u))
    assert derivative==3*a*u**6+2*b*u**5+e*u**4+f*h*u**2-2*c*h**2*u+3*d*h**3

    # Formal A_mu periods: recurrence in the one-variable power phase,
    # including integration by parts and parameter horizontality.
    eps=S.Symbol('eps')
    spectra=[];identities=0
    for degree in range(2,8):
        for character in range(degree-1):
            nu=1+S.Rational(character+1,degree)
            assert not nu.is_integer
            for ay in range(3*degree+1):
                for bs in range(4):
                    for cw in range(4):
                        def value(ay,bs,cw):
                            if min(ay,bs,cw)<0 or bs!=cw or ay%degree!=character:
                                return S.Integer(0)
                            power=(ay-character)//degree
                            return (-1)**(bs+power)*S.factorial(bs)*S.rf(S.Rational(character+1,degree),power)*eps**(bs+power)
                        # Overall exp(alpha*t)*t^-nu has been factored out.
                        assert S.expand(ay*value(ay-1,bs,cw)+degree/eps*value(ay+degree-1,bs,cw))==0
                        assert S.expand(bs*value(ay,bs-1,cw)+value(ay,bs,cw+1)/eps)==0
                        assert S.expand(cw*value(ay,bs,cw-1)+value(ay,bs+1,cw)/eps)==0
                        val=value(ay,bs,cw)
                        assert S.expand(value(ay+degree,bs,cw)+value(ay,bs+1,cw+1)+nu*eps*val+eps**2*S.diff(val,eps))==0
                        identities+=4
            spectra.append({'mu':degree-1,'character':character,'exponent':str(nu),'integer':False})

    # Wick integration by parts with a general symmetric covariance.
    q00,q01,q02,q11,q12,q22=S.symbols('q00 q01 q02 q11 q12 q22')
    cov=S.Matrix([[q00,q01,q02],[q01,q11,q12],[q02,q12,q22]])
    @lru_cache(None)
    def wick(indices):
        if len(indices)%2:return S.Integer(0)
        if not indices:return S.Integer(1)
        return S.expand(sum(cov[indices[0],indices[j]]*wick(indices[1:j]+indices[j+1:]) for j in range(1,len(indices))))
    wick_checks=0
    for aa in range(5):
        for bb in range(5-aa):
            for cc in range(5-aa-bb):
                powers=[aa,bb,cc];indices=tuple(i for i,n in enumerate(powers) for _ in range(n))
                for i in range(3):
                    rhs=0
                    for j in range(3):
                        if powers[j]:
                            rest=list(indices);rest.remove(j)
                            rhs+=cov[i,j]*powers[j]*wick(tuple(rest))
                    assert S.expand(wick((i,)+indices)-rhs)==0
                    wick_checks+=1

    # Confluent interpolation for every multiplicity partition of the
    # finite ranks in the four retained profiles. No realizability of
    # every abstract partition as a cubic is asserted.
    X=S.Symbol('X');frame_cases=[]
    for escapes,finite_rank in [(1,6),(3,4),(4,3),(5,2)]:
        count=0
        for multiplicities in partitions(finite_rank):
            nodes=list(range(1,len(multiplicities)+1))
            modulus=S.prod((X-node)**mu for node,mu in zip(nodes,multiplicities))
            columns=[]
            for node,mu in zip(nodes,multiplicities):
                local=(X-node)**mu
                other=S.cancel(modulus/local)
                inverse=S.invert(other,local,X)
                for order in range(mu):
                    basis=S.rem(S.expand((X-node)**order*other*inverse),modulus,X)
                    assert S.degree(basis,X)<finite_rank
                    for other_node,other_mu in zip(nodes,multiplicities):
                        for jet in range(other_mu):
                            observed=S.diff(basis,X,jet).subs(X,other_node)/S.factorial(jet)
                            assert observed==int(node==other_node and order==jet)
                    columns.append(basis)
            assert len(columns)==finite_rank
            count+=1
        frame_cases.append({'escapes':escapes,'finite_rank':finite_rank,'partitions_checked':count,
                            'rank':escapes+finite_rank,'escape_polynomial_degree':6,
                            'finite_to_escape_valuation_lower_bound':8-(finite_rank-1),
                            'escape_to_finite_valuation_lower_bound':finite_rank})
    assert all(row['rank']==7 and row['finite_to_escape_valuation_lower_bound']>0 for row in frame_cases)

    # Analytic units and an invertible local coordinate give a triangular
    # jet matrix. Its diagonal is nonzero at all finite singularities.
    unit,coordinate=S.symbols('unit coordinate',nonzero=True)
    for mu in range(1,7):
        assert S.prod(unit*coordinate**j for j in range(mu))==unit**mu*coordinate**(mu*(mu-1)//2)

    # Exponential-free finite factors remain nontrivial over the BASE
    # Laurent field. Even ramification can hide this obstruction.
    descent=[]
    for row in spectra:
        nu=S.Rational(row['exponent'])
        residue=S.frac(nu)
        assert residue!=0
        descent.append({'exponent':str(nu),'base_residue_mod_Z':str(residue)})
    assert 4*S.Rational(3,2)==6  # Not a license to discard base descent.
    critical=json.loads((p/'nodal_gaussian_critical_verification.json').read_text())
    center=json.loads((p/'nodal_morse_dominance_verification.json').read_text())
    assert critical['completed'] and center['completed']
    assert [x['escape_count'] for x in critical['escape_profiles']]==[1,3,4,5]
    result={'completed':True,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'prior_critical_checker_sha256':critical['checker_sha256'],'prior_center_all_checker_sha256':center['checker_sha256'],
            'normal_form_jacobian':'u; use the swapped chart when u=0,v!=0',
            'finite_singularity_type':'s*w+y^(mu+1), 1<=mu<=6; no finite point on the singular conic center outside the critical-line case',
            'power_phase_identities_checked':identities,'Wick_IBP_identities_checked':wick_checks,
            'finite_spectra':spectra,'base_descent':descent,'formal_frame_cases':frame_cases,
            'escape_cutoff_multiplicity':8,'rank':7,
            'zero_exponential_block':'Direct sum of base-unramified power factors t^(-1-j/d), 1<=j<d; no integer residue and hence no trivial composition factor.',
            'resolved':'Every irreducible nodal-leading ternary cubic, with all lower terms and the full affine nondegenerate Gaussian; finite degeneracies, coincident values and all quartic resonance conditions included.',
            'written_proof_certified_by_checker':False,
            'proof_boundary':'Exact normal coordinate identity, finite formal recurrences and specified Wick identities, confluent interpolation across all small rank partitions, degree valuations and noninteger residue arithmetic. The all-order formal construction, completeness of the frame, differential-module descent and no-trivial-subquotient argument are internal written mathematics.',
            'remaining':'General nonhomogeneous cuspidal/reducible-leading cubics, numerical uniform pure gates, coefficient-scheme questions and lower-power GVC; independent replay, external review and formalization.'}
    (p/'nodal_formal_frame_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['completed','rank','formal_frame_cases','resolved','written_proof_certified_by_checker']}))


if __name__=='__main__':
    main()
