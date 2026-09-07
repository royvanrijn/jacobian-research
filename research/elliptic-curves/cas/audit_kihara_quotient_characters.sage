#!/usr/bin/env sage-python
"""Eight fixed good-prime quotient traces on each of five retained parents."""
import hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,prime_range,kronecker
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/kihara-quotient-characters-v1'
def main():
    protocol=json.loads((D/'protocol.json').read_text())
    for name,h in protocol['sources'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h
    parents=json.loads((ART/'kihara_five_parent_distinctness_v1.json').read_text())['rows']
    involution=json.loads((ART/'kihara_section_involution_v1.json').read_text())
    action=matrix(QQ,involution['involution_action'])[:11,:11];assert action**2==matrix.identity(QQ,11) and 11-(action-matrix.identity(QQ,11)).rank()==6
    rows=[]
    for parent in parents:
        R=PolynomialRing(QQ,'T');T=R.gen();K=R.fraction_field();A=R(parent['raw_A']);B=R(parent['raw_B']);E=EllipticCurve(K,[A,B])
        P=[E([K(c) for c in v]) for v in parent['generic_sections'][:11]]
        for i in range(11):assert E([c(-T) for c in P[i].xy()])==sum((ZZ(action[j,i])*P[j] for j in range(11)),E(0))
        a=R([A[2*i] for i in range(5)]);b=R([B[2*i] for i in range(7)]);disc=-16*(4*a**3+27*b**2)
        assert A(-T)==A and B(-T)==B and (a.degree(),b.degree(),disc.degree())==(4,6,10) and disc.gcd(disc.derivative())==1 and disc.gcd(a)==1
        counts=[];rejected=[]
        for prime in prime_range(31,200):
            p=int(prime);F=GF(p);S=PolynomialRing(F,'s')
            try:aa,bb=S(a),S(b)
            except (ValueError,ZeroDivisionError):rejected.append(p);continue
            dd=-16*(4*aa**3+27*bb**2)
            if (aa.degree(),bb.degree(),dd.degree())!=(4,6,10) or dd.gcd(dd.derivative())!=1 or dd.gcd(aa)!=1:rejected.append(p);continue
            node=-3*bb[6]/(2*aa[4])
            if not (3*node).is_square():rejected.append(p);continue
            chars=[0]+[-1]*(p-1)
            for x in range(1,p):chars[x*x%p]=1
            fibres=[]
            for i in range(p+1):
                av,bv=(aa[4],bb[6]) if i==p else (aa(F(i)),bb(F(i)))
                n=p+1+sum(chars[(x*x*x+int(av)*x+int(bv))%p] for x in range(p))
                if 4*av**3+27*bv**2:assert n==EllipticCurve(F,[av,bv]).cardinality()
                fibres.append(n)
            total=sum(fibres)+p;sign=QQ(total-1-p*p-9*p)/p;assert sign in (-1,1)
            counts.append({'prime':p,'fibre_counts':fibres,'infinity_resolution_correction':p,'total':total,'extra_divisor_sign':int(sign),'minus_three_symbol':int(kronecker(-3,p))})
            if len(counts)==8:break
        assert len(counts)==8
        row={'path_parameter':parent['path_parameter'],'quotient_A':list(map(str,a.list())),'quotient_B':list(map(str,b.list())),
             'invariant_visible_rank':6,'geometric_MW_rank':7,'counts':counts,'rejected_primes':rejected,
             'rational_MW_rank':6 if any(c['extra_divisor_sign']==-1 for c in counts) else 'UNKNOWN_6_OR7',
             'quadratic_field':'UNKNOWN','minus_three_agrees_at_all_eight':all(c['extra_divisor_sign']==c['minus_three_symbol'] for c in counts)}
        rows.append(row)
        print('PARENT',parent['path_parameter'],'signs',[(c['prime'],c['extra_divisor_sign']) for c in counts],'minus3',row['minus_three_agrees_at_all_eight'],flush=True)
    out={'schema':'kihara-quotient-characters.v1','status':'PASS','rows':rows,'sources':protocol['sources'],
         'scope':'Five retained Kihara parents, no new parameter. Exact visible-section involution supplies six rational quotient directions; rational elliptic quotient geometry gives geometric rank7. A negative trace excludes rationality of the remaining divisor. Finite character agreement is not identification of its quadratic field or a theorem over the parent parameter. No point search or new curve.'}
    with (D/'result.json').open('x') as f:json.dump(out,f,indent=2);f.write('\n')
if __name__=='__main__':main()
