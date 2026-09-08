"""Exact geometric Picard19 via two-prime point counts and Artin--Tate.

One worker,120seconds. Verify19 rational divisor directions directly from
the17 section heights, count at149,149^2,151,151^2, and reconstruct the full
Frobenius polynomials by orthogonality of the residual dimension3 space.
No controlled-reduction backend is executed or required for the proof.
"""
import argparse,json,runpy,signal
from hashlib import sha256
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,vector

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'artifacts/generated-results/elliptic-curves'
PARENT=BASE/'curve302_recovered_mw17_parent_v1.json'
LOADER=ROOT/'elliptic-curves/cas/load_curve302_recovered_parent.sage'
OUT=BASE/'curve302_parent_geometric_picard19_v1.json'


def build():
    data=json.loads(PARENT.read_text())
    E,points,_=runpy.run_path(str(LOADER))['load_curve302_recovered_parent']()
    K=E.base_ring();R=K.ring();A,B=R(-E.c4()/48),R(-E.c6()/864)
    delta=-16*(4*A**3+27*B**2)
    assert (A.degree(),B.degree(),delta.degree())==(8,12,24)
    assert delta.is_squarefree() and A.gcd(delta).degree()==0
    J=EllipticCurve(K,[A,B]);b2=E.b2()
    points=[J(P[0]+b2/12,P[1]+(E.a1()*P[0]+E.a3())/2) for P in points]
    def h(x):return max(x.numerator().degree(),x.denominator().degree()+4)
    H=matrix(QQ,17)
    for i,P in enumerate(points):
        H[i,i]=h(P[0])
        for j in range(i):
            Q=points[j];assert P[0]!=Q[0]
            xd=((P[1]+Q[1])/(P[0]-Q[0]))**2-P[0]-Q[0]
            H[i,j]=H[j,i]=(H[i,i]+H[j,j]-h(xd))/2
    assert H==matrix(QQ,data['generic_height_gram']) and H.is_positive_definite() and H.det()==1092
    # A rank17 positive section lattice plus U gives19 rational divisor
    # classes. All finite geometric fibres are I1, so their span is primitive
    # after the only possible index2 enlargement is excluded below.
    parity=H.change_ring(GF(2)).right_kernel_matrix();assert parity.nrows()==1
    half=vector(QQ,list(map(ZZ,parity.row(0))))/2
    assert half*H*half==45
    records=[]
    for p in [149,151]:
        Rp=PolynomialRing(GF(p),'t');a,b=Rp(A),Rp(B);d=-16*(4*a**3+27*b**2)
        assert d.degree()==24 and d.is_squarefree() and a.gcd(d).degree()==0
        # Semistable24-I1 model and smooth infinity give good K3 reduction.
        counts=[]
        for n in [1,2]:
            F=GF(p**n,'v');q=F.order();Rn=PolynomialRing(F,'t');an,bn=Rn(A),Rn(B)
            total=ZZ(0);singular=0
            for t in list(F)+[None]:
                aa,bb=(an(t),bn(t)) if t is not None else (an[8],bn[12])
                if 4*aa**3+27*bb**2:
                    total+=q+1-EllipticCurve(F,[aa,bb]).cardinality()
                else:
                    assert aa
                    node=-3*bb/(2*aa)
                    total+=1 if (3*node).is_square() else -1
                    singular+=1
            N=q*q+1+2*q-total
            counts.append({'extension_degree':n,'surface_count':int(N),'singular_fibres':singular,
                           'field_modulus':list(map(int,F.modulus().list()))})
            print('COUNT',p,n,int(N),flush=True)
        s1=ZZ(counts[0]['surface_count'])-1-p*p-19*p
        s2=ZZ(counts[1]['surface_count'])-1-p**4-19*p*p
        # The residual orthogonal dimension3 representation has eigenvalues
        # epsilon*p, alpha,beta with alpha*beta=p^2. Hence
        # s2=s1^2-2*epsilon*p*s1. This determines the entire residual cubic.
        signs=[e for e in [-1,1] if s2==s1*s1-2*e*p*s1]
        assert len(signs)==1
        epsilon=signs[0];trace=s1-epsilon*p
        assert epsilon==-1 and abs(trace)<2*p and QQ(trace)/p not in ZZ
        S=PolynomialRing(QQ,'T');T=S.gen()
        full=(T-p)**19*(T-epsilon*p)*(T*T-trace*T+p*p)
        assert full.degree()==22
        if p==149:
            assert full==(T-p)**2*S(data['frobenius']['primitive_polynomial'])
        # Over Fp^2 all20 root-of-unity eigenvalues are p^2. Artin--Tate
        # gives disc(NS) modulo rational squares as trace^2-4p^2.
        disc=trace*trace-4*p*p
        records.append({'prime':p,'counts':counts,'residual_trace_1':int(s1),'residual_trace_2':int(s2),
                        'epsilon':epsilon,'nonalgebraic_pair_trace':int(trace),
                        'full_H2_polynomial':list(map(str,full.list())),
                        'geometric_reduction_Picard_rank':20,'Artin_Tate_discriminant_representative':int(disc)})
    assert [r['nonalgebraic_pair_trace'] for r in records]==[-248,-244]
    ratio=QQ(records[1]['Artin_Tate_discriminant_representative'])/records[0]['Artin_Tate_discriminant_representative']
    assert ratio==QQ(29)/25 and not ratio.is_square()
    # If geometric Picard rank were20, specialization would have full rank
    # in each reduction, forcing the two discriminants to agree mod squares.
    paths=[Path(__file__),PARENT,LOADER]
    return {'schema':'curve302.geometric-picard19.v1','status':'PASS',
            'sources':{str(path.relative_to(ROOT)):sha256(path.read_bytes()).hexdigest() for path in paths},
            'directly_verified_rational_divisor_rank':19,'reductions':records,
            'Artin_Tate_discriminant_ratio':str(ratio),'ratio_is_square':False,
            'geometric_Picard_rank':19,'arithmetic_Picard_rank':19,
            'geometric_generic_MW_rank':17,'arithmetic_generic_MW_rank':17,
            'full_geometric_basis_is_displayed_rational_basis':True,
            'absolute_full_geometric_NS_determinant':1092,'transcendental_rank':3,
            'rank19_moduli_CM_point':False,
            'external_controlled_reduction_required':False,
            'boundary':'Exact finite fibre cardinalities use Sage/PARI. Theorems used: good-reduction divisor specialization, Weil/Poincare duality, Tate and Artin--Tate for K3 surfaces in characteristic>=5, and Shioda--Tate. No BSD/GRH, new rational-fibre rank or original-provenance claim.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--build',action='store_true');args=parser.parse_args();signal.alarm(120)
    result=build()
    if args.build:
        assert not OUT.exists();OUT.write_text(json.dumps(result,indent=2)+'\n')
    assert result==json.loads(OUT.read_text())
    print('PASS geometric Picard19, full geometric MW17, discriminant1092; Artin--Tate ratio29/25')
