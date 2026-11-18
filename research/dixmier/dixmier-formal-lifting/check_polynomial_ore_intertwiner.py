#!/usr/bin/env python3
"""Exact finite controls for the polynomial Ore-intertwiner admission theorem."""
from fractions import Fraction as Q
from pathlib import Path
from math import comb
import argparse, hashlib, json, platform, tempfile

HERE=Path(__file__).resolve().parent
ONE={(0,0):Q(1)};X={(1,0):Q(1)};D={(0,1):Q(1)}
def clean(p):return {e:c for e,c in p.items() if c}
def add(*polys):
    out={}
    for p in polys:
        for e,c in p.items():out[e]=out.get(e,Q(0))+c
    return clean(out)
def scale(p,c):return clean({e:a*c for e,a in p.items()})
def falling(n,j):
    out=1
    for i in range(j):out*=n-i
    return out
def ore(a,b):
    out={}
    for (i,j),c in a.items():
        for (k,l),d in b.items():
            for n in range(min(j,k)+1):
                e=(i+k-n,j+l-n);out[e]=out.get(e,Q(0))+c*d*comb(j,n)*falling(k,n)
    return clean(out)
def comm(a,b):return add(ore(a,b),scale(ore(b,a),-1))
def rank(matrix):
    a=[list(row) for row in matrix];nr=len(a);nc=len(a[0]) if a else 0;r=0
    for j in range(nc):
        pivot=next((i for i in range(r,nr) if a[i][j]),None)
        if pivot is None:continue
        a[r],a[pivot]=a[pivot],a[r];u=a[r][j]
        a[r]=[c/u for c in a[r]]
        for i in range(r+1,nr):
            if a[i][j]:
                u=a[i][j];a[i]=[c-u*d for c,d in zip(a[i],a[r])]
        r+=1
        if r==nr:break
    return r
def shift(S):return max(i-j for i,j in S)
def sigma(S,n):
    m=shift(S)
    return sum((c*falling(n,j) for (i,j),c in S.items() if i-j==m),Q(0))
def finite_matrix(S,N):
    m=shift(S);out=[[Q(0)]*(N+1) for _ in range(N+m+1)]
    for n in range(N+1):
        for (i,j),c in S.items():
            if n>=j:out[n+i-j][n]+=c*falling(n,j)
    return out
def ore_eps(a,b):
    out={}
    for (i,j,e),c in a.items():
        for (k,l,f),d in b.items():
            for n in range(min(j,k)+1):
                key=(i+k-n,j+l-n,e+f)
                out[key]=out.get(key,Q(0))+c*d*comb(j,n)*falling(k,n)
    return clean(out)

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output-dir',type=Path);args=parser.parse_args()
    H=ore(X,D)
    assert comm(D,X)==ONE
    operators=[('D',D,-1,1,0),('x',X,1,0,1),
               ('xD-2',add(H,scale(ONE,-2)),0,1,1),
               ('D2+xD+1',add(ore(D,D),H,ONE),0,0,0),
               ('xD+1',add(H,ONE),0,0,0)]
    controls=[]
    for name,S,m,ker,coker in operators:
        assert shift(S)==m
        for N in [6,10]:
            r=rank(finite_matrix(S,N))
            assert (N+1-r,N+m+1-r)==(ker,coker)
            assert ker-coker==-m
        controls.append({'operator':name,'m':m,'cutoffs':[6,10],'kernel_dimension':ker,'cokernel_dimension':coker,'index':-m})
    for n in range(20):
        assert sigma(add(H,ONE),n)==n+1
        assert sigma(add(ore(D,D),H,ONE),n)==n+1
    assert comm(X,add(H,ONE))==scale(X,-1)
    # A true scalar intertwiner for the actual standard pair.
    scalar=scale(ONE,2)
    assert ore(X,scalar)==ore(scalar,X) and ore(D,scalar)==ore(scalar,D)
    # Exact one-operator Darboux control; not a pair construction.
    S=add(H,scale(ONE,-2));P={(4,2):Q(1)}
    Pprime={(4,2):Q(1),(3,1):Q(-2),(2,0):Q(2)}
    assert ore(P,S)==ore(S,Pprime)
    assert shift(S)==0 and sigma(S,2)==0
    r=max(j for i,j in S);sr={i:c for (i,j),c in S.items() if j==r};c=max(sr)
    assert (r,c)==(1,1)
    assert Pprime.get((3,1),Q(0))-P.get((3,1),Q(0))==2*c-4*r==-2
    # Same-zero-subprincipal single-operator control, excluded by m!=0.
    Ssquare=add({(2,1):Q(1)},X);Psquare=ore(Ssquare,Ssquare)
    assert Psquare=={(4,2):Q(1),(3,1):Q(4),(2,0):Q(2)}
    assert Psquare[(3,1)]-Q(4*2,2)==0
    assert ore(Psquare,Ssquare)==ore(Ssquare,Psquare)
    assert shift(Ssquare)==1
    assert rank(finite_matrix(Ssquare,8))==9
    # Leading coefficient with a forbidden nonzero root: the necessary
    # equation would make x+1 divide 2x^4, impossible at x=-1.
    Sbad=add(ore(add(X,ONE),D),ONE)
    assert shift(Sbad)==0
    top_residual={i:c for (i,j),c in add(ore(P,Sbad),scale(ore(Sbad,P),-1)).items() if j==2}
    assert top_residual=={4:Q(-2),3:Q(-4)}
    assert sum((c*(-1)**i for i,c in top_residual.items()),Q(0))==2
    # Independent auxiliary parameter: [D,x]=1, epsilon is central.
    Seps={(0,0,0):Q(1),(1,1,1):Q(1)}
    Punit={(4,2,0):Q(1)}
    first_residual=add(ore_eps(Punit,Seps),scale(ore_eps(Seps,Punit),-1))
    assert first_residual=={(4,2,1):Q(-2)}
    # A vanishing-corner single-operator control, outside the unit
    # hypothesis, satisfies exact intertwining and equal residues.
    Psmall={(3,1,1):Q(1),(2,0,0):Q(1),(2,0,1):Q(2)}
    Psmallprime={(3,1,1):Q(1),(2,0,0):Q(1)}
    assert ore_eps(Psmall,Seps)==ore_eps(Seps,Psmallprime)
    assert {key:c for key,c in Psmall.items() if key[2]==0}=={key:c for key,c in Psmallprime.items() if key[2]==0}

    def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
    receipt={'status':'PASS','claim':'Finite controls for the actual-global polynomial Ore-intertwiner Fredholm and subprincipal descent theorem','environment':{'python':platform.python_version(),'arithmetic':'Independent Fraction Ore products and finite exact Gaussian ranks'},'Fredholm_controls':controls,'written_universal_index':'dimker-dimcoker=-max(i-j)','actual_pair_consequence':'kernel=cokernel=0,m=0,allnonnegativeintegerdiagonalvaluesnonzero','exact_leading_coefficient_relation':"pprime_(b-1)=p_(b-1)+b*x^a*s_rprime/s_r-a*r*x^(a-1)",'leading_constraint':'s_r=lambda*x^c,c<=r','subprincipal_shift':'A_Wprime-A_W=b*c-a*r; negativeintegerfornonconstantSanda>b','Darboux_control':{'first_intertwining':True,'shift':-2,'full_pair_admission':False,'reason':'nonzerofinitekernel'},'same_zero_subprincipal_single_operator':{'first_intertwining':True,'m':1,'full_pair_admission':False,'reason':'nonzerofinitecokernel'},'bijective_insufficient_control':'xD+1 bijectiveonpolynomialsbutxS!=Sx;no polynomialstandardPprime','forbidden_nonzero_root_control':'S=(x+1)D+1, P=x4D2: required(x+1)Z=-2x4-4x3 hasresidueobstruction2atx=-1','auxiliary_parameter_corollary':'Fixed[D,x]=1, finiteuniformPBWsupport, S=1modepsilon, unitfixedcenteredcorner implySscalar; normalizedcornerdifferenceisnonzerointegerbutmustvanishmodepsilon','unit_corner_control':'S=1+epsilon*xD,P=x4D2 hasPS-SP=-2epsilon*x4D2','vanishing_corner_single_operator_control':'P=epsilon*x3D+(1+2epsilon)x2,Pprime=epsilon*x3D+x2 satisfyPS=SPprime; corneralphaepsilonisnonunit; nofullpairclaim','auxiliary_parameter_limits':['physicalhbarWeylrelation','positive-parameter-ordercorner','unboundedPBWsupport','generalrationalOreintertwiner'],'hashes':{'POLYNOMIAL_ORE_INTERTWINER_AND_SUBPRINCIPAL_DESCENT.md':sha(HERE/'POLYNOMIAL_ORE_INTERTWINER_AND_SUBPRINCIPAL_DESCENT.md'),'check_polynomial_ore_intertwiner.py':sha(Path(__file__))},'scope':'Scalar polynomial Ore intertwiners between two actual polynomial Weyl pairs. Finite examples are not extrapolated to universal Fredholm theorem. No general rational Weyl-division, matrix, multivariable or formal-gauge classification; noabsoluteA_Wzero theorem or DC1 resolution.','formal_verification':False,'external_review':False}
    out=args.output_dir or Path(tempfile.mkdtemp(prefix='polynomial-ore-intertwiner-'));out.mkdir(parents=True,exist_ok=True);dest=out/'polynomial_ore_intertwiner_controls.json';dest.write_text(json.dumps(receipt,indent=2)+'\n')
    print('PASS exact finite-index controls, one-operator Darboux shift and admission-failure controls');print('Receipt:',dest)

if __name__=='__main__':main()
