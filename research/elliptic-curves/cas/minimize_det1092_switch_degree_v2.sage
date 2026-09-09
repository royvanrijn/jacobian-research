#!/usr/bin/env sage-python
"""Two exact, bounded CVPs over all original generic MW translations.

The degree ceiling is supplied by the already available untranslated 2B curve.
No target fibre or exceptional point enters. One --sign per 25-second process.
"""
import argparse,hashlib,json
from fractions import Fraction
from math import isqrt
from pathlib import Path
from sage.all import QQ,ZZ,matrix,vector,identity_matrix,block_diagonal_matrix
from visibility_lattice_v2 import ExactParity
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_two_fibration_action_v1';PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,data):
    p=OUT/name
    if p.exists():assert read(p)==data
    else:
        with p.open('x') as s:json.dump(data,s,indent=2,sort_keys=True);s.write('\n')
parser=argparse.ArgumentParser();parser.add_argument('--sign',type=int,choices=[-1,1],required=True);args=parser.parse_args()
source=read(OUT/'action.json');p=read(PARENT);G=matrix(QQ,p['generic_height_gram'])
NS=block_diagonal_matrix(matrix(QQ,[[-2,1],[1,0]]),-G);F=vector(QQ,[0,1,*([0]*17)])
tau=matrix(QQ,source['alternate_B_translation']);C=vector(QQ,source['source']['NS'])
bound=int((tau*C)*NS*F);assert bound==63
paths=[PARENT,OUT/'action.json',Path(__file__).with_name('visibility_lattice_v2.py'),Path(__file__)]
save('degree-minima-protocol-v2.json',{'classification':'exact generic-only low-degree translation-orbit calculation',
    'rule':'From the unique degree13 previous rational curve C, minimize the original projection degree after one arbitrary original generic translation and then either sign of alternate B. Enumerate the full closed ellipsoid through the known nonsection witness tau_B(C). No coefficient box or basis prefix truncation.',
    'degree_ceiling':bound,'signs':[-1,1],'source_case':source['source']['case'],
    'limits':{'seconds_per_sign':25,'nodes_per_sign':200000,'parameter_evaluations':0,'point_searches':0,'rational_map_expansions':0},
    'inputs':{str(q.relative_to(ROOT)):sha(q) for q in paths}})
K=(tau**(-args.sign))*F;a,b=C[:2];v=vector(QQ,C[2:]);c,e=K[:2];r=vector(QQ,K[2:])
assert a>0 and c>0
coefficient=a*c/2;center=r/c-v/a
constant=C*NS*K-coefficient*(center*G*center)
radius=(bound-constant)/coefficient;assert radius>0
U=matrix(ZZ,G).LLL_gram();assert abs(U.det())==1
Gr=U.transpose()*G*U;ct=U.inverse()*center
metric=ExactParity(Gr.rows());mu,d=metric.mu,metric.d
ct=[Fraction(int(q.numerator()),int(q.denominator())) for q in ct];limit=Fraction(int(radius.numerator()),int(radius.denominator()));word=[0]*17;nodes=0;answers=[]
def visit(i,used):
    global nodes
    nodes+=1
    if nodes>200000:raise RuntimeError('NODE_CAP_EXHAUSTED_NO_COMPLETENESS_CLAIM')
    if i<0:
        q=U*vector(ZZ,word);degree=constant+coefficient*((q-center)*G*(q-center))
        assert degree in ZZ and degree>=1 and degree<=bound
        shifted=vector(QQ,[a,b+v*G*q+a*(q*G*q)/2,*(v+a*q)])
        out=(tau**args.sign)*shifted
        assert out*NS*F==degree and out*NS*out==-2
        answers.append({'word':list(map(int,q)),'degree':int(degree),'NS':list(map(int,out))})
        return
    remaining=limit-used
    if remaining<0:return
    shift=-ct[i]+sum((mu[j][i]*(word[j]-ct[j]) for j in range(i+1,17)),Fraction(0))
    den=shift.denominator;num=shift.numerator
    sq=remaining/d[i]*den*den;rad=isqrt(sq.numerator//sq.denominator)
    lo=-((rad+num)//den);hi=(rad-num)//den
    for val in sorted(range(lo,hi+1),key=lambda val:abs(den*val+num)):
        word[i]=val;visit(i-1,used+d[i]*(Fraction(val)+shift)**2)
visit(16,Fraction(0))
answers.sort(key=lambda row:(row['degree'],row['word']))
nonsections=[r for r in answers if r['degree']>1]
report={'status':'COMPLETE_EXACT_DEGREE_ELLIPSOID','alternate_sign':args.sign,
    'degree_ceiling':bound,'nodes':nodes,'center':list(map(str,center)),
    'constant':str(constant),'quadratic_coefficient':str(coefficient),'radius':str(radius),
    'LLL_columns':[[int(x) for x in row] for row in U.rows()],
    'answers':answers,'minimum_nonsection_degree':min((r['degree'] for r in nonsections),default=None),
    'degree_histogram':{str(dd):sum(r['degree']==dd for r in answers) for dd in sorted(set(r['degree'] for r in answers))},
    'boundary':'Complete only through the stated degree ceiling; this certifies a global nonsection minimum only when one is present. No incidence or specialized rank claim.'}
save('degree-minima-v2-%s.json'%('minus' if args.sign<0 else 'plus'),report)
print(args.sign,'nodes',nodes,'answers',len(answers),'histogram',report['degree_histogram'],flush=True)
