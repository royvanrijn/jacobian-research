#!/usr/bin/env sage-python
"""Read-only exact replay from portable JSON, without producer checkpoints."""
import hashlib
import json
import zipfile
import argparse
from pathlib import Path
from sage.all import *

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--artifacts',type=Path,default=ART)
parser.add_argument('--class-index',type=int,default=3)
args=parser.parse_args(); ART=args.artifacts
assert args.class_index not in (1,6)
P='x1092_class'+str(args.class_index)+'_realization_'
def read(name):return json.loads((ART/(P+name+'_v1.json')).read_text())
def mat(x):return matrix(ZZ,x)
source=json.loads((ART/'curve302_recovered_mw17_parent_v1.json').read_text())
mark=read('marking'); tr=read('trace'); rr=read('rr'); eq=read('equation'); sec=read('sections')
census=json.loads((ART/'det1092_pruned_rootless_j2_census_v1.json').read_text())
for name,digest in mark['inputs'].items():
    assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest
G=mat(source['generic_height_gram']); J=matrix(ZZ,[[0,1],[1,0]])
N=block_diagonal_matrix(J,-G); T=mat(mark['transport_rows_D_D_plus_O_complement'])
H=mat(mark['frame_gram']);D=vector(ZZ,mark['fibre_D']);O=vector(ZZ,mark['rational_zero_O'])
assert mat(mark['source_NS_gram'])==N
assert T*N*T.T==block_diagonal_matrix(J,-H) and abs(T.det())==1
assert D==T[0] and D+O==T[1] and O==vector(ZZ,[-1,1]+[0]*17)
assert D*N*O==1 and D*N*D==0 and O*N*O==-2
iso=mat(mark['representative_columns_in_complement'])
target=mat(census['rootless_classes'][args.class_index-1]['gram'])
assert abs(iso.det())==1 and iso.T*H*iso==target and H.det()==1092
with zipfile.ZipFile(ART/'det1092_pruned_anchor_packets_v1.zip') as z:
    requested=mat(json.loads(z.read('anchor-16.json'))['embeddings'][mark['requested_shared_core_witness']['new_embedding_index']]['gram'])
req=mat(mark['requested_embedding_frame_rows_in_source_NS'])
assert req*N*req.T==-requested and T[:2]*N*req.T==0
assert abs(T[:2].stack(req).det())==1
assert pari(H).qfminim(2)[0]==0 and pari(G).qfminim(2)[0]==0
assert not pari(G).qfisom(pari(H))
assert not pari(mat(census['rootless_classes'][0]['gram'])).qfisom(pari(H))
assert mark['class_index']==args.class_index and list(D[:2])==[3,2]
assert list(D[2:])==mark['generic_trace_word']
witness=mark['requested_shared_core_witness']
with zipfile.ZipFile(ART/'det1092_pruned_anchor_packets_v1.zip') as z:
    packet=json.loads(z.read('anchor-16.json'))
core=mat(witness['common_core_basis_in_niemeier'])
assert mat(witness['common_core_gram']).det()==witness['common_core_determinant']
for index in (witness['known_embedding_index'],witness['new_embedding_index']):
    embedding=packet['embeddings'][index]
    assert embedding['sixth_index']==witness['shared_sixth_index']
    frame=mat(embedding['complement_basis_in_ambient'])
    coordinates=mat(frame.T.solve_right(core.T).T)
    assert coordinates*mat(embedding['gram'])*coordinates.T==mat(witness['common_core_gram'])
    assert coordinates.row_module().saturation()==coordinates.row_module()
# All potentially negative effective components have old degree <= D.F=2.
# Enumerate the homogenized distance ellipsoid and check actual intersections.
w=vector(ZZ,D[2:]); wall_counts=[]
for degree in (1,2):
    Q=block_matrix(ZZ,[[4*G,(-2*degree*G*w).column()],
       [(-2*degree*w*G).row(),matrix(ZZ,[[degree*degree*(w*G*w)+1]])]])
    raw=pari(Q).qfminim(7); wall_counts.append(int(raw[0]))
    for col in matrix(ZZ,raw[2]).columns():
        if abs(col[-1])!=1:continue
        v=vector(ZZ,col[:-1])*col[-1]; height=v*G*v
        if (height-2)%(2*degree):continue
        C=vector(ZZ,[(height-2)//(2*degree),degree]+list(v))
        assert C*N*C==-2 and D*N*C>=0
assert wall_counts==[r['full_short_vector_count'] for r in mark['nef_certificate']['complete_horizontal_wall_shells']]
print('PASS full NS transport, target-class isometry, distinct frame, nef fibre and rational zero',flush=True)

R=PolynomialRing(QQ,'t'); F=R.fraction_field();t=R.gen()
def dec(r,K=F):
    R0=K.ring();return K(R0(r['numerator']))/R0(r['denominator'])
Esource=EllipticCurve(F,[dec(a) for a in source['a_invariants']])
basis=[Esource([dec(c) for c in p]) for p in source['basis_weierstrass_coordinates']]
Eold=EllipticCurve(F,[-Esource.c4()/48,-Esource.c6()/864])
oldpts=[Eold(p[0]+Esource.b2()/12,p[1]+(Esource.a1()*p[0]+Esource.a3())/2) for p in basis]
A,B,h,nx,ny=[R(tr[k]) for k in ('A','B','h','Nx','Ny')]
assert Eold.a4()==A and Eold.a6()==B
trace=sum((ZZ(c)*p for c,p in reversed(list(zip(mark['generic_trace_word'],oldpts)))),Eold(0))
assert trace==Eold(F(nx/h**2),F(ny/h**3)) and h.degree()==4
a0,b0,a1,b1=[R(rr[k]) for k in ('a0','b0','a1','b1')]
cols=[(t**i*nx)%(h*h) for i in range(8)]+[(-t**i*ny)%(h*h) for i in range(2)]
M=matrix(QQ,8,10,lambda i,j:cols[j][i]); kr=matrix(QQ,rr['kernel'])
assert M.rank()==8 and kr.rank()==2 and M*kr.T==0
assert all((a*nx-b*ny)%(h*h)==0 for a,b in ((a0,b0),(a1,b1)))
Ru=PolynomialRing(QQ,'u'); K=Ru.fraction_field();u=Ru.gen()
S=PolynomialRing(K,'t');F2=S.fraction_field();tt=S.gen()
lift=lambda p:S(list(p))
q=S([dec(c,K) for c in eq['quartic_coefficients']]);t0,q0=[dec(c,K) for c in eq['quartic_zero']]
sf=S([dec(c,K) for c in eq['radical_square_factor']]);ds=S([dec(c,K) for c in eq['radical_denominator_sqrt']])
m=F2((lift(a1)-u*lift(a0))/((u*lift(b0)-lift(b1))*lift(h)))
xp,yp=F2(lift(nx)/lift(h)**2),F2(lift(ny)/lift(h)**3)
assert m**4-6*xp*m*m-8*yp*m-3*xp*xp-4*lift(A)==q*(sf/ds)**2
assert q(t0)==q0*q0 and q.gcd(q.derivative()).degree()==0
V=PolynomialRing(F2,'v');L=F2.extension(V.gen()**2-F2(q),'v');v=L.gen()
z=L(tt-t0);e,d,c,b,a=[K(q(tt+t0)[i]) for i in range(5)]
X=(2*q0*(v+q0)+d*z)/z**2
Y=((X*X-4*e*a)*z-d*X-2*e*b)/(2*q0)
E=EllipticCurve(K,[dec(r,K) for r in eq['a_invariants']])
assert E.a_invariants()==(0,c,0,b*d-4*e*a,e*b*b+a*d*d-4*e*a*c)
assert Y*Y==X**3+E.a2()*X*X+E.a4()*X+E.a6()
zi=(2*q0*Y+d*X+2*e*b)/(X*X-4*e*a)
assert zi==z and (X*zi*zi-d*zi-2*e)/(2*q0)==v
xx=(L(m*m-xp)+L(sf/ds)*v)/2; yy=L(m)*(xx-L(xp))-L(yp)
assert yy*yy==xx**3+L(lift(A))*xx+L(lift(B))
print('PASS exact equation and mutually inverse rational maps',flush=True)

child_rows=[]
for r in sec['sections']:
    xnew,ynew,tnew=[dec(r[k],K) for k in ('X','Y','t')]
    E(xnew,ynew)
    znew=(2*q0*ynew+d*xnew+2*e*b)/(xnew*xnew-4*e*a)
    vnew=(xnew*znew*znew-d*znew-2*e)/(2*q0)
    assert znew+t0==tnew and vnew*vnew==q(tnew)
    xs=K(((m*m-xp).numerator()(tnew)/(m*m-xp).denominator()(tnew)+sf(tnew)/ds(tnew)*vnew)/2)
    ms=m.numerator()(tnew)/m.denominator()(tnew)
    ys=ms*(xs-xp.numerator()(tnew)/xp.denominator()(tnew))-yp.numerator()(tnew)/yp.denominator()(tnew)
    word=vector(ZZ,r['source_word'])
    if r['source_kind']=='old_section':
        original=sum((ZZ(c0)*p for c0,p in zip(word,oldpts)),Eold(0))
        assert xs==original[0].numerator()(tnew)/original[0].denominator()(tnew)
        assert ys==original[1].numerator()(tnew)/original[1].denominator()(tnew)
        div=vector(ZZ,[(word*G*word-2)//2,1]+list(word))
    else:
        glue=read('glue'); assert glue['source_word']==list(word)
        x0,x1,y0,y1,branch=[lift(R(glue[k])) for k in ('x0','x1','y0','y1','branch')]
        svalue=(xs-x0(tnew))/x1(tnew)
        assert svalue*svalue==branch(tnew) and ys==y0(tnew)+y1(tnew)*svalue
        assert word*G*word==10
        div=vector(ZZ,[2,2]+list(word))
    child=div*T.inverse();assert child[1]==1
    assert list(child[2:])==r['child_frame_coordinates'];child_rows.append(child[2:])
C=matrix(ZZ,child_rows); hg=C*H*C.T
assert abs(C.det())==1 and hg==mat(sec['height_gram']) and hg.det()==1092
assert hg.is_positive_definite() and sec['rank_lower']==sec['rank_upper']==17
print('PASS 17 exact rational sections, saturated height Gram determinant1092, generic rank17',flush=True)

parent=read('parent');compact=read('compact_parent')
Ep=EllipticCurve(K,[dec(r,K) for r in parent['a_invariants']]);g=dec(parent['gauge_from_pointed_cubic'],K)
assert Ep.a4()==(-E.c4()/48)*g**4 and Ep.a6()==(-E.c6()/864)*g**6
ub=dec(compact['u_of_s'],K); den=ub.denominator();sc=QQ(compact['scale'])
def subst(f):return f.numerator()(ub)/f.denominator()(ub)
Ec=EllipticCurve(K,[dec(r,K) for r in compact['a_invariants']])
assert Ec.a4()==subst(Ep.a4())*den**8/sc**4 and Ec.a6()==subst(Ep.a6())*den**12/sc**6
for i,r in enumerate(sec['sections']):
    px,py=[dec(k,K) for k in parent['basis_weierstrass_coordinates'][i]]
    assert px==g*g*(dec(r['X'],K)+E.b2()/12) and py==g**3*dec(r['Y'],K)
    cx,cy=[dec(k,K) for k in compact['basis_weierstrass_coordinates'][i]]
    assert cx==subst(px)*den**4/sc**2 and cy==subst(py)*den**6/sc**3
    Ec(cx,cy)
print('PASS normalized and compact parent transports; no specialization used',flush=True)

print('PASS separate exact replay; nonisometric to class1 and curve302; arithmetic strict classes UNKNOWN',flush=True)
