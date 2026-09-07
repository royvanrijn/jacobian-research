#!/usr/bin/env sage-python
"""Exact single-outlier recovery of a rank17 form from302, with lattice gates.

No MILP replay is needed: linear algebra identifies the unique removable row
and reconstructs the Gram. Rational LDL enumerates norms2 and4. Verify the
primitive K3-lattice embedding and the local L_546 conditions. The published
Shimura-curve identification is cited, not independently rederived here.
One worker,120seconds, two shell enumerations capped5million nodes each.
No elliptic K3 equation, generic section, or302 parent is certified.
"""
import argparse,json,runpy,signal
from hashlib import sha256
from pathlib import Path
from sage.all import QQ,ZZ,GF,EllipticCurve,QuaternionAlgebra,matrix,vector,block_diagonal_matrix,lcm,gcd,kronecker
ROOT=Path(__file__).resolve().parents[2]
INPUT=ROOT/'artifacts/generated-results/elliptic-curves/curve302_integral_shell_inputs_v1.json'
PUBLIC=ROOT/'elliptic-curves/cas/icarm_curve302.py'
ENUM=ROOT/'elkies-k3/scripts/certify_curve302_anchor6_triangle_gate.sage'
OUT=ROOT/'artifacts/generated-results/elliptic-curves/curve302_det1092_recovery_v1.json'
def rows(A):return [list(map(str,r)) for r in A.rows()]
def canon(v):return min(tuple(v),tuple(-a for a in v))
def build():
    data=json.loads(INPUT.read_text())['target_core'];words=data['integral_words'];n=17
    pairs=[(i,j) for i in range(n) for j in range(i,n)]
    A=matrix(QQ,[[v[i]*v[j]*(1 if i==j else 2) for i,j in pairs]+[-1] for v in words])
    assert A.nrows()==401 and A.rank()==154
    left=A.left_kernel_matrix();deletable=[i for i in range(A.nrows()) if not any(left.column(i))]
    assert deletable==[356]
    retained=A.matrix_from_rows([i for i in range(len(words)) if i not in deletable]);K=retained.right_kernel_matrix()
    assert retained.rank()==153 and K.nrows()==1 and K[0,-1]
    G=matrix(QQ,n)
    for a,(i,j) in zip(4*K[0]/K[0,-1],pairs):G[i,j]=G[j,i]=a
    assert all(a in ZZ for a in G.list());G=matrix(ZZ,G)
    assert G.is_positive_definite() and all(a==4 for a in G.diagonal()) and G.det()==1092
    heights=[vector(ZZ,v)*G*vector(ZZ,v) for v in words]
    assert heights.count(4)==400 and heights[356]==6
    enum=runpy.run_path(str(ENUM))['exact_shell'];roots,nodes2=enum(G,2);shell,nodes4=enum(G,4)
    assert not roots and len(shell)==2436
    smith=list(G.smith_form()[0].diagonal());assert smith==[1]*16+[1092]
    inv=G.inverse();g=inv.column(0);assert lcm([a.denominator() for a in g])==1092
    assert g*G*g==QQ(4997)/1092 and g*G*g-QQ(629)/1092 in 2*ZZ
    # The sole possible prime index of a proper even overlattice is2.
    # Its unique order-two discriminant class has odd norm, not zero mod2.
    assert (546*g)*G*(546*g)==4997*273 and (4997*273)%2==1
    B=matrix(ZZ,data['public_embedding']);assert B.nrows()==31 and B.ncols()==17
    assert list(B.smith_form()[0].diagonal())==[1]*17
    pub=runpy.run_path(str(PUBLIC));E=EllipticCurve(QQ,pub['GENERAL_WEIERSTRASS_COEFFICIENTS'])
    public=[E(QQ(x),QQ(y)) for x,y in pub['POINTS']]
    pts=[sum((a*P for a,P in zip(w,public) if a),E(0)) for w in B.columns()]
    def point(w):return sum((a*P for a,P in zip(w,pts) if a),E(0))
    assert all(point(w)[0].denominator()==1 for w in words)
    fitted={canon(w) for w in words};rays=sorted({canon(w) for w in shell});integral=[];new=[]
    assert len(rays)==1218
    for w in rays:
        P=point(w)
        if P[0].denominator()!=1:continue
        integral.append(w)
        if w not in fitted:new.append(dict(word=list(map(int,w)),x=str(P[0]),y=str(P[1])))
    assert len(integral)==404 and len(new)==4
    # Explicit anti-isotropic graph glue of N=U+(-G) and a compatible T.
    N=block_diagonal_matrix(matrix(ZZ,[[0,1],[1,0]]),-G)
    T=matrix(ZZ,[[-2,1,0],[1,2,2],[0,2,220]])
    assert T.det()==-1092
    d=[QQ(-2),QQ(5)/2,QQ(1092)/5]
    L=matrix(QQ,[[1,0,0],[-QQ(1)/2,1,0],[0,QQ(4)/5,1]])
    assert L*matrix.diagonal(QQ,d)*L.transpose()==T
    t=T.inverse().column(2);assert lcm([a.denominator() for a in t])==1092 and t*T*t==QQ(5)/1092
    unit=-25;assert (629*unit**2-5)%(2*1092)==0
    glue=vector(QQ,list(unit*N.inverse().column(2))+list(t))
    ambient=block_diagonal_matrix(N,T);assert glue*ambient*glue in 2*ZZ
    generators=matrix(ZZ,1092*matrix.identity(ZZ,22)).stack(matrix(ZZ,[1092*glue]))
    change=matrix(QQ,generators.row_module().basis_matrix())/1092
    unimod=change*ambient*change.transpose()
    assert all(a in ZZ for a in unimod.list()) and all(a%2==0 for a in unimod.diagonal()) and unimod.det()==-1
    Ni=matrix(ZZ,change.inverse()[:19,:]);Ti=matrix(ZZ,change.inverse()[19:,:])
    assert list(Ni.smith_form()[0].diagonal())==[1]*19 and list(Ti.smith_form()[0].diagonal())==[1]*3
    assert Ni*unimod*Ni.transpose()==N and Ti*unimod*Ti.transpose()==T and Ni*unimod*Ti.transpose()==0
    local=[]
    for p in [3,7,13]:
        c=ZZ(4997*(1092//p));actual=kronecker(c,p);required=-kronecker(1092//p,p)
        assert actual==required
        local.append(dict(prime=p,positive_frame_dual_numerator=int(c),legendre=int(actual),required=int(required)))
    quaternion=QuaternionAlgebra(QQ,5,QQ(2184)/5);assert quaternion.discriminant()==546
    # Published model of X(546)/<w546> (Rotger, thesis Table6.4).
    C=EllipticCurve(QQ,[1,0,1,-137,380]);P=C(-9,-26)
    counts=[]
    for p in [5,11,17]:
        assert C.discriminant()%p
        counts.append(dict(prime=p,order=int(C.change_ring(GF(p)).cardinality())))
    torsion_bound=gcd([r['order'] for r in counts]);assert torsion_bound==4 and 4*P!=C(0)
    return dict(schema='curve302.det1092-recovery.v1',status='EXACT_ROOTLESS_FORM_AND_PRIMITIVE_K3_LATTICE_EMBEDDING_NO_PARENT',
        input_sha256={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),INPUT,PUBLIC,ENUM]},
        recovery=dict(input_words=401,quadratic_constraint_rank=154,unique_deletable_row=356,retained_constraint_rank=153,
                      outlier_word=words[356],outlier_height=6,gram=rows(G),determinant=1092,minimum=4,smith_factors=list(map(int,smith)),
                      exact_norm2_nodes=nodes2,exact_norm4_nodes=nodes4,signed_norm4_count=len(shell),no_proper_even_overlattice=True,
                      discriminant_generator=list(map(str,g)),generator_norm='4997/1092',discriminant_norm_mod_2Z='629/1092',primitive_embedding_in_displayed_D=True),
        point_diagnostic=dict(norm4_rays=len(rays),fit_norm4_rays=400,unused_norm4_rays=818,integral_norm4_rays=404,
                              unused_integral_points=new,ray_words_sha256=sha256(json.dumps(rays,separators=(',',':')).encode()).hexdigest()),
        primitive_K3_embedding=dict(N_gram=rows(N),T_gram=rows(T),T_signature=[2,1],T_discriminant_norm='5/1092',matching_unit=unit,
                                    glue_vector=list(map(str,glue)),unimodular_basis=rows(change),unimodular_gram=rows(unimod),
                                    ambient_signature=[3,19],N_primitive=True,T_primitive=True),
        arithmetic_lead=dict(local_L546_conditions=local,even_Clifford_parameters=['5','2184/5'],quaternion_discriminant=546,
            cited_moduli_curve='X(546)/<w546>',cited_elliptic_model=[1,0,1,-137,380],rational_point=['-9','-26'],
            good_reduction_counts=counts,torsion_order_upper_bound=int(torsion_bound),four_times_point=list(map(str,(4*P).xy())),
            rational_point_is_nontorsion=True,individual_point_non_CM_status='UNKNOWN',
            references=['https://arxiv.org/pdf/0802.1301 (section2, pp6-7)',
                        'https://web.mat.upc.edu/victor.rotger/docs/Tesi.pdf (Table6.4)'],
            identification_boundary='Local L546 tests and the elliptic arithmetic are exact. The moduli interpretation and its published elliptic model are literature inputs, not independently reconstructed period maps. No explicit rational non-CM marked K3 or universal K3 equation is supplied.'),
        boundary='The normalized quadratic form is recovered from a finite target point set, not proved to be a generic height pairing. Its primitive complex K3-lattice embedding and local arithmetic type do not identify a surface through302. No generic MW basis, specialization parameter, discoverer provenance, or alternative302 parent is certified.')
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--build',action='store_true');args=parser.parse_args();signal.alarm(120)
    result=build();serialized=json.dumps(result,indent=2,sort_keys=True,default=lambda x:int(x))+'\n'
    if args.build:assert not OUT.exists();OUT.write_text(serialized)
    assert json.loads(serialized)==json.loads(OUT.read_text());print(result['status'])
