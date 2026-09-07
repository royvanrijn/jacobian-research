#!/usr/bin/env sage-python
"""Standalone exact replay of Kihara geometry, seeds and parent separation.
Requires only Sage and the exported JSON bundle. No repository imports.
"""
import argparse,json
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,lcm,gcd,prod

def finite_rank(model,points,primes,ell):
    E=EllipticCurve(QQ,model);P=[E([QQ(x),QQ(y)]) for x,y in points];projective=[]
    for point in P:
        den=lcm([c.denominator() for c in point]);v=[ZZ(c*den) for c in point];g=gcd(v);projective.append([c//g for c in v])
    rows=[];records=[]
    for prime in primes:
        F=GF(prime);e=EllipticCurve(F,[F(c) for c in model])
        if not e.discriminant():raise ArithmeticError('good finite specialization required')
        key=lambda point:tuple(int(c) for c in point)
        elements=e.points();multiples={key(ell*P):ell*P for P in elements};mask={key(P):0 for P in multiples.values()};reps=[e(0)]
        while len(mask)<len(elements):
            P=next(P for P in elements if key(P) not in mask);old=list(reps);size=len(old)
            for digit in range(1,ell):
                for i,R in enumerate(old):
                    rep=R+digit*P;reps.append(rep)
                    for T in multiples.values():
                        k=key(rep+T)
                        if k in mask:raise ArithmeticError('quotient cosets overlap')
                        mask[k]=i+size*digit
        dimension=ZZ(len(reps)).valuation(ell)
        if ell**dimension!=len(reps) or dimension>2:raise ArithmeticError('elliptic quotient dimension differs')
        reduced=[e([F(c) for c in P]) for P in projective]
        for j in range(dimension):rows.append([(mask[key(P)]//ell**j)%ell for P in reduced])
        records.append({'prime':prime,'group_order':len(elements),'ell_multiple_subgroup_order':len(multiples),'quotient_dimension':int(dimension)})
    rank=int(matrix(GF(ell),rows).rank())
    if rank!=14:raise ArithmeticError('selected generic seed is not injective in finite ell quotients')
    return {'modulus':ell,'rank':rank,'groups':records}

def check_geometry(rawA,rawB,g):
 R=rawA.parent();A=R(g['minimal_A']);B=R(g['minimal_B']);s=R(g['removed_scale'])
 assert rawA==A*s**4 and rawB==B*s**6
 for f,_ in A.gcd(B).factor():assert min(A.valuation(f)//4,B.valuation(f)//6)==0
 D=-16*(4*A**3+27*B**2);assert D.gcd(A).degree()==0
 factors=[(R(r['coefficients']),r['multiplicity']) for r in g['finite_factors']]
 assert prod(f**m for f,m in factors).monic()==D.monic()
 for i,(f,m) in enumerate(factors):
  assert f.gcd(f.derivative()).degree()==0
  for h,_ in factors[:i]:assert f.gcd(h).degree()==0
 chi=max((A.degree()+3)//4,(B.degree()+5)//6);n=12*chi-D.degree()
 assert A.degree()==4*chi and B.degree()==6*chi and n>0
 assert [int(A.degree()),int(B.degree()),int(D.degree())]==g['minimal_degrees']
 assert g['chi']==chi and g['infinity_type']=='I'+str(n)
 root=n-1+sum(f.degree()*(m-1) for f,m in factors)
 assert root==g['geometric_root_rank'] and g['h11']==10*chi and g['geometric_MW_Hodge_upper_bound']==10*chi-2-root
 return A,B

def count_parent(A,B,row):
 p=row['prime'];F=GF(p);R=PolynomialRing(F,'T');a,b=R(A),R(B);D=-16*(4*a**3+27*b**2)
 assert (a.degree(),b.degree(),D.degree())==(8,12,20) and D.gcd(D.derivative())==1 and D.gcd(a)==1
 node=-3*b[12]/(2*a[8]);assert (3*node).is_square()
 assert 3*node**2+a[8]==0 and node**3+a[8]*node+b[12]==0
 counts=[]
 for t in range(p+1):
  av,bv=(a[8],b[12]) if t==p else (a(F(t)),b(F(t)))
  if 4*av**3+27*bv**2:n=EllipticCurve(F,[av,bv]).cardinality()
  else:n=1+sum(1 if v==0 else 2 if v.is_square() else 0 for v in [x**3+av*x+bv for x in F])
  counts.append(int(n))
 assert counts==row['weierstrass_fibre_counts'] and sum(counts)+3*p==row['surface_point_count']

def main(path):
 d=json.loads(path.read_text());c=d['certificates'];R=PolynomialRing(QQ,'t')
 e,h,v,b,a=[R(f) for f in c['curve302_inverse_kihara_and_rank14_16_intake_v1.json']['quartic_scaled_coefficients_low_to_high']]
 I=12*a*e-3*b*h+v*v;J=72*a*v*e+9*b*v*h-27*a*h*h-27*b*b*e-2*v**3
 g=c['kihara_global_parent_geometry_v1.json']['geometry'];check_geometry(-27*I,-27*J,g);assert g['chi']==36
 print('PASS minimal global surface chi36',flush=True)
 parents=c['kihara_five_parent_distinctness_v1.json'];verified_counts=0
 for row,seed in zip(parents['rows'],d['rows']):
  assert row['path_parameter']==seed['parameter'];t=QQ(seed['parameter']);p=t*t*(8+3*t*t);q=-6*(2+t*t)*(4+t*t)
  roots=[0,(2*p*p+p*q+2*q*q)**2,2*(p+q)**2*(2*p*p+p*q+q*q),q*q*(4*p*p-p*q+4*q*q),p*(2*p-q)*(2*p*p+4*p*q+5*q*q),4*p**4+8*p**3*q+9*p*p*q*q-2*p*q**3+2*q**4]
  normalizer=roots[1];roots=[r/normalizer for r in roots];assert roots==list(map(QQ,row['normalized_roots']))
  R=PolynomialRing(QQ,'T');T=R.gen();K=R.fraction_field();S=PolynomialRing(K,'x');x=S.gen()
  F=prod(x-r-sign*T for sign in [1,-1] for r in roots);G=x**6
  for j in range(5,-1,-1):G+=(F[6+j]-(G*G)[6+j])/2*x**j
  quartic=(G*G-F)/T**2;assert quartic==S([K(R(f)) for f in row['quartic_coefficients']])
  e,h,v,b,a=quartic.list();I=12*a*e-3*b*h+v*v;J=72*a*v*e+9*b*v*h-27*a*h*h-27*b*b*e-2*v**3
  A=R(row['raw_A']);B=R(row['raw_B']);assert A==-27*I and B==-27*J
  check_geometry(A,B,row['geometry']);assert row['geometry']['chi']==2
  delta=-16*(4*A**3+27*B**2);assert delta.degree()==20 and delta.gcd(delta.derivative())==1
  node=-3*B[12]/(2*A[8]);assert (3*node).is_square()
  Q=[(K(X),K(Y)) for X,Y in row['quartic_points']];assert len(Q)==13 and all(Y*Y==quartic(X) for X,Y in Q)
  assert [X for X,Y in Q[:12]]==[r+sign*T for sign in [1,-1] for r in roots]
  assert all(Y==G(X)/T for X,Y in Q[:12])
  E=EllipticCurve(K,[A,B]);P=[E([K(X),K(Y)]) for X,Y in row['generic_sections']];assert len(P)==12
  X0,Y0=Q[0];translated=quartic(x+X0);ee,dd,cc,bb,aa=translated.list()
  for (X,Y),point in zip(Q[1:],P):
   z=X-X0;wx=(2*Y0*(Y+Y0)+dd*z)/(z*z);wy=(2*(wx*wx-4*Y0*Y0*aa)*z-2*dd*wx-4*Y0*Y0*bb)/(4*Y0)
   assert point==E([9*wx+3*cc,27*wy])
  Tc=QQ(row['control_T']);oldE=EllipticCurve(QQ,seed['old_model']);old=[oldE([QQ(X),QQ(Y)]) for X,Y in seed['old_points']]
  Et=EllipticCurve(QQ,[A(Tc),B(Tc)]);special=[Et([z(Tc) for z in point.xy()]) for point in P]
  isos=[iso for iso in Et.isomorphisms(oldE) if all(iso(point)==v-old[0] for point,v in zip(special,old[1:13]))]
  assert len(isos)==1 and list(isos[0].tuple())==list(map(QQ,row['control_isomorphism']))
  if t==2:
   cert=c['kihara_rank14_t2_v1.json']['independence_certificate'];ell=cert['relation_prime'];assert ell==5
   finite_rank(seed['old_model'],seed['old_points'],[r['prime'] for r in cert['rows']],ell)
   witness=cert['torsion_witness'];n=EllipticCurve(GF(witness['prime']),seed['old_model']).cardinality();assert n==witness['group_order'] and n%5
   fixed=c['kihara_fixed_parent_geometry_v1.json'];assert fixed['raw_A']==row['raw_A'] and fixed['raw_B']==row['raw_B'] and fixed['generic_sections']==row['generic_sections']
  else:
   s=seed['seed'];En=EllipticCurve(QQ,s['curve']);oldimages=[En([QQ(X),QQ(Y)]) for X,Y in s['original_points']];new=[En([QQ(X),QQ(Y)]) for X,Y in s['points']]
   assert len(new)==14 and oldimages[1:]==new[1:] and 6*new[0]==sum(oldimages[:12],En(0))
   assert any(all(iso(P)==Q for P,Q in zip(old,oldimages)) for iso in oldE.isomorphisms(En))
   oldintake=seed['old_intake'];E1=EllipticCurve(QQ,oldintake['curve']);P1=[E1([QQ(X),QQ(Y)]) for X,Y in oldintake['points']]
   assert any(all(iso(P)==Q for P,Q in zip(P1,oldimages)) for iso in E1.isomorphisms(En))
   finite=seed['finite'];assert finite['curve']==s['curve'] and finite['points']==s['points']
   finite_rank(s['curve'],s['points'],[r['prime'] for r in finite['signatures']],2)
   prime=finite['rank_certificate']['no_rational_2_torsion_prime'];F2=GF(prime);R2=PolynomialRing(F2,'z');z=R2.gen()
   assert (z**3+F2(En.a4())*z+F2(En.a6())).is_irreducible()
   assert s['old_subgroup_index_in_new']==6
  for count in row['surface_counts']:
   if count['status']=='PASS':count_parent(A,B,count);verified_counts+=1
  print('PASS parent',t,'generic12 and control rank14',flush=True)
 old=c['mestre_parent_portfolio_intake_v1.json'];allrows=parents['rows']+old['rows']+[old['reference_parent']]
 labels=['kihara-path-'+r['path_parameter'] for r in parents['rows']]+[r['id'] for r in old['rows']+[old['reference_parent']]]
 counts=[{r['prime']:r['surface_point_count'] for r in row['surface_counts'] if r.get('status','PASS')=='PASS'} for row in allrows];expected=[]
 for i in range(5):
  for j in range(i+1,len(allrows)):
   ps=[p for p in [131,239,251] if p in counts[i] and p in counts[j] and counts[i][p]!=counts[j][p]]
   expected.append({'left':labels[i],'right':labels[j],'separating_primes':ps,'status':'PROVED_Q_DISTINCT' if ps else 'UNKNOWN'})
 assert expected==parents['pairwise_separations'] and sum(r['status']=='PROVED_Q_DISTINCT' for r in expected)==44
 assert [r for r in expected if r['status']=='UNKNOWN']==[{'left':'kihara-path-2','right':'kihara-path-5/2','separating_primes':[],'status':'UNKNOWN'}]
 print('PASS four new Q-distinct parents; index6 seeds; 44/45 comparisons;',verified_counts,'new good-prime counts',flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);a=p.parse_args();main(a.input)
