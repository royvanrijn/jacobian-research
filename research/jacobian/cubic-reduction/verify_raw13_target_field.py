#!/usr/bin/env python3
"""Independent exact cubic identities and companion-matrix specialization minor."""
from pathlib import Path
from fractions import Fraction as Q
from math import isqrt
import hashlib,json,time,ast
from verify import decode,variable,add,scale,mul,compose,derivative,determinant3,evaluate
B=Path(__file__).resolve().parent;O=B/'raw13-target-field';ROOT=B.parents[1]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def con(c,n):return {(0,)*n:Q(c)} if c else {}
def total(*ps):
 z={}
 for p in ps:z=add(z,p)
 return z

def pw(p,k,n):
 z=con(1,n)
 for _ in range(k):z=mul(z,p)
 return z

def madd(A,B,p):return [[(a+b)%p for a,b in zip(r,s)] for r,s in zip(A,B)]
def ms(A,c,p):return [[x*c%p for x in r] for r in A]
def mm(A,B,p):return [[sum(A[i][k]*B[k][j] for k in range(3))%p for j in range(3)] for i in range(3)]
def eye():return [[int(i==j) for j in range(3)] for i in range(3)]
def det(A,p):return (A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])-A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])+A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]))%p
def inv_adj(A,p):
 d=det(A,p);assert d;di=pow(d,p-2,p);Z=[]
 for i in range(3):
  row=[]
  for j in range(3):
   rs=[k for k in range(3) if k!=j];cs=[k for k in range(3) if k!=i];co=(A[rs[0]][cs[0]]*A[rs[1]][cs[1]]-A[rs[0]][cs[1]]*A[rs[1]][cs[0]])*((-1)**(i+j));row.append(co*di%p)
  Z.append(row)
 assert mm(A,Z,p)==mm(Z,A,p)==eye();return Z,d

def field_rank_certificate(D,basis,F,p,point_count,minor_size,extra_known=None):
 extra_known=[] if extra_known is None else extra_known
 # Every selected field functional is recomputed using companion matrices;
 # the producer used polynomial triples,division by the cubic and a3x3 solve.
 allrows=[];point_dets=[];positions={e:i for i,e in enumerate(basis)};known=[[Q(0)]*len(basis) for _ in range(4+len(extra_known))];known[0][0]=Q(1)
 for row,z in enumerate(F+extra_known,1):
  for e,co in z.items():known[row][positions[e]]=co
 for point in D['target_points']:
  ss,tt,uu=point;assert all(isinstance(x,int) for x in point) and ss%p
  I=eye();Z=[[0,0,2*ss*ss*uu%p],[1,0,-ss*tt%p],[0,1,1]];Z2=mm(Z,Z,p);Z3=mm(Z2,Z,p)
  assert madd(madd(Z3,ms(Z2,-1,p),p),madd(ms(Z,ss*tt,p),ms(I,-2*ss*ss*uu,p),p),p)==[[0]*3 for _ in range(3)]
  EE=madd(madd(ms(Z2,3,p),ms(Z,-2,p),p),ms(I,ss*tt,p),p);Ei,ed=inv_adj(EE,p);point_dets.append(ed);is_=pow(ss,p-2,p);i2=pow(2,p-2,p);AA=ms(Ei,ss,p);BB=ms(madd(Z,ms(EE,-1,p),p),i2*is_%p,p);E2=mm(EE,EE,p);E3=mm(E2,EE,p);CC=ms(madd(madd(ms(E2,5,p),ms(mm(Z,EE,p),-3,p),p),ms(E3,-2,p),p),i2*is_*is_%p,p)
  assert mm(AA,BB,p)==mm(BB,AA,p) and mm(AA,CC,p)==mm(CC,AA,p) and mm(BB,CC,p)==mm(CC,BB,p)
  powers=[]
  for z in [AA,BB,CC]:
   stack=[I]
   for j in range(max(map(sum,basis))):stack.append(mm(stack[-1],z,p))
   powers.append(stack)
  matrices=[mm(mm(powers[0][i],powers[1][j],p),powers[2][k],p) for i,j,k in basis]
  for row,z in enumerate(F):
   actual=[[0]*3 for _ in range(3)]
   for e,co in z.items():actual=madd(actual,ms(matrices[positions[e]],int(co),p),p)
   assert actual==ms(I,point[row],p)
  rows=[[z[j][0] for z in matrices] for j in [1,2]]
  assert all(sum(int(x)*y for x,y in zip(vv,row))%p==0 for vv in known for row in rows)
  allrows+=rows
 assert len(D['target_points'])==point_count and len(allrows)==minor_size and len(set(map(tuple,D['target_points'])))==point_count and not D['excluded_singular_points']
 assert len(D['selected_rows'])==len(D['minor_columns'])==len(set(D['minor_columns']))==minor_size
 selected=[]
 for z in D['selected_rows']:
  assert z['component'] in [1,2] and z['row_index']==2*z['point']+z['component']-1;selected.append([allrows[z['row_index']][j] for j in D['minor_columns']])
 # Independent column elimination of the square minor,not producer row-basis selection.
 A=[row[:] for row in selected];determinant=1
 for i in range(minor_size):
  j=next((j for j in range(i,minor_size) if A[i][j]),None);assert j is not None
  if j!=i:
   for row in A:row[i],row[j]=row[j],row[i]
   determinant=-determinant
  pivot=A[i][i];determinant=determinant*pivot%p;inverse=pow(pivot,p-2,p)
  for j in range(i+1,minor_size):
   fac=A[i][j]*inverse%p
   if fac:
    for k in range(i,minor_size):A[k][j]=(A[k][j]-fac*A[k][i])%p
 assert determinant
 return {"determinant_mod_prime":determinant,"matrix_sha256":hashlib.sha256(json.dumps(selected,separators=(",",":")).encode()).hexdigest(),"E_determinants":point_dets,"points":point_count,"minor_size":minor_size,"columns":len(basis)}

def main():
 started=time.monotonic();cp=O/'checkpoint.json';C=json.loads(cp.read_text());ip=O/'input.json';D=json.loads(ip.read_text());raw=B/'raw13-specialization/input.json';R=json.loads(raw.read_text());assert C['status']=='complete' and C['source_sha256']==D['source_sha256']==sha(B/'raw13_target_field.py') and C['input_sha256']==sha(ip) and C['source_input_sha256']==D['source_input_sha256']==sha(raw)
 p=D['modulus'];assert p==1000000007 and p>2 and all(p%d for d in range(2,isqrt(p)+1))
 basis=sorted([(i,j,k) for i in range(8) for j in range(8) for k in range(8) if i+j+k<=7],key=lambda e:(sum(e),-e[0],-e[1]));assert list(map(list,basis))==D['source_basis'] and len(basis)==120
 n=3;a,b,c=[variable(j,n) for j in range(3)];F=[a,b,c]
 for atom in R['atoms']:F[atom['row']]=add(F[atom['row']],{tuple(atom['exponent']):Q(atom['coefficient'])})
 assert determinant3([[derivative(z,j) for j in range(3)] for z in F],n)==con(1,n) and [{e:co for e,co in z.items() if sum(e)<=1} for z in F]==[a,b,c]
 s,t,u=F;g={tuple([e[0]-1,e[1],e[2]]):co for e,co in s.items()};assert all(e[0]>=1 for e in s);v=add(con(1,n),scale(mul(a,b),2));W=mul(g,v);E=total(scale(pw(W,2,n),3),scale(W,-2),mul(s,t))
 assert E==g and total(pw(W,3,n),scale(pw(W,2,n),-1),mul(mul(s,t),W),scale(mul(pw(s,2,n),u),-2))=={}
 assert mul(a,E)==s and scale(mul(s,b),2)==add(W,E,Q(-1))
 assert scale(mul(pw(s,2,n),c),2)==total(scale(pw(E,2,n),5),scale(mul(W,E),-3),scale(pw(E,3,n),-2))
 rank7=field_rank_certificate(D,basis,F,p,58,116);determinant=rank7['determinant_mod_prime'];selected_hash=rank7['matrix_sha256'];point_dets=rank7['E_determinants']
 # A nonlinear target shear attains source degree8;it is a control,not a simplification.
 n=3;Q8=pw(F[0],2,n);assert max(sum(e) for e in Q8)==8;G=[F[0],add(F[1],Q8),F[2]];assert max(sum(e) for z in G for e in z)==8 and [sum(sum(e)>=2 for e in z) for z in G]==[2,11,6]
 assert determinant3([[derivative(z,j) for j in range(3)] for z in G],n)==con(1,n);points=[[0,1,-8],[-1,1,4]];assert [[evaluate(z,x) for z in G] for x in points]==[[Q(0),Q(1),Q(0)]]*2
 # Complete degree8 source-field intersection,with Q,S,T as three independent new classes.
 Q8=pw(F[0],2,3);S8=total(pw(F[1],3,3),scale(mul(F[0],pw(F[2],2,3)),27),scale(mul(F[1],F[2]),-9));T8=add(pw(F[1],2,3),scale(mul(F[0],S8),-4));extra=[Q8,S8,T8]
 saved=json.loads((O/'known-degree8.json').read_text());assert extra==[decode(z,3) for z in saved['source_polynomials']]
 heads=[{e:c for e,c in z.items() if sum(e)==8} for z in extra];assert [list(z) for z in heads]==[[(6,0,2)],[(3,4,1)],[(4,2,2)]]
 p8=O/'degree8';cp8=p8/'checkpoint.json';C8=json.loads(cp8.read_text());ip8=p8/'input.json';D8=json.loads(ip8.read_text());assert C8['status']=='complete' and C8['source_sha256']==D8['source_sha256']==sha(B/'raw13_target_field8.py') and C8['input_sha256']==sha(ip8) and C8['source_input_sha256']==sha(raw)
 basis8=sorted([(i,j,k) for i in range(9) for j in range(9) for k in range(9) if i+j+k<=8],key=lambda e:(sum(e),-e[0],-e[1]));assert len(basis8)==165 and list(map(list,basis8))==D8['source_basis']
 rank8=field_rank_certificate(D8,basis8,F,p,79,158,extra)
 # Reconstruct the complete tangent-identity target determinant,then multiply all seven linear memberships.
 loc=O/'target-locus';LI=json.loads((loc/'input.json').read_text());names=LI['parameters'];nn=12;xx,yy,zz=[variable(i,nn) for i in range(3)];vv={name:variable(3+i,nn) for i,name in enumerate(names)};qq=pw(xx,2,nn);ss=total(pw(yy,3,nn),scale(mul(xx,pw(zz,2,nn)),27),scale(mul(yy,zz),-9));tt=add(pw(yy,2,nn),scale(mul(xx,ss),-4));RR=[total(z,mul(vv['u'+str(i)],qq),mul(vv['v'+str(i)],ss),mul(vv['w'+str(i)],tt)) for i,z in enumerate([xx,yy,zz])];J=add(determinant3([[derivative(z,j) for j in range(3)] for z in RR],nn),con(1,nn),Q(-1));mons=sorted({e[:3] for e in J},reverse=True);eqs=[{e[3:]:co for e,co in J.items() if e[:3]==mon} for mon in mons];assert len(eqs)==27 and eqs==[decode(z,9) for z in LI['equations']]
 def parse(text):
  def visit(t):
   if isinstance(t,ast.Constant):return con(t.value,9)
   if isinstance(t,ast.Name):return variable(names.index(t.id),9)
   if isinstance(t,ast.UnaryOp) and isinstance(t.op,ast.USub):return scale(visit(t.operand),-1)
   assert isinstance(t,ast.BinOp)
   if isinstance(t.op,ast.Pow):return pw(visit(t.left),ast.literal_eval(t.right),9)
   aa,bb=visit(t.left),visit(t.right)
   if isinstance(t.op,ast.Add):return add(aa,bb)
   if isinstance(t.op,ast.Sub):return add(aa,bb,Q(-1))
   if isinstance(t.op,ast.Mult):return mul(aa,bb)
   assert isinstance(t.op,ast.Div) and len(bb)==1 and (0,)*9 in bb;return scale(aa,1/bb[(0,)*9])
  return visit(ast.parse(text.replace('^','**'),mode='eval').body)
 cert_cp=json.loads((loc/'certificate-checkpoint.json').read_text());assert cert_cp['status']=='complete-primary-certificates' and cert_cp['input_sha256']==sha(loc/'input.json') and cert_cp['script_sha256']==sha(loc/'certify.sing') and cert_cp['log_sha256']==sha(loc/'certificate.log')
 text=(loc/'certificate.log').read_text();targets=[];certificate_terms=0;max_certificate_degree=0
 for block in text.split('CERT_BEGIN|')[1:]:
  header,tail=block.split('\n',1);lines=tail.split('CERT_END')[0].strip().splitlines();assert len(lines)==27;result={}
  for i,(line,eq) in enumerate(zip(lines,eqs),1):
   index,m=line.split('|',1);assert int(index)==i;z=parse(m);certificate_terms+=len(z);max_certificate_degree=max(max_certificate_degree,max((sum(e) for e in z),default=0));result=add(result,mul(z,eq))
  assert result==variable(names.index(header),9);targets.append(header)
 assert targets==['u0','v0','v1','v2','w0','w1','w2']
 # Generic two-parameter shear and its two-sided target inverse.
 pp=[variable(3+i,nn) for i in range(9)];plus=[xx,add(yy,mul(vv['u1'],qq)),add(zz,mul(vv['u2'],qq))];minus=[xx,add(yy,mul(vv['u1'],qq),Q(-1)),add(zz,mul(vv['u2'],qq),Q(-1))]
 assert [compose(z,minus+pp,nn) for z in plus]==[xx,yy,zz] and [compose(z,plus+pp,nn) for z in minus]==[xx,yy,zz]
 target_locus={'parameters':9,'native_determinant_generators':27,'linear_memberships':7,'multiplier_terms':certificate_terms,'max_multiplier_degree':max_certificate_degree,'all_memberships_multiplied_exactly':True,'surviving_parameters':['u1','u2'],'normal_form':'Y+(0,u1,u2)*Y1^2','generic_two_sided_inverse_checked':True,'fixed_source_nonlinear_postcomposition_minimum':19,'fixed_source_bound_uses13_old_and6_disjoint_new_source_columns':True,'input_sha256':sha(loc/'input.json'),'certificate_checkpoint_sha256':sha(loc/'certificate-checkpoint.json'),'certificate_log_sha256':sha(loc/'certificate.log')}
 old_columns={e for z in F for e in z if sum(e)>=2};assert len(old_columns)==13 and not (old_columns&set(Q8)) and len(Q8)==6
 assert determinant3([[derivative(z,j) for j in range(3)] for z in plus],nn)==con(1,nn)
 ledger=json.loads((ROOT/'MATH_STATUS.json').read_text());entries={e['id']:e for e in ledger['entries']};e=entries['KCUB41'];assert e['state']=='proved';pin={'KCUB41':{'scope_fingerprint':hashlib.sha256((e['state']+'\0'+e['title']+'\0'+e['scope']).encode()).hexdigest()[:16],'proof_sha256':sha(ROOT/e['canonical_source'])}}
 receipt={'status':'pass','claim':'KCUB42','checker_sha256':sha(Path(__file__)),'producer_sha256':C['source_sha256'],'arithmetic_sha256':sha(B/'verify.py'),'input_sha256':sha(ip),'checkpoint_sha256':sha(cp),'discovery_sha256':sha(O/'discovery.json'),'full_generic_cubic_and_all_source_recovery_identities_checked':True,'irreducibility_by_written_primitive_linear_target_argument':True,'source_degree_cap':7,'source_degree_caps':[7,8],'degree8_field_rank':rank8,'degree8_kernel_basis':['1','H1','H2','H3','H1^2','S','T'],'degree8_kernel_dimension':7,'degree8_target_locus':target_locus,'degree8_input_sha256':sha(ip8),'degree8_checkpoint_sha256':sha(cp8),'known_degree8_sha256':sha(O/'known-degree8.json'),'source_monomial_columns':120,'target_degree_cap':None,'target_denominator_cap':None,'rational_field_remainder_dimensions':2,'exact_rational_target_specialization_points':58,'selected_minor_size':116,'independent_companion_matrix_and_adjugate_reconstruction':True,'independent_column_elimination':True,'modulus':p,'prime_checked_by_trial_division':True,'nonzero_minor_determinant_mod_prime':determinant,'selected_minor_matrix_sha256':selected_hash,'all_specialization_denominators_nonzero_mod_prime':True,'specialization_E_determinants':point_dets,'kernel_dimension_exact':4,'kernel_basis':['1','H1','H2','H3'],'all_rational_targets_with_source_degree_le7_are_affine':True,'degree7_rigidity_persists_in_arbitrary_affine_source_target_frames':True,'all_target_only_compositions_preserving_degree7_are_affine':True,'minimum_source_degree_for_nonaffine_rational_target_with_polynomial_pullback':8,'matching_nonlinear_target_shear_control':{'target':'(Y1,Y2+Y1^2,Y3)','degree':8,'nonlinear_counts':[2,11,6],'full_determinant_one':True,'collision_points':points,'common_image':[0,1,0],'global_improvement':False},'imports':pin,'goal_complete':False,'no_global_positive_bound_or_unrestricted_profile_change':True,'trust_boundary':'Generic cubic irreducibility,normal-form field membership,specialization-to-characteristic-zero rank and affine-frame consequences are written,not formalized or externally reviewed. Generic identities,both116/158-square minors,the complete27-generator target model and all seven17-term memberships are independently checked. Point evaluations certify exact rank minors,not an empirical candidate grid. Unknown target denominators are not specialized;only fixed source normal-form denominators are. KCUB41 supplies the affine13 count. The nonlinear degree8 minimum19 keeps the original source coordinates;arbitrary source-affine frames of the resulting two-parameter family remain outside that count.','seconds':time.monotonic()-started}
 proof=B/'RAW13_TARGET_FUNCTION_FIELD_RIGIDITY.md'
 if proof.exists():receipt['proof_sha256']=sha(proof)
 (O/'validation.json').write_text(json.dumps(receipt,indent=2)+'\n');print('PASS full raw13 target field:kernel4/120 through7,kernel7/165 through8;independent116/158 minors;complete9-parameter target locus,7 exact memberships/'+str(certificate_terms)+' terms;only two-parameter quadratic shears;fixed-source nonlinear minimum19;global goal open',flush=True)
if __name__=='__main__':main()
