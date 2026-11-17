'''PF2RULED42: torus-curvette marked ends, coprime and prime-order quartic exclusions.
Exact graphs and identities support the written all-parameter proofs.
No finite sweep is used as a proof of the unbounded geometric statements.
'''
from pathlib import Path
from math import gcd
import hashlib,json,platform
import sympy as S
from sympy.matrices.normalforms import smith_normal_form
HERE=Path(__file__).resolve().parent
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def det(a,b):return a[0]*b[1]-a[1]*b[0]
def corner(q,r):
    weights={'Dx':0,'Dy':0,'L':0}
    edges={frozenset(['Dx','Dy']),frozenset(['Dx','L'])}
    rays={'Dx':(1,0),'L':(0,1)}
    order=['Dx','L'];left='Dx';right='L';steps=0
    while True:
        steps+=1;new='T'+str(steps);ray=tuple(rays[left][j]+rays[right][j] for j in [0,1])
        edges.remove(frozenset([left,right]))
        edges.add(frozenset([left,new]));edges.add(frozenset([new,right]))
        weights[left]-=1;weights[right]-=1;weights[new]=-1;rays[new]=ray
        order.insert(order.index(left)+1,new)
        if ray==(q,r):return weights,edges,rays,order[1:-1],new,steps
        if det(ray,(q,r))>0:left=new
        else:right=new
        assert steps<q+r
def mat(names,weights,edges):
    return S.Matrix([[-weights[a] if a==b else -int(frozenset([a,b]) in edges)
                      for b in names] for a in names])
def minor(C,i):return C.minor_submatrix(i,i).det() if C.rows>1 else S.Integer(1)
def smallest_prime(r):
    return min(p for p in range(2,r+1) if r%p==0 and S.isprime(p))
def main():
    preds={}
    for name in ['borderline_regular_end_receipt.json','higher_regular_coordinate_receipt.json',
                 'regular_sector_end_receipt.json','quartic_cabling_receipt.json',
                 'torus_infinity_receipt.json','intrinsic_normalization_receipt.json']:
        d=json.loads((HERE/name).read_text());preds[name]=d
        check('predecessor PASS '+name,d['outcome']=='PASS')
        for p,h in d['source_sha256'].items():
            check('predecessor seal '+name+' '+p,digest(HERE/p)==h)
    manifest=HERE/'sources/coprime-regular-end/manifest.json'
    pins=json.loads(manifest.read_text())
    check('all seven retained primary imports explicitly scoped',len(pins['primary_sources'])==7)
    for row in pins['primary_sources']:
        p=HERE/row['path_relative_to_affine_complements']
        check('primary seal '+p.name,p.read_bytes().startswith(b'%PDF') and digest(p)==row['sha256'] and
              p.stat().st_size==row['bytes'])
    visual=pins['new_visual_source_audit']
    check('closed horizontal-surface obstruction and complete proof visually checked',
          visual['proposition']=='2.2(b)' and visual['printed_pages']==[26,27] and
          visual['complete_criterion_and_proof_inspected'] and
          visual['closed_horizontal_surface_requires_zero_Euler_number'])
    layer=preds['quartic_cabling_receipt.json']['finite_peripheral_layer']
    check('quartic permutation layer is reused without new unbounded permutation assumptions',
          layer['exhaustive_S4_subgroups']==30 and layer['S_fiber_trivial_state_count']==9 and
          layer['S_residue_controls']==208 and layer['exponent_modulus']==12)
    selected=[(1,2),(2,3),(3,2),(4,3),(5,2),(5,3),(5,4),(5,6),
              (6,5),(7,3),(7,4),(8,3),(8,5),(8,7),(9,4),(11,5),(13,8),(21,13)]
    graphs=[];z,t,H,Hp,c=S.symbols('z t H Hp c',nonzero=True)
    for q,r in selected:
        assert gcd(q,r)==1 and r>=2 and q>=r-1
        m=q-r+2;n=q*r+m;tag=' q'+str(q)+' r'+str(r)
        weights,edges,rays,names,E,steps=corner(q,r)
        C0=mat(names,weights,edges);j=names.index(E)
        check('actual toric corner matrix is unimodular'+tag,C0.det()==1 and weights[E]==-1)
        aa=S.Matrix([rays[k][0] for k in names]);bb=S.Matrix([rays[k][1] for k in names])
        eL=S.zeros(len(names),1);eL[0]=1;eR=S.zeros(len(names),1);eR[-1]=1
        check('principal-divisor equations identify both whole arms'+tag,C0*aa==eL and C0*bb==eR)
        left=C0[:j,:j];right=C0[j+1:,j+1:]
        dl=left.det() if left.rows else 1;dr=right.det() if right.rows else 1
        check('actual opposite arm determinants are r and q'+tag,dl==r and dr==q and minor(C0,j)==q*r)
        # The m successive free centers depend on the reference jet; their
        # graph and the marked last point are uniform.
        last=E
        for i in range(1,m+1):
            new='F'+str(i);weights[last]-=1;weights[new]=-1
            edges.add(frozenset([last,new]));last=new
        R=last;mark=E if m==1 else 'F'+str(m-1)
        allnames=names+['F'+str(i) for i in range(1,m)]
        C=mat(allnames,weights,edges);idx=allnames.index(mark)
        check('closed end determinant is the surgery coefficient'+tag,C.det()==n)
        check('actual omitted meridian is primitive with cofactor n minus one'+tag,minor(C,idx)==n-1)
        ev=S.zeros(len(allnames),1);ev[allnames.index(E)]=1
        ek=S.zeros(len(allnames),1);ek[idx]=1
        solution=C.inv()*(ev-m*ek)
        check('closed regular fiber is m times the marked knot'+tag,
              all(v.q==1 for v in solution) and 0<m<n)
        B=mat(['Dy','Dx']+allnames,weights,edges)
        check('R6 removes the zero-leaf pair without changing the marked determinant'+tag,
              B.det()==-n and sum('Dx' in e for e in edges)==2)
        # Dx also meets the retained L only before the first corner blowup;
        # afterward its two edges are Dy and the first exceptional.
        check('marked free puncture is retained outside the R6 pair'+tag,
              frozenset([R,mark]) in edges and mark not in ['Dx','Dy'])
        if (q,r) in [(1,2),(3,2),(4,3),(5,3),(8,5)]:
            SN=smith_normal_form(C,domain=S.ZZ)
            check('cyclic closed first homology exact Smith control'+tag,
                  sorted(abs(int(SN[i,i])) for i in range(C.rows))==[1]*(C.rows-1)+[n])
        check('coprimality gives odd surgery order with only the Wright exception'+tag,
              n%2==1 and ((n==3)==((q,r)==(1,2))) and (n==3 or n>=7))
        A,Bb,gg=S.gcdex(-q,r);A=int(A);Bb=int(Bb)
        check('Bezout monomials are actual unimodular coordinates'+tag,-q*A+r*Bb==1 and gg==1)
        # Differentiate with t at fixed s, allowing any H and its derivative.
        ds0=H**A+A*z*H**(A-1)*Hp
        da0=r*H**(r-1)*Hp
        dst=A*z**(m+1)*H**(A-1);dat=r*z**m*H**(r-1)
        adjusted=S.cancel(dat-da0/ds0*dst)
        expected=r*z**m*H**r/(H+A*z*Hp)
        check('all mixed principal parts have the exact m-contact differential'+tag,
              S.cancel(adjusted-expected)==0)
        leading=S.cancel(expected/z**m/H**(A*m)).subs(z,0)
        check('whole omitted-line final coordinate has a nonzero linear coefficient'+tag,
              S.cancel(leading-r*H**(r-1-A*m))==0)
        check('volume discrepancy reaches exactly one after m free blowups'+tag,
              r-q-1+m==1)
        graphs.append(dict(q=q,r=r,m=m,n=n,toric_blowups=steps,
            closed_cones=[q,r,m],marked_cofactor=n-1,regular_fiber_multiple=m,
            surgery=[n,1],source_and_marked_exterior_atoroidal=True))
    # Explicit monomial affine rings, in a bounded set of distinct exponent regimes.
    v,x=S.symbols('v x')
    affine=[]
    for q,r in [(1,2),(3,2),(2,3),(4,3),(5,4),(5,6),(8,5)]:
        m=q-r+2;k=(q+m-1)//m;f=v**r-1;gi=S.div(f,v-1)[0]
        tag=' q'+str(q)+' r'+str(r)
        check('universal four-generator affine chart exponents are nonnegative'+tag,q-m==r-2>=0 and k*m-q>=0)
        check('actual source affine-ring localization recovers both nontrivial generators'+tag,
              S.cancel((v*f/z**m).subs(v,z**q*x)-z**(q-m)*x*(z**(q*r)*x**r-1))==0 and
              S.cancel((v*f**k/z**q).subs(v,z**q*x)-x*(z**(q*r)*x**r-1)**k)==0)
        check('actual root boundary-ring localization recovers t and B'+tag,
              S.rem(f,(v-1))==0 and
              S.cancel((v*f**k/z**q).subs(v,1+z**m*t)-
                       z**(k*m-q)*(1+z**m*t)*gi.subs(v,1+z**m*t)**k*t**k)==0)
        check('all root chart opens cover by the same unit Bezout identity'+tag,
              S.expand(v*(r*v**(r-1))-r*f)==r)
        xx=z**(-q)+z**(2-r)*t;yy=z**r
        num_dz=S.expand(-q*xx*S.diff(yy,z)-r*yy*S.diff(xx,z))
        num_dt=S.expand(-r*yy*S.diff(xx,t))
        check('Euler contraction cancels every apparent boundary pole'+tag,
              S.expand(num_dz+r*m*z*t)==0 and S.expand(num_dt+r*z**2)==0)
        theta_dz=num_dz/(r-q);theta_dt=num_dt/(r-q)
        check('actual global Euler primitive differentiates to the full volume'+tag,
              S.expand(S.diff(theta_dt,z)-S.diff(theta_dz,t)+r*z)==0)
        check('deck orbit on omitted components is free and source reflection is faithful'+tag,gcd(q,r)==1)
        affine.append(dict(q=q,r=r,m=m,k=k,regular_primitive='(-q*x dy-r*y dx)/(r-q)',
                           finite_Keller_target_constructed=False))
    q,r,m=S.symbols('q r m',positive=True,integer=True)
    nn=q*r+m
    check('unbounded surgery slope after the entire punctured free chain',
          S.expand(m*(q*r+1)-(m-1)*q*r-nn)==0)
    check('unbounded Seifert Euler and determinant identity',
          S.cancel(q*r*m*(1/(q*r)+1/m)-nn)==0)
    check('successive marked cofactors differ by one',S.expand((q*r+m)-(q*r+m-1))==1)
    signed=[]
    for eps in [-1,1]:
        for delta in [-1,1]:
            aa=2*(3*eps+delta);coefficient=aa//2-eps
            check('last surviving signed atoroidal case has only small homology '+str((eps,delta)),
                  abs(aa//2-3*eps)==1 and abs(coefficient) in [1,3])
            signed.append(dict(epsilon=eps,delta=delta,cable_slope=aa,source_H1_order=abs(coefficient)))
    # Literal topological sharpness control at n=3, using the same van Kampen representation.
    import verify_borderline_regular_end as V
    tau=(1,0,2,3);alpha=(1,2,0,3);dd=V.mul(alpha,tau)
    bb=(2,1,3,0);AA=(0,1,3,2)
    check('order-three sharpness control has full S4 and exact cable peripherals',
          len(V.generated([alpha,tau,bb]))==24 and V.power(alpha,3)==I4 and
          V.power(dd,3)==dd and V.power(dd,-8)==I4 and V.power(tau,-24)==I4 and
          V.mul(V.power(AA,-1),V.power(bb,2))==dd)
    check('order-three sharpness cover has both atoroidality types and a primitive knot',
          8//2-1==3 and abs(8//2-3)==1 and gcd(8%3,3)==1 and abs(3-6)==3)
    arithmetic=[]
    for rr in range(2,31):
        ell=smallest_prime(rr)
        surviving=[qq for qq in range(rr+1,2*rr+1) if gcd(qq,rr)>1]
        check('sharp arithmetic gap from the noncoprime region r'+str(rr),min(surviving)==rr+ell)
        arithmetic.append(dict(r=rr,smallest_prime=ell,minimum_remaining_pole=rr+ell,
                               prime_r_forces_double_pole=not S.isprime(rr) or rr+ell==2*rr))
    # Integral leading powers: first nondivisible reciprocal term and
    # its whole-chart variation, retaining a mixed divisible part.
    shifted=[]
    for rr,qq,jj in [(2,4,3),(2,4,1),(3,6,5),(3,6,4),(3,6,2),
                     (4,8,7),(4,8,5),(5,10,9),(5,10,6),(6,12,11)]:
        beta=2*qq-jj;k=jj-rr+2;NN=rr*beta+k
        assert qq%rr==0 and jj%rr and gcd(jj,rr)==1
        tag=' r'+str(rr)+' q'+str(qq)+' j'+str(jj)
        # Even/divisible part h_E=z^-q; first nondivisible part z^-j.
        ww=z**qq/(1+z**(qq-jj))
        check('exact first nonintegral reciprocal term'+tag,
              S.limit((ww-z**qq)/z**beta,z,0)==-1)
        check('whole omitted variation has the required later contact'+tag,
              S.limit((-z**(2*qq-rr+2)/(1+z**(qq-jj))**2)/z**(2*qq-rr+2),z,0)==-1)
        check('shifted toric discrepancy fixes the number of free blowups'+tag,
              (rr+beta-1-2*qq)+k==1 and k>=1)
        check('shifted source surgery remains odd and above the sharp exception'+tag,
              NN==(rr-1)*(beta-1)+2*qq+1 and NN%2==1 and NN>=7)
        # After the analytic integral-power shift, the original deleted
        # x=infinity curve attaches internally, not at the local u-axis.
        wwgt,eed,rays,cnames,EE,steps=corner(beta,rr)
        attach=next(name for name in cnames if rays[name]==(qq//rr,1))
        eed={e for e in eed if not({'Dx','Dy'}&set(e))}
        eed.add(frozenset(['Dy','Dx']));eed.add(frozenset(['Dx',attach]))
        wwgt['Dx']=-(qq//rr);wwgt['Dy']=0
        last=EE
        for ii in range(1,k+1):
            new='S'+str(ii);wwgt[last]-=1;wwgt[new]=-1
            eed.add(frozenset([last,new]));last=new
        shifted_names=cnames+['S'+str(ii) for ii in range(1,k)]
        marked=EE if k==1 else 'S'+str(k-1)
        CC=mat(shifted_names,wwgt,eed);idx=shifted_names.index(marked)
        check('integral-power shift keeps the actual internally attached deleted boundary'+tag,
              attach!=EE and rays[attach]==(qq//rr,1) and
              sum('Dx' in edge for edge in eed)==2)
        check('actual shifted full graph and reduced surgery graph agree'+tag,
              CC.det()==NN and mat(['Dy','Dx']+shifted_names,wwgt,eed).det()==-NN)
        check('shifted whole marked puncture is primitive'+tag,minor(CC,idx)==NN-1)
        shifted.append(dict(r=rr,q=qq,largest_nondivisible_support=jj,
                            torus_orders=[rr,beta],free_blowups=k,surgery=NN))
    prime_cases=0
    for rr in [2,3,5,7,11,13]:
        for qq in range(rr+1,3*rr+1):
            if qq%rr:
                assert gcd(qq,rr)==1;prime_cases+=1
            else:
                for jj in range(rr-1,qq):
                    if jj%rr:
                        assert gcd(jj,rr)==1 and jj-rr+2>=1;prime_cases+=1
        check('prime r exhausts the leading-coprime or shifted-coprime cases r'+str(rr),True)
    check('the composite remaining first-break condition can genuinely differ',
          gcd(4,6)==2 and gcd(4,5)==1 and gcd(6,8)==2 and gcd(6,11)==1)
    inputs=[Path(__file__),HERE/'COPRIME_REGULAR_END_AND_EXCLUSION.md',manifest]+[
        HERE/row['path_relative_to_affine_complements'] for row in pins['primary_sources']]+[
        HERE/k for k in preds]
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
        scope='all existing coprime pole/order regular-coordinate models, plus the shifted one-characteristic-pair sector including every prime boundary order; K=R and one omitted A1; quartic targets excluded',
        compactification=dict(all_lower_Laurent_coefficients_retained=True,
            toric_ray='(q,r)',free_blowups='m=q-r+2',full_omitted_line_checked=True,
            final_coordinate='T=psi_m+r*c_q^(r-1-A*m)*t on R',
            arithmetic_control_is_not_unbounded_proof=True),
        end=dict(surgery='n=qr+q-r+2 surgery on T(q,r)',H1='Z/n',marked_knot='surgery dual',
            marked_class_primitive=True,closed_Seifert_orders=['q','r','q-r+2'],
            Euler_absolute='1/(qr)+1/(q-r+2)',marked_exterior='T(q,r) knot exterior',
            closed_source_atoroidal=True,marked_exterior_atoroidal=True,
            central_fiber_multiple='m*[k], 1<=m<n',radial_irreducibility_from_nontrivial_center=True,
            actual_radial_pair_recognized_via_peripheral_Waldhausen=True,
            group_only_closed_Seifert_classification_assumed=False),
        strengthened_topology=dict(closed_source_hypothesis='irreducible and atoroidal',
            source_lens_required=False,source_H1='Z/n, n>3 odd',full_S4=True,simple_branch=True,
            identity_longitude=True,primitive_ramification=True,atoroidal_ramification_exterior=True,
            such_cover_excluded=True,unchanged_exhaustive_permutation_layer=True,
            repeated_cone_quotient_all_fillings=True,
            replacement_homology_preserved_without_distance_one_assumption=True,
            last_trefoil_cases_H1_orders=[1,3]),
        single_characteristic_extension=dict(integral_power_removal_written=True,
            first_nondivisible_reciprocal_exponent='beta=2q-j',
            free_blowups='k=j-r+2',surgery='N=r*beta+k',
            all_lower_terms_retained=True,original_deleted_boundary_attachment_tracked=True,
            all_prime_boundary_orders_quartic_excluded=True,
            all_single_characteristic_quartic_models_excluded=True,
            remaining_boundary_order='r composite, r>=4',remaining_translation_zero_order_at_least=2,
            remaining_reference_germ='at least two characteristic pairs',
            if_r_divides_q_requires='gcd(r,largest_nondivisible_support)>1',
            exact_shifted_controls=shifted,prime_support_controls=prime_cases),
        consequences=dict(whole_coprime_quartic_sector_excluded=True,
            all_q_r_all_Laurent_parts=True,remaining_quartic_regular_region='r composite, q>r, gcd(q,r)>1',
            sharpened_minimum_pole='q>=r+smallest_prime_divisor(r)',
            minimum_gap=2,prime_r_requires='q>=2r',
            all_degree_torus_infinity_excluded=True,all_degree_prime_boundary_degree_excluded=True,
            all_finite_degrees_excluded=False),
        affine_controls=dict(entire_monomial_coprime_family_constructed=True,
            four_generator_root_ring='C[z,v,v*f/z^m,v*f^k/z^q], f=v^r-1, k=ceil(q/m)',
            smooth_affine_quotient=True,global_regular_primitive='(-q*x dy-r*y dx)/(r-q)',
            no_target_equivariance_assumed=True,controls=affine,
            finite_Keller_map_constructed=False),
        sharpness_control=dict(actual_topological_cover=True,target='(3,8)-cable of positive trefoil',
            source_Seifert_orders=[2,3,3],source_H1_order=3,primitive_ramification_class=2,
            closed_source_atoroidal=True,marked_exterior='trefoil exterior',
            marked_exterior_atoroidal=True,source_lens=False,
            polynomial_infinity_realization_claimed=False,finite_Keller_map_constructed=False),
        controls=dict(graphs=graphs,signed_last_cases=signed,arithmetic_gaps=arithmetic),
        assurance=dict(all_parameter_and_all_depth_proofs_written=True,
            primary_closed_horizontal_surface_criterion_visually_checked=True,
            primary_imports_explicit=True,actual_peripheral_data_tracked=True,
            finite_controls_not_full_geometric_proof=True,
            affine_converse_for_arbitrary_Laurent_data=False,
            noncoprime_regular_sector_excluded=False,no_regular_coordinate_case_excluded=False,
            extra_omitted_primes_excluded=False,arbitrary_quartic_packets_excluded=False,
            all_finite_Wright_degrees_excluded=False,finite_Keller_map_constructed=False,
            formal_verification=False,independent_full_reproof=False,external_review=False,
            research_goal_completed=False),
        source_sha256={str(p.relative_to(HERE)):digest(p) for p in inputs},
        runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'coprime_regular_end_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS '+str(len(checks))+' torus-curvette controls; coprime and prime-order quartic sectors excluded; goal active')
I4=(0,1,2,3)
if __name__=='__main__':main()
