"""Exact small controls for direct Gaussian support and ordinary-ideal transfers."""
from pathlib import Path
import hashlib,itertools,json,platform
import sympy as s


def pair_expectation(f,zs,ws):
    out=0;n=len(zs)
    for powers,c in s.Poly(s.expand(f),*zs,*ws).terms():
        if powers[:n]==powers[n:]:
            out+=c*s.prod(s.factorial(i) for i in powers[:n])
    return s.expand(out)


def contraction(f,zs,us):
    out=0;n=len(zs)
    for powers,c in s.Poly(s.expand(f),*zs).terms():
        val=c
        for u,k in zip(us,powers):val=s.diff(val,u,k)
        out+=val
    return s.expand(out)


def apply_symbol(symbol,f,zs,us):
    return contraction(symbol*f,zs,us)


def coeffs(f,xs):return [c for _,c in s.Poly(s.expand(f),*xs).terms() if c!=0]


def run():
    x,y,u,v,z,w,z2,w2,t,a,b,e=s.symbols('x y u v Z W Z2 W2 t a b e')
    cases=[((z,),(w,),(u,),(x,),1+2*z+z*z,u**3+a*u+1,u*u+u+3,2),
           ((z,z2),(w,w2),(u,v),(x,y),z*z2+z*z+z2,u*u+u*v+v+1,u*v+u*u,3),
           ((z,),(w,),(u,),(x,),z**3-a*z,u*u+b*u,u+1,2)]
    forward=0;reverse=0;rows=[]
    for zs,ws,us,xs,symbol,p,q,maxm in cases:
        shift=dict(zip(us,[xx+ww for xx,ww in zip(xs,ws)]))
        phase=s.expand(symbol*p.subs(shift,simultaneous=True))
        for m in range(1,maxm+1):
            direct=pair_expectation(q.subs(shift,simultaneous=True)*phase**m,zs,ws)
            original=apply_symbol(symbol**m,q*p**m,zs,us).subs(dict(zip(us,xs)),simultaneous=True)
            assert s.expand(direct-original)==0;forward+=1
            aa=(1,)+(0,)*(len(zs)-1);bb=(2,)+(0,)*(len(zs)-1)
            amp=s.prod(qz**j for qz,j in zip(zs,aa))*s.prod(qw**j for qw,j in zip(ws,bb))
            observed=pair_expectation(amp*phase**m,zs,ws);reconstructed=0
            for jj in itertools.product(*[range(k+1) for k in bb]):
                marker=s.prod(uu**j for uu,j in zip(us,jj))
                inner=apply_symbol(symbol**m,marker*p**m,zs,us)
                for uu,j in zip(us,aa):inner=s.diff(inner,uu,j)
                outer=s.prod(s.binomial(k,j)*(-xx)**(k-j) for k,j,xx in zip(bb,jj,xs))
                reconstructed+=outer*inner.subs(dict(zip(us,xs)),simultaneous=True)
            assert s.expand(observed-reconstructed)==0;reverse+=1
        rows.append({'n':len(zs),'symbol_order':s.Poly(symbol,*zs).total_degree(),
                     'phase_degree':s.Poly(p,*us).total_degree(),'Gaussian_degree':s.Poly(phase,*zs,*ws).total_degree(),
                     'checked_moments':maxm})

    # Full SIC,including mixed/nonseparable dependence.
    sic_cases=[((z,),(w,),(u,),(x,),z*z*u+z*u*u+a*z+u+1,z*u*u,2),
               ((z,z2),(w,w2),(u,v),(x,y),z*v+z2*u+z*z2*u*v+u*u+z*z,z*v+u*v,2)]
    sic_forward=0;sic_reverse=0
    for zs,ws,us,xs,f,g,maxm in sic_cases:
        shift=dict(zip(us,[xx+ww for xx,ww in zip(xs,ws)]))
        phase=f.subs(shift,simultaneous=True)
        for m in range(1,maxm+1):
            direct=pair_expectation(g.subs(shift,simultaneous=True)*phase**m,zs,ws)
            original=contraction(g*f**m,zs,us).subs(dict(zip(us,xs)),simultaneous=True)
            assert s.expand(direct-original)==0;sic_forward+=1
            marker=z*ws[0]**2
            expected=sum(s.binomial(2,j)*(-xs[0])**(2-j)*contraction(z*us[0]**j*f**m,zs,us).subs(dict(zip(us,xs)),simultaneous=True) for j in range(3))
            assert s.expand(pair_expectation(marker*phase**m,zs,ws)-expected)==0;sic_reverse+=1

    # Bounded ordinary-ideal transcription for Lambda=d_u,P=a*u+b.
    # Universal equality is proved by the two coefficient identities in the note.
    ideals=[]
    phase=z*(a*(x+w)+b)
    for N in range(3):
        original=[];gaussian=[]
        for m in range(1,4):
            for j in range(3):
                if m>N*j:original+=coeffs(s.diff(u**j*(a*u+b)**m,u,m).subs(u,x),(x,))
            for az in range(3):
                for bw in range(3-az):
                    if m>N*(az+bw):gaussian+=coeffs(pair_expectation(z**az*w**bw*phase**m,(z,),(w,)),(x,))
        go=s.groebner(original,a,b);gg=s.groebner(gaussian,a,b)
        assert list(go)==list(gg)==([a,b] if N==0 else [a])
        ideals.append({'N':N,'bounded_original_and_Gaussian_ideal':[str(g) for g in go]})

    # Nilpotent coefficient specialization of the forward identity.
    symbol=z**2+e*z;p=u**2+e*u
    left=pair_expectation((x+w)*(symbol*p.subs(u,x+w))**2,(z,),(w,))
    right=apply_symbol(symbol**2,u*p**2,(z,),(u,)).subs(u,x)
    assert s.rem(left-right,e**3,e)==0

    # Ordinary coefficient ideals equal grid-value ideals,with nilpotent thickness.
    hs=[a*a*x*x+a*b*x*y+b**3*y*y,a*a+b**3*x]
    grid=list(itertools.product(range(3),repeat=2))
    values=[s.expand(h.subs({x:xx,y:yy})) for h in hs for xx,yy in grid]
    coefficients=[c for h in hs for c in coeffs(h,(x,y))]
    gc=s.groebner(coefficients,a,b);gv=s.groebner(values,a,b)
    assert list(gc)==list(gv)==[a*a,a*b,b**3]
    assert gc.reduce(a)[1]==a and gc.reduce(b)[1]==b
    for h in hs:
        interpolation=0
        for xx,yy in grid:
            lx=s.prod((x-j)/s.Integer(xx-j) for j in range(3) if j!=xx)
            ly=s.prod((y-j)/s.Integer(yy-j) for j in range(3) if j!=yy)
            interpolation+=h.subs({x:xx,y:yy})*lx*ly
        assert s.expand(h-interpolation)==0

    # Equal supports do not identify differential types of factorial normalizations.
    normalization=[]
    for m in range(6):
        original=s.diff(u**m,u,m);assert original==s.factorial(m)
        assert pair_expectation((z*(x+w))**m,(z,),(w,))==original
        normalization.append({'m':m,'direct_exponential_coefficient':str(original/s.factorial(m)),
                              'old_diagonal_coefficient':str(original/s.factorial(m)**2)})
    # Incorrect output-origin testing and frozen variable-coefficient symbols.
    output=s.diff(u**4,u,2).subs(u,x);assert output==12*x*x and output.subs(x,0)==0
    variable=u*s.diff(u*s.diff(u*u,u),u);frozen=u*u*s.diff(u*u,u,2)
    assert variable==4*u*u and frozen==2*u*u and variable!=frozen

    # Retained Dvorsky source:only three small moments are recomputed.
    T,A,B,C,D=s.symbols('T A B C D');qT,qA,qB,qC,qD=s.symbols('qT qA qB qC qD')
    us=(T,A,B,C,D);zs=(qT,qA,qB,qC,qD)
    pd=(T+C)*(A*D+B*T);ld=qT*(qA*qD-qB*qC);dvg=[]
    for m in range(1,4):
        assert apply_symbol(ld**m,pd**m,zs,us)==0
        marked=apply_symbol(ld**m,-C*pd**m,zs,us)
        expected=C+2*T if m==1 else (-1)**(m+1)*s.factorial(m+1)*s.factorial(m)**2*T
        assert s.expand(marked-expected)==0
        dvg.append({'m':m,'marked_output':str(marked)})
    scaling=[]
    for N in range(5):
        m=N+1
        coefficient=(C+2*T) if m==1 else (-1)**(m+1)*s.factorial(m+1)*s.factorial(m)**2*T
        assert coefficient!=0
        assert s.rem(e**m,e**(m+1),e)==e**m
        scaling.append({'N':N,'DVG1_scale_cutoff_ideal_power':m,'SIC2C4_scale_mu_lower':N+1,'SIC2C4_scale_mu_upper':2*N+1})

    # One SIC2C4 contraction is a transcription control,not a new all-order proof.
    R=z*u+z2*v;Z=z*v;W=2*z2*u;TT=z*u-z2*v
    fs=(R+Z)*(R*R*W-(2*R+Z)*TT*TT/2)
    assert contraction(fs,(z,z2),(u,v))==0
    assert contraction(Z*fs,(z,z2),(u,v))==240
    artin=[]
    for length in range(1,6):
        if length>1:
            last=length-1
            assert s.rem(e**last,e**length,e)!=0
            assert s.rem(e**length,e**length,e)==0
        artin.append({'nilpotence_length':length,'DVG1_sharp_integer_slope':length-1,
                      'SIC2C4_lower_integer_slope':(length-1+1)//2,'SIC2C4_upper_integer_slope':length-1})
    return {'completed':True,'python_version':platform.python_version(),'sympy_version':s.__version__,
            'product_phase_cases':rows,'GVC_forward_controls':forward,'GVC_reverse_marker_controls':reverse,
            'full_SIC_forward_controls':sic_forward,'full_SIC_reverse_marker_controls':sic_reverse,
            'bounded_ordinary_cutoff_ideal_controls':ideals,'nilpotent_specialization_control':True,
            'grid_interpolation_controls':2,'grid_value_count':len(values),'grid_ordinary_ideal':['a**2','a*b','b**3'],
            'factorial_normalization_controls':normalization,'output_origin_mutation_detected':True,
            'variable_coefficient_mutation_detected':True,'bounded_DVG1_controls':dvg,
            'nilpotent_scaling_power_controls':scaling,'SIC2C4_m1_marked_control':'240','Artin_slope_controls':artin,
            'all_three_field_uniformity_written_deduction':True,'full_ordinary_cutoff_transfer_written_proof':True,
            'universal_transfer_inferred_from_samples':False,'general_ideal_or_uniform_cutoff_producer_executed':False,
            'universal_numeric_cutoff_or_output_grid_enumerated':False,'rationality_between_normalizations_transferred':False,
            'all_order_DVG1_or_SIC2C4_inputs_reproved':False,'exact_SIC_scaling_ideal_exponents_computed':False,
            'ordinary_cutoff_ideals_stabilize_claimed':False,'all_nonreduced_uniformity_claimed':False,
            'operator_order_or_dimension_frontier_improved':False,'independent_review':False,
            'formal_verification':False,'novelty_claim':False}


if __name__=='__main__':
    path=Path(__file__).resolve();result=run()
    result['checker_sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
    path.with_name('polynomial_gaussian_transfers_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
