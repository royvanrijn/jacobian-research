"""Physical auxiliary Riordan structure and uniform non-TN obstruction.

The proofs apply to every integer gap d>=1. Finite replay controls below
check definitions, inverse identities, and explicit negative minors; they
are not used to extrapolate the uniform theorem.
"""
from pathlib import Path
from math import factorial
import hashlib,json,platform
import sympy as s
from higher_gap_adjoint_rank import auxiliary, gap_three_W


def inverse_polynomial(d,u):
    if not isinstance(d,int) or d<1:
        raise ValueError('gap must be a positive integer')
    return s.expand(((1-u)-(1-u)**(d+1))/d)


def inverse_prefactor(d,beta,u):
    return s.diff(inverse_polynomial(d,u),u)*(1-u)**(s.sympify(beta)-2)


def physical_array(N,d,beta,z):
    return s.Matrix(N+1,N+1,lambda n,k: auxiliary(n,n+beta-1,d,z).nth(k) if k<=n else 0)


def inverse_array(N,d,beta,u):
    P=inverse_polynomial(d,u)
    g=inverse_prefactor(d,beta,u)
    columns=[s.Poly(s.series(g*P**k,u,0,N+1).removeO().expand(),u) for k in range(N+1)]
    return s.Matrix(N+1,N+1,lambda n,k:s.expand(s.factorial(n)*columns[k].nth(n)/s.factorial(k)))


def production_data(d,beta,u):
    A=1/s.diff(inverse_polynomial(d,u),u)
    Z=s.diff(A,u)+(s.sympify(beta)-2)*A/(1-u)
    return s.cancel(A),s.cancel(Z)


def production_array(N,d,beta,u):
    A,Z=production_data(d,beta,u)
    aa=s.Poly(s.series(A,u,0,N+2).removeO().expand(),u)
    zz=s.Poly(s.series(Z,u,0,N+1).removeO().expand(),u)
    at=lambda row,k:row.nth(k) if k>=0 else s.Integer(0)
    return s.Matrix(N+1,N+2,lambda i,j:s.expand(s.factorial(i)*(at(zz,i-j)/s.factorial(j)
             +(at(aa,i-j+1)/s.factorial(j-1) if j>=1 else 0))) if i-j>=-1 else 0)


def reverted_series(N,d,u,t):
    P=inverse_polynomial(d,u)
    return sum(s.Rational(1,n)*s.expand(s.series((u/P)**n,u,0,n).removeO()).coeff(u,n-1)*t**n
               for n in range(1,N+1))


def fixed_gamma_array(N,d,gamma,z):
    return s.Matrix(N+1,N+1,lambda n,k:auxiliary(n,gamma,d,z).nth(k) if k<=n else 0)


def fixed_gamma_production_factors(N,d,gamma):
    c=s.sympify(gamma)+d-1
    T=s.Matrix(N+1,N+1,lambda i,j:s.binomial(i-j+d-1,d-1) if i>=j else 0)
    U=s.Matrix(N+1,N+1,lambda i,j:c-j if i==j else (j if i==j-1 else 0))
    D=s.diag(*(s.factorial(j) for j in range(N+1)))
    return D,T,U


def cyclic_row(polynomial,z,u):
    p=s.Poly(polynomial,z)
    return [s.Poly(sum(p.nth(k)*u**((k-r)//3) for k in range(r,p.degree()+1,3)),u)
            for r in range(3)]


def cyclic_plucker_sections(B,C,z,u):
    """Cyclic sections of the rotated divided difference, in their actual variable u."""
    A=cyclic_row(B,z,u);E=cyclic_row(C,z,u)
    return [s.Poly(s.expand((A[i].as_expr()*E[j].as_expr()-A[j].as_expr()*E[i].as_expr()).subs(u,-u)),u)
            for i,j in [(0,1),(0,2),(1,2)]]


def check_cyclic_plucker_identity(B,C,W,z,u):
    actual=cyclic_row(W,z,u)
    expected=cyclic_plucker_sections(B,C,z,u)
    return all(s.expand(a.as_expr()-b.as_expr())==0 for a,b in zip(actual,expected))


def sine_weight_polynomial(B,C,z):
    B,C=s.Poly(B,z),s.Poly(C,z)
    N=max(B.degree(),C.degree())
    chi=(0,1,1,0,-1,-1)
    return s.Poly(s.expand(sum(chi[(k-j)%6]*(B.nth(j)*C.nth(k)-B.nth(k)*C.nth(j))*z**(j+k-1)
                              for j in range(N+1) for k in range(j+1,N+1))),z)


def cyclic_hurwitz_minor(index,descending_coefficients):
    if index==0:
        return s.Integer(1)
    if index<0:
        raise ValueError('minor index must be nonnegative')
    size=(index+1)//2
    first_row=2 if index%2 else 1
    at=lambda k:descending_coefficients[k] if 0<=k<len(descending_coefficients) else s.Integer(0)
    return s.expand(s.Matrix(size,size,lambda i,j:at(3*(j+1)-(first_row+i))).det(method='domain-ge'))


def cyclic_tail_certificate(B,C,z,t):
    """Identity-only certificate; permits zero constants and zero determinants."""
    B,C=s.Poly(B,z),s.Poly(C,z)
    n=B.degree()
    if n<1 or C.degree()!=n+1 or B.LC()!=1 or C.LC()!=1:
        raise ValueError('expected monic consecutive degrees n,n+1 with n>=1')
    W=sine_weight_polynomial(B.as_expr(),C.as_expr(),z)
    w=W.all_coeffs()
    reverse=s.Poly(sum(w[j]*t**j for j in range(2*n+1)),t)
    reverseB=s.Poly(sum((-1)**(n+k)*B.nth(k)*t**(n-k) for k in range(n+1)),t)
    product=reverse*reverseB
    rows=[s.Matrix([[w[3*j-i] if 0<=3*j-i<len(w) else 0 for j in range(1,n+1)]])
          for i in range(2,n+3)]
    relation=rows[0]+sum((reverseB.nth(i-2)*rows[i-2] for i in range(3,n+3)),s.zeros(1,n))
    previous=cyclic_hurwitz_minor(2*n-2,w)
    odd=cyclic_hurwitz_minor(2*n-1,w)
    even=cyclic_hurwitz_minor(2*n,w)
    return {'missing_residue':all(s.expand(product.nth(3*j-2))==0 for j in range(1,n+1)),
            'row_relation':all(s.expand(v)==0 for v in relation),
            'odd_identity':s.expand(odd-B.nth(0)*previous)==0,
            'even_identity':s.expand(even-W.nth(0)*previous)==0,
            'previous':previous,'odd':odd,'even':even,'B_constant':B.nth(0),'W_constant':W.nth(0)}


def main():
    HERE=Path(__file__).resolve().parent
    counts={}
    def check(label,condition):
        assert condition,label
        counts[label]=counts.get(label,0)+1
    u,t,z,beta,gamma,x=s.symbols('u t z beta gamma x')
    profiles=[]
    for d in range(1,6):
        R=physical_array(5,d,beta,z)
        S=inverse_array(5,d,beta,u)
        check('symbolic_shape_Riordan_inverse_replayed',all(s.expand(v)==0 for v in R*S-s.eye(6)))
        f=reverted_series(5,d,u,t)
        check('formal_reversion_input',s.series(inverse_polynomial(d,u).subs(u,f)-t,t,0,6).removeO().expand()==0)
        A,Z=production_data(d,beta,u)
        P=inverse_polynomial(d,u)
        check('production_A_identity',s.cancel(A*s.diff(P,u)-1)==0)
        check('inverse_prefactor_elementary_identity',s.expand(s.diff(P,u)
              -((d+1)*(1-u)**d-1)/d)==0)
    for d in range(1,6):
        R=fixed_gamma_array(5,d,gamma,z)
        D,T,U=fixed_gamma_production_factors(4,d,gamma)
        Pi=D*T*U*D.inv()
        direct=R[:5,:5].inv()*R[1:6,:5]
        check('fixed_gamma_production_factorization_identity',all(s.expand(v)==0 for v in direct-Pi))
        G=s.Matrix(5,5,lambda i,j:1 if i>=j else 0)
        check('fixed_gamma_Toeplitz_factor_is_geometric_power',T==G**d)
        Ds,Ts,Us=fixed_gamma_production_factors(4,d,5-d+x)
        check('fixed_gamma_bidiagonal_nonnegative_on_required_range',
              all(all(c>=0 for c in s.Poly(v,x).all_coeffs()) for v in Us))
    R=fixed_gamma_array(5,3,gamma,z)
    D,T,U=fixed_gamma_production_factors(5,3,gamma)
    Pi=D*T*U*D.inv()
    Q=s.Matrix.vstack(s.Matrix([[1,0,0,0,0,0]]),Pi[:5,:])
    check('fixed_gamma_recursive_TN_factorization_identity',
          all(s.expand(v)==0 for v in s.diag(s.Integer(1),R[:5,:5])*Q-R))
    for n in range(1,6):
        B=auxiliary(n,gamma,3,z);C=auxiliary(n+1,gamma,3,z)
        W=gap_three_W(B.as_expr(),z)
        check('actual_adjacent_auxiliary_raising_identity',
              s.expand(C.as_expr()-(z+3*n+gamma+2)*B.as_expr()+z*B.diff().as_expr())==0)
        check('uniform_shape_cyclic_Plucker_transfer_identity',
              check_cyclic_plucker_identity(B.as_expr(),C.as_expr(),W.as_expr(),z,u))
        check('uniform_shape_sine_weight_transfer_identity',
              s.expand(sine_weight_polynomial(B.as_expr(),C.as_expr(),z).as_expr()-W.as_expr())==0)
        if n==3:
            check('wrong_adjacent_row_transfer_rejected',
                  not check_cyclic_plucker_identity(B.as_expr(),C.as_expr()+1,W.as_expr(),z,u))
    B=1+2*z+3*z**2+5*z**3+7*z**4
    C=11+13*z+17*z**2+19*z**3+23*z**4+29*z**5
    alpha=(1+s.sqrt(-3))/2;alpha_bar=(1-s.sqrt(-3))/2
    divided=s.Poly(s.cancel(s.expand(C.subs(z,alpha*z)*B.subs(z,alpha_bar*z)
                      -B.subs(z,alpha*z)*C.subs(z,alpha_bar*z))/(s.sqrt(-3)*z)),z)
    check('generic_divided_difference_cyclic_Plucker_identity',
          check_cyclic_plucker_identity(B,C,divided.as_expr(),z,u))
    check('generic_divided_difference_sine_weight_identity',
          s.expand(sine_weight_polynomial(B,C,z).as_expr()-divided.as_expr())==0)
    # Independent symbolic pairs: these need not satisfy an auxiliary recurrence.
    for n in [1,2,3]:
        bb=s.symbols('bb0:'+str(n));cc=s.symbols('cc0:'+str(n+1))
        B=z**n+sum(bb[k]*z**k for k in range(n))
        C=z**(n+1)+sum(cc[k]*z**k for k in range(n+1))
        cert=cyclic_tail_certificate(B,C,z,t)
        for label in ['missing_residue','row_relation','odd_identity','even_identity']:
            check('generic_symbolic_tail_'+label,cert[label])
    for label,B,C in [
          ('zero_B_constant',z+2*z**2+z**3,1+3*z+5*z**2+7*z**3+z**4),
          ('zero_previous_minor',1+2*z+z**2,3+4*z+2*z**2+z**3),
          ('monomial_degeneracy',z**2,z**3)]:
        cert=cyclic_tail_certificate(B,C,z,t)
        check('degenerate_tail_identities_'+label,all(cert[k] for k in
              ['missing_residue','row_relation','odd_identity','even_identity']))
        if label=='zero_B_constant':check('zero_B_constant_is_retained',cert['B_constant']==0)
        if label=='zero_previous_minor':check('zero_previous_minor_is_retained',cert['previous']==0)
    # TN auxiliary blocks alone do not force TN of the cyclic Wronskian.
    polys=[s.Poly(1,z)]
    for j in range(9):
        B=polys[-1].as_expr()
        polys.append(s.Poly(s.expand((z+9)*B-z*s.diff(B,z)),z))
    constP=s.Matrix(10,10,lambda i,j:9-i if i==j else (1 if j==i+1 else 0))
    constR=s.Matrix(10,10,lambda i,j:polys[i].nth(j))
    check('off_family_TN_bidiagonal_production',all(c>=0 for c in constP))
    check('off_family_bidiagonal_generates_array',constR[1:,:]==constR[:9,:]*constP)
    W=gap_three_W(polys[8].as_expr(),z);C=W.all_coeffs()
    minor=s.Matrix(3,3,lambda i,j:C[3*(j+1)-(i+1)]).det()
    check('off_family_TN_auxiliary_does_not_transfer_to_W_H3',minor==-60115263786123)
    # Independent direct derivative-array replay of the production formula.
    for d in [1,2,3,5]:
        R=physical_array(5,d,beta,z)
        S=inverse_array(4,d,beta,u)
        direct=S*R[1:6,:6]
        P=production_array(4,d,beta,u)
        check('symbolic_shape_production_formula_replayed',all(s.expand(v)==0 for v in direct-P))
    # A low-order direct EGF reconstruction, rather than only an inverse test.
    f=reverted_series(4,3,u,t)
    G=s.series(s.diff(f,t)*(1-f)*s.exp(z*f),t,0,4).removeO().expand()
    for n in range(4):
        check('physical_EGF_direct_reconstruction',s.expand(s.factorial(n)*G.coeff(t,n)-auxiliary(n,n,3,z).as_expr())==0)
    # Infinite-family proof is in the note; these are bounded exact controls.
    for d in range(1,9):
        m=d if d%2==0 else d+1
        for order in [m,m+2]:
            R=physical_array(order,d,s.Integer(1),z)
            minor=R[1:order+1,:order].det(method='domain-ge')
            expected=-s.Rational(factorial(order),d)
            check('uniform_negative_minor_bounded_controls',minor==expected and minor<0)
            profiles.append({'d':d,'beta':1,'rows':[1,order],'columns':[0,order-1],
                             'minor':str(minor),'uniform_formula':'-m!/d for every even m>=d'})
    R=physical_array(4,3,beta,z)
    minor=s.factor(R[1:5,:4].det(method='domain-ge'))
    check('symbolic_shape_first_gap_three_negative_minor',s.expand(minor-(beta-2)*(beta**3+4*beta**2-17*beta+20))==0)
    R1=physical_array(4,3,s.Integer(1),z)
    expected=s.Matrix([[3,1,0,0],[28,10,1,0],[440,160,21,1],[9720,3560,520,36]])
    check('actual_gap_three_negative_array_minor',R1[1:5,:4]==expected and expected.det()==-8)
    P=production_array(3,3,s.Integer(1),u)
    block=P[[2,3],[0,1]]
    check('actual_gap_three_negative_production_minor',block==s.Matrix([[166,62],[1874,698]]) and block.det()==-320)
    changed=expected.copy();changed[3,0]+=1
    check('corrupted_negative_minor_control_rejected',changed.det()!=-8)
    try:
        inverse_polynomial(0,u)
        rejected=False
    except ValueError:
        rejected=True
    check('invalid_gap_rejected',rejected)
    receipt={'status':'pass','canonical_proof':'gap_three_hurwitz_structure.tex',
      'command':'python research/factorial/image-rigidity/gap_three_hurwitz_structure.py',
      'python':platform.python_version(),'sympy':s.__version__,
      'uniform_physical_auxiliary_Riordan_EGF_proved':True,
      'uniform_Riordan_inverse_proved':True,
      'fixed_gamma_finite_positive_production_factorization_proved':True,
      'exact_cyclic_Plucker_Wronskian_transfer_proved':True,
      'exact_sine_weight_Wronskian_transfer_proved':True,
      'uniform_odd_and_even_tail_special_minor_identities_proved':True,
      'tail_scope':'every monic B,C of degrees n,n+1, n>=1: Delta_(2n-1)=B(0)Delta_(2n-2), Delta_(2n)=W(0)Delta_(2n-2)',
      'tail_proof_uses_division_by_constant_or_minor':False,
      'fixed_gamma_array_TN_scope':'leading rows and columns0..N whenever gamma+d-1>=N',
      'TN_auxiliary_block_alone_transfers_to_Wronskian_H3_refuted':True,
      'off_family_Wronskian_negative_H3_minor':-60115263786123,
      'uniform_negative_auxiliary_array_minors_proved':True,
      'negative_minor_scope':'every integer d>=1, beta=1, every even m>=d: det R[1..m,0..m-1]=-m!/d',
      'actual_gap_three_production_TN_refuted':True,
      'H3_Wronskian_TN_refuted':False,'all_index_gap_three_SFC_closed':False,
      'unrestricted_binomial_SFC_closed':False,
      'bounded_replay_profiles':profiles,'negative_gap_three_array_minor':-8,
      'negative_gap_three_production_minor':-320,
      'finite_extrapolation':False,'independent_replay':False,'formal_verification':False,
      'external_review':False,'novelty_claim':False,'total':sum(counts.values()),'counts':counts,
      'files':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
        ['gap_three_hurwitz_structure.tex','gap_three_hurwitz_structure.py','higher_gap_adjoint_rank.py']}}
    (HERE/'gap_three_hurwitz_structure_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()
