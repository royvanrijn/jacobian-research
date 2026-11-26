"""Verify all normal charts at the w=0 endpoint of L3 with SageMath 10.9."""
from math import comb
from sage.all import GF, PolynomialRing, QQ, matrix, vector

from verify_cm_newton_full_width_reconstruction_certificates import (
    reduced_residual_rows, weighted_exponents,
)


def weight(exponent, weights):
    return sum(a*b for a,b in zip(exponent, weights, strict=True))


def local_rows(rows):
    names = ("r", "L", "H", "w", "tau")
    ring = PolynomialRing(QQ, names=names)
    r, l, h, w, tau = ring.gens()
    result = []
    for shift, row in enumerate(rows, 7):
        result.append(ring.sum(
            coefficient*r**e[0]*(l+8)**e[1]*(h+17*w)**e[2]*w**e[4]
            * tau**((shift-1-weight(e, (2,1,2,1,2)))//2)
            for e, coefficient in row.dict().items()
        ))
    return result


def initial_forms(rows, weights, residue):
    field = residue.base_ring()
    basis = {}
    initial = []
    for row in rows:
        terms = dict(row.dict())
        while terms:
            leading = min(terms, key=lambda e: (weight(e, weights), e))
            coefficient = terms[leading]
            if leading not in basis:
                assert field(coefficient)
                normalized = {e: c/coefficient for e,c in terms.items()}
                basis[leading] = normalized
                degree = weight(leading, weights)
                initial.append(residue({e:field(c) for e,c in normalized.items()
                                        if weight(e, weights)==degree}))
                break
            for e,c in basis[leading].items():
                value = terms.get(e, 0)-coefficient*c
                if value:
                    terms[e] = value
                else:
                    terms.pop(e, None)
    return initial


def verify_radial_weight_eight(rows, residue):
    weights = (2,1,1,1,1)
    cap = 8
    field = residue.base_ring()
    monomials = sorted((e for e in weighted_exponents(weights, cap)
                        if weight(e, weights)<cap or not e[4]),
                       key=lambda e: (weight(e, weights), e))
    positions = {e:j for j,e in enumerate(monomials)}
    low_count = sum(weight(e, weights)<cap for e in monomials)
    truncated = [{tuple(e):c for e,c in row.dict().items() if weight(e,weights)<=cap}
                 for row in rows]
    labels = []
    entries = {}
    for index, terms in enumerate(truncated):
        minimum = min(weight(e, weights) for e in terms)
        for shift in weighted_exponents(weights, cap-minimum):
            current = len(labels)
            labels.append((index, shift))
            for exponent, coefficient in terms.items():
                displacement = tuple(a+b for a,b in zip(exponent,shift,strict=True))
                if displacement in positions:
                    entries[current,positions[displacement]] = field(coefficient)
    coefficient_matrix = matrix(field,len(labels),len(monomials),entries)
    assert coefficient_matrix.dimensions() == (1162,591)
    target = vector(field,[int(e==(4,0,0,0,0)) for e in monomials])
    modular_solution = coefficient_matrix.solve_left(target)
    assert modular_solution*coefficient_matrix == target
    low_matrix = coefficient_matrix[:,:low_count]
    pivot_rows, pivot_columns = low_matrix.pivot_rows(), low_matrix.pivots()
    assert len(pivot_rows)==len(pivot_columns)==238

    def exact_lows(row):
        index,shift = labels[row]
        return vector(QQ,[
            truncated[index].get(tuple(a-b for a,b in zip(e,shift,strict=True)),0)
            for e in monomials[:low_count]
        ])

    seed = sum(QQ(int(modular_solution[i]))*exact_lows(i)
               for i in modular_solution.nonzero_positions())
    exact = matrix(QQ,[exact_lows(i) for i in pivot_rows])
    correction = exact.matrix_from_columns(pivot_columns).solve_left(
        vector(QQ,[-seed[j] for j in pivot_columns]))
    assert correction*exact == -seed
    assert all(not field(c) for c in correction)


def second_normalization_forms(rows, residue):
    # Scalar row operations alone leave points in the r=1 and tau=1
    # charts for the short prefix. Here use their monomial multiples.
    weights = (3, 2, 1, 1, 2)
    cap = 9
    field = residue.base_ring()
    monomials = sorted(
        weighted_exponents(weights, cap),
        key=lambda exponent: (weight(exponent, weights), exponent),
    )
    positions = {exponent: index for index, exponent in enumerate(monomials)}
    degrees = [weight(exponent, weights) for exponent in monomials]
    entries = {}
    row_count = 0
    for row in rows:
        terms = {e: c for e, c in row.dict().items() if weight(e, weights) <= cap}
        minimum = min(weight(e, weights) for e in terms)
        for multiplier in weighted_exponents(weights, cap - minimum):
            for exponent, coefficient in terms.items():
                shifted = tuple(a + b for a, b in zip(exponent, multiplier, strict=True))
                position = positions.get(shifted)
                if position is not None:
                    entries[row_count, position] = coefficient
            row_count += 1
    exact = matrix(QQ, row_count, len(monomials), entries)
    assert exact.dimensions() == (231, 342)
    modular = exact.change_ring(field)
    pivot_rows = modular.transpose().pivots()
    pivot_columns = modular.pivots()
    assert len(pivot_rows) == len(pivot_columns) == 62
    selected = exact.matrix_from_rows(pivot_rows)
    pivot = selected.matrix_from_columns(pivot_columns)
    echelon = pivot.solve_right(selected)
    assert pivot * echelon == selected
    for row, column in enumerate(pivot_columns):
        assert all(not echelon[row, earlier] for earlier in range(column))
        assert echelon[row, column] == 1
    reduced = echelon.change_ring(field)
    return [
        residue({
            monomials[k]: reduced[i, k] for k in range(j, len(monomials))
            if degrees[k] == degrees[j] and reduced[i, k]
        })
        for i, j in enumerate(pivot_columns)
    ]


def verify_first_normalizations(rows):
    base = rows[0].parent()
    field = GF(331)
    residue = PolynomialRing(field,names=base.variable_names())
    rho,x,y,z,u = residue.gens()
    _,l,h,_,_ = base.gens()

    first = initial_forms(rows,(2,1,1,1,1),residue)
    assert residue.ideal([*first,u-1]).groebner_basis()==[1]
    verify_radial_weight_eight(rows,residue)
    direction = residue.ideal([*first,rho,u])
    assert not ((x+124*y)**3).reduce(direction.groebner_basis())

    shifted = [f.subs({l:l+h/8}) for f in rows]
    second = initial_forms(shifted,(3,2,1,1,2),residue)
    refined = second_normalization_forms(shifted, residue)
    assert residue.ideal([*refined,rho-1]).groebner_basis()==[1]
    assert residue.ideal([*refined,u-1]).groebner_basis()==[1]
    assert residue.ideal([*second,rho,u,y,z,x-1]).groebner_basis()==[1]
    assert residue.ideal([*second,rho,u,z,y-1]).groebner_basis()==[1]
    final = residue.ideal([*second,rho,u,z-1])
    assert list(final.radical().groebner_basis())==[
        y*y+17*y, rho, x+95*y, z-1, u
    ]
    print("Exact checks for both normalization steps at the w=0 endpoint passed.")


def rescale(rows):
    ring = rows[0].parent()
    _, ell, h, _, _ = ring.gens()
    weights = (3, 2, 1, 1, 2)
    basis = {}
    for row in rows:
        terms = dict(row.subs({ell:ell+h/8}).dict())
        while terms:
            leading = min(terms, key=lambda e: (weight(e,weights),e))
            scalar = terms[leading]
            if leading not in basis:
                assert GF(331)(scalar)
                degree = weight(leading,weights)
                basis[leading] = {e:c/scalar for e,c in terms.items()}
                yield ring({(e[0],e[1],e[2],weight(e,weights)-degree,e[4]):c/scalar
                            for e,c in terms.items()})
                break
            for e,c in basis[leading].items():
                value = terms.get(e,0)-scalar*c
                if value: terms[e] = value
                else: terms.pop(e,None)


def homogeneous_part(poly,d):
    return poly.parent()({e:c for e,c in poly.dict().items() if sum(e)==d})


def branch_jets(rows):
    exact = PolynomialRing(QQ,"w").fraction_field()
    w = exact.gen()
    ring = PolynomialRing(exact,names=("X","Y","T"))
    x,y,t = ring.gens()
    N=7

    def truncate(poly,degree):
        return ring({e:c for e,c in poly.dict().items() if sum(e)<=degree})

    def evaluate(blocks,radial,degree):
        value,power = ring.zero(), ring.one()
        for f in blocks:
            value += truncate(f*power,degree)
            power = truncate(power*radial,degree)
            if not power: break
        return value

    expansions=[]
    for row in rows:
        blocks=[ring.zero() for _ in range(N+1)]
        for (r,i,j,k,l),c in row.dict().items():
            if r+i+j+l<=N:
                blocks[r]+=c*w**k*x**i*y**j*t**l
        expansions.append(blocks)

    pivot=expansions[0][1].constant_coefficient()
    assert pivot==QQ(7860239)/288028122 and GF(331)(pivot)
    radial=ring.zero()
    for degree in range(1,N+1):
        correction=-evaluate(expansions[0],radial,degree)/pivot
        assert correction==homogeneous_part(correction,degree)
        radial+=correction
    jets=[evaluate(block,radial,N) for block in expansions]
    assert not jets[0] and all(not truncate(f,2) for f in jets[1:])
    return exact,jets[1:]


def specialize331(coefficient):
    numerator = GF(331)(coefficient.numerator()(0))
    denominator = GF(331)(coefficient.denominator()(0))
    assert denominator
    return numerator/denominator


def witness(leading,jets,position,expected_degree,expected_norm):
    field=leading[0].base_ring()
    exact=PolynomialRing(field,names=("A","B"),order="lex")
    A,B=exact.gens()
    dehomogenize=leading[0].parent().hom(position(A,B,exact.one()),exact)
    leading=list(map(dehomogenize,leading))
    jets=list(map(dehomogenize,jets))

    linear,univariate=leading[0].parent().ideal(leading).groebner_basis()
    assert linear.degree(A)==1 and linear.monomial_coefficient(A)==1
    assert univariate.degree(A)==0
    coordinate=PolynomialRing(field,"b")
    b=coordinate.gen()
    monic=exact.hom([0,b],coordinate)(univariate)
    assert monic.degree()==expected_degree and monic.is_monic()
    exact_algebra=coordinate.quotient(monic,"b")
    beta=exact_algebra.gen()
    alpha=-exact.hom([0,beta],exact_algebra)(linear)
    exact_evaluate=exact.hom([alpha,beta],exact_algebra)
    assert all(not exact_evaluate(f) for f in leading)

    residue=PolynomialRing(GF(331),names=exact.variable_names())
    a,b=residue.gens()
    specialize=lambda f:residue({e:specialize331(c) for e,c in f.dict().items()})
    reduced_leading=list(map(specialize,leading))
    reduced_jets=list(map(specialize,jets))
    poly=PolynomialRing(GF(331),"b")
    m=poly({i:specialize331(c) for i,c in monic.dict().items()})
    assert m.is_monic() and m.is_squarefree()
    algebra=poly.quotient(m,"b")
    beta=algebra.gen()
    alpha=algebra(poly({i:specialize331(c) for i,c in alpha.lift().dict().items()}))
    ev=residue.hom([alpha,beta],algebra)
    assert all(not ev(f) for f in reduced_leading)
    residue_basis=residue.ideal(reduced_leading).groebner_basis()
    assert not specialize(linear).reduce(residue_basis)
    assert not specialize(univariate).reduce(residue_basis)
    selected=(0,2)
    pivot=matrix(algebra,[[ev(reduced_leading[i].derivative(z)) for z in (a,b)] for i in selected])
    assert pivot.det().is_unit()
    correction=pivot.inverse()*vector(algebra,[-ev(reduced_jets[i]) for i in selected])
    f=reduced_leading[4]
    remainder=ev(reduced_jets[4])+sum(ev(f.derivative(z))*c for z,c in zip((a,b),correction,strict=True))
    assert m.resultant(remainder.lift())==expected_norm


def verify_L3_direction(rows):
    field,jets=branch_jets(rows)
    cubic=[homogeneous_part(f,3) for f in jets]
    quartic=[homogeneous_part(f,4) for f in jets]
    mod=PolynomialRing(GF(331),names=jets[0].parent().variable_names())
    X,Y,T=mod.gens()
    specialize=lambda f:mod({e:specialize331(c) for e,c in f.dict().items()})
    assert mod.ideal([*map(specialize,cubic),T,X-1]).groebner_basis()==[1]
    witness(cubic,quartic,lambda a,b,s:[b,a,s],3,GF(331)(324))

    chart=PolynomialRing(field,names=("A","B","s"))
    a,b,s=chart.gens()
    rows=[chart.sum(c*a**e[0]*b**e[2]*s**(sum(e)-3) for e,c in f.dict().items()) for f in jets]
    assert all(not homogeneous_part(f,d) for f in rows for d in range(2))
    pivot=homogeneous_part(rows[0],2)
    scale=pivot.monomial_coefficient(a*a)
    assert specialize331(scale)
    residuals=[f-homogeneous_part(f,2).monomial_coefficient(a*a)/scale*rows[0] for f in rows[1:]]
    assert all(not homogeneous_part(f,2) for f in residuals)
    lead=[pivot, *(homogeneous_part(f,3) for f in residuals)]
    nxt=[homogeneous_part(rows[0],3),*(homogeneous_part(f,4) for f in residuals)]
    mod=PolynomialRing(GF(331),names=chart.variable_names())
    A,B,S=mod.gens()
    leading=[mod({e:specialize331(c) for e,c in f.dict().items()}) for f in lead]
    for z in (A,B):
        assert mod.ideal([*leading,S,z-1]).groebner_basis()==[1]
    witness(lead,nxt,lambda a,b,s:[a,b,s],6,GF(331)(285))
    print("Both subcharts of the L3-tangent direction at w=0 are excluded.")


def verify_L2_direction(rows):
    prime=GF(331)
    field=PolynomialRing(QQ,"a").fraction_field()
    a=field.gen()
    w_exact=a*(8-a)/17
    x_exact=-QQ(289)/(8*a*a)
    ring=PolynomialRing(field,names=("Y","T"))
    Y,T=ring.gens()

    def reduce(coefficient):
        coefficient=field(coefficient)
        numerator=prime(coefficient.numerator()(8))
        denominator=prime(coefficient.denominator()(8))
        assert denominator
        return numerator/denominator

    assert not reduce(w_exact) and reduce(x_exact)==prime(-QQ(289)/512)
    truncation=6
    blocks=[]
    for row in rows:
        coefficients={}
        for (r,x,y,w,t),c in row.dict().items():
            if r+2*t>truncation:
                continue
            constant=c*w_exact**w*T**t
            for i in range(min(x,truncation-r-2*t)+1):
                for j in range(min(y,truncation-r-2*t-i)+1):
                    coefficients[r,i]=coefficients.get((r,i),ring.zero()) + (
                        constant*comb(x,i)*x_exact**(x-i)*comb(y,j)*(-17)**(y-j)*Y**j
                    )
        blocks.append({e:c for e,c in coefficients.items() if c})

    def truncate(f,bound):
        return ring({e:c for e,c in f.dict().items() if e[0]+2*e[1]<=bound})

    def part(f,degree):
        return ring({e:c for e,c in f.dict().items() if e[0]+2*e[1]==degree})

    def evaluate(block,values,bound):
        powers=[[ring.one()],[ring.one()]]
        for degree in range(1,bound+1):
            for j in range(2):
                powers[j].append(truncate(powers[j][-1]*values[j],bound))
        return truncate(ring.sum(
            truncate(f*powers[0][i]*powers[1][j],bound)
            for (i,j),f in block.items() if i+j<=bound),bound)

    pivot=matrix(field,[[block[1,0].constant_coefficient(),block[0,1].constant_coefficient()]
                        for block in blocks[:2]])
    assert reduce(pivot.det())==255
    inverse=pivot.inverse()
    values=vector(ring,[0,0])
    for order in range(1,truncation+1):
        error=vector(ring,[evaluate(block,values,order) for block in blocks[:2]])
        assert all(part(f,order)==f for f in error)
        values-=inverse*error
    for value in values:
        for coefficient in value.dict().values():
            reduce(coefficient)

    jets=[evaluate(block,values,truncation) for block in blocks]
    assert not any(jets[:2])
    assert all(not truncate(f,4) for f in jets[2:])
    head=part(jets[3],5)
    winning=head.monomial_coefficient(Y*T*T)
    assert reduce(winning)==37
    cancelled=[]
    for row in jets[2:]:
        multiple=part(row,5).monomial_coefficient(Y*T*T)/winning
        reduce(multiple)
        f=row-multiple*jets[3]
        assert not truncate(f,5)
        cancelled.append(f)
    residue=PolynomialRing(prime,names=ring.variable_names())
    y,t=residue.gens()
    mod=lambda f:residue({e:reduce(c) for e,c in f.dict().items()})
    head_mod=mod(head)
    obstruction=mod(part(cancelled[3],6))
    assert head_mod==37*y*t*t+165*y**3*t+142*y**5
    assert obstruction==119*t**3
    polynomial=PolynomialRing(prime,"u")
    u=polynomial.gen()
    chart=residue.hom([1,u],polynomial)
    assert chart(head_mod).resultant(chart(obstruction))==325
    print("The L2-tangent direction at w=0 is excluded by the exact weighted-six jets.")


def verify_wzero_endpoint(rows):
    assert len(rows) == 7
    original = local_rows(rows)
    verify_first_normalizations(original)
    scaled = list(rescale(original))
    verify_L3_direction(scaled)
    verify_L2_direction(scaled)
    print("At 331: the w=0 endpoint of L3 is excluded using S7 through S13.")


if __name__ == "__main__":
    verify_wzero_endpoint(reduced_residual_rows(5, 13, even_linear=True, primitive=False))
