"""Exact retained group-law DAGs over finite rational-function fields.

Discovery uses normalized Sage elliptic-curve arithmetic. Replay derives
secant/tangent numerators and checks coordinate identities by polynomial cross multiplication, without elliptic
addition, rational-function normalization, inversion, or a point search.
"""

class RationalPair:
    def __init__(self,n,d):
        if not d:raise ArithmeticError('zero rational-function denominator')
        self.n,self.d=n,d
    def __add__(self,other):return RationalPair(self.n*other.d+other.n*self.d,self.d*other.d)
    def __neg__(self):return RationalPair(-self.n,self.d)
    def __sub__(self,other):return self+-other
    def __mul__(self,other):return RationalPair(self.n*other.n,self.d*other.d)
    def __pow__(self,n):return RationalPair(self.n**n,self.d**n)
    def __eq__(self,other):return self.n*other.d==other.n*self.d
    def inverted_parameter(self,weight):
        R=self.n.parent();t=R.gen()
        def reverse(f):return R(list(reversed(f.list())))
        if not self.n:return RationalPair(R(0),R(1))
        exponent=weight+self.d.degree()-self.n.degree()
        return RationalPair(reverse(self.n)*t**max(0,exponent),reverse(self.d)*t**max(0,-exponent))


def poly_record(p):
    # Compact coefficients keep the retained graph from duplicating millions
    # of Python integer objects beside Sage's normalized discovery points.
    return ','.join(map(str,map(int,p.list()))) or '0'
def polynomial(R,row):return R([int(c) for c in row.split(',')]) if isinstance(row,str) else R(row)
def rational_record(value):return [poly_record(value.numerator()),poly_record(value.denominator())]
def point_record(point):
    if point.is_zero():return None
    X,Y=point[0],point[1];h=X.denominator().sqrt()
    if h.leading_coefficient()!=1:h=-h
    if h*h!=X.denominator() or h**3!=Y.denominator():
        raise ArithmeticError('nonintegral short-model denominator frame')
    return {'jacobian':[poly_record(X.numerator()),poly_record(Y.numerator()),poly_record(h)]}
def read_rational(R,row):
    if len(row)!=2:raise ValueError('invalid rational-function witness')
    return RationalPair(polynomial(R,row[0]),polynomial(R,row[1]))
def read_point(R,row):
    if row is None:return None
    if isinstance(row,dict):
        x,y,z=(polynomial(R,c) for c in row['jacobian'])
        return RationalPair(x,z*z),RationalPair(y,z**3)
    return tuple(read_rational(R,r) for r in row)


class GroupLawDiscovery:
    def __init__(self,curve,basis,*,retained=None):
        self.curve,self.basis=curve,tuple(basis)
        self.dimension=len(basis);self.nodes=[];self.points=[];self.words=[];self.index={}
        if retained is not None:
            R=curve.base_ring().ring()
            replay_group_law(R,R(curve.a4()),R(curve.a6()),basis,retained)
            for row in retained:
                word=tuple(row['word']);point=row['point']
                if point is None:P=curve(0)
                else:
                    K=curve.base_ring();P=curve(*[K(r.n/r.d) for r in read_point(R,point)])
                self.index[word]=len(self.nodes);self.nodes.append(row);self.points.append(P);self.words.append(word)
        if not self.nodes:
            self._store((0,)*self.dimension,curve(0),{'kind':'identity'})
            for i,P in enumerate(basis):
                self._store(tuple(int(j==i) for j in range(self.dimension)),P,{'kind':'basis','basis_index':i})

    def _store(self,word,point,evidence):
        index=len(self.nodes);word=tuple(word)
        if word in self.index:return self.index[word]
        self.index[word]=index;self.points.append(point);self.words.append(word)
        self.nodes.append({'word':list(word),'point':point_record(point),**evidence})
        return index

    def negate(self,index):
        word=tuple(-c for c in self.words[index])
        return self._store(word,-self.points[index],{'kind':'negative','source':index})

    def add(self,left,right):
        word=tuple(a+b for a,b in zip(self.words[left],self.words[right]))
        if word in self.index:return self.index[word]
        P,Q=self.points[left],self.points[right];result=P+Q
        return self._store(word,result,{'kind':'sum','left':left,'right':right})

    def trace(self,word):
        word=tuple(map(int,word))
        if len(word)!=self.dimension:raise ValueError('wrong trace-word dimension')
        if word in self.index:return self.index[word]
        support=[i for i,c in enumerate(word) if c]
        i=support[-1];c=word[i]
        if len(support)==1:
            if c<0:
                positive=tuple(-x for x in word);return self.negate(self.trace(positive))
            left=list(word);left[i]=c//2
            right=list(word);right[i]=c-c//2
            return self.add(self.trace(left),self.trace(right))
        prefix=list(word);prefix[i]=0
        last=[0]*self.dimension;last[i]=c
        return self.add(self.trace(prefix),self.trace(last))

    def trace_in_basis(self,coefficients,basis_nodes):
        coefficients=tuple(map(int,coefficients))
        if len(coefficients)!=len(basis_nodes):raise ValueError('wrong transported basis dimension')
        word=tuple(sum(c*self.words[j][i] for c,j in zip(coefficients,basis_nodes)) for i in range(self.dimension))
        if word in self.index:return self.index[word]
        support=[i for i,c in enumerate(coefficients) if c]
        i=support[-1];c=coefficients[i]
        if len(support)==1:
            if c==1:return basis_nodes[i]
            if c<0:return self.negate(self.trace_in_basis(tuple(-v for v in coefficients),basis_nodes))
            left=list(coefficients);left[i]=c//2
            right=list(coefficients);right[i]=c-c//2
        else:
            left=list(coefficients);left[i]=0
            right=[0]*len(coefficients);right[i]=c
        return self.add(self.trace_in_basis(left,basis_nodes),self.trace_in_basis(right,basis_nodes))


def replay_group_law(R,A,B,basis,nodes):
    """Verify the labelled DAG; roots can then be referenced by node number."""
    dimension=len(basis)
    constants={i:RationalPair(R(i),R(1)) for i in (2,3)}
    AA=RationalPair(A,R(1));BB=RationalPair(B,R(1))
    expected=[read_point(R,point_record(P)) for P in basis]
    points=[];words=[];seen=set()
    for index,row in enumerate(nodes):
        word=tuple(row['word']);point=read_point(R,row['point']);kind=row['kind']
        if len(word)!=dimension or any(type(c) is not int for c in word) or word in seen:
            raise ArithmeticError('invalid or duplicate labelled group-law node')
        seen.add(word)
        def parent(name):
            j=row[name]
            if type(j) is not int or not 0<=j<index:raise ArithmeticError('non-topological group-law witness')
            return points[j],words[j]
        if kind=='identity':
            if any(word) or point is not None:raise ArithmeticError('invalid identity node')
        elif kind=='basis':
            i=row['basis_index']
            if type(i) is not int or not 0<=i<dimension or word!=tuple(int(j==i) for j in range(dimension)) or point!=expected[i]:
                raise ArithmeticError('basis point differs from exact model input')
            if point is None or point[1]**2!=point[0]**3+AA*point[0]+BB:raise ArithmeticError('base section is off curve')
        elif kind=='negative':
            P,v=parent('source')
            if word!=tuple(-c for c in v) or point!=(None if P is None else (P[0],-P[1])):
                raise ArithmeticError('invalid inverse witness')
        elif kind=='sum':
            P,v=parent('left');Q,w=parent('right')
            if word!=tuple(a+b for a,b in zip(v,w)):raise ArithmeticError('sum has a different trace word')
            if P is None or Q is None:
                if point!=(Q if P is None else P):raise ArithmeticError('invalid identity addition')
            elif point is None:
                if P[0]!=Q[0] or P[1]!=-Q[1]:raise ArithmeticError('invalid vertical chord')
            else:
                if P[0]==Q[0]:
                    if P[1]!=Q[1] or not P[1].n:raise ArithmeticError('invalid tangent')
                    numerator=constants[3]*P[0]**2+AA;denominator=constants[2]*P[1]
                else:numerator=Q[1]-P[1];denominator=Q[0]-P[0]
                slope=RationalPair(numerator.n*denominator.d,numerator.d*denominator.n)
                if 'slope' in row and slope!=read_rational(R,row['slope']):raise ArithmeticError('invalid retained slope')
                if point[0]!=slope**2-P[0]-Q[0] or point[1]!=slope*(P[0]-point[0])-P[1]:
                    raise ArithmeticError('group-law coordinate identity failed')
        else:raise ValueError('unknown group-law witness kind')
        points.append(point);words.append(word)
    return points,words


def chord_witness(trace,A,R,K,helper):
    X,Y=K(trace[0]),K(trace[1]);frame=helper['trace_chord_frame'](X,Y,R);inverted=False
    if frame['h'].degree()!=3:
        X=helper['invert_rational'](X,4,R,K);Y=helper['invert_rational'](Y,6,R,K)
        frame=helper['trace_chord_frame'](X,Y,R);inverted=True
    if frame['h'].degree()!=3:raise ArithmeticError('trace poles do not expose the chord')
    localA=helper['reciprocal_with_bound'](A,8,R) if inverted else A
    h,Nx,Ny,M=(frame[k] for k in ('h','Nx','Ny','M0'))
    q,remainder=(M**4-6*M*M*Nx-8*M*Ny-3*Nx*Nx-4*localA*h**4).quo_rem(h**6)
    if remainder:raise ArithmeticError('invalid chord quotient')
    if inverted:q=helper['reciprocal_with_bound'](q,2,R)
    return {'inverted':inverted,'h':poly_record(h),'Nx':poly_record(Nx),'Ny':poly_record(Ny),
            'M':poly_record(M),'q':poly_record(q)}


def replay_chord(R,A,delta,point,witness):
    if point is None:raise ArithmeticError('a chord cannot be attached to the identity')
    X,Y=point;inverted=witness['inverted']
    if type(inverted) is not bool:raise ValueError('invalid affine chart flag')
    def reciprocal(p,bound):
        if p.degree()>bound:raise ArithmeticError('reciprocal degree bound failed')
        return R([p[bound-i] for i in range(bound+1)])
    if inverted:X,Y=X.inverted_parameter(4),Y.inverted_parameter(6)
    h,Nx,Ny,M,q=(polynomial(R,witness[k]) for k in ('h','Nx','Ny','M','q'))
    if h.degree()!=3 or h.leading_coefficient()!=1 or X!=RationalPair(Nx,h*h) or Y!=RationalPair(Ny,h**3):
        raise ArithmeticError('chord denominator frame differs from its trace')
    modulus=h*h
    if (Nx.gcd(modulus)!=1 or M.degree()>=6 or (M*Nx+Ny)%modulus or (M*M-Nx)%modulus or M.gcd(h).degree()!=0):
        raise ArithmeticError('invalid regular chord slope witness')
    localq=reciprocal(q,2) if inverted else q
    localA=reciprocal(A,8) if inverted else A
    if M**4-6*M*M*Nx-8*M*Ny-3*Nx*Nx-4*localA*h**4!=h**6*localq:
        raise ArithmeticError('chord quotient identity failed')
    if q.degree()!=2 or q.gcd(q.derivative()).degree()!=0 or q.gcd(delta).degree()!=0:
        raise ArithmeticError('invalid finite bisection validation')
    return q
