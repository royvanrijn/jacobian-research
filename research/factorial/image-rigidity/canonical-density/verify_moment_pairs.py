#!/usr/bin/env python3
"""Independent standard-library replay of the explicit moment/Bezout witnesses."""
import argparse,json
from hashlib import sha256
from fractions import Fraction as F
from math import comb,factorial,isqrt
from pathlib import Path
HERE=Path(__file__).resolve().parent
PROOF_NAME='moment_pair_obstructions.md'
PROOF_SHA256='2f40b43023e38b83f160420c44983305f810225964a29696ae6340711da91322'
CERTIFICATE_SHA256='0f3e8f68902d838a3038e7913407ef76c2621d72dbb3e59d44cfb45513bec378'

def require(condition,message):
    if not condition:raise ArithmeticError(message)

def beta_moment(a,b,j):
    return F(factorial(a+j-1)*factorial(a+b-1),
             factorial(a-1)*factorial(a+b+j-1))

def coefficients(k,moments):
    return [(-1)**j*comb(k,j)*moments[j] for j in range(k+1)]

def expected_pairs():
    e=F(1,10000);w=F(701,1344)
    nu=[(1-e)*(w*beta_moment(111,143,j)+(1-w)*beta_moment(132,122,j))+e*beta_moment(111,122,j) for j in range(12)]
    mu=[(1-e)*beta_moment(122,134,j)+e*beta_moment(112,123,j) for j in range(11)]
    first=(coefficients(11,nu),coefficients(10,mu))
    h=[sum(F(55,2)**j/F(factorial(j))*comb(2*j,i)*2**i*(-F(112,123))**(2*j-i)
           for j in range((i+1)//2,6)) for i in range(11)]
    norm=sum(v*beta_moment(56,67,i) for i,v in enumerate(h))
    nu=[sum(v*beta_moment(56,67,i+j) for i,v in enumerate(h))/norm for j in range(12)]
    mu=[beta_moment(62,74,j) for j in range(11)]
    second=(coefficients(11,nu),coefficients(10,mu))
    nu=[F(9,20)*beta_moment(4,4,j)+F(11,20)*beta_moment(5,5,j) for j in range(5)]
    mu=[beta_moment(2,2,j) for j in range(4)]
    third=(coefficients(4,nu),coefficients(3,mu))
    return [first,second,third]

def verify_row(row,expected):
    upper=[F(v) for v in row['upper_coefficients_descending']]
    lower=[F(v) for v in row['lower_coefficients_descending']]
    require((upper,lower)==expected,'independent factorial-moment reconstruction')
    p=row['prime'];require(p>2 and all(p%d for d in range(2,isqrt(p)+1)),'prime modulus')
    require(all(v.denominator%p for v in upper+lower),'unit denominators')
    a=[v.numerator*pow(v.denominator,-1,p)%p for v in upper[::-1]]
    b=[v.numerator*pow(v.denominator,-1,p)%p for v in lower[::-1]]
    require(a[-1]==b[-1]==1,'degree preservation')
    require(a==row['upper_mod_p_ascending'] and b==row['lower_mod_p_ascending'],'reduced coefficients')
    u=row['bezout_upper_multiplier_ascending'];v=row['bezout_lower_multiplier_ascending']
    require(len(u)<=len(b)-1 and len(v)<=len(a)-1,'Bezout degree budget')
    out=[0]*max(len(u)+len(a)-1,len(v)+len(b)-1)
    for left,right in [(u,a),(v,b)]:
        for i,c in enumerate(left):
            for j,d in enumerate(right):out[i+j]=(out[i+j]+c*d)%p
    require(out[0]==1 and not any(out[1:]),'modular Bezout identity')
    aa,bb=upper,lower;pivots=[]
    for iteration in range(2):
        k=len(bb)-1;q0=aa[1]-bb[1]
        remainder=[x-y-q0*z for x,y,z in zip(aa,bb+[F(0)],[F(0)]+bb)]
        require(remainder[:2]==[0,0],'rational leading-term cancellation')
        remainder=remainder[2:];beta=remainder[0]
        pivots.append({'k':k,'pivot':str(beta)})
        if beta==0:break
        aa,bb=bb,[v/beta for v in remainder]
    require(pivots==row['first_pivots'],'independent first-pivot reconstruction')
    require(row['rational_gcd_degree']==0,'claimed gcd degree')
    return {'name':row['name'],'upper_degree':len(a)-1,'lower_degree':len(b)-1,
            'prime':p,'rational_gcd_degree':0}

def determinant(matrix):
    a=[row[:] for row in matrix];result=F(1)
    for i in range(len(a)):
        pivot=next((j for j in range(i,len(a)) if a[j][i]),None)
        if pivot is None:return F(0)
        if pivot!=i:a[i],a[pivot]=a[pivot],a[i];result=-result
        value=a[i][i];result*=value
        for j in range(i+1,len(a)):
            scale=a[j][i]/value
            for k in range(i+1,len(a)):a[j][k]-=scale*a[i][k]
            a[j][i]=F(0)
    return result

def sylvester(upper,lower):
    m,n=len(upper)-1,len(lower)-1
    rows=[]
    for shift in range(n):rows.append([F(0)]*shift+upper+[F(0)]*(n-1-shift))
    for shift in range(m):rows.append([F(0)]*shift+lower+[F(0)]*(m-1-shift))
    return determinant(rows)

def centered(coefficients):
    degree=len(coefficients)-1;out=[F(0)]*(degree+1)
    for i,v in enumerate(coefficients):
        for j in range(degree-i+1):out[j]+=v*comb(degree-i,j)*F(1,2)**(degree-i-j)
    return out

def verify_family(data):
    expected_fields={'current_beta_shapes':[2,2],
        'previous_beta_mixture':'theta*Beta(4,4)+(1-theta)*Beta(5,5)',
        'theta_interval':['0','1'],'represented_degrees':[4,3],
        'center':'1/2','A':'3/20','b_constant':'-3/220','b_linear':'1/33',
        'c_constant':'3/2288','c_linear':'1/1716','resultant_formula':'c*(c-A*b)^2',
        'zero_nominal_pivot_theta':'9/20','zero_pivot_constant_remainder':'9/5720',
        'common_root_theta':'72/85','common_quadratic':'y^2+3/20',
        'upper_other_factor':'y^2+9/748','common_quadratic_gcd_degree':2}
    require(all(data.get(k)==v for k,v in expected_fields.items()),'family coefficient fixture')
    A=F(3,20);points=[F(0),F(1,3),F(2,3),F(1)];controls=[]
    P=coefficients(3,[beta_moment(2,2,j) for j in range(4)])
    require(centered(P)==[0,A,0,1],'centered current polynomial')
    for theta in points+[F(9,20),F(72,85)]:
        Q=coefficients(4,[theta*beta_moment(4,4,j)+(1-theta)*beta_moment(5,5,j) for j in range(5)])
        b=(20*theta-9)/660;c=F(3,2288)+theta/1716
        require(centered(Q)==[c,0,A+b,0,1],'independent centered family')
        actual=sylvester(Q,P);expected=c*(c-A*b)**2
        require(actual==expected,'independent Sylvester determinant')
        controls.append({'theta':str(theta),'resultant':str(actual),'nominal_pivot':str(b)})
        if theta==F(9,20):
            require(b==0 and c==F(9,5720) and actual>0,'harmless degree-drop point')
        if theta==F(72,85):
            B=F(9,748)
            require(b==B and c==A*B and actual==0,'genuine generic collision')
            require(centered(Q)==[A*B,0,A+B,0,1],'common quadratic factor')
            require(A!=B,'gcd exactly quadratic')
    require(data['zero_nominal_pivot_theta']=='9/20' and data['common_root_theta']=='72/85',
            'family special parameters')
    return {'controls':controls,'uniform_identity_certificate':'Both sides have theta degree at most3; four distinct exact Sylvester evaluations certify equality.',
            'collision_gcd_degree':2,'scope':'Generic moment pair; not a factorial-profile counterexample.'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--certificate',type=Path,default=Path(__file__).with_name('moment_pair_certificates.json'))
    ap.add_argument('--output',type=Path,default=Path(__file__).with_name('moment_pair_verification.json'));args=ap.parse_args()
    require(sha256((HERE/PROOF_NAME).read_bytes()).hexdigest()==PROOF_SHA256,'proof authentication')
    require(sha256(args.certificate.read_bytes()).hexdigest()==CERTIFICATE_SHA256,'certificate authentication')
    data=json.loads(args.certificate.read_text());rows=data['generic_counterpairs']+[data['zero_pivot_control']]
    require(len(rows)==3,'complete case set')
    require([r['name'] for r in rows]==['strict_Bernstein_and_transport_counterpair',
            'three_logconcave_densities_counterpair','zero_pivot_with_constant_nonzero_remainder'],
            'case identity')
    done=[verify_row(row,ex) for row,ex in zip(rows,expected_pairs())]
    result={'status':'PASS','proof_sha256':PROOF_SHA256,'certificate_sha256':CERTIFICATE_SHA256,'checker_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'cases':done,'parameter_family':verify_family(data['parameter_family']),
            'arithmetic':'Python standard-library exact integers/Fraction; no sibling imports',
            'scope':'The two explicit obstruction moment pairs and the zero-pivot control have no common root over any characteristic-zero field.'}
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result))
if __name__=='__main__':main()
