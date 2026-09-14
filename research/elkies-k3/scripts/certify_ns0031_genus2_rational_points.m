// Exact rational-point proof for v^2 = u^6 - 4*u^4 + 4.
// Run with Magma >= 2.29 (unconditional class/unit group defaults).
// No GRH, analytic-rank assumption, height cutoff, or unsaturated MW group.
SetColumns(0); SetSeed(1);
major,minor,patch:=GetVersion();
assert major gt 2 or (major eq 2 and minor ge 29);
print "MAGMA_VERSION",major,minor,patch;
Q:=Rationals(); P<t>:=PolynomialRing(Q);
E0:=EllipticCurve([Q!0,0,1,-1,0]);
A0,m0,r0,f0:=MordellWeilGroup(E0);
assert r0 and f0 and Invariants(A0) eq [0];
assert m0(A0.1) in {E0![0,0,1],E0![0,-1,1]};
print "BASE_MW_Z_GENERATED_BY_0_0_PROVED";

// theta is a root of 4*theta^3 - 4*theta + 1.  For a square
// x-coordinate s^2 on E0, the Kummer class of s^2-theta is
// either 1 or -theta, because E0(Q)/2E0(Q) = {0,P}.
K<a>:=NumberField(t^3-4*t+2); th:=a/2;
assert 4*th^3-4*th+1 eq 0;
OK:=MaximalOrder(K);
assert Discriminant(OK) eq 148;
Cl:=ClassGroup(OK); assert #Cl eq 1;
UnitGroup(OK);
print "CUBIC_FIELD_AND_CLASS_GROUP_PROVED";
R<s>:=PolynomialRing(K);
for case_id in [0,1] do
  delta:=case_id eq 0 select K!1 else -th;
  C:=HyperellipticCurve((s^4+th*s^2+th^2-1)/delta);
  pt:=case_id eq 0 select C![1,th/(2*th-1),1] else C![0,1/(2*th),1];
  E,mp:=EllipticCurve(C,pt);
  lo,hi:=RankBounds(E); assert lo eq 1 and hi eq 1;
  A,mw,rank_proved,full_group_proved:=MordellWeilGroup(E);
  assert rank_proved and full_group_proved;
  assert Invariants(A) eq [2,0];
  print "CASE",case_id,"RANK_BOUNDS",lo,hi,"FULL_MW_GROUP_PROVED",full_group_proved;
  print "ELLIPTIC_AINVS",aInvariants(E);
  print "MW_GENERATORS",[mw(A.i):i in [1..Ngens(A)]];
  pols:=DefiningPolynomials(Inverse(mp));
  common:=GCD(pols[1],pols[3]);
  cov:=map<E -> ProjectiveSpace(Q,1) | [pols[1] div common,pols[3] div common]>;
  mwc:=map<A -> E | z :-> mw(z)>;
  N,V,index_condition,cosets:=Chabauty(mwc,cov,5 : Aux:={3,7,11,13,17,19});
  assert N eq 6 and #V eq 6 and IsEmpty(cosets[2]);
  // The index condition is satisfied because the full MW group is proved.
  images:=[cov(mw(z)):z in V];
  P1:=Codomain(cov);
  expected:=case_id eq 0 select {P1![1,1],P1![-1,1],P1![1,0]}
                         else {P1![0,1],P1![1,2],P1![-1,2]};
  assert Seqset(images) eq expected;
  assert forall{z:z in expected | #[w:w in images | w eq z] eq 2};
  print "CASE",case_id,"CHABAUTY_EXACT_BOUND",N,"INDEX_CONDITION",index_condition,"RESIDUAL_COSETS_EMPTY",IsEmpty(cosets[2]);
  print "RATIONAL_IMAGES",images;
  print "CASE",case_id,"PASS";
end for;
// s=1/u, Y=v/u^3 identifies the original genus-two curve with
// Y^2=4*s^6-4*s^2+1; each rational s in {0,+/-1/2,+/-1,infinity}
// gives exactly two points, including the two smooth points at infinity.
print "PASS_NS0031_GENUS2_RATIONAL_POINTS_EXACT_12";
