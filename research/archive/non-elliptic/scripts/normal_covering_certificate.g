# Exact GAP front end for faithful normal-cover certificates.
#
# Load with:
#   Read("scripts/normal_covering_certificate.g");
# and call:
#   NormalCoveringCertificate(G, [basepoint_1, ..., basepoint_r]);
#
# The group G must already act faithfully on the disjoint union of the
# component root orbits.  The dependency-free Python replay in
# jcsearch/normal_covering.py verifies the checked-in small certificates.

NormalConjugateUnion := function(G, H)
    local result, conjugate;
    result := [];
    for conjugate in AsList(ConjugacyClassSubgroups(G, H)) do
        UniteSet(result, AsList(conjugate));
    od;
    return Set(result);
end;

ExactNormalCoveringNumber := function(G)
    local classes, covers, full, size, selected, union, cover;
    if IsCyclic(G) then
        return fail;
    fi;
    classes := ConjugacyClassesMaximalSubgroups(G);
    covers := List(
        classes,
        class -> NormalConjugateUnion(G, Representative(class))
    );
    full := Set(AsList(G));
    for size in [1..Length(covers)] do
        for selected in Combinations(covers, size) do
            union := [];
            for cover in selected do
                UniteSet(union, cover);
            od;
            if Set(union) = full then
                return size;
            fi;
        od;
    od;
    Error("finite noncyclic group has no computed normal covering");
end;

NormalCoveringCertificate := function(G, componentPoints)
    local stabilizers, orbits, indices, cores, covered, H, commonCore;
    stabilizers := List(componentPoints, point -> Stabilizer(G, point));
    orbits := List(componentPoints, point -> Set(Orbit(G, point)));
    indices := List(stabilizers, H -> Index(G, H));
    cores := List(stabilizers, H -> Core(G, H));
    covered := [];
    for H in stabilizers do
        UniteSet(covered, NormalConjugateUnion(G, H));
    od;
    commonCore := Intersection(cores);
    return rec(
        groupOrder := Size(G),
        componentCount := Length(componentPoints),
        factorizationShape := SortedList(indices),
        stabilizerOrders := List(stabilizers, Size),
        normalCover := Set(covered) = Set(AsList(G)),
        commonCoreOrder := Size(commonCore),
        faithful := Size(commonCore) = 1,
        normalCoveringNumber := ExactNormalCoveringNumber(G),
        indexSum := Sum(indices),
        componentOrbits := orbits
    );
end;
