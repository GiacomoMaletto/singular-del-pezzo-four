def is_subset_id(I, J):
    return all([f in J for f in I.groebner_basis()])

def is_eq_id(I, J):
    return is_subset_id(I, J) and is_subset_id(J, I)

def is_subset_sch(X, Y):
    return is_subset_id(Y.defining_ideal(), X.defining_ideal())

def is_eq_sch(X, Y):
    return is_subset_sch(X, Y) and is_subset_sch(Y, X)

print("TYPE IX, Z=P1xP1:")

al = 3
be = 5
a = al^2
b = be^2

A2.<x, y> = AffineSpace(2, QQ)
P4.<z0, z1, z2, z3, z4> = ProjectiveSpace(4, QQ)

I1 = ideal([x^2 - a, y^2 - b])

Z0 = x^2 - a
Z1 = y*(x^2 - a)
Z2 = y^2*(x^2 - a)
Z3 = y^2 - b
Z4 = x*(y^2 - b)

I = [Z0, Z1, Z2, Z3, Z4]
assert(is_subset_id(ideal(I), I1))

X = P4.subscheme([z1^2 - z0*z2, b*z0*z3 - z2*z3 - a*z3^2 + z4^2])
print(X)

f = A2.hom(I, P4)
assert(is_eq_sch(X, f.image()))

assert(x == Z4/Z3)
assert(y == Z1/Z0)

print("\nSINGULARITIES:")
sings = P4.subscheme(X.defining_ideal() + X.Jacobian_matrix().minors(2)).rational_points()
print(sings)
assert(P4(0, 0, 0, 1, al) in sings)
assert(P4(0, 0, 0, 1, -al) in sings)
assert(P4(1, -be, be^2, 0, 0) in sings)
assert(P4(1, be, be^2, 0, 0) in sings)

print("\nLINES:")
L = [P4.subscheme([z1 - be*z0, z2 - b*z0, z4 - al*z3]),
     P4.subscheme([z1 + be*z0, z2 - b*z0, z4 - al*z3]),
     P4.subscheme([z1 - be*z0, z2 - b*z0, z4 + al*z3]),
     P4.subscheme([z1 + be*z0, z2 - b*z0, z4 + al*z3])]

for l in L:
    assert(is_subset_sch(l, X))
    print(l)
