def is_subset_id(I, J):
    return all([f in J for f in I.groebner_basis()])

def is_eq_id(I, J):
    return is_subset_id(I, J) and is_subset_id(J, I)

def is_subset_sch(X, Y):
    return is_subset_id(Y.defining_ideal(), X.defining_ideal())

def is_eq_sch(X, Y):
    return is_subset_sch(X, Y) and is_subset_sch(Y, X)

print("TYPE IV, Z=P1xP1:")

al = 4
be = 5
a = al^2
b = be^2

K.<x, y> = PolynomialRing(QQ, 2)
I1 = ideal([x^2 - a, y])
I2 = ideal([x, y^2 - b])
I3 = I1.intersection(I2)
assert(is_eq_id(I3, ideal([x*y, y^3 - b*y, b*x^2 + a*y^2 - a*b])))

A2.<x, y> = AffineSpace(2, QQ)
P4.<z0, z1, z2, z3, z4> = ProjectiveSpace(4, QQ)

Z0 = x*y
Z1 = x^2*y
Z2 = x*y^2
Z3 = x^2*y^2
Z4 = b*x^2 + a*y^2 - a*b

I = [Z0, Z1, Z2, Z3, Z4]
assert(is_subset_id(ideal(I), I3))

X = P4.subscheme([z0*z3 - z1*z2, z3*z4 - b*z1^2 - a*z2^2 + a*b*z0^2])
print(X)

f = A2.hom(I, P4)
assert(is_eq_sch(X, f.image()))

assert(x == Z1/Z0)
assert(y == Z2/Z0)

print("\nSINGULARITIES:")
sings = P4.subscheme(X.defining_ideal() + X.Jacobian_matrix().minors(2)).rational_points()
print(sings)
assert(P4(0, 0, 0, 0, 1) in X)

print("\nLINES:")
L = [P4.subscheme([z1 - al*z0, z2, z3]),
     P4.subscheme([z1 + al*z0, z2, z3]),
     P4.subscheme([z1, z2 - be*z0, z3]),
     P4.subscheme([z1, z2 + be*z0, z3]),
     P4.subscheme([z1 - al*z0, z3 - al*z2, z4 - al*z2]),
     P4.subscheme([z1 + al*z0, z3 + al*z2, z4 + al*z2]),
     P4.subscheme([z2 - be*z0, z3 - be*z1, z4 - be*z1]),
     P4.subscheme([z2 + be*z0, z3 + be*z1, z4 + be*z1])]

for l in L:
    assert(is_subset_sch(l, X))
    print(l)