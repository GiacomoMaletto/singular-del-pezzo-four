def is_subset_id(I, J):
    return all([f in J for f in I.groebner_basis()])

def is_eq_id(I, J):
    return is_subset_id(I, J) and is_subset_id(J, I)

def is_subset_sch(X, Y):
    return is_subset_id(Y.defining_ideal(), X.defining_ideal())

def is_eq_sch(X, Y):
    return is_subset_sch(X, Y) and is_subset_sch(Y, X)

print("TYPE XII:")

A2.<x, y> = AffineSpace(2, QQ)
P4.<z0, z1, z2, z3, z4> = ProjectiveSpace(4, QQ)

I1 = ideal([y - x^2, x^3])
I2 = ideal([x - 1, y^2])
I3 = I1.intersection(I2)
assert(is_eq_id(I3, ideal(y^2, x^2*y - x*y, x^3 - x^2 - x*y + y)))

Z0 = y^2
Z1 = x*y^2
Z2 = y^3
Z3 = x^2*y - x*y
Z4 = x^3 - x^2 - x*y + y

I = [Z0, Z1, Z2, Z3, Z4]
assert(is_subset_id(ideal(I), I3))

X = P4.subscheme([z0*z1 - z1^2 + z2*z3, z0^2 - z0*z1 + z1*z3 - z2*z4])
print(X)

f = A2.hom(I, P4)
assert(is_eq_sch(X, f.image()))

assert(x == Z1/Z0)
assert(y == Z2/Z0)

print("\nSINGULARITIES:")
sings = P4.subscheme(X.defining_ideal() + X.Jacobian_matrix().minors(2)).rational_points()
print(sings)
assert(P4(0, 0, 0, 0, 1) in sings)

print("\nLINES:")
L = [P4.subscheme([z0, z1, z2]),
     P4.subscheme([z0 - z1, z3, z4]),
     P4.subscheme([z0 - z2, z1, z2])]

for l in L:
    assert(is_subset_sch(l, X))
    print(l)