def is_subset_id(I, J):
    return all([f in J for f in I.groebner_basis()])

def is_eq_id(I, J):
    return is_subset_id(I, J) and is_subset_id(J, I)

def is_subset_sch(X, Y):
    return is_subset_id(Y.defining_ideal(), X.defining_ideal())

def is_eq_sch(X, Y):
    return is_subset_sch(X, Y) and is_subset_sch(Y, X)

print("TYPE IX, Z=Q_alpha:")

al = 56
be = 23
a = al^2
b = be^2

A2.<x, y> = AffineSpace(2, QQ)
P4.<z0, z1, z2, z3, z4> = ProjectiveSpace(4, QQ)

I1 = ideal([x, y^2 - b])
I2 = ideal([x^2 - b/a, y])
I3 = I1.intersection(I2)
assert(is_eq_id(I3, ideal([x*y, y^2 + a*x^2 - b, y^3 - b*y])))

Z0 = x*y
Z1 = y^2 + a*x^2 - b
Z2 = x*(y^2 - a*x^2 + b)
Z3 = y*(y^2 - a*x^2 - b)
Z4 = (y^2 - a*x^2)^2 - b^2

I = [Z0, Z1, Z2, Z3, Z4]
assert(is_subset_id(ideal(I), I3))

X = P4.subscheme([z2*z3 - z0*z4, 4*a*b*z0^2 - b*z1^2 - a*z2^2 - z3^2 + z1*z4])
print(X)

f = A2.hom(I, P4)
assert(is_eq_sch(X, f.image()))

assert(y - al*x == (al*Z2 + Z3)/(2*al*Z0 + Z1))
assert(y + al*x == (al*Z2 - Z3)/(2*al*Z0 - Z1))

print("\nSINGULARITIES:")
sings = P4.subscheme(X.defining_ideal() + X.Jacobian_matrix().minors(2)).rational_points()
print(sings)
assert(P4(1, -2*al, 2*be, -2*al*be, -4*al*b) in sings)
assert(P4(1, 2*al, 2*be, 2*al*be, 4*al*b) in sings)
assert(P4(1, -2*al, -2*be, 2*al*be, -4*al*b) in sings)
assert(P4(1, 2*al, -2*be, -2*al*be, 4*al*b) in sings)

print("\nLINES:")
L = [P4.subscheme([z2 - 2*be*z0, z3 - be*z1, z4 - 2*be*z3]),
     P4.subscheme([z2 + 2*be*z0, z3 + be*z1, z4 + 2*be*z3]),
     P4.subscheme([z2 + be/al*z1, z3 + 2*al*be*z0, z4 - 2*b*z1]),
     P4.subscheme([z2 - be/al*z1, z3 - 2*al*be*z0, z4 - 2*b*z1])]

for l in L:
    assert(is_subset_sch(l, X))
    print(l)