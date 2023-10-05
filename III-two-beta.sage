def is_subset_id(I, J):
    return all([f in J for f in I.groebner_basis()])

def is_eq_id(I, J):
    return is_subset_id(I, J) and is_subset_id(J, I)

def is_subset_sch(X, Y):
    return is_subset_id(Y.defining_ideal(), X.defining_ideal())

def is_eq_sch(X, Y):
    return is_subset_sch(X, Y) and is_subset_sch(Y, X)

print("TYPE III, X ISKOVSKIH, TWO (-1)-CURVES CAN BE IMMEDIATELY BLOWN DOWN, k(alpha)!=k(beta):")

a = 3
be = 5
b = be^2
K.<al> = NumberField([x^2 - a])

A = 7
B = 11

A2.<x, y> = AffineSpace(K, 2)
P4.<z0, z1, z2, z3, z4> = ProjectiveSpace(K, 4)

I1 = ideal([y + B/2*x + be*(A/(2*b)*x - 1), x^2])
I2 = ideal([y + B/2*x - be*(A/(2*b)*x - 1), x^2])
I3 = I1.intersection(I2)
assert(is_eq_id(I3, ideal(y^2 + B*x*y + A*x - b, x^2)))

Z0 = x^2
Z1 = y^2 - a*x^2 - b + A*x + B*x*y
Z2 = x*(y^2 - a*x^2 - b)
Z3 = y*(y^2 - a*x^2 - b) + A*x*y + B*b*x
Z4 = (y^2 - a*x^2 - b)^2

I = [Z0, Z1, Z2, Z3, Z4]
assert(is_subset_id(ideal(I), I3))

X = P4.subscheme([z2^2 - z0*z4,
                  a*(b*B^2 - A^2)*z0^2 - b*z1^2 - 2*a*A*z0*z2 - A*z1*z2 + B*z2*z3 + z3^2 - a*z0*z4 - z1*z4])
print(X)

f = A2.hom(I, P4)
assert(is_eq_sch(X, f.image()))

assert(y - al*x == (Z3 - al*A*Z0 + (B - al)*Z2)/(Z1 + al*B*Z0) == y - al*x)
assert(y + al*x == (Z3 + al*A*Z0 + (B + al)*Z2)/(Z1 - al*B*Z0))

print("\nSINGULARITIES:")
sings = P4.subscheme(X.defining_ideal() + X.Jacobian_matrix().minors(2)).rational_points()
print(sings)
assert(P4(0, 1, 0, be, 0))
assert(P4(0, 1, 0, -be, 0))