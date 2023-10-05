def is_subset_id(I, J):
    return all([f in J for f in I.groebner_basis()])

def is_eq_id(I, J):
    return is_subset_id(I, J) and is_subset_id(J, I)

def is_subset_sch(X, Y):
    return is_subset_id(Y.defining_ideal(), X.defining_ideal())

def is_eq_sch(X, Y):
    return is_subset_sch(X, Y) and is_subset_sch(Y, X)

print("TYPE V, k(alpha)=k(beta):")

a = 3
K.<al> = NumberField([x^2 - a])

A = 7

A2.<x, y> = AffineSpace(K, 2)
P4.<z0, z1, z2, z3, z4> = ProjectiveSpace(K, 4)

I1 = ideal([al*x - a - A/(2*a)*y, y^2])
I2 = ideal([al*x + a + A/(2*a)*y, y^2])
I3 = I1.intersection(I2)
assert(is_eq_id(I3, ideal(y^2 - a*x^2 + a^2 + A*y, y^2)))

Z0 = y^2
Z1 = y^2 - a*x^2 + a^2 + A*y
Z2 = x*(y^2 - a*x^2 + a^2 + A*y)
Z3 = y*(y^2 - a*x^2 + a^2 + A*y)
Z4 = (y^2 - a*x^2 + a^2 + A*y)^2

I = [Z0, Z1, Z2, Z3, Z4]
assert(is_subset_id(ideal(I), I3))

X = P4.subscheme([z3^2 - z0*z4, a^2*z1^2 - a*z2^2 + A*z1*z3 + z0*z4 - z1*z4])
print(X)

f = A2.hom(I, P4)
assert(is_eq_sch(X, f.image()))

assert(x == Z2/Z1)
assert(y == Z3/Z1)

print("\nSINGULARITIES:")
sings = P4.subscheme(X.defining_ideal() + X.Jacobian_matrix().minors(2)).rational_points()
print(sings)
assert(P4(0, 1, 0, al, 0))
assert(P4(0, 1, 0, -al, 0))
assert(P4(1, 0, 0, 0, 0))