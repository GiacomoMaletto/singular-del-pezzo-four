def is_subset_id(I, J):
    return all([f in J for f in I.groebner_basis()])

def is_eq_id(I, J):
    return is_subset_id(I, J) and is_subset_id(J, I)

def is_subset_sch(X, Y):
    return is_subset_id(Y.defining_ideal(), X.defining_ideal())

def is_eq_sch(X, Y):
    return is_subset_sch(X, Y) and is_subset_sch(Y, X)

print("TYPE III, X ISKOVSKIH, FOUR (-1)-CURVES CAN BE IMMEDIATELY BLOWN DOWN:")

a = 3
K.<al> = NumberField([x^2 - a])
s = K.galois_group().gen()

A = 23 + 12*al
B = 35 + 55*al

s0 = A*B + s(A*B)
r0 = al*(A*B - s(A*B))
s1 = -(A + B + s(A + B))
r1 = -al*(A + B - s(A + B))

A2.<x, y> = AffineSpace(K, 2)
P4.<z0, z1, z2, z3, z4> = ProjectiveSpace(K, 4)

I1 = ideal([x - A, y - al])
I2 = ideal([x - B, y - al])
I3 = ideal([x - s(A), y + al])
I4 = ideal([x - s(B), y + al])
I5 = I1.intersection(I2, I3, I4)
assert(is_eq_id(ideal([y^2 - a,
                       2*x^2*y + s1*x*y + s0*y + r1*x + r0,
                       r1*x*y + r0*y + 2*a*x^2 + a*s1*x + a*s0,
                       (x - A)*(x - B)*(x - s(A))*(x - s(B))]), I5))

Z0 = y^2 - a
Z1 = x*(y^2 - a)
Z2 = x^2*(y^2 - a)
Z3 = 2*x^2*y + s1*x*y + s0*y + r1*x + r0
Z4 = r1*x*y + r0*y + 2*a*x^2 + a*s1*x + a*s0

I = [Z0, Z1, Z2, Z3, Z4]
assert(is_subset_id(ideal(I), I5))

X = P4.subscheme([z1^2 - z0*z2,
                  r0*z0*z3 + r1*z1*z3 + a*z3^2 - s0*z0*z4 - s1*z1*z4 - 2*z2*z4 - z4^2])
print(X)

f = A2.hom(I, P4)
assert(is_eq_sch(X, f.image()))

assert(x == Z1/Z0)
assert(y == (Z4 + 2*Z2 + s1*Z1 + s0*Z0)/Z3)

print("\nSINGULARITIES:")
sings = P4.subscheme(X.defining_ideal() + X.Jacobian_matrix().minors(2)).rational_points()
print(sings)
assert(P4(0, 0, 0, 1, al) in sings)
assert(P4(0, 0, 0, 1, -al) in sings)