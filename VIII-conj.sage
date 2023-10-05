def is_subset_id(I, J):
    return all([f in J for f in I.groebner_basis()])

def is_eq_id(I, J):
    return is_subset_id(I, J) and is_subset_id(J, I)

def is_subset_sch(X, Y):
    return is_subset_id(Y.defining_ideal(), X.defining_ideal())

def is_eq_sch(X, Y):
    return is_subset_sch(X, Y) and is_subset_sch(Y, X)

print("TYPE VIII, C1, C2 CONJUGATE OVER k(alpha)/k:")

a = 3
K.<al> = NumberField([x^2 - a])
s = K.galois_group().gen()

A = 7 + 9*al
B = 11 + 13*al

s0 = A*B + s(A*B)
r0 = al*(A*B - s(A*B))
s1 = -(A + B + s(A + B))
r1 = -al*(A + B - s(A + B))

A2.<x, y> = AffineSpace(K, 2)
P4.<z0, z1, z2, z3, z4> = ProjectiveSpace(K, 4)

I1 = ideal([x - al, y - A])
I2 = ideal([x - al, y - B])
I3 = ideal([x + al, y - s(A)])
I4 = ideal([x + al, y - s(B)])
I5 = I1.intersection(I2, I3, I4)
assert(is_eq_id(ideal([x^2 - a,
                       2*x*y^2 + s1*x*y + s0*x + r1*y + r0,
                       r1*x*y + r0*x + 2*a*y^2 + a*s1*y + a*s0,
                       (y - A)*(y - B)*(y - s(A))*(y - s(B))]), I5))

Z0 = y*(x^2 - a)
Z1 = y^2*(x^2 - a)
Z2 = x*y^2*(x^2 - a)
Z3 = x^2*y^2*(x^2 - a)
Z4 = 2*r0*x*y^2 - 2*a*s0*y^2 + (s1*r0 - s0*r1)*x*y + (r0*r1 - a*s0*s1)*y + (r0^2 - a*s0^2)

I = [Z0, Z1, Z2, Z3, Z4]
assert(is_subset_id(ideal(I), I5))

X = P4.subscheme([z2^2 - z1*z3,
                  (r0^2 - a*s0^2)*z0^2 + (r0*r1 - a*s0*s1)*z0*z1 - 2*a*s0*z1^2 + (s1*r0 - s0*r1)*z0*z2 + 2*r0*z1*z2 + a*z1*z4 - z3*z4])
print(X)

f = A2.hom(I, P4)
assert(is_eq_sch(X, f.image()))

assert(x == Z2/Z1)
assert(y == Z1/Z0)

print("\nSINGULARITIES:")
sings = P4.subscheme(X.defining_ideal() + X.Jacobian_matrix().minors(2)).rational_points()
print(sings)
assert(P4(0, 0, 0, 0, 1) in sings)