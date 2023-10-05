def is_subset_id(I, J):
    return all([f in J for f in I.groebner_basis()])

def is_eq_id(I, J):
    return is_subset_id(I, J) and is_subset_id(J, I)

def is_subset_sch(X, Y):
    return is_subset_id(Y.defining_ideal(), X.defining_ideal())

def is_eq_sch(X, Y):
    return is_subset_sch(X, Y) and is_subset_sch(Y, X)

print("TYPE VIII, C1, C2 DEFINED OVER k:")

A = 1
B = 2
C = 3
D = 4

s0 = A*B
s1 = -(A + B)
r0 = C*D
r1 = -(C + D)

A2.<x, y> = AffineSpace(2, QQ)
P4.<z0, z1, z2, z3, z4> = ProjectiveSpace(4, QQ)

I1 = ideal([x - 1, y - A])
I2 = ideal([x - 1, y - B])
I3 = ideal([x + 1, y - C])
I4 = ideal([x + 1, y - D])
I5 = I1.intersection(I2, I3, I4)
assert(is_eq_id(I5, ideal([x^2 - 1,
                           (x - 1)*(y^2 + r1*y + r0),
                           (x + 1)*(y^2 + s1*y + s0)])))

Z0 = y*(x^2 - 1)
Z1 = y^2*(x^2 - 1)
Z2 = x*y^2*(x^2 - 1)
Z3 = x^2*y^2*(x^2 - 1)
Z4 = r0*(x + 1)*(y^2 + s1*y + s0) - s0*(x - 1)*(y^2 + r1*y + r0)

I = [Z0, Z1, Z2, Z3, Z4]
assert(is_subset_id(ideal(I), I5))

X = P4.subscheme([z2^2 - z1*z3,
                  (r0 - s0)*z1*z2 + (r0*s1 - s0*r1)*z0*z2 + (r0 + s0)*z1^2 + (r0*s1 + s0*r1)*z0*z1 + 2*r0*s0*z0^2 - z4*(z3 - z1)])
print(X)

f = A2.hom(I, P4)
assert(is_eq_sch(X, f.image()))

assert(x == Z2/Z1)
assert(y == Z1/Z0)

print("\nSINGULARITIES:")
sings = P4.subscheme(X.defining_ideal() + X.Jacobian_matrix().minors(2)).rational_points()
print(sings)
assert(P4(0, 0, 0, 0, 1) in sings)

print("\nLINES:")
L = [P4.subscheme([z1 - A*z0, z2 - z1, z3 - z2]),
     P4.subscheme([z1 - B*z0, z2 - z1, z3 - z2]),
     P4.subscheme([z1 - C*z0, z2 + z1, z3 + z2]),
     P4.subscheme([z1 - D*z0, z2 + z1, z3 + z2])]

for l in L:
    assert(is_subset_sch(l, X))
    print(l)