def is_subset_id(I, J):
    return all([f in J for f in I.groebner_basis()])

def is_eq_id(I, J):
    return is_subset_id(I, J) and is_subset_id(J, I)

def is_subset_sch(X, Y):
    return is_subset_id(Y.defining_ideal(), X.defining_ideal())

def is_eq_sch(X, Y):
    return is_subset_sch(X, Y) and is_subset_sch(Y, X)

print("TYPE I, Z=P1xP1:")

A = 1
B = 2
C = 3
D = 4

s0 = A*B*C*D
s1 = -(A*B*C + A*B*D + A*C*D + B*C*D)
s2 = (A*B + A*C + A*D + B*C + B*D + C*D)
s3 = -(A + B + C + D)

A2.<x, y> = AffineSpace(QQ, 2)
P4.<z0, z1, z2, z3, z4> = ProjectiveSpace(QQ, 4)

I1 = ideal([x - A, y - A])
I2 = ideal([x - B, y - B])
I3 = ideal([x - C, y - C])
I4 = ideal([x - D, y - D])
I5 = ideal(I1).intersection(I2, I3, I4)
assert(is_eq_id(I5, ideal([x - y, x^2*y^2 + s3*x*y^2 + s2*x*y + s1*x + s0])))

Z0 = x - y
Z1 = x*(x - y)
Z2 = y*(x - y)
Z3 = x*y*(x - y)
Z4 = x^2*y^2 + s3*x*y^2 + s2*x*y + s1*x + s0

I = [Z0, Z1, Z2, Z3, Z4]
assert(is_subset_id(ideal(I), I5))

X = P4.subscheme([z0*z3 - z1*z2, z3^2 + s3*z2*z3 + s2*z1*z2 + s1*z0*z1 + s0*z0^2 - z4*(z1 - z2)])
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
L = [P4.subscheme([z1 - A*z0, z2 - A*z0, z3 - A^2*z0]),
     P4.subscheme([z1 - B*z0, z2 - B*z0, z3 - B^2*z0]),
     P4.subscheme([z1 - C*z0, z2 - C*z0, z3 - C^2*z0]),
     P4.subscheme([z1 - D*z0, z2 - D*z0, z3 - D^2*z0]),
     P4.subscheme([z1 - A*z0, z3 - A*z2, z4 + (A + s3)*z3 + (s2+A*s3+A^2)*z1]),
     P4.subscheme([z1 - B*z0, z3 - B*z2, z4 + (B + s3)*z3 + (s2+B*s3+B^2)*z1]),
     P4.subscheme([z1 - C*z0, z3 - C*z2, z4 + (C + s3)*z3 + (s2+C*s3+C^2)*z1]),
     P4.subscheme([z1 - D*z0, z3 - D*z2, z4 + (D + s3)*z3 + (s2+D*s3+D^2)*z1]),
     P4.subscheme([z2 - A*z0, z3 - A*z1, z4 - A*z3 + s0/A*z0]),
     P4.subscheme([z2 - B*z0, z3 - B*z1, z4 - B*z3 + s0/B*z0]),
     P4.subscheme([z2 - C*z0, z3 - C*z1, z4 - C*z3 + s0/C*z0]),
     P4.subscheme([z2 - D*z0, z3 - D*z1, z4 - D*z3 + s0/D*z0])]

for l in L:
    assert(is_subset_sch(l, X))
    print(l)