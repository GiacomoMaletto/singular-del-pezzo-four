def is_subset_id(I, J):
    return all([f in J for f in I.groebner_basis()])

def is_eq_id(I, J):
    return is_subset_id(I, J) and is_subset_id(J, I)

def is_subset_sch(X, Y):
    return is_subset_id(Y.defining_ideal(), X.defining_ideal())

def is_eq_sch(X, Y):
    return is_subset_sch(X, Y) and is_subset_sch(Y, X)

print("TYPE III, X NOT ISKOVSKIH:")

a = 6
b = 3
c = 2

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

I1 = ideal([1 + a*y + b*x*y + c*x^2*y, x - A])
I2 = ideal([1 + a*y + b*x*y + c*x^2*y, x - B])
I3 = ideal([1 + a*y + b*x*y + c*x^2*y, x - C])
I4 = ideal([1 + a*y + b*x*y + c*x^2*y, x - D])
I5 = I1.intersection(I2, I3, I4)
assert(is_eq_id(I5, ideal([1 + a*y + b*x*y + c*x^2*y, x^4 + s3*x^3 + s2*x^2 + s1*x + s0])))

Z0 = 1 + a*y + b*x*y + c*x^2*y
Z1 = y*(1 + a*y + b*x*y + c*x^2*y)
Z2 = x*y*(1 + a*y + b*x*y + c*x^2*y)
Z3 = x^2*y*(1 + a*y + b*x*y + c*x^2*y)
Z4 = y^2*(x^4 + s3*x^3 + s2*x^2 + s1*x + s0)

I = [Z0, Z1, Z2, Z3, Z4]
assert(is_subset_id(ideal(I), I5))

X = P4.subscheme([z2^2 - z1*z3, z3^2 + s3*z2*z3 + s2*z2^2 + s1*z1*z2 + s0*z1^2 - z4*(z0 + a*z1 + b*z2 + c*z3)])
print(X)

f = A2.hom(I, P4)
assert(is_eq_sch(X, f.image()))

assert(x == Z2/Z1)
assert(y == Z1/Z0)

print("\nSINGULARITIES:")
sings = P4.subscheme(X.defining_ideal() + X.Jacobian_matrix().minors(2)).rational_points()
print(sings)
assert(P4(1, 0, 0, 0, 0) in sings)
assert(P4(0, 0, 0, 0, 1) in sings)

print("\nLINES:")
L = [P4.subscheme([z2 - A*z1, z3 - A^2*z1, z4]),
     P4.subscheme([z2 - B*z1, z3 - B^2*z1, z4]),
     P4.subscheme([z2 - C*z1, z3 - C^2*z1, z4]),
     P4.subscheme([z2 - D*z1, z3 - D^2*z1, z4]),
     P4.subscheme([z2 - A*z1, z3 - A^2*z1, z0 + a*z1 + b*A*z1 + c*A^2*z1]),
     P4.subscheme([z2 - B*z1, z3 - B^2*z1, z0 + a*z1 + b*B*z1 + c*B^2*z1]),
     P4.subscheme([z2 - C*z1, z3 - C^2*z1, z0 + a*z1 + b*C*z1 + c*C^2*z1]),
     P4.subscheme([z2 - D*z1, z3 - D^2*z1, z0 + a*z1 + b*D*z1 + c*D^2*z1])]

for l in L:
    assert(is_subset_sch(l, X))
    print(l)