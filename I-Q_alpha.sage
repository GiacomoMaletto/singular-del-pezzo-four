def is_subset_id(I, J):
    return all([f in J for f in I.groebner_basis()])

def is_eq_id(I, J):
    return is_subset_id(I, J) and is_subset_id(J, I)

def is_subset_sch(X, Y):
    return is_subset_id(Y.defining_ideal(), X.defining_ideal())

def is_eq_sch(X, Y):
    return is_subset_sch(X, Y) and is_subset_sch(Y, X)

print("TYPE I, Z=Q_alpha:")

al = 2
be = 3
a = al^2
b = be^2

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

I1 = ideal([x - (b/A - A)/(2*al), y - (b/A + A)/2])
I2 = ideal([x - (b/B - B)/(2*al), y - (b/B + B)/2])
I3 = ideal([x - (b/C - C)/(2*al), y - (b/C + C)/2])
I4 = ideal([x - (b/D - D)/(2*al), y - (b/D + D)/2])
I5 = I1.intersection(I2, I3, I4)
assert(is_eq_id(I5, ideal([y^2 - a*x^2 - b, (y - al*x)^4 + s3*(y - al*x)^3 + s2*(y - al*x)^2 + s1*(y - al*x) + s0])))

Z0 = y^2 - a*x^2 - b
Z1 = x*(y^2 - a*x^2 - b)
Z2 = y*(y^2 - a*x^2 - b)
Z3 = (y^2 - a*x^2)^2 - b^2
Z4 = (y - al*x)^2*b^2 + s3*(y - al*x)*b^2 + s2*b^2 + s1*(y + al*x)*b + s0*(y + al*x)^2

I = [Z0, Z1, Z2, Z3, Z4]
assert(is_subset_id(ideal(I), I5))

X = P4.subscheme([z0*z3 - z2^2 + a*z1^2 - b*z0^2,
                  b^2*(z2 - al*z1)^2 + s3*b^2*(z2 - al*z1)*z0 + s2*b^2*z0^2 + s1*b*(z2 + al*z1)*z0 + s0*(z2 + al*z1)^2 - (z3 - 2*b*z0)*z4])
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
L = [P4.subscheme([2*al*z1 - (b/A - A)*z0, 2*z2 - (b/A + A)*z0, z3 - 2*b*z0]),
     P4.subscheme([2*al*z1 - (b/B - B)*z0, 2*z2 - (b/B + B)*z0, z3 - 2*b*z0]),
     P4.subscheme([2*al*z1 - (b/C - C)*z0, 2*z2 - (b/C + C)*z0, z3 - 2*b*z0]),
     P4.subscheme([2*al*z1 - (b/D - D)*z0, 2*z2 - (b/D + D)*z0, z3 - 2*b*z0]),
     P4.subscheme([z2 - al*z1 - A*z0, z3 - 2*A*al*z1 - (A^2 + b)*z0, A^2*z4 - s0*z3 - A*b*s1*z0]),
     P4.subscheme([z2 - al*z1 - B*z0, z3 - 2*B*al*z1 - (B^2 + b)*z0, B^2*z4 - s0*z3 - B*b*s1*z0]),
     P4.subscheme([z2 - al*z1 - C*z0, z3 - 2*C*al*z1 - (C^2 + b)*z0, C^2*z4 - s0*z3 - C*b*s1*z0]),
     P4.subscheme([z2 - al*z1 - D*z0, z3 - 2*D*al*z1 - (D^2 + b)*z0, D^2*z4 - s0*z3 - D*b*s1*z0]),
     P4.subscheme([z2 + al*z1 - (b/A)*z0, z3 + 2*(b/A)*al*z1 - ((b/A)^2 + b)*z0, z4 - A^2*z3 - s3*b*A*z0]),
     P4.subscheme([z2 + al*z1 - (b/B)*z0, z3 + 2*(b/B)*al*z1 - ((b/B)^2 + b)*z0, z4 - B^2*z3 - s3*b*B*z0]),
     P4.subscheme([z2 + al*z1 - (b/C)*z0, z3 + 2*(b/C)*al*z1 - ((b/C)^2 + b)*z0, z4 - C^2*z3 - s3*b*C*z0]),
     P4.subscheme([z2 + al*z1 - (b/D)*z0, z3 + 2*(b/D)*al*z1 - ((b/D)^2 + b)*z0, z4 - D^2*z3 - s3*b*D*z0])]

for l in L:
    assert(is_subset_sch(l, X))
    print(l)