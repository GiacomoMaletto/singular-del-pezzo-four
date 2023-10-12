var('m00, m10, m01, m20, m11, m02, m30, m21, m12, m03, Ax, Ay, Bx, By, Cx, Cy, Dx, Dy, Ex, Ey, Fx, Fy')

# coordinates of the points P1,...,P6.
Ax = 0
Ay = 0

Bx = 1
By = 1

Cx = -1
Cy = 1

Dx = 2
Dy = 4

Ex = -2
Ey = 4

Fx = 0
Fy = 3

eqns = [m00 + m10*Ax + m01*Ay + m20*Ax^2 + m11*Ax*Ay + m02*Ay^2 + m30*Ax^3 + m21*Ax^2*Ay + m12*Ax*Ay^2 + m03*Ay^3,
        m00 + m10*Bx + m01*By + m20*Bx^2 + m11*Bx*By + m02*By^2 + m30*Bx^3 + m21*Bx^2*By + m12*Bx*By^2 + m03*By^3,
        m00 + m10*Cx + m01*Cy + m20*Cx^2 + m11*Cx*Cy + m02*Cy^2 + m30*Cx^3 + m21*Cx^2*Cy + m12*Cx*Cy^2 + m03*Cy^3,
        m00 + m10*Dx + m01*Dy + m20*Dx^2 + m11*Dx*Dy + m02*Dy^2 + m30*Dx^3 + m21*Dx^2*Dy + m12*Dx*Dy^2 + m03*Dy^3,
        m00 + m10*Ex + m01*Ey + m20*Ex^2 + m11*Ex*Ey + m02*Ey^2 + m30*Ex^3 + m21*Ex^2*Ey + m12*Ex*Ey^2 + m03*Ey^3,
        m00 + m10*Fx + m01*Fy + m20*Fx^2 + m11*Fx*Fy + m02*Fy^2 + m30*Fx^3 + m21*Fx^2*Fy + m12*Fx*Fy^2 + m03*Fy^3]

sol = solve(eqns, m00, m10, m01, m20, m11, m02, m30, m21, m12, m03)[0]

for s in sol: print(s)

# Choose the rows (starting from 0) which have just r1, r2, r3, r4 and nothing else.
r1 = sol[9].rhs().variables()[0]
r2 = sol[8].rhs().variables()[0]
r3 = sol[7].rhs().variables()[0]
r4 = sol[6].rhs().variables()[0]

A2.<x, y> = AffineSpace(2, QQ)
P3.<z0, z1, z2, z3> = ProjectiveSpace(3, QQ)

AiAj = [1, x, y, x^2, x*y, y^2, x^3, x^2*y, x*y^2, y^3]

I = [sum([sol[i].rhs().subs(r1 == 1, r2 == 0, r3 == 0, r4 == 0)*AiAj[i] for i in range(0, 10)]),
     sum([sol[i].rhs().subs(r1 == 0, r2 == 1, r3 == 0, r4 == 0)*AiAj[i] for i in range(0, 10)]),
     sum([sol[i].rhs().subs(r1 == 0, r2 == 0, r3 == 1, r4 == 0)*AiAj[i] for i in range(0, 10)]),
     sum([sol[i].rhs().subs(r1 == 0, r2 == 0, r3 == 0, r4 == 1)*AiAj[i] for i in range(0, 10)])]

print("BASE OF THE VECTOR SPACE OF CUBIC POLYNOMIALS VANISHING IN P1,...,P6:")

for i in I: print(i)

print("CUBIC SURFACE:")

f = A2.hom(I, P3)
X = f.image().irreducible_components()[0]
print(X)