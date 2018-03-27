using ApproxFun, Plots
import ApproxFun: setdomain

#####
# This solves
# u'(x) + u(x-1) = 0
# u(0) = 2.3
# on 0..3 by breaking the problem into three pieces
#####

S = Chebyshev(0..1)
D = Derivative(S)
L = [D            0I      0I;
     I             D      0I;
     0I           I        D]


B = [ldirichlet(S)   0                 0;
     rdirichlet(S)  -ldirichlet(S)     0;
     0              rdirichlet(S)      -ldirichlet(S)]

v = [B; L] \ [2.3; zeros(2); zeros(3)]

u = v[1] + setdomain(v[2], Domain(1..2)) + setdomain(v[3], Domain(2..3))

plot(u)

u'(1.3) + u(1.3-1)  # ≈ 0
