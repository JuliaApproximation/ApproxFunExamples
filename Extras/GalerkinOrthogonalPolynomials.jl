using ApproxFun
using Plots
using LinearAlgebra

# In this example, we numerically compute with Galerkin orthogonal polynomials.
# If H is a separable Hilbert space and B contains ν linear functional constraints,
# then the Galerkin orthogonal polynomials are those orthogonal polynomials with
# respect to H_B := \{ u \in H : Bu = 0 \}.
# For example, if H = L^2([-1,1]) and we wish to enforce that u(-1) = u'(-1) =
# u'(1) = 0, then the 2-normalized orthonormal polynomials may be created by
# taking the QR factorization of the connection coefficients between the
# normalized Legendre polynomials \sqrt{n+1/2}P_n(x) and a degree-graded
# polynomial basis that satisfies the constraints.

N = 501

NS = NormalizedLegendre()
B = [Evaluation(NS, -1, 0); Evaluation(NS, -1, 1); Evaluation(NS, 1, 1)]
ν = size(B, 1)
QS = QuotientSpace(B)
A = Conversion(QS, NS)
Q, R = qr(A[1:N+ν,1:N])

p = Vector{Fun}(undef, N)
for n = 1:N
    en = [zeros(n-1);1.0;zeros(N+ν-n)]
    p[n] = Conversion(NS, Legendre())*Fun(NS, (Q*en)[1:n+ν])
    pad!(p[n], 3n+50)
end

# The mass matrix with entries M_{i,j} = \int_{-1}^1 p_i(x) p_j(x) dx,
# is numerically the identity, as we expect.

M = [innerproduct(p[i], p[j]) for i = 1:N, j = 1:N]
@show norm(M-I) ≤ 4*norm(M)*eps()
@show opnorm(M-I) ≤ 4*sqrt(N)*opnorm(M)*eps()

# Although the polynomials now satisfy deg(p_n) = n+ν instead of deg(p_n) = n,
# the index appears to determine the number of roots of the Galerkin orthogonal
# polynomials. In particular, p_n has exactly n roots in (-1,1).

pl = plot(p[1]; legend = false)
for n in 2:10
    plot!(p[n])
end
plot(pl)

# This need not be the case! If the second constraint is modified to impose
# u'(0) = 0, then the odd-index polynomials have one more root than their index
# and the even-index polynomials have one fewer.

B = [Evaluation(NS, -1, 0); Evaluation(NS, 0, 1); Evaluation(NS, 1, 1)]
QS = QuotientSpace(B)
A = Conversion(QS, NS)
Q, R = qr(A[1:N+ν,1:N])

q = Vector{Fun}(undef, N)
for n = 1:N
    en = [zeros(n-1);1.0;zeros(N+ν-n)]
    q[n] = Conversion(NS, Legendre())*Fun(NS, (Q*en)[1:n+ν])
    pad!(q[n], 3n+50)
end

pl = plot(q[1]; legend = false)
for n in 2:10
    plot!(q[n])
end
plot(pl)

# More information on Galerkin orthogonal polynomials is available in
#
#   P. W. Livermore. Galerkin orthogonal polynomials, J. Comp. Phys., 229:2046–2060, 2010.
#
# The banded QR factorization of the connection coefficients is a linear
# complexity algorithm to represent the Galerkin orthogonal polynomials in H_B
# in terms of the polynomial basis for H. This is described in
#
#   J. L. Aurentz and R. M. Slevinsky. On symmetrizing the ultraspherical spectral method for self-adjoint problems, arXiv:1903.08538, 2019.
#
