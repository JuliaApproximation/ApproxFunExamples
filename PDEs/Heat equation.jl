####
# The following solves heat equation with mixed boundary conditions.
# we first formulate the construction as
#
# B*u = r
# C*u_t = L*u
#
# where B represents the mixed boundary conditions, r is their values
# C is a conversion matrix and L is the Laplacian.
####

using ApproxFun, Plots, LinearAlgebra

S = Chebyshev()
D = Derivative() : S

C = I : domainspace(D) → rangespace(D)

n = 10
#
# Here we construct the boundary conditions. it is beneficial if the principle
# sub block is identity
#
B = [ldirichlet(S); rneumann(S)]
    B = B[1:2,1:n]
    B = B[1:2,1:2]\B



L = D^2

#Differentiating with respect
# to time we see
#
# B*u_t = 0
#
# We can then use this to eliminate dependence on the first two columns to get
#
# M*u_t = L*u
#
# where M is zero in the first two columns. As L is also 0 we can drop two columns
# to get
#
# M[:,3:end]*u_t = L[:,3:end]*u
#
# Now we multiply through to get
#
# u_t = A * u
#
# for A = M[:,3:end] \ L[:,3:end]
# This a standard linear ODE and the solution is u(t) = exp(A*t)*u(0).



M = C[1:n-2,1:n]
    M = M - M[:,1:2]*B # remove degrees of freedom
    A = M[:,3:end]\L[1:n-2,3:n]


u₀ = Fun(exp, S)
c₀ = pad(u₀.coefficients,n)
r = B*c₀ # boundary condition values, assume they stay fixed here

p = plot(u₀; legend=false, ylims=1:6)
@gif for t in 0.0:0.1:10
    c̃_t = exp(A*t)*c₀[3:n]

    c_t =
    [r - B[:,3:end]*c̃_t;
    c̃_t]

    plot(Fun(S, c_t); legend=false, ylims=(0,6))
end
