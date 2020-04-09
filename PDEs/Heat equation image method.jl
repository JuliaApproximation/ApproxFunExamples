using ApproxFun, LinearAlgebra, Plots

# In this example, we solve the heat equation with left Neumann and
# right Dirichlet boundary conditions.
# To allow for a solution through matrix exponentiation, we approximate
# the initial condition with a basis that automatically constructs
# even and odd extensions on the left and right boundaries respectively. (Method of images)

n = 100
xlim = 10
xarray = (chebyshevpoints(n) .+ 1) .* xlim/2 # Chebyshev points on 0..xlim
r = 1.0 # right boundary value (left is assumed to be zero)

S = CosSpace(0..4*xlim)

u₀(x) = atan(x)*(π + atan(x-9) + atan(5-x))/π^2 - r # subtract by `r` to homogenise the right BC
   plot(xarray, u₀)
v = u₀.(xarray) # values to use in fitting

m = div(n,2) # m << n for regularisation

# Create a Vandermonde matrix by evaluating the basis at the grid:
# We ensure the basis satisfies the boundary conditions by using
# the entries corresponding to {cos(x), cos(3x), cos(5x), ...}.
V = zeros(n, m)
for k = 1:m
   V[:,k] = Fun(S,[zeros(2*k-1);1]).(xarray)
end

c₀ = zeros(2*m)
c₀[2:2:end] = V\v # do the fit

f = Fun(S, c₀) # Likely not a good approximation due to boundary discontinuities, but high frequencies will disappear quickly anyway
   plot([xarray.+10; xarray; xarray.-10], f, xlims=(-5, 15))
   scatter!(xarray, u₀, markersize=0.2)

L = Derivative(S, 2)
   L = Matrix(L[1:2m,1:2m])

@manipulate for T=[zeros(5); 0.0:0.01:1]
   uₜ = Fun(S, exp(T*L)*c₀) + r # calculate solution and add back homogenisation term
   plot(-0:0.05:10, uₜ, ylims=(0,1), legend=:none)
end

@gif for T=[zeros(5); 0.0:0.5:50]
   uₜ = Fun(S, exp(T*L)*c₀) + r # calculate solution and add back homogenisation term
   plot(-0:0.05:10, uₜ, ylims=(0,1), legend=:none)
end
