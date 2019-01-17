using ApproxFun

###
# Fractional Derivative
###

S = Legendre()

n = 50; # number of data points
p = range(-1,stop=1,length=n);   # a non-default grid

v = exp.(p);           # values at the non-default grid
m = 20   # number of basis functions
V = Array{Float64}(undef,n,m); # Create a Vandermonde matrix by evaluating the basis at the grid
for k = 1:size(V,2)
   V[:,k] = Fun(S,[zeros(k-1);1]).(p)
end
f = Fun(S,V\v); # least squares fit to data

f_rl = LeftDerivative(0.5)*f # RL left half-derivative
f_cap = LeftIntegral(0.5)*Fun(f', S) # Caputo left half-derivative


using Plots
plot(f_rl;ylims=(0,3),legend=:bottomright,label="RL")
plot!(f_cap; label="Caputo")
savefig("RLvCaputo.png")
