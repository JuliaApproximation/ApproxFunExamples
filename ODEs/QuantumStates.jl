using ApproxFun, Plots

x = Fun(-20 .. 20)
V = abs(x)
D = Derivative()
S = space(V)
B = [Dirichlet(S); continuity(S,0:1)]
λ, v = eigs(B, -D^2 + V, 500,tolerance=1E-10)

p = plot(V; legend=false)
    for k=1:20
        plot!(real(v[k]/norm(v[k]) + λ[k]))
    end
    p
