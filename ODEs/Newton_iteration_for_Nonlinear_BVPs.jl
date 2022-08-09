using ApproxFun, Plots

x=Fun()
u0=0.0x
N=u->[u(-1.),u(1.),0.001u''+6*(1-x^2)*u'+u^2-1.]
u=newton(N,u0)


x=Fun()
u0=0.0x
N=u->[u(-1.),u(1.),0.1u''+sin(u)-1.]
u=newton(N,u0)


##
# Heaven and hell for PI
##

M = 10.0
x = Fun(0..M)
N = u -> [u(0), u(M)-sqrt(M), u'' - u^2 + x]
u0 = 0.0x
@time heaven = newton(N, u0)

u0 = heaven-10x * exp(-x)
@time hell = newton(N, u0)
plot(heaven; label="heaven", linewidth=1.5, legend=:bottomright, title="Painleve I: Heaven and Hell")
plot!(hell; label="hell", linewidth=1.5, linestyle=:dash)
savefig("heavenhell.pdf")

scatter(abs.(heaven.coefficients); yscale=:log10,  label="heaven", markersize=1.5, title="Painleve I: Heaven and Hell Coefficients")
scatter!(abs.(hell.coefficients); yscale=:log10,  label="hell", markersize=1.5, marker=:square)
savefig("heavenhellcoefficients.pdf")