######
# This example shows how to do function approximation in 2D using Padua points
# and in 3D by hand
######

using ApproxFun
    import ApproxFun: transform
f = (x,y) -> exp(x*y)

S = Chebyshev(0..1)^2

N = 10^2; p = points(S, N) # 105 padua points
f̃ = (xy) -> f(xy[1], xy[2])
F = f̃.(p)
F̌ = transform(S, F)  # 105 Chebyshev^2 coefficients ordered by polynomial degree
f̌ = Fun(S, F̌)
f̌(0.1,0.2) ≈ f(0.1,0.2)


S = Chebyshev(0..1)
p_x = points(S, 10)
p_y = points(S, 10)'
p_z = reshape(points(S, 10), 1, 1, 10)

f = (x,y,z) -> exp(x*y+z)

F = f.(p_x, p_y, p_z)

for k = 1:size(F,2),j=1:size(F,3)
    F[:,k,j] = transform(S,F[:,k,j])
end
for k = 1:size(F,1),j=1:size(F,3)
    F[k,:,j] = transform(S,F[k,:,j])
end
for k = 1:size(F,1),j=1:size(F,2)
    F[k,j,:] = transform(S,F[k,j,:])
end

f̌ = function(x,y,z)
    ret1 = Array{Float64}(undef, size(F,2), size(F,3))
    for j=1:size(F,2),l=1:size(F,3)
        ret1[j,l] = Fun(S,F[:,j,l])(x)
    end
    ret2 = Array{Float64}(undef, size(F,3))
    for l=1:size(F,3)
        ret2[l] = Fun(S,ret1[:,l])(y)
    end
    Fun(S, ret2)(z)
end

f̌(0.1,0.2,0.3) ≈ f(0.1,0.2,0.3)
