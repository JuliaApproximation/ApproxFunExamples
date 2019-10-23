using ApproxFun, LinearAlgebra, DualNumbers
import ApproxFun: bandwidth

#
# Consider the following Sturm--Liouville problem with polynomial p, q, and w.
#
# [(-𝒟)(p𝒟) + q] u = λ w u,  (u±bu')(±1) = 0.
#
# For which value of b, the Neumann part of the Robin boundary conditions,
# is λ₀ equal to 1?
#
# To answer this question, we use a dual number to observe the sensitivity in λ₀
# with respect to perturbations in b. The `realpart` of the objective function
#
#    f(λ(b)) = λ(b) - 1,
#
# extracts the difference, whereas the `dualpart` evaluates to ∂f/∂b (λ(b)).
#

n = 100
λ = dual(0.0, 0.0)
b = dual(0.0, 1.0)
v = zeros(Dual{Float64}, n)
v[1] = 1
v .+= 0.003.*randn(Float64, n)
for _ in 1:8
    d = Segment(-1..1)
    S = Ultraspherical(0.5, d)
    NS = NormalizedPolynomialSpace(S)
    D = Derivative(S)
    p = Fun(x->1+x^2, S)
    q = Fun(x->x^3, S)
    L = -D*(p*D)+q
    w = Fun(x->1+x/4, S)
    M = Multiplication(w, S)
    C = Conversion(domainspace(L), rangespace(L))
    B = [Evaluation(S, -1, 0) - b*Evaluation(S, -1, 1); Evaluation(S, 1, 0) + b*Evaluation(S, 1, 1)]
    QS = QuotientSpace(B)
    Q = Conversion(QS, S)
    D1 = Conversion(S, NS)
    D2 = Conversion(NS, S)
    R = D1*Q
    P = cache(PartialInverseOperator(C, (0, bandwidth(L, 1) + bandwidth(R, 1) + bandwidth(C, 2))))
    A = R'D1*P*L*D2*R
    B = R'D1*M*D2*R
    SA = Symmetric(A[1:n,1:n], :L)
    SB = Symmetric(B[1:n,1:n], :L)
    for _ = 1:7
        global λ = dot(v, SA*v)/dot(v, SB*v)
        ldiv!(ldlt!(SA-λ*SB), v)
        global v = v/norm(v)
        println("\tThis is λ: ", λ)
    end
    global b -= (realpart(λ)-1)/dualpart(λ)
    println("This is b: ", b)
end
