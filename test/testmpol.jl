x = monomial(MPol{Float64}, :x)
y = monomial(MPol{Float64}, :y)

@test 1 + x == 1.0 + x
@test eltype(x + im) == Complex{Float64}
@test (x+y)^2 == x^2 + 2x*y + y^2
@test vars(x+y+one(MPol{Float64}, [:u, :v, :w])) == [:x, :y, :u, :v, :w]

p = x^2 * y + 3x + 9

@test p[0,0] == 9.0
@test vars(p) == [:x, :y]

r = zero(MPol{Float64}, [:x, :y])
r += 10y
r += 42x

@test r[[1,0]] == 42.0