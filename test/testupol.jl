pNULL = UPol(Float32[])
p0 = UPol([0])
p1 = UPol([0,0,0,0,0,0,0,0,0,0,0,0,0,1])
p2 = UPol([0,0,1,1])
p3 = UPol([0,0,0,0,1,2,1])
p4 = UPol([0,1,3,3,1])
p5 = UPol([0,0,0,0,0,0,0,0,0,0,0,0,0,1,4,6,4,1])
pN = UPol([0,24,15,87,3,276])
pR = UPol([1//1, -2//1, 3//4])
p1000 = UPol(randn(1000))

@test length(pNULL) == 0
@test length(p1000) == 1000
sprint(show, p1000)
sprint(show, pNULL)

@test p3 == UPol([1,2,1])
@test pN*10 == UPol([240, 150, 870, 30, 2760])
@test pN/10 == UPol([2.4, 1.5, 8.7, 0.3, 27.6])
@test 10*pNULL + pN == pN
@test 10*p0 + pN == pN
@test p5 + 2*p1 == UPol([1,4,6,4,3])
@test 10*pNULL - pN == -pN
@test p0 - pN == -pN
@test p5 - 2*p1 == UPol([1,4,6,4,-1])
@test p2*p2*p2 == p4
@test p2^4 == p5
@test pNULL^3 == pNULL
@test pNULL*pNULL == pNULL

@test polyval(pN, -.125) == 276.9609375
@test polyval(pNULL, 10) == 0
@test polyval(p0, -10) == 0
@test polyval(poly([1//2, 3//2]), 1//2) == 0//1
@test polyder(polyint(pN)) == pN
@test polyder(pR) == UPol([2//1, -2//1])
@test polyint(pNULL,1) == p1
@test polyint(UPol(Rational[3, 2, 1])) == UPol(Rational[1, 1, 1, 0])
@test polyder(p3) == UPol([2,2])
@test polyder(p1) == polyder(p0) == polyder(pNULL) == pNULL

@test poly([-1,-1]) == p3
@test roots(p0)==roots(p1)==roots(pNULL)==[] 
@test roots(p2) == [-1]
a_roots = copy(pN.a)
@test all(abs(sort(roots(poly(a_roots))) - sort(a_roots)) .< 1e6)
@test length(roots(p5)) == 4
@test roots(pNULL) == []
@test roots(pR) == [1//2, 3//2]

@test pNULL + 2 == p0 + 2 == 2 + p0 == UPol([2])
@test p2 - 2 == -2 + p2 == UPol([1,-1])
@test 2 - p2 == UPol([-1,1])

p0 = UPol([0])
p1 = UPol([1])
p2 = UPol([4, 2, -3, 6, 5])
p3 = UPol([6, 2, -3, 7])
p4 = p2 * p3
@test divrem(p4, p2) == (p3, zero(p3))
@test p3%p2 == p3
@test all((abs((p2/p3 - UPol([2/3,1/9])).a)) .< eps())
@test divrem(p0,p1) == (p0,p0)
@test divrem(p1,p1) == (p1,p0)
@test divrem(p2,p2) == (p1,p0)
@test divrem(pR, pR) == (one(pR), zero(pR))
@test_throws p1/p0
@test_throws divrem(p0,p0)

#Tests for multivariable support
pX = UPol([1, 2, 3, 4, 5])
#pS1 = UPol([1, 2, 3, 4, 5], "s")
#pS2 = UPol([1, 2, 3, 4, 5], 's')
#pS3 = UPol([1, 2, 3, 4, 5], :s)
#@test pX != pS1
#@test pS1 == pS2
#@test pS1 == pS3
#@test_throws pS1 + pX
#@test_throws pS1 - pX
#@test_throws pS1 * pX
#@test_throws pS1 / pX
#@test_throws pS1 % pX
