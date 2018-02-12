using LogarithmicGreensFunctions
using Base.Test

using PyPlot

symset!(A,i,j,v) = (A[i,j] = A[j,i] = v)


# One interval

n = 100
H = full(Tridiagonal(-ones(n-1),zeros(n),-ones(n-1)))
symset!(H,n,1,-1)
d = [minimum(abs(i-j+s*n) for s = -1:1) for i = 1:n, j = 1:n]

figure()
for (i,z) in enumerate((3,im,0.1-0.1*im))
    semilogy(vec(d), vec(abs.(inv(H - z*I))), "C$(i-1)o", ms=2,mew=0);
    semilogy(0:n÷2, exp.(-greens(domain(-2,2),z)*(0:n÷2)), "C$(i-1)--")
end
show()


# Two intervals

H = full(Tridiagonal(-ones(n-1),repmat([-1,1],n÷2),-ones(n-1)))
symset!(H,n,1,-1)
d = [minimum(abs(i-j+s*n) for s = -1:1) for i = 1:n, j = 1:n]

figure()
for (i,z) in enumerate((0,0.5,3,2+0.1im))
    semilogy(vec(d), vec(abs.(inv(H - z*I))), "C$(i-1)o", ms=2,mew=0);
    semilogy(0:n÷2, exp.(-greens(domain(-sqrt(5),-1,1,sqrt(5)),z)*(0:n÷2)), "C$(i-1)--")
end
show()
