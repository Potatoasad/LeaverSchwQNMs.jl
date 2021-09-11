## Whether the derivative function works
ψ = RadialMode(2,2,0)
ψ2 = ∂r(ψ)

ϵ = 0.001
r = 1.01:ϵ:3.01
ds = ψ.(r);
dst = ψ2.(r)
dsf = [(ds[i+1]-ds[i])/ϵ for i ∈ 1:(length(ds)-1)]

using Plots
p1 = plot(r[1:(end-1)],real.(dsf))
p2 = plot!(p1,r[1:(end-1)],real.(dst[1:(end-1)]))
display(p2)
