function ∂r(ψᵣ::RadialMode)
    α = ψ.α; β = ψ.β; γ = ψ.γ;
    ρ = ψ.ρ; aₙ = ψ.aₙ;
    s = ψ.s; l = ψ.l; n = ψ.n;
    """Add a sum over different copies of qnm with some
    change"""
    ψα = RadialMode(s,l,n,α-1,β,γ,ρ,aₙ)
    ψβ = RadialMode(s,l,n,α,β-1,γ,ρ,aₙ)
    aₙshift = circshift(aₙ,-1)
    aₙshift[end] = 0
    nn = 1:length(aₙshift)
    aₙshift = aₙshift .*nn
    ψaₙ = RadialMode(s,l,n,α,β-2,γ,ρ,aₙshift)
    α*ψα + β*ψβ - γ*ψ + ψaₙ
end
