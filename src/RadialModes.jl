struct RadialMode <: CallableAtom
    s::Int64; l::Int64; n::Int64;
    α::ComplexF64; β::ComplexF64; γ::ComplexF64
    ρ::Complex{Float64};
    aₙ::Vector{Complex{Float64}};
end

function (ψ::RadialMode)(r)
    asymptoticpart = ((im*(r-1))^ψ.α)*((im*r)^(ψ.β))*(exp(-ψ.γ*(r-1)))
    z = (r-1)/r
    s=Complex(0.0)
    #sum(ψ.aₙ[n]*z^(n-1) for n ∈ 1:length(ψ.aₙ))
    for n ∈ 1:length(ψ.aₙ)
        s += ψ.aₙ[n]*z^(n-1)
    end
    asymptoticpart*s
end

function RadialMode(s,l,n)
    ρ = -im*GetFreq(s,l,n);
    α = ρ;
    β = -2*ρ;
    γ = ρ;

    D₀ = 1 + 2*ρ;
    D₁ = -8*ρ-4;
    D₂ = 4*ρ + 3;
    D₃ = (1-s^2) - 8*ρ^2 - 4*ρ - l*(l+1);
    D₄ = 4*ρ - (1-s^2) - 1 + 4*ρ^2;

    aₙ = RadialCoefficients(D₀, D₁, D₂, D₃, D₄; N = 250)

    RadialMode(s,l,n,α,β,γ,ρ,aₙ)
end

import Base.show

function Int2Sub(num,converter)
    if num < 0
       strnum = string(num)[2]
       return "₋"*converter[strnum]
    else
       strnum = string(num)[1]
       return converter[strnum]
    end
end

function Base.show(io::IO, ψ::QuasinormalModeFunction)
    BigDigits = "+-0123456789"
    SmallDigits = "₊₋₀₁₂₃₄₅₆₇₈₉"
    Subdict = Dict(zip(BigDigits,SmallDigits))
    s = Int2Sub(ψ.s,Subdict)
    l = Int2Sub(ψ.l,Subdict)
    n = Int2Sub(ψ.n,Subdict)
    print(io,s*"Ψ"*l*n)
end
