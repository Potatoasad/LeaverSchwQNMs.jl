__precompile__()
module LeaverSchwQNM

include("LinearCombinations.jl")

## Call the qnm package for the Schwarzschild modes
using PyCall
const qnm = PyNULL()

function __init__()
    copy!(qnm, pyimport("qnm"))
    qnm.download_data()
end
__init__()

function GetFreq(s,l,n; qnm=qnm)
    grav_freq = qnm.modes_cache(s=s,l=l,m=0,n=n)
    ω, Alm, Cllʼ = grav_freq(a=0.0)
    ω
end

## Series Calculation
function rNCoeffs(D0,D1,D2,D3,D4)
    fp = sqrt(-D0 - D1 - D2)
    u1 = -fp
    u2 = (1/2)*(-4 - 2*D0 - D1)
    u3 = (fp)*(8 + 16*D0 + 8*(D0^2) + 12*D1 + 8*D0*D1 + (D1^2) + 8*D2 +
    4*D0*D2 - 4*D3 - 4*D4)/(8*(D0 + D1 + D2))
    u4 = (1/2)*(4 + 4*D0 + 2*D0^2 + D1 + D0*D1 - D3)
    u1,u2,u3,u4
end

function ComputeSeriesFromab(an::Function,bn::Function; N=250, rN = 0.0*im, PreN = 40)
    ##Initialize rn and fn vectors
    rₙ = zeros(Complex{Float64},N+1)
    fₙ = zeros(Complex{Float64},N+1)

    rold = rN
    ##Startup Pass for rₙ
    for n = (N+PreN):-1:(N+1)
        rnew = -bn(n)/(an(n) + rold)
        rold = rnew
    end
    rₙ[N+1] = rold;
    ##Pass for rn
    for n= N:-1:1
        rₙ[n] = -bn(n)/(an(n) + rₙ[n+1])
    end

    fₙ[1] = 1;
    ##Pass for fn
    for n = 2:(N+1)
        fₙ[n] = rₙ[n-1]*fₙ[n-1]
    end
    fₙ
end

function RadialCoefficients(D₀, D₁, D₂, D₃, D₄; N = 250)
    αₙ(n) = (n+1)*(n+D₀)
    βₙ(n) = -2*n^2 + (D₁+2)*n + D₃
    γₙ(n) = (n-1)*(n+D₂-2) + D₄

    an(n) = βₙ(n)/αₙ(n) ; bn(n) = γₙ(n)/αₙ(n);

    ##Set largest rN value
    u1,u2,u3,u4 = rNCoeffs(D₀,D₁,D₂,D₃,D₄)
    rN = 1 + u1*N^(-0.5) + u2*N^(-1) + u3*N^(1.5) + u4*N^(-2)

    ComputeSeriesFromab(an,bn; N=N, rN=rN, PreN=0)
end

## Define the Radial Modes
struct RadialMode <: CallableAtom
    s::Int64; l::Int64; n::Int64;
    α::ComplexF64; β::ComplexF64; γ::ComplexF64
    ρ::Complex{Float64};
    aₙ::Vector{Complex{Float64}};
end

function (ψ::RadialMode)(r)
    asymptoticpart = ((r-1)^ψ.α)*(r^(ψ.β))*(exp(-ψ.γ*(r-1)))
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

## Define the derivative of the function
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

export ∂r, RadialMode, GetFreq

end # module
