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
