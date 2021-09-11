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
    ω/2
end

## Series Calculation
include("SeriesCalculation.jl")

## Define the Radial Modes
include("RadialModes.jl")

## Define the derivative of the function
include("Derivatives.jl")

export ∂r, RadialMode, GetFreq

end # module
