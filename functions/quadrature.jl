#########################################################################################################################################
############################## quadrature.jl contains basic univariate quadrature routines ##############################################
#########################################################################################################################################
#
# Types:
#       QuadRule               - object that stores n quadrature points and weights on the interval [-1,1]
#
# Functions:
#       dim()                  - determine the order of the quadrature rule
#
#
#########################################################################################################################################
#
# File part of the Julia IsoGeometric Analysis toolbox
# R.R.Hiemstra
# 19-04-2016
#
#########################################################################################################################################


"container for quadrature formula"
struct QuadRule
    points  :: Vector{Float64}
    weights :: Vector{Float64}
end

"order of the quadrature rule - number of quadrature points"
dim(qr::QuadRule) = length(qr.points)

"""
compute an 'n'-point Gausian integration rule of type 'method' on the interval [-1,1]

## input:
    n       :: Int             - number of quadrature points
    method  :: String          - choose among the following methods: "legendre", "lobatto", "chebyshev", "hermite", "jacobi", #                                 "radau". If left blank then the standard is "legendre"

## output:
    u       :: Vector{Float64} - quadrature points
    w       :: Vector{Float64} - quadrature weights
"""
function QuadRule(n::Int,method="legendre")

    if method=="legendre"
        u, w = gausslegendre(n)
    elseif method=="lobatto"
        u, w = gausslobatto(n)
    elseif method=="chebyshev"
        u, w = gausschebyshev(n)
    elseif method=="hermite"
        u, w = gausshermite(n)
    elseif method=="jacobi"
        u, w = gaussjacobi(n)
    elseif method=="radau"
        u, w = gaussradau(n)
    else
        error("Not a valid quadrature method")
    end

    return QuadRule(u,w)
end
