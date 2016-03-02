#########################################################################################################################################
################################ quadrature.jl features different Gauss integration rules for use in FEA ################################
#########################################################################################################################################
#
# Functions:
#       gaussquadrule()     - provides different n-point quadrature rules for use in FEA
#                           - flexibility in the choice of interval
#
#########################################################################################################################################
#
# File part of the Julia IsoGeometric Analysis toolbox
# R.R.Hiemstra
# 14-09-2015
#
#########################################################################################################################################

include("quadrature/gausslegendre.jl")
include("quadrature/gausschebyshev.jl")
include("quadrature/gausshermite.jl")
include("quadrature/gaussjacobi.jl")
include("quadrature/gausslobatto.jl")
include("quadrature/gaussradau.jl")

# compute an 'n'-point Gausian integration rule of type 'method' on the interval [-1,1]
function gaussquadrule(n::Int,method="legendre")
# input:
#    n       :: Int             - number of quadrature points
#    method  :: String          - choose among the following methods: "legendre", "lobatto", "chebyshev", "hermite", "jacobi", #                                 "radau". If left blank then the standard is "legendre"
# output:
#    u       :: Vector{Float64} - quadrature points
#    w       :: Vector{Float64} - quadrature weights

    if method=="legendre"
        return gausslegendre(n)
    elseif method=="lobatto"
        return gausslobatto(n)
    elseif method=="chebyshev"
        return gausschebyshev(n)
    elseif method=="hermite"
        return gausshermite(n)
    elseif method=="jacobi"
        return gaussjacobi(n)
    elseif method=="radau"
        return gaussradau(n)
    else
        error("Not a valid quadrature method")
    end
end

# compute an 'n'-point Gausian integration rule of type 'method' on the interval [a,b]
function gaussquadrule(a::Float64,b::Float64,n::Integer,method="legendre")
# input:
#    a       :: Float64         - left boundary of integration interval
#    b       :: Float64         - right boundary of integration interval
#    n       :: Int             - number of quadrature points
#    method  :: String          - choose among the following methods: "legendre", "lobatto", "chebyshev", "hermite", "jacobi", 
#                                 "radau". If left blank then the standard is "legendre"
# output:
#    u       :: Vector{Float64} - quadrature points
#    w       :: Vector{Float64} - quadrature weights

    # gauss quadrature points and weights on [-1,1]
    (xx, ww) = gaussquadrule(n,method)

    # compute weights and nodes B-spline parametric space
    center = 0.5 * (a + b)
    det    = 0.5 * (b - a)
    for i in 1:n
      xx[i] = center + xx[i] * det
      ww[i] = ww[i] * det
    end

    return xx, ww
end

# compute an 'n'-point Gaussian integration rule of type 'method' on the the partition of elements 'unknots'
function gaussquadrule(unknots::Array{Float64,1},n::Integer,method="legendre")
# input:
#    unknots :: Vector{Float64} - vector of increasing numbers representing a univariate partition
#    n       :: Int             - number of quadrature points
#    method  :: String          - choose among the following methods: "legendre", "lobatto", "chebyshev", "hermite", "jacobi", 
#                                 "radau". If left blank then the standard is "legendre"
# output:
#    u       :: Matrix{Float64} - Matrix with quadrature points:  u[:,k] are the quadrature points in the kth element of       #                                 'unknots'
#    w       :: Matrix{Float64} - Matrix with quadrature weights: w[:,k] are the quadrature weights in the kth element of       #                                 'unknots'

    # number of elements
    m = length(unknots)-1

    # compute weights and nodes B-spline parametric space
    u = zeros(Float64,n,m)
    w = zeros(Float64,n,m)
    for j in 1:m
        u[:,j], w[:,j] = gaussquadrule(unknots[j],unknots[j+1],n,method)
    end
    return u, w
end

# compute optimal quadrature rule on a uniform mesh of 'e' elements
function optimalquadrule(p,r,e,a,b)
# this function computes the quadrature rule for the target splinspace 
# S^p_1 on a uniform partition with open knotvector.
# input:
#    p :: Int     - polynomial degree
#    r :: Int     - internal regularity in between elements
#    m :: Int     - number of elements
#    a :: Float64 - left boundary of parametric domain (first knot)
#    b :: Float64 - right boundary of parametric domain (last knot)

    # initialize quadrature rules
    n = (p+1)*2 + (e-1)*(p-r) - p-1   # dimension of the spline space
    m = int64(ceil(n/2))              # dimension of the quadrature rule
    nodes   = zeros(Float64,m)        # allocate space for nodes
    weights = zeros(Float64,m)        # allocate space for weights

    # load data & initialize
    if p==4 && r==0
        if iseven(e)
            qr = matread(string("quadrule_p$p","_r$r","/quadrule_even_r$r","_p$p.mat"))
        else
            qr = matread(string("quadrule_p$p","_r$r","/quadrule_odd_r$r","_p$p.mat"))
        end
    else
        if iseven(m) 
            qr = matread(string("quadrule_p$p","_r$r","/quadrule_even_r$r","_p$p.mat"))
        else
            println("hoi")
            qr = matread(string("quadrule_p$p","_r$r","/quadrule_odd_r$r","_p$p.mat"))
        end
    end
    bnodes = qr["bnodes"]; bweights = qr["bweights"]; m1 = length(bnodes)
    inodes = qr["inodes"]; iweights = qr["iweights"]; m2 = length(inodes)
    mnodes = qr["mnodes"]; mweights = qr["mweights"]; m3 = length(mnodes)


    # rule for even number of nodes
    if 2m==n

        # rule with repetetive structure
        if m>2m1

            # construct quadrature rule in the uniform parameter domain with unit spacing
            z = ceil(0.5m)
            nodes[1:m1]   = bnodes
            weights[1:m1] = bweights
            left  = ceil(bnodes[end])
            right = ceil(inodes[end]) 
            ii = m1
            while ii<z
                nodes[ii+1:ii+m2]   = inodes + left
                weights[ii+1:ii+m2] = iweights
                left                = left + right
                ii+=m2
            end
            nodes[m-z+1:m] = e - flipud(nodes[1:z])
            weights[m-z+1:m] = flipud(weights[1:z])

        # load precomputed rule    
        else
            qr = matread(string("quadrule_p$p","_r$r","/quadrule_e$e","_r$r","_p$p.mat"))
            nodes, weights = qr["nodes"], qr["weights"]
        end

    # rule for odd number of nodes
    else

        # rule with repetetive structure
        if m>2*(m1+m3)

            # construct quadrature rule in the uniform parameter domain with unit spacing
            z = ceil(0.5m)
            nodes[1:m1]   = bnodes
            weights[1:m1] = bweights
            left  = ceil(bnodes[end])
            right = ceil(inodes[end]) 
            ii = m1
            while ii<z
                nodes[ii+1:ii+m2]   = inodes + left
                weights[ii+1:ii+m2] = iweights
                left              = left + right
                ii+=m2
            end
            nodes[z-m3+1:z]   = 0.5e + mnodes 
            weights[z-m3+1:z] = mweights 

            nodes[m-z+1:m] = e - flipud(nodes[1:z])
            weights[m-z+1:m] = flipud(weights[1:z])

        # load precomputed rule
        else
            qr = matread(string("quadrule_p$p","_r$r","/quadrule_e$e","_r$r","_p$p.mat"))
            nodes, weights = qr["nodes"], qr["weights"]
        end

    end

    # rescale rule for elementsize h
    h = (b-a) / e
    nodes = a + h*nodes
    weights = h*weights

    return nodes, weights
end