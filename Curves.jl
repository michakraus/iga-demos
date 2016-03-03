module Curves

# export types
export Curve, Degree, KnotVector, BsplineBasisFun


import PyPlot.plot
include("functions/BaseExtension.jl")
include("functions/knotvector.jl")
include("functions/algebraicoperators.jl")
include("functions/bsplinebasisfuns.jl")
include("functions/knotrefinement.jl")
include("functions/projection.jl")
include("functions/quadrature.jl")

##################################################################################################
##################################### Type BsplineBasisFun #######################################
##################################################################################################

# typedef single univariate basis function
type BsplineBasisFun <: AbstractVector{Float64}
    knots  :: Vector{Float64}
    vals   :: Vector{Float64}
end
BsplineBasisFun(knots::Vector{Float64}) = BsplineBasisFun(knots,Float64[])
BsplineBasisFun() = BsplineBasisFun([0.0,1.0])

# set and get methods
getindex(kts::BsplineBasisFun, I::RangeIndex) = kts.knots[I]
setindex!(kts::BsplineBasisFun, X, I::RangeIndex) = setindex!(kts.knots, X, I)

# SubMatrix functionality
size(kts::BsplineBasisFun) = size(kts.knots)
length(kts::BsplineBasisFun) = length(kts.knots)

# plot b-spline basis function
function plot(b::BsplineBasisFun,gran::Int64,linetype::String="-b")
    # evaluate b-spline bais function
    ukts = unique(b.knots)
    uu   = linspace(ukts[1],ukts[end],gran*(length(ukts)-1))
    ff   = onebasisfuneval(b,uu)
    
    # call matplotlib
    plot(uu,ff,linetype)
end
plot(b::BsplineBasisFun,linetype::String="-b") = plot(b,10,linetype)
plot(b::Vector{BsplineBasisFun},gran::Int64,linetype::String="-b") = [plot(b[i],gran,linetype) for i in 1:length(b)]
plot(b::Vector{BsplineBasisFun},linetype::String="-b") = [plot(b[i],linetype) for i in 1:length(b)]

# evaluate a single basis function
onebasisfuneval(b::BsplineBasisFun,u::Vector{Float64}) = onebasisfuneval(length(b.knots)-2,b.knots,u)
function onebasisfuneval!(b::BsplineBasisFun,u::Vector{Float64})
    b.vals = onebasisfuneval(length(b.knots)-2,b.knots,u)
end



##################################################################################################
##################################### Type Curve #################################################
##################################################################################################

type Curve
    degree  :: Degree
    knots   :: KnotVector
    control :: Matrix{Float64}
    weights :: Vector{Float64}
end
function Curve(p,kts) 
	n = dimsplinespace(p,kts)
	P = zeros(Float64,n,2); P[:,1] = grevillepoints(p,kts)
	Curve(p,kts, P, ones(Float64,n))
end


function lsquarederror(C::Curve,func::Function,nquad::Int)
    # initialize
    p   = C.degree
    kts = C.knots
    CP  = C.control
    CW  = C.weights
    n = length(CP)
    l2e = 0.0
    
    # compute gauss nodes on [-1,1]
    u,w =  Base.gauss(Float64,nquad)
    
    # loop over elements
    for k in p+1:length(kts)-p-1

        if kts[k+1]!=kts[k]
            
            # detrminant of mapping from [-1,1] to spline parameter space
            detJ = 0.5 * (kts[k+1]-kts[k])
            uu   =  mean([kts[k+1],kts[k]]) + u * detJ 

            for i in 1:nquad
                # evaluate basis functions
                Nu = bsplinebasisfuns(p, kts, k, uu[i])
                W = dot(Nu, CW[k-p:k])
                f = dot(Nu, (CW[k-p:k].*CP[k-p:k,end]) / W)

                # compute contribution of kth quadrature point to Grammian matrix 
                l2e += (f - func(uu[i])).^2 * (w[i] * detJ)    
            end
        end
    end
    return sqrt(l2e)
end
lsquarederror(C::Curve,func::Function) = lsquarederror(C,func,C.degree+1)

function plot(C::Curve, n::Int=5,cps::String="off",knots::String="off")

    # initialize variables
    p   = C.degree
    kts = C.knots
    CP  = C.control
    CW  = C.weights
    dim = size(CP,2)
    # CP = broadcast(*,CW, CP)

    # loop over elements
    for k in p+1:length(kts)-p-1
        if kts[k+1]!=kts[k]
            u = linspace(kts[k],kts[k+1], n)
            Nu = dersbsplinebasisfuns(p, kts, fill(k,n), u, 1)[1]
            W = Nu * CW[k-p:k]
            f = (Nu * broadcast(*,(CW[k-p:k]), CP[k-p:k,:]))
            f = broadcast(/,f,W)

            # plot element k
    		if dim==1
	            plot(u, f,"-b")
                if knots=="on"
	               plot(u[1],   f[1],"or")
	               plot(u[end], f[end],"or")
                end
	        elseif dim==2
	            plot(f[:,1], f[:,2],"-b")
                if knots=="on"
	               plot(f[1,1], f[1,2],".r")
	               plot(f[end,1], f[end,2],".r")
                end
	        elseif dim==3
	        	plot3d(f[:,1],f[:,2],f[:,3])
                if knots=="on"
	               plot3d(f[1,1], f[1,2],f[1,3],"or")
	               plot3d(f[end,1], f[end,2],f[end,3],"or")
                end
	        end
        end
    end

    # plot controlpoints
     if cps=="on"
        if dim==1
            plot(grevillepoints(p,kts),CP,"-oy") 
        elseif dim==2
            plot(CP[:,1],CP[:,2],"-oy") 
        elseif dim==3
            plot3d(CP[:,1], CP[:,2], CP[:,3],"or")
        end
    end

end

end