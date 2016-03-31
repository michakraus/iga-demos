module IGA

# import dependent packages
using FastGaussQuadrature

# export types
export Degree, KnotVector, BezierPatch, NURBS, QuadRule

# export functionality from other libs
export ndgrid, allcomb, linearspace, boundary
export buildvector, uniqueknots, globalrefine, findspan
export dimsplinespace, grevillepoints, onebasisfuneval, bsplinebasiseval, dersbsplinebasisfuns, dersbsplinebasisfunsinterpolationmatrix
export knotinsertionoperator!, knotinsertionoperator, degreeelevationoperator, bezierdecompoperator
export dim

# export libs from NURBS
export nsdim, nbezierpatches, ienarray
export refine!, globalrefine!, degreeelevate!
export boundary, ∂

include("functions/baseextension.jl")
include("functions/knotvector.jl")
include("functions/bsplinebasisfuns.jl")
include("functions/knotrefinement.jl")

# quadrature
type QuadRule
    points  :: Vector{Float64}
    weights :: Vector{Float64}
end

dim(qr::QuadRule) = length(qr.points)

# compute an 'n'-point Gausian integration rule of type 'method' on the interval [-1,1]
function QuadRule(n::Int,method="legendre")
# input:
#    n       :: Int             - number of quadrature points
#    method  :: String          - choose among the following methods: "legendre", "lobatto", "chebyshev", "hermite", "jacobi", #                                 "radau". If left blank then the standard is "legendre"
# output:
#    u       :: Vector{Float64} - quadrature points
#    w       :: Vector{Float64} - quadrature weights

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


###################################################################################################################
############################################ Type definition BezierPatch ##########################################
###################################################################################################################

typealias Point Vector{Float64}

# type definition of Bezier patch
type BezierPatch{N}
    deg    :: NTuple{N,Degree}
    cpts   :: Array{Float64}
    wpts   :: Vector{Float64}
    dec    :: NTuple{N,Matrix{Float64}}
end

# contruct a NURBS patch in parametric space
function BezierPatch{N}(p::NTuple{N,Degree})
    cpts = allcomb(ntuple(k -> grevillepoints(p[k],buildvector([-1.0;1.0], [p[k]+1,p[k]+1])), N)...)
    wpts = ones(Float64,size(cpts,1))
    C    = ntuple(k -> eye(p[k]+1), N)
    return BezierPatch(p, cpts, wpts, C)
end

############################################ Elementary operations ################################################

nsdims(S::BezierPatch) = length(S.deg)

function boundary(S::BezierPatch{1}, k::Int)
    if k==1
        return vec(S.cpts[1,:])
    elseif k==2
        return vec(S.cpts[end,:])
    else
        error("Out of bounds")
    end
end

# boundary of a 2d Bezier Patch
function boundary{N}(S::BezierPatch{N}, dir::Int, k::Int)
    
    # initialize
    dims = S.deg+1  
    
    # local numbering
    I = squeeze(slicedim(reshape(1:prod(dims),dims), dir, k),dir)[:]
    J = 1:N.!=dir

    # return boundary k
    return BezierPatch(S.deg[J], S.cpts[I,:], S.wpts[I], S.dec[J])
end

# entire boundary of 
boundary(S::BezierPatch{1})    = [∂(S,1) ∂(S,2)]
boundary{N}(S::BezierPatch{N}) = BezierPatch{N-1}[∂(S,dir,k) for k in 1:2, dir in 1:N][:]

# short name
∂(S::BezierPatch, dir::Int, k::Int) = boundary(S, dir, k)
∂(S::BezierPatch, k::Int) = boundary(S, k)
∂(S::BezierPatch) = boundary(S)


###################################################################################################################
############################################ Type definition NURBS ###############################################
###################################################################################################################

# type definition of NURBS patch
type NURBS{N}
    deg    :: Vector{Degree}
    kts    :: Vector{KnotVector}
    cpts   :: Matrix{Float64}
    wpts   :: Vector{Float64}
    dec    :: NTuple{N,Array{Float64}}
    ien    :: Array{Int64}
    id     :: Array{Int64}
end
NURBS(p, kts, cpts, wpts) = NURBS(p, kts, cpts, wpts, ntuple(k->Float64[], length(p)),Int64[], Int64[])

# contruct a NURBS patch in parametric space
function NURBS(p, kts)
    cpts = allcomb(ntuple(k -> grevillepoints(p[k],kts[k]), length(p))...)
    wpts = ones(Float64,size(cpts,1))
    return NURBS(p, kts, cpts, wpts)
end

# dimension of patch
nsdim(S::NURBS) = length(S.deg)
dimsplinespace(S::NURBS) = ntuple(k -> dimsplinespace(S.deg[k],S.kts[k]), nsdim(S))
nbezierpatches(S::NURBS) = ntuple(k -> length(unique(S.kts[k]))-1, nsdim(S))

# bezier decomposition operators
bezierdecompoperator(S::NURBS) = ntuple(k -> bezierdecompoperator(S.deg[k],S.kts[k]), nsdim(S))


# construct IEN array for single NURBS patch S
function ienarray(S::NURBS)
    # initialize
    p     = S.deg
    nsd   = length(p)
    ndofs = dimsplinespace(S)
    ukts, umult, inda, indc = uniqueknots(S.kts)
    nelms = ntuple(k -> length(ukts[k])-1, nsd)
    
    
    # allocate space for IEN array
    m = prod(nelms)    # determine number of elements
    n = prod(S.deg+1)  # number of Bezier degrees of freedom
    
    # construct IEN-array
    IEN = zeros(Int64,n,m)
    for k in 1:m
        subs = ind2sub(nelms,k)
        span = ntuple(i -> indc[i][subs[i]+1]-1, nsd)
        ind  = allcomb(ntuple(i -> collect(span[i]-p[i]:span[i]), nsd)...)
        IEN[:,k] = sub2ind(ndofs, ntuple(i->ind[:,i], nsd)...)
    end
    return IEN 
end

# refine NURBS patch
function refine!(S::NURBS, dir::Int, u::Vector{Float64})

    # initialize
    p     = S.deg
    nsd   = nsdim(S)
    dims  = dimsplinespace(S)
    perm  = [dir:nsd; 1:dir-1]
    
    # compute refinement operator
    C, nkts = knotinsertionoperator(S.deg[dir],S.kts[dir],u)
    
    # compute new control points
    m = size(S.cpts,2)
    n = prod(dims[1:nsd.!=dir])*dimsplinespace(p[dir],nkts)
    P = zeros(n,m); [P[:,k] = ipermutedims(C * permutedims(reshape(S.cpts[:,k],dims),perm), perm)[:] for k in 1:m]
    W = ipermutedims(C * permutedims(reshape(S.wpts,dims),perm), perm)[:]
    
    # update data
    S.kts[dir] = nkts
    S.cpts = P
    S.wpts = W;
end

# global knot insertion in direction dir
function globalrefine!(S::NURBS, dir::Int, k::Int)      
    u = linearspace(unique(S.kts[dir]),k+2)[2:end-1,:][:]
    refine!(S, dir, u)
end

# global knot insertion in direction dir
function globalrefine!(S::NURBS, k::Int)    
    for i in 1:nsdim(S)
        u = linearspace(unique(S.kts[i]),k+2)[2:end-1,:][:]
        refine!(S, i, u)
    end
end

# global knot insertion in direction dir
function degreeelevate!(S::NURBS, dir::Int, t::Int)      

    # initialize
    nsd   = nsdim(S)
    dims  = dimsplinespace(S)
    perm  = [dir:nsd; 1:dir-1]
    
    # compute refinement operator
    C, nkts, p = degreeelevationoperator(S.deg[dir],S.kts[dir],t)
    
    # compute new control points
    m = size(S.cpts,2)
    n = prod(dims[1:nsd.!=dir])*dimsplinespace(p,nkts)
    P = zeros(n,m); [P[:,k] = ipermutedims(C * permutedims(reshape(S.cpts[:,k],dims),perm), perm)[:] for k in 1:m]
    W = ipermutedims(C * permutedims(reshape(S.wpts,dims),perm), perm)[:]
    
    # update data
    S.deg[dir] = p
    S.kts[dir] = nkts
    S.cpts = P;
    S.wpts = W;
    
end
degreeelevate!(S::NURBS,t::Int) = [degreeelevate!(S,dir,t) for dir in 1:nsdim(S)]


############################################ Boundary routines ##############################################

# boundary of a 2d Bezier Patch
function boundary{N}(S::NURBS{N}, dir::Int, k::Int)

    # initialize
    dims = dimsplinespace(S) 
    
    # local numbering
    I = squeeze(slicedim(reshape(1:prod(dims),dims), dir, k==1 ? k : dims[dir]),dir)[:]
    J = 1:length(S.deg).!=dir
    
    # set degree, controlpoints and weights of boundary (dir,k)
    p    = S.deg[J]
    kts  = S.kts[J]
    cpts = S.cpts[I,:]
    wpts = S.wpts[I]
    dec  = S.dec[J]
    
    # set IEN-array if non-empty
    ien = Int64[]
    if isempty(S.ien)==false
        
        # determine boundary dofs
        dims2 = S.deg+1
        I     = squeeze(slicedim(reshape(1:prod(dims2),dims2...), dir, k==1 ? k : dims2[dir]),dir)[:]
        
        # determine boundary elements
        dims3 = nbezierpatches(S) 
        J = squeeze(slicedim(reshape(1:prod(dims3),dims3...), dir, k==1 ? k : dims3[dir]),dir)[:]
        
        ien = S.ien[I,J]
    end
    
    # return boundary k
    return NURBS(p, kts, cpts, wpts, dec, ien, Int64[])
end

function boundary(S::NURBS{1}, dir::Int, k::Int)
    if k==1
        return vec(S.cpts[1,:])
    elseif k==2
        return vec(S.cpts[end,:])
    end
end
boundary(S::NURBS{1}) = [∂(S,1,1) ∂(S,1,2)]
boundary{N}(S::NURBS{N}) = NURBS{N-1}[∂(S,dir,k) for k in 1:2, dir in 1:N][:]

∂(S::NURBS, dir::Int, k::Int) = boundary(S, dir, k)
∂(S::NURBS, k::Int) = boundary(S, k)
∂(S::NURBS) = boundary(S)


end