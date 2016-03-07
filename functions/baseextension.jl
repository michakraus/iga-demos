#########################################################################################################################################
############################# baseextension.jl features some functionality that ought to be in Julia.Base. ##############################
#########################################################################################################################################
#
# Functions:
#       +()             -  addition for tuples
#       allcomb()       -  take all possible combinations of the elements of a set of n vectors
#       ndgrid()        -  construct an n-dimensional cartesian grid using n input vectors with coordinates (ndim-generalization of meshgrid in MatLab)
#       blockmat()      -  construct an abstract matrix of type T with N abstract matrices of type T on its diagonal 
#
#########################################################################################################################################
#
# File part of the Julia IsoGeometric Analysis toolbox
# R.R.Hiemstra
# 14-09-2015
#
#########################################################################################################################################

import Base.getindex, Base.setindex!, Base.size, Base.length, Base.similar, Base.*, Base.+, Base.-


# multiplication by a scalar for tuples
*{N,S<:Real,T<:Real}(a::S,q::NTuple{N,T}) = ntuple(i -> a*q[i], N)
*{N,S<:Real,T<:Real}(q::NTuple{N,T},a::S) = ntuple(i -> a*p[i], N)

# addition for tuples
+{N,T<:Real}(p::NTuple{N,T},q::NTuple{N,T}) = ntuple(i -> p[i]+q[i], N)
-{N,T<:Real}(p::NTuple{N,T},q::NTuple{N,T}) = ntuple(i -> p[i]-q[i], N)
# input:
#    p :: Tuple{N,T}  - ntuple of real numbers of type T
#    q :: Tuple{N,T}  - ntuple of real numbers of type T
# output:
#    r :: Tuple{N,T}  - ntuple of real numbers of type T

# adding a real number to a tuple
+{N,S<:Real,T<:Real}(p::NTuple{N,S},q::T)   = ntuple(i -> p[i]+q, N)
-{N,S<:Real,T<:Real}(p::NTuple{N,S},q::T)   = ntuple(i -> p[i]-q, N)
+{N,S<:Real,T<:Real}(p::NTuple{N,S},q::T)   = ntuple(i -> p[i]+q, N)
-{N,S<:Real,T<:Real}(p::NTuple{N,S},q::T)   = ntuple(i -> p[i]-q, N)
+{N,S<:Real,T<:Real}(p::S,q::NTuple{N,T})   = ntuple(i -> p[i]+q, N)
-{N,S<:Real,T<:Real}(p::S,q::NTuple{N,T})   = ntuple(i -> p[i]-q, N)

# create linearly spaced vector
linearspace(a,b,n) = collect(linspace(a,b,n))

# create linearly spaced Matrix
function linearspace(v::Vector{Float64},k)
    m = length(v)-1
    V = zeros(k,m) 
    [V[:,i] = linearspace(v[i],v[i+1],k) for i in 1:m]
    return V
end

# Take all possible combinations of the elements of a set of n vectors
allcomb(v::AbstractVector) = copy(v)
allcomb{T}(v1::AbstractVector{T}, v2::AbstractVector{T}) = [repmat(v1, length(v2), 1) repmat(reshape(v2,1,length(v2)), size(v1,1), 1)[:]]
allcomb{T}(a::AbstractMatrix{T}, v::AbstractVector{T}) = [repmat(a, length(v), 1) repmat(reshape(v,1,length(v)), size(a,1), 1)[:]]

function allcomb{T}(vs::AbstractVector{T}...)
# input:
#    vs :: AbstractVector{T}...  -  set of n abstract vectors of type T
# output:
#    out :: AbstractMatrix{T}    -  one output of abstract matrix of type T which contains all possible combinations of the 
#                                   elements of the n input vectors
    n = length(vs)
    out = allcomb(vs[1])
    for i in 2:n
        out = allcomb(out,vs[i])
    end
    return out
end

# Construct an n-dimensional cartesian grid using n input vectors with coordinates (ndim-generalization of meshgrid in MatLab)
function ndgrid{T}(vs::AbstractVector{T}...)
# input:
#    vs :: AbstractVector{T}...  -  set of n abstract vectors of type T representing univariate coordinates
# output:
#    out :: AbstractMatrix{T}...    -  set of n abstract matrices of type T representing an n-dimensional grid
    n = length(vs)
    sz = map(length, vs)
    out = ntuple(n, i->Array(T, sz))
    s = 1
    for i=1:n
        a = out[i]::Array
        v = vs[i]
        snext = s*size(a,i)
        ndgrid_fill(a, v, s, snext)
        s = snext
    end
    out
end

# Construct an abstract matrix of type T with N abstract matrices of type T on its diagonal
function blockmat{T}(A::AbstractMatrix{T}...)
# input:
#    A :: AbstractMatrix{T}...  -  set of n abstract matrices of type T
# output:
#    B :: AbstractMatrix{T}     -  abstract matrix of type T with the N abstract matrices on its diagonal

    # determine size of matrices
    n = [0,cumsum([size(A[k],1) for k in 1:length(A)])]
    m = [0,cumsum([size(A[k],2) for k in 1:length(A)])]

    # allocate space for blockmatrix
    B = similar(A[1],n[end],m[end]); B[:] = 0

    # fill blockmatrix
    for k in 1:length(A)
        B[n[k]+1:n[k+1],m[k]+1:m[k+1]] = A[k]
    end
    return B
end