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
ndgrid(v::AbstractVector) = copy(v)

function ndgrid{T}(v1::AbstractVector{T}, v2::AbstractVector{T})
    m, n = length(v1), length(v2)
    v1 = reshape(v1, m, 1)
    v2 = reshape(v2, 1, n)
    (repmat(v1, 1, n), repmat(v2, m, 1))
end

function ndgrid_fill(a, v, s, snext)
    for j = 1:length(a)
        a[j] = v[div(rem(j-1, snext), s)+1]
    end
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

function rref{T}(A::Matrix{T})
    nr, nc = size(A)
    U = copy!(similar(A, T <: Complex ? Complex128 : Float64), A)
    e = eps(norm(U,Inf))
    i = j = 1
    while i <= nr && j <= nc
        (m, mi) = findmax(abs(U[i:nr,j]))
        mi = mi+i - 1
        if m <= e
            U[i:nr,j] = 0
            j += 1
        else
            for k=j:nc
                U[i, k], U[mi, k] = U[mi, k], U[i, k]
            end
            d = U[i,j]
            for k = j:nc
                U[i,k] /= d
            end
            for k = 1:nr
                if k != i
                    d = U[k,j]
                    for l = j:nc
                        U[k,l] -= d*U[i,l]
                    end
                end
            end
            i += 1
            j += 1
        end
    end
    return U
end
rref(x::Number) = one(x)


# type allias for frequently used union
# typealias RangeIndex Union{Int, Vector{Int}, Range{Int}, UnitRange{Int}, AbstractVector{Int}}

# ################# constructor of type SubMatrix ##################
# type SubMatrix{T,I<:RangeIndex,J<:RangeIndex} <: AbstractMatrix{T}
#     data :: Matrix{T}
#     rows :: I
#     cols :: J
#     dims :: (Int,Int) 
# end
# function SubMatrix(data, rows, cols)
#     rows = isa(rows, Colon) ? (1:size(data,1)) : rows
#     cols = isa(cols, Colon) ? (1:size(data,2)) : cols
#     return SubMatrix(data, rows, cols,(length(rows),length(cols)))
# end

# # get methods
# getindex(S::SubMatrix, i::Int, j::Int) = S.data[S.rows[i],S.cols[j]]
# getindex(S::SubMatrix, k::Int) = getindex(S,ind2sub(S.dims, k)...)
# getindex(S::SubMatrix, I::RangeIndex, J::RangeIndex) = S.data[S.rows[I],S.cols[J]]
# getindex(S::SubMatrix, K::RangeIndex) = [getindex(S,ind2sub(S.dims,k)...) for k in K]

# # set methods
# function setindex!(S::SubMatrix, X, I::RangeIndex, J::RangeIndex)
#     S.data[S.rows[I],S.cols[J]] = X
# end
# setindex!(S::SubMatrix, X, K::RangeIndex) = [setindex!(S,X[k],ind2sub(S.dims, k)...) for k in K]

# # SubMatrix functionality
# size(S::SubMatrix) = S.dims
# size(S::SubMatrix,i::Int) = S.dims[i]
# length(S::SubMatrix) = prod(S.dims)

# (*){T,S}(A::SubMatrix{T}, B::DenseArray{S}) = A.data[A.rows,A.cols] * B
# (*){T,S}(A::DenseArray{T}, B::SubMatrix{S}) = A * B.data[B.rows,B.cols]

#  ############# constructor of type SubVector #################
# type SubVector{T,I<:RangeIndex} <: AbstractVector{T}
#     data :: Vector{T}
#     rows :: I
#     dims :: Int 
# end
# function SubVector(data, rows)
#     rows = isa(rows, Colon) ? (1:size(data,1)) : rows
#     return SubVector(data, rows,length(rows))
# end

# # set and get methods
# getindex(V::SubVector, I::RangeIndex) = V.data[V.rows[I]]
# setindex!(V::SubVector, X, I::RangeIndex) = setindex!(V.data,X,V.rows[I])

# # SubMatrix functionality
# size(V::SubVector) = (V.dims,1)
# size(V::SubVector,i::Int) = size(V)[i]
# length(V::SubVector) = V.dims
# similar{T}(A::SubVector{T}) = SubVector(A.data,copy(A.rows),copy(A.dims))