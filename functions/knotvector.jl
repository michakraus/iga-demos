#########################################################################################################################################
######### knotvector.jl contains all basic functionality in order to construct and manipulate B-spline knotvectors ######################
#########################################################################################################################################
#
# Functions:
#       buildvector()           -  construct a knotvector using a vector with increasing values of type T and a vector #                                          specifying the knot multiplicity
#       uniqueknots()           -  decompose a knotvector into an increasing vector of Float64, the knotmultiplicity, and vectors specifying indexing into the knotvector
#       globalrefine()          -  put k new knots inside the non-zero knotspans of the knotvector
#       findspan()              -  determine the knot span index of a point 'u' in the knotvector 'kts'
#       findnonzerospan()       -  find span of all non-zero span elements in knotvector kts
#       getelementknotvectors() -  determine the local knotvectors of each non-zero-span element in the knovector
#
#########################################################################################################################################
#
# File part of the Julia IsoGeometric Analysis toolbox
# R.R.Hiemstra
# 14-09-2015
#
########################################################################################################################################

#########################################################################################################################################

# typealiases
typealias Degree Int
typealias KnotVector Vector{Float64}

# construct a knotvector using a vector with increasing values of type T and a vector specifying the knot multiplicity
function buildvector{T<:Number}(uvals::Array{T,1},mult::Array{Int64,1})
# input:
#    uvals :: Vector{T}   - vector with increasing values of type T
#    mult  :: Vector{Int} - vector specifying the knot multiplicity
# output:
#    kts   :: Vector{T}   - knotvector with knot values of type T  

  # allocate space
  v = zeros(T,sum(mult))

  # set values
  ind = 1
  for i in 1:length(uvals)
      v[ind:ind+mult[i]-1] = fill(uvals[i],mult[i])
      ind = ind + mult[i]
  end

  return v
end

# decompose a knotvector into an increasing vector of Float64, the knotmultiplicity, and vectors specifying indexing into the knotvector
function uniqueknots(knots::Array{Float64,1})
# input:
#    kts :: Vector{Float64}  - knotvector
# output:
#    ukts   :: Vector{Float64} - unique knot values in increasing order
#    mult   :: Vector{Int}     - corresponding knot multiplicity
#    inda   :: Vector{Int}     - indexing such that kts==ukts[inda]
#    indb   :: Vector{Int}     - indexing such that kts[indb]==ukts

    n = length(knots)
    mult = ones(Int64,n)
    uknots = zeros(Float64,n)
    inda = zeros(Int64,n)
    indc = zeros(Int64,n)

    j = 1; uknots[1] = knots[1]; inda[1] = 1; indc[1] = 1
    for i in 2:n
        if knots[i]==knots[i-1]
            mult[j] += 1
        else
            j+=1
            uknots[j] = knots[i]
            indc[j] = i
        end
        inda[i] = j
    end
    return uknots[1:j], mult[1:j], inda, indc[1:j]
end

function uniqueknots(kts)
  n = length(kts)
  temp = ntuple(i -> uniqueknots(kts[i]), n)
  return ukts, umult, inda, indc = ntuple(j -> ntuple(i -> temp[i][j], n), 4)
end


function globalrefine(kts::Vector{Float64},k::Int)
# global refine knotvector: put k new knots in non-zero knotspans
#
# input:
#    kts  :: Vector{Float64}    - knotvector
#    k    :: Int                - integer prescribing the number of new knots to put in current non-zero knot-spans
# output:
#    newkts :: Vector{Float64}  - refined knotvector    

    # get unique knots and multiplicity
    ukts, umult, inda, indc = uniqueknots(kts)
    n = length(kts); m = length(ukts)
    
    # compute new knotvector
    newkts = zeros(n + k*(m-1)); j = 1
    for i in 1:m-1
        newkts[j:j+umult[i]-1] = ukts[i]
        j += umult[i]-1
        newkts[j:j+k+1] = linspace(ukts[i],ukts[i+1],k+2)
        j += k+1
    end
    newkts[j:j+umult[m]-1] = ukts[m]
    return newkts
end

# Determine the knot span index of a point 'u' in th knotvector 'kts'
function findspan(p::Integer, kts::Vector{Float64}, u::Float64)
# Determine the knot span index of a point 'u' in th knotvector 'kts'
#
# input:
#    p      :: Int              - polynomial degree
#    kts    :: Vector{Float64}  - knotvector
#    u      :: Float64          - evaluation point
# output:
#    span   :: Int              - knot span index of a point 'u' in th knotvector 'kts'

    if u<kts[1] || u>kts[end]
        error("Boundserror findspan(p::Integer, kts::Vector{Float64}, u::Float64)")
    end

    if u==kts[end]
        span = length(kts)-1
        while kts[span]==kts[span+1]
            span+=-1
        end
    else
        low = 0
        high = length(kts)
        mid = round(Int64,(low + high)/2)

        while u < kts[mid] || u >= kts[mid + 1]
          if u < kts[mid]
            high = mid
          else
            low = mid
          end
          mid = round(Int64,(low + high)/2)
        end
        span = round(Int64,(low + high)/2)
    end
    return span
end

# Determine the knot span index of a vector of points 'u' in the knotvector 'kts'
function findspan(p::Integer, kts::Vector{Float64}, u::Vector{Float64})
# input:
#    p      :: Int              - polynomial degree
#    kts    :: Vector{Float64}  - knotvector
#    u      :: Vector{Float64}  - evaluation point
# output:
#    span   :: Vector{Int}      - vector with knot span index of each point u[i] in the knotvector 'kts'

    # spline space dimension
    n = length(kts)-p-1

    # allocate space for span
    m = length(u)
    span = zeros(Int64,m)
    span[1] = findspan(p,kts,u[1])
    span[m] = findspan(p,kts,u[m])

    # find span at u[i]
    for i in 2:m-1
        span[i] = span[i-1]
        while u[i] > kts[span[i]+1]
            span[i] += 1
        end
    end

    return span
end

function findfunswithcompactsupp(p::Int,kts::Vector{Float64},a::Float64,b::Float64)
    n = length(kts); k = 1
    while a!=kts[k] && k<n
        k+=1
    end
    a = k
    while b>=kts[k] && k<=n
        k+=1
    end
    b = k-2-p
    return a:b
end


# find span of all non-zero span elements in knotvector 'kts'
function findnonzerospan(p::Int,kts::Vector{Float64})
# input:
#    p      :: Int              - polynomial degree
#    kts    :: Vector{Float64}  - knotvector
#
# output:
#    span   :: Vector{Int}      - vector with knot span index of each non-zero-span elements in the knotvector 'kts'

    m = length(kts)
    span = zeros(Int64,m)
    
    count = 1
    for j in p+1:m-p-1
        if kts[j]!=kts[j+1]
            span[count] = j
            count+=1
        end
    end
    
    return span[1:count-1]
end
findnonzerospan{N}(kts::NTuple{N,Vector{Float64}}) = ntuple(length(kts),(j)->findnonzerospan(kts[j]))

# find the elements in support of basisfunctions
function findsupp(p::Int, kts::Vector{Float64})
    span = findnonzerospan(p,kts)
    n    = dimsplinespace(p,kts)
    m    = length(span)
    
    graph = zeros(Bool,m,n)
    for k in 1:m
        graph[k,span[k]-p:span[k]] = true
    end
    
    return [[1:m][graph[:,k]] for k in 1:n]
end

# get element knot vectors from global knot partition
getelementknotvectors(p::Int, kts::Vector{Float64}, span::Vector{Int}) = [kts[span[i]-p:span[i]+1+p] for i in 1:length(span)]
# input:
#    p      :: Int              - polynomial degree
#    kts    :: Vector{Float64}  - knotvector
#    span   :: Vector{Int}      - span-index of each non-zero knot-span in the knotvector
#
# output:
#    ekts   :: Vector{Vector{Float64}} - vector knotvectors local to each element in knotvector kts


# get element knot vectors from global knot partition
getelementknotvectors(p::Int, kts::Vector{Float64}) = getelementknotvectors(p,kts,findnonzerospan(p,kts))
# input:
#    p      :: Int              - polynomial degree
#    kts    :: Vector{Float64}  - knotvector
#
# output:
#    vkts   :: Vector{Vector{Float64}} - vector knotvectors local to each element in knotvector kts