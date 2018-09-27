
#########################################################################################################################################
######### bsplinebaisfuns.jl contains all basic functionality in order to compute univariate B-splines and their derivative. ############
#########################################################################################################################################
#
# Functions:
# 		dimsplinespace()       -  compute the dimension of the univariate splinespace
#       grevillepoints()       -  compute the Greville Absiscae of the univariate splinespace
#		onebasisfuneval()      -  evaluate a single basisfunction at a single or a vector of points
#		bsplinebasiseval()     -  evaluate univariate B-pline basisfunction at a point using the Cox-DeBoor algorithm
#		dersbsplinebasisfuns() -  evaluate univariate B-spline basisfunctions and its derivatives at a single point or a vector of points using the Cox-DeBoor algorithm
#		dersbasisfunsinterpolationmatrix()  -  compute the collocation matrices of the space of univariate B-splines up to its nth order derivatives evaluated at a set of points
#
#########################################################################################################################################
#
# File part of the Julia IsoGeometric Analysis toolbox
# R.R.Hiemstra
# 14-09-2015
#
#########################################################################################################################################

"""
compute the dimension of the spline space defined by degree 'p' and knotsvector 'kts'

## input:
    p   :: Int             - polynomial degree
    kts :: Vector{Float64} - knotvector

## output:
    n   :: Int             - dimension of the splinespace
"""
dimsplinespace(p::Integer,kts::Vector{Float64}) = length(kts)-p-1


"""
compute the Greville Absiscae of the spline space defined by degree 'p' and knotvector 'kts'

## input:
    p   :: Int             - polynomial degree
    kts :: Vector{Float64} - knotvector

## output:
    y   :: Vector{Float64} - vector with Greville Absiscae
"""
grevillepoints(p::Integer,kts::Array{Float64,1}) = [sum(kts[j+1:j+p]) / p for j in 1:length(kts)-p-1]


"""
compute value of a single B-spline basis-function of degree 'p' with local knotvector 'knots' evaluated at point 'u'

## input:
    p     :: Int             - polynomial degree
    knots :: Vector{Float64} - knotvector
    u     :: Float64         - evaluation point

## output:
    b     :: Float64         - basisfunction defined by degree 'p' and local knotvector 'knots' evaluated at point 'u'
"""
function onebasisfuneval(p::Integer,kts::Vector{Float64},u::Float64)
    # initialize
    funs = zeros(Float64,p+1)

    if u<kts[1] || u > kts[p+2]
      funs = 0.0
    elseif sum(u.==kts)==p+1
    	funs = 1.0
    else
      # initialise zeroth degree basisfunction
      for j in 1:p+1
        funs[j] = 0.0
        if (u>=kts[j] && u<kts[j+1])
          funs[j] = 1.0
        end
      end

      # compute triangular table
      for k in 1:p
        if funs[1]==0.0
          saved = 0.0
        else
          saved = ((u-kts[1])*funs[1]) / (kts[k+1]-kts[1])
        end

        for j in 1:p+1-k
          ktsleft  = kts[j+1]
          ktsright = kts[j+k+1]

          if funs[j+1]==0.0
            funs[j] = saved
            saved = 0.0
          else
            temp = funs[j+1] / (ktsright - ktsleft)
            funs[j] = saved + (ktsright - u) * temp
            saved = (u - ktsleft) * temp
          end
        end
      end
    end
  return funs[1]
end


"""
compute value of a single B-spline basis-function of degree 'p' with local knotvector 'knots' evaluated at points 'u'

## input:
    p     :: Int             - polynomial degree
    knots :: Vector{Float64} - knotvector
    u     :: Vector{Float64} - vector with evaluation points

## output:
    b     :: Vector{Float64} - basisfunction defined by degree 'p' and local knotvector 'knots' evaluated at points 'u'
"""
onebasisfuneval(p::Integer,kts::Vector{Float64},u::Vector{Float64}) = Float64[onebasisfuneval(p,kts,u[i]) for i in 1:length(u)]


"""
Compute triangular table of the nonvanishing b-spline basisfunctions up to degree 'p' at the site 'u' and their corresponding support

## input:
    p     :: Int             - polynomial degree
    kts   :: Vector{Float64} - knotvector
    span  :: Int             - index of point u in knotvector
    u     :: Float64         - evaluation point

## output:
    Nu    :: Matrix{Float64} - Triangular table of basis functions - Nu[:,k] denotes the kth order basis
    supp  :: Matrix{Float64} - Triangular table of the support of the basis functions - supp[:,k] denotes the support of the kth order basis
"""
function bsplinebasiseval(p::Integer, knots::Vector{Float64}, span::Integer, u::Float64)
    # initialize
    supp = zeros(Float64,p+1,p+1); supp[1,1] = knots[span+1] - knots[span]
    funs = zeros(Float64,p+1,p+1); funs[1,1] = 1

    # Calculate truncated table of pth basis functions
    for i in 1:p
        for j in 1:i

            # determine support
            supp[j,i+1] = knots[span+j] - knots[span-1+j-i]
            supp[j+1,i+1] = knots[span+j+1] - knots[span+j-i]

            # determine factor of affine map
            alfa = (u - knots[span+j-i]) / supp[j,i]

            # calculate basis function as multi affine map
            funs[j,i+1]   = funs[j,i+1]   + (1-alfa) * funs[j,i]
            funs[j+1,i+1] = funs[j+1,i+1] +    alfa  * funs[j,i]
        end
    end
    return funs, supp
end

"""
Compute triangular table of the nonvanishing b-spline basisfunctions up to degree 'p' at the site 'u' and their corresponding support

## input:
    p     :: Int             - polynomial degree
    knots :: Vector{Float64} - knotvector
    u     :: Float64         - evaluation point

## output:
    Nu    :: Matrix{Float64} - Triangular table of basis functions - Nu[:,k] denotes the kth order basis
    supp  :: Matrix{Float64} - Triangular table of the support of the basis functions - supp[:,k] denotes the support of the
                               kth order basis
"""
bsplinebasiseval(p::Integer,kts::Vector{Float64},u::Float64) = bsplinebasiseval(p,kts,findspan(p,kts,u),u)

"""
compute the nonvanishing B-spline basis-functions of degree 'p' and its '1,2,...,nout-1'th order derivatives at the site 'u'

## input:
    p     :: Int             - polynomial degree
    kts   :: Vector{Float64} - knotvector
    i     :: Int             - index of point u in knotvector
    u     :: Float64         - evaluation point
    nout  :: Int             - number of outputs

## output:
    ders  :: Matrix{Float64} - [B(u); B'(u); B''(u) ... ] where B(u) denote the vector of B-spline basis functions evaluated
                               at u
"""
function dersbsplinebasisfuns(p::Degree, kts::KnotVector, i::Int, u::Float64, nout::Int)
    # initialize
    n = dimsplinespace(p,kts)

    # compute triangular table of basis functions
    ndu   = zeros(Float64,p+1,p+1)
    left  = zeros(Float64,p+1)
    right = zeros(Float64,p+1)
    ndu[1,1] = 1
    for j in 1:p
        left[j+1]  = u - kts[i+1-j]
        right[j+1] = kts[i+j] - u
        saved = 0.0

        for r in 0:j-1
            # lower trinagle
            ndu[j+1,r+1] = right[r+2]+ left[j-r+1]
            temp = ndu[r+1,j] / ndu[j+1,r+1]

            # upper triangle
            ndu[r+1,j+1] = saved + right[r+2] * temp;
            saved = left[j-r+1]*temp;
        end
        ndu[j+1,j+1] = saved;
    end

    # load the basisfunctions
    ders = zeros(Float64,nout,p+1)
    ders[1,:] = ndu[:,end]'

    # This section computes the derivatives
    a = zeros(nout,p+1)
    r = 0
    for r in 0:p
        s1=0; s2=1; a[1,1] = 1

        # loop to compute kth derivative
        for k in 1:nout-1
            d=0
            rk=r-k; pk = p-k
            if (r >= k)
                a[s2+1,1] = a[s1+1,1] / ndu[pk+2,rk+1]
                d = a[s2+1,1] * ndu[rk+1,pk+1]
            end

            if (rk >= -1)
                j1 = 1
            else
                j1 = -rk
            end

            if (r-1 <= pk)
                j2 = k-1
            else
                j2 = p-r
            end

            for j in j1:j2
                a[s2+1,j+1] = (a[s1+1,j+1] - a[s1+1,j]) / ndu[pk+2,rk+j+1]
                d = d + a[s2+1,j+1] * ndu[rk+j+1,pk+1]
            end

            if (r <= pk)
                a[s2+1,k+1] = - a[s1+1,k] / ndu[pk+2,r+1]
                d = d + a[s2+1,k+1] * ndu[r+1,pk+1]
            end
            ders[k+1,r+1] = d
            j = s1; s1 = s2; s2 = j
        end
    end

    # multiply by the correct factors
    r = p
    for k in 1:nout-1
        for j=0:p
            ders[k+1,j+1] = ders[k+1,j+1] * r
        end
        r = r*(p-k);
    end

    return ders
end

"""
compute B-spline basis-functions of degree 'p' and its '1,2,...,nout-1'th order derivatives at the site 'u'

## input:
    p     :: Int             - polynomial degree
    kts   :: Vector{Float64} - knotvector
    u     :: Float64         - evaluation point
    nout  :: Int             - number of outputs

## output:
    ders  :: Matrix{Float64} - [B(u); B'(u); B''(u) ... ] where B(u) denote the vector of B-spline basis functions evaluated
                               at u
"""
dersbsplinebasisfuns(p::Degree,kts::KnotVector,u::Float64,nout::Int) = dersbsplinebasisfuns(p, kts, findspan(p,kts,u), u, nout)

"""
compute B-spline basis-functions of degree 'p' and its '1,2,...,nout-1'th order derivatives at the vector of sites 'u'

## input:
    p     :: Int             - polynomial degree
    kts   :: Vector{Float64} - knotvector
    span  :: Vector{Int}     - vector knot-span-index of points u in knotvector
    u     :: Vector{Float64} - vector with evaluation points
    nout  :: Int             - number of outputs

## output:
    ders  :: Vector{Matrix{Float64}} - [B(u); B'(u); B''(u) ... ] where B(u) denote the vector of B-spline basis functions
                                       evaluated at the vector of sites u
"""
function dersbsplinebasisfuns(p::Int, kts::Vector{Float64}, span::Vector{Int64}, u::Vector{Float64}, nout::Int)
    # dimension Spline space and number of collocation points
    n = dimsplinespace(p,kts)
    m = length(u)

    # allocate space for matrices
    Nu = [zeros(Float64,m,p+1) for k in 1:nout]

    for i in 1:m
        ders = dersbsplinebasisfuns(p, kts, span[i], u[i], nout)
        [Nu[k][i,:] = ders[k,:] for k in 1:nout]
    end

    return Nu
end

"""
compute B-spline basis-functions of degree 'p' and its '1,2,...,nout-1'th order derivatives at the vector of sites 'u'

## input:
    p     :: Int             - polynomial degree
    kts   :: Vector{Float64} - knotvector
    u     :: Vector{Float64} - vector with evaluation points
    nout  :: Int             - number of outputs

## output:
    ders  :: Vector{Matrix{Float64}} - [B(u); B'(u); B''(u) ... ] where B(u) denote the vector of B-spline basis functions
                                       evaluated at the vector of sites u
"""
dersbsplinebasisfuns(p::Int, kts::Vector{Float64}, u::Vector{Float64}, nout::Int) = dersbsplinebasisfuns(p, kts, findspan(p,kts,u), u, nout)

"""
compute the collocation matrix and the matrices with 1st, 2nd,..., derivatives for the space of B-splines at a set of collocation points

## input:
    p     :: Int             - polynomial degree
    kts   :: Vector{Float64} - knotvector
    u     :: Vector{Float64} - vector with evaluation points
    nout  :: Int             - number of outputs

## output:
    ders  :: Vector{Matrix{Float64}} - [B(u); B'(u); B''(u) ... ] where B(u) denotes the vector of B-spline basis functions
                                       evaluated at the vector of sites u
"""
function dersbsplinebasisfunsinterpolationmatrix(p::Int, kts::Vector{Float64}, u::Vector{Float64}, nout::Int)
    # initialize
    n = dimsplinespace(p,kts)
    m = length(u)

    # find span and compute basis functions and derivatives
    span = findspan(p, kts, u)
    Nu = dersbsplinebasisfuns(p, kts, span, u, nout)

    # indices sparse matrix format
    I = repmat(collect(1:m)',p+1,1)
    J = Int64[span[i]-p-1+j for j in 1:p+1, i in 1:m]

    # return matrixes
    return [sparse(I[:],J[:],Nu[k]'[:],m,n) for k in 1:nout]
end
