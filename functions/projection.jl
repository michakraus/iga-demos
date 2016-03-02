# functions associated with B-spline local or global projection
require("knotvector.jl")

# compute Greville absisscae of Bspline of degree p and knotvector kts
grevillepoints(kts::Array{Float64,1}) = sum(kts[2:end-1]) / (length(kts)-2)
grevillepoints(p::Integer,kts::Array{Float64,1}) = [sum(kts[j+1:j+p]) / p for j in 1:length(kts)-p-1]
grevillepoints{N}(p::NTuple{N,Integer},kts::NTuple{N,Vector{Float64}}) = ntuple(N,(i) -> grevillepoints(p[i],kts[i]))

# determine conversion rules from B-splines to Schoenmaker B-splines
converttoedge(p::Integer,kts::Array{Float64,1}) = [(p+1) / (kts[i+p+1] - kts[i]) for i in 1:length(kts)-p-1]
converttonode(p::Integer,kts::Array{Float64,1}) = [(kts[i+p+1] - kts[i])/(p+1) for i in 1:length(kts)-p-1]

# Compute collocation matrix associated with global b-spline interpolation
function globalinterpolationmatrix(p::Integer,knots::Array{Float64,1},u::Array{Float64,1})

    # dimension Spline space
    n = length(knots)-p-1;
    m = length(u);

    # determine knotspan of collocation points
    spanu = findspan(p,knots,u)

    # reset interpolation matrix
    Nu  = zeros(Float64,m,n)

    # Calculate entries of matrix
    for i in 1:m
        # compute triangular table of B-spline basis functions and support
        (funs,supp) = bsplinebasiseval(p,knots,u[i],spanu[i])

        # B-spline basis functions of degree p
        Nu[i,spanu[i]-p:spanu[i]] = funs[:,p+1]
    end

    return sparse(Nu)
end


# Compute histopolation matrix associated with global b-spline interpolation
function globalhistopolationmatrix(p::Integer,knots::Array{Float64,1},u::Array{Float64,1})

    # dimension Spline space
    n = length(knots)-p-1;
    m = length(u) - 1;

    # initialize degree and knotvector of primitive B-spline function
    p = p+1
    knots  = [knots[1], knots, knots[end]]

    # determine knotspan of collocation points primitive b-spline
    spanu = findspan(p,knots,u)

    # reset histopolation matrix
    Mu  = zeros(Float64,m,n)

    # initialize iteration
    (funs,supp) = bsplinebasiseval(p,knots,u[1],spanu[1])
    zfill = zeros(Float64,1,spanu[2]-spanu[1])
    oldfuns = funs[:,p+1]

    # Calculate entries of matrix
    for i in 2:m+1
        # compute triangular table of B-spline basis functions and support
        (funs,supp) = bsplinebasiseval(p,knots,u[i],spanu[i])

        # compute difference of basis functions
        zfill = zeros(Float64,spanu[i]-spanu[i-1])
        newfuns = funs[:,p+1]
        dNu = [oldfuns, zfill] - [zfill, newfuns]

        # compute histopolation matrix coefficients
        Mu[i-1,spanu[i-1]-p] = dNu[1]
        for j in 2:length(dNu)-1
            Mu[i-1,spanu[i-1]-p+j-1] = Mu[i-1,spanu[i-1]-p+j-2] + dNu[j]
        end
        oldfuns = newfuns

    end

    return sparse(Mu)
end

# compute devided difference table
function devideddifference(k1::Integer,x::Array{Float64,1})
    # initialize the divided difference table
    DD = zeros(Float64,k1,k1)
    DD[1,1] = 1.0

    # Calculate rest of table
    for i in 2:k1
        for j in 1:i
            # calculate expanded form of devided differences
            temp = 1.0
            for k in 1:i
                if x[k]!=x[j]
                    temp = temp / (x[j]-x[k])
                end
            end

            # substitute in table
            DD[i,j] = temp
        end
    end
    return DD
end

function bsplinedualfunctional(kts::Vector{Float64},x::Vector{Float64})
# this function computes the dual functional to the B-spline function
# defined by knots kts
# input:
# p - polynomial degree
# x - active interpolation sites
# kts - active knots

    # local approximation order k1 \leq p+1
    k1 = length(x);
    t = sum(kts[2:end-1]) / (length(kts)-2.0)
    γ = zeros(Float64,k1)

    # calculate newton polynomials and derivatives in truncated table
    ϕ = zeros(Float64,k1,k1); ϕ[1,1] = 1
    for p in 1:k1-1
        # Newton polynomial of degree p
        ϕ[p+1,1] = ϕ[p,1] * (t - kts[p+1])

        # ith derivative of degree p Newton polynomial
        for i in 1:p
            ϕ[p+1,i+1] =  ϕ[p,i] +  (t - kts[p+1]) * ϕ[p,i+1]
        end
    end

    # calculate temporary coefficients
    for m in k1:-1:1
        γ[m] = (-1)^(m-1) * (prod(1:m-1)/prod(k1-m+1:k1-1)) * ϕ[k1,k1-m+1]
    end

    # calculate devided difference table
    DD = devideddifference(k1,x)

    # calculate another set of temporary coefficients
    μ = zeros(Float64,k1)
    μ[1] = 1.0
    for m in 2:k1
        dd = DD * [(x-t).^(m-1)]
        temp = γ[m]
        for j in 1:m-1
            temp = temp - dd[j] * μ[j]
        end
        μ[m] = temp / dd[m]
    end

    # fill local interpolation matrix
    return DD' * μ
end

# compute approximate inverse of collocation matrix
function localinterpolationmatrix(p::Integer,kts::Array{Float64,1},y::Array{Float64,1},k1::Int)
# This function calculates the approximate inverse A ~ B^{-1}, where
# B is the consistent B-spline interpolation matrix. The approximation is of
# order p+1
    
    # initialize
    n = dimsplinespace(p,kts)
    A = zeros(n,n); A[1,1] = 1.0; A[n,n] = 1.0

    for i in 2:n-1
        # determining the active interpolating sites
        if i < int(0.5*k1)+1
            index1 = 1:k1
            index2 = 1:k1
        elseif i>n-int(0.5*k1)
            index1 = n-k1+1:n
            index2 = n-k1+1:n
        else
            if iseven(k1)
                index1 = i-int(0.5*k1):i-int(0.5*k1)+k1-1
                index2 = i-int(0.5*k1)+1:i-int(0.5*k1)+k1
            else
                index1 = i-int(0.5*k1)+1:i-int(0.5*k1)+k1
                index2 = i-int(0.5*k1)+1:i-int(0.5*k1)+k1
            end
        end
        
        A[i,index1] = 0.5*bsplinedualfunctional(kts[i:i+p+1],y[index1])'
        A[i,index2] += 0.5*bsplinedualfunctional(kts[i:i+p+1],y[index2])'
    end
    return A
end

# Approximate inverse of the hisopolation matrix
function localhistopolationmatrix(p::Integer,kts::Array{Float64,1},y::Array{Float64,1})
    # # input: KnotVector, degree, function (optional)

    # compute local interpolation matrix
    (A,y) = localinterpolationmatrix(p+1,[kts[1], kts[:], kts[end]],y)
    A = full(A)

    # dimension spline space
    n = size(A,1)-1
    m = size(A,2)-1

    # calculate Histopolation matrix
    B  = zeros(Float64,n,m)
    dA = zeros(Float64,n,m+1)

    for i in 1:n
        dA[i,:] = A[i,:] - A[i+1,:]
    end

    B[:,1] = dA[:,1]
    for j in 2:m
        B[:,j] = B[:,j-1] + dA[:,j]
    end

    return sparse(B), y
end

localhistopolationmatrix(p::Integer,kts::Array{Float64,1}) = localhistopolationmatrix(p,kts,grevillepoints(p,kts))


function ndarraylocalinterp(Data::Array{Float64,4},dir::Array{Int64,1},deg::(Int64,Int64,Int64),kts::(Array{Float64,1},Array{Float64,1},Array{Float64,1}),limiter::Int64)
# This function performs succesive 1D least squares interpolation on all the respective
# directions of the Data array. The order in which the interpolation is performed is
# given the vector 'index'
# input:
    # Data   - 4D array of voxel data
    # index  - vector specifying the order in which least squares projection is done
    # degree - polynomial degree for al four directions
    # parametric direction
    # knotvector - knotvector defining the spline space in each parametric
    # direction

    # initialize datasize
    datasize = [size(Data)...]

    # perform succesive 1D interpolation
    for k in 1:3

        # Direction to be interpolated
        i = dir[1]
        if k==1
            p = deg[i]
            X = kts[i]
        else
            p = deg[i]-1
            X = kts[i][2:end-1]
        end

        # reshape the Data array such that rowspace corresponds to
        # direction id
        Data = reshape(permutedims(Data,dir),datasize[i],prod(datasize[1:length(datasize).!=i]))

        # get local interpolant
        (A,y) = localhistopolationmatrix(p,X,X[ceil(p/2)+1:1:end-floor(p/2)])
        if k==1
            A = spdiagm(converttoedge(p,X)) * A
        end

        # solve 1d local interpolation
        Data = [repmat(Data[1,:],int64(ceil(p/2)),1), Data,  repmat(Data[end,:],int64(floor(p/2)),1)]
        temp = full(A) * Data

        # update size of Data and
        datasize[i] = size(temp,1)

        # Schoenberg variation diminishing spline interpolation for
        # critical regions with wiggles
        if limiter>0
            D = difference(datasize[i])
            index = abs(D'*D * temp) .> limiter
            if round(p/2)==(p/2)
                temp[index] = Data[index]
            else
                temp[index] = 0.5 * (Data[index] + Data[index[[1,1:datasize[i]-1],:]]);
            end
        end
        Data[2:end-1,:] = temp[2:end-1,:]

        # permute data back to original order and reshape
        # to 4D array of new size
        Data = ipermutedims(reshape(Data,tuple(datasize[dir]...)),dir)

        # initialize next loop to interpolate next direction
        dir[1:3] = dir[[2, 3, 1]];

    end

    return Data
end

# do local lsquared projection of a function func
function locallsquaredprojection(p::Degree,kts::KnotVector,func::Function,nquad::Int)
    n = dimsplinespace(p,kts)        # dimension splinespace
    ϕ = zeros(Float64,n)             # allocate space for result
    u,w =  Base.gauss(Float64,nquad) # compute gauss nodes on [-1,1]

    # loop over elements
    for k in p+1:length(kts)-p-1
    
        # non-zero element spans
        if kts[k+1]!=kts[k]
            # determinant of mapping from [-1,1] to spline parameter space
            detJ = 0.5 * (kts[k+1]-kts[k])
            uu = mean([kts[k+1],kts[k]]) + u * detJ
            ww = w * detJ

            # compute contribution of kth quadrature point to Grammian matrix 
            f = zeros(Float64,p+1)
            Nu = bsplinebasisfuns(p, kts, fill(k,nquad), uu)
            for i in 1:nquad
                f += Nu[:,i] * (ww[i] * func(uu[i]))
            end
            P = elementlocallsquaredprojection(p,kts[k-p:k+1+p],Nu,ww)
            
            # Add contribution to approximate inverse 
            ϕ[k-p:k] += P * f
        end
    end
    return ϕ
end
locallsquaredprojection(p::Degree,kts::KnotVector,func::Function) = locallsquaredprojection(p,kts,func,p+1)

# compute element lsquared projection matrix
function elementlocallsquaredprojection(p::Degree, kts::KnotVector, Nu::Matrix{Float64}, ww::Vector{Float64})
    # initialize
    a = zeros(Float64,p+1,p+1)
    b = zeros(Float64,p+1)

    # compute contribution of kth quadrature point to Grammian matrix 
    for i in 1:length(ww)
        temp = Nu[:,i] * ww[i]
        a += temp * Nu[:,i]'
        b += temp
    end

    # compute inverse of Grammian matrix and apply weighing scheme
    b = Float64[b[j] * ((p+1) / (kts[p+j+1] - kts[j])) for j in 1:p+1]

    # return element local l2 inverse
    return broadcast(*,b,inv(a))
end