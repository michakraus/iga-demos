#########################################################################################################################################
############### knotrefinement.jl features univariate knotrefinement, degree elevation and Bezier extraction  ###########################
#########################################################################################################################################
#
# Functions:
#       knotinsertionoperator()    -  compute matrix operator from coarse space to h-refined space
#       degreeelevationoperator()  -  compute matrix operator from coarse space to p-refined space
#       bezierdecompoperator()     -  compute bezier extraction operator
#       
#########################################################################################################################################
#
# File part of the Julia IsoGeometric Analysis toolbox
# R.R.Hiemstra
# 14-09-2015
#
#########################################################################################################################################


# insert a single knot into knotvector kts and output the transformation operator from the coarse to the refined space
function knotinsertionoperator!(n::Integer,p::Integer,kts::Vector{Float64},span::Integer,u::Float64)
# input:
#    n     :: Int              -  dimension of the spline space
#    p     :: Int              -  polynomial degree
#    kts   :: Vector{Float64}  -  the initial knotvector
#    span  :: Int              -  knot span index of new knot 'u' in 'kts'
#    u     :: Float64          -  knot to insert
#
# output:
#    C     :: Matrix{Float64}  -  transformation matrix from dofs coarse space to dofs refined space
#                                 (knotvector is updated in-place)

    # allocate space for temporary refinement operator
    C = zeros(Float64,n+1,n)

    # Initialize temporary refinement operator
    inda = 1:span-p; indb = span-p+1:n
    C[inda,inda]   = eye(Float64,length(inda))
    C[indb+1,indb] = eye(Float64,length(indb))

    for j in 1:p
       # determine factor of affine map
       alpha = (u - kts[span-p+j]) / (kts[span+j]-kts[span-p+j])
       if alpha!=0
           C[span-p+j,span-p+j-1] = 1.0 - alpha
           C[span-p+j,span-p+j]   = alpha
       end
    end

    insert!(kts,span+1,u)

    return C
end

# same as above function, however, returns a copy of the updated knotvector
function knotinsertionoperator(n::Integer,p::Integer,kts::Vector{Float64},span::Integer,u::Float64)
    newkts = copy(kts)
    C = knotinsertionoperator!(n, p, newkts, span, u)
    return C, newkts
end


# insert a multiple knots into knotvector kts and output the transformation operator from the coarse to the refined space
function knotinsertionoperator!(p::Int,kts::Vector{Float64},u::Vector{Float64})
# input:
#    p     :: Int              -  polynomial degree
#    kts   :: Vector{Float64}  -  the initial knotvector
#    u     :: Vector{Float64}  -  vector of knots to insert
#
# output:
#    C     :: Matrix{Float64}  -  transformation matrix from dofs coarse space to dofs refined space
#                                 (knotvector is updated in-place)

    # dimension spline space
    n = dimsplinespace(p,kts)
    r = length(u)

    # initialize knot insertion matrix
    C = zeros(Float64,n+r,n)
    C[1:n,1:n] = eye(n,n)

    # loop over refinement vector u
    for i in 1:r

        # Find active span of new knots
        span = findspan(p,kts,u[i])

        # compute next refinement operator
        temp = knotinsertionoperator!(n,p,kts,span,u[i])

        # update knot refinement matrix
        C[1:n+1,:] = temp * C[1:n,:]

        n += 1
    end

    return C
end

# same as above function, however, returns a copy of the updated knotvector
function knotinsertionoperator(p::Int,kts::Vector{Float64},u::Vector{Float64})
    newkts = copy(kts)
    C = knotinsertionoperator!(p, newkts, u)
    return C, newkts
end



function bezierdecompoperator(p::Int,kts::Vector{Float64})
# input:
#    p     :: Int              -  polynomial degree
#    kts   :: Vector{Float64}  -  the initial knotvector
#
# output:
#    C     :: Array{Float64,3}  -  3D array; C[:,:,k] denotes the kth Bezier decomp operator
    
    # Initialization:
    a = p+1
    b = a+1
    nb = 1
    m = length(kts)
    alphas = zeros(p-1)
    C = zeros(p+1,p+1,m-2*p)
    C[:,:,1] = eye(p+1)

    while(b<m)
        C[:,:,nb+1] = eye(p+1)  # Initialize the next extraction operator.
        i = b

        # compute knot multiplicity
        while ((b<m) && (kts[b+1]==kts[b])) b += 1; end
        mult = b-i+1

        if mult < p

            # numerator of alpa
            numer = kts[b] - kts[a]

            # compute and store the alfas
            for j in p:-1:mult+1
                alphas[j-mult] = numer / (kts[a+j]-kts[a])
            end

            # Update the matrix coefficients for r new knots
            r = p - mult
            for j in 1:r
                save = r-j+1
                s = mult+j

                for k in p+1:-1:s+1
                    alpha = alphas[k-s]
                    C[k,:,nb] = alpha * C[k,:,nb] + (1-alpha)*C[k-1,:,nb]
                end

                if b<m
                    C[save,save:j+save,nb+1] = C[p+1,p-j+1:p+1,nb]
                end
            end
        end
        
        # initialize next operator
        if b<m
            a = b
            b += 1
            nb+=1 # finished with current operator
        end
    end

    return C[:,:,1:nb]
end


# compute degree elevation operator from a basis with degree 'p' to a basis with degree 'p+t'
function degreeelevationoperator(p::Integer,knots::Array{Float64,1},t::Integer)
# input:
#    p     :: Int              -  polynomial degree
#    kts   :: Vector{Float64}  -  the initial knotvector
#    t     :: Int              -  raise the degree with 't'
# output:
#    C     :: Matrix{Float64}  -  3D array; C[:,:,k] denotes the kth Bezier decomp operator
#    newknots :: Vector{Float64} - the updated knotvector
#    ph       :: Int             - the updated degree 'ph = p + t'

    # initialize
    (knotsun,multun,inda,indc) = uniqueknots(knots)

    # append knots at begin and end of knotvector
    xa = p+1-multun[1]; xb = p+1-multun[end]
    multun[1] = multun[end] = p+1                   # update knot multiplicity
    knots = buildvector(knotsun,multun)             # new knotvector
    ph = p+t; ph2 = int(floor(ph/2))                # new degree
    newknots = buildvector(knotsun,multun+t)

    # B-spline decomposition into Bezier elements
    numbez = length(knotsun) - 1
    alpha = p - multun; alpha[1] = alpha[end] = 0
    insertknots = buildvector(knotsun,alpha)
    (A,Ubez) = knotinsertionoperator(p,knots,insertknots)

    # Degree elevation of Bezier segments
    bezalfs = zeros(Float64,ph+1,p+1)
    bezalfs[ph+1,p+1] = 1.0
    bezalfs[1,1] = 1.0

    for i in 1:ph2
       inv = 1.0 / binomial(ph,i)
       mpi = minimum([p i])
       for j in maximum([0 i-t]):mpi
           bezalfs[i+1,j+1] = inv * binomial(p,j) * binomial(t,i-j)
       end
    end

    for i in ph2+1:ph-1
       mpi = minimum([p,i])
       for j in maximum([0 i-t]):mpi
           bezalfs[i+1,j+1] = bezalfs[ph-i+1,p-j+1]
       end
    end

    # Compose Bezier degree elevation matrix
    B = spzeros(numbez*ph+1,numbez*p+1)
    for k in 1:numbez
       row = [1:ph+1] + (k-1)*ph
       col = [1:p+1]  + (k-1)*p
       B[row,col] = bezalfs
    end

    # B-spline composition from Bezier elements
    (Z,Ubez) = knotinsertionoperator(ph,newknots,insertknots)
    L = cholfact(Z'*Z); C = spzeros(size(L,2),size(Z,1))
    for i in 1:size(Z,1)
        C[:,i] = (L \ full(Z[i,:]'))
    end

    # Degree elevation matrix
    temp = C * B * A; (row,col) = findnz(temp .> 10e-15)
    D = zeros(size(temp,1),size(temp,2))
    [D[row[i],col[i]] = temp[row[i],col[i]] for i in 1:length(row)]
    index = xa+1:length(newknots)-ph-xb

    return D, newknots, ph
end
