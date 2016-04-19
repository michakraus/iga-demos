module IGA

# import dependent packages
using FastGaussQuadrature, PyCall, PyPlot

# export types
export Degree, KnotVector, BezierPatch, NURBS, QuadRule

# export functionality from other libs
export ndgrid, allcomb, linearspace, boundary
export buildvector, uniqueknots, globalrefine, findspan
export dimsplinespace, grevillepoints, onebasisfuneval, bsplinebasiseval, dersbsplinebasisfuns, dersbsplinebasisfunsinterpolationmatrix
export knotinsertionoperator!, knotinsertionoperator, degreeelevationoperator, bezierdecompoperator
export dim

# export functionality from nurbs.jl
export nsdim, nbezierpatches, ienarray
export refine!, globalrefine!, degreeelevate!
export boundary, ∂

# export functionality from this module
export dersbezierbasisfuns, assembly, plot, plot_stress_contours, L2residual

include("functions/baseextension.jl")
include("functions/knotvector.jl")
include("functions/bsplinebasisfuns.jl")
include("functions/knotrefinement.jl")
include("functions/nurbs.jl")
include("functions/quadrature.jl")


######################################################################################################
####################################### Shape function routines ######################################
######################################################################################################

# bernstein basis functions on [-1,1] and their derivatives
dersbezierbasisfuns(p::Int, u::Float64, nout::Int) = dersbsplinebasisfuns(p, buildvector([-1.0;1.0],[p+1;p+1]), p+1, u, nout)'

# bernstein basis functions on [-1,1] and their derivatives at a set of points u
function dersbezierbasisfuns(p::Int, u::Vector{Float64}, nout::Int)
    ϕ = zeros(p+1,length(u),nout) 
    for k in 1:length(u)
        ϕ[:,k,:] = dersbezierbasisfuns(p, u[k], nout)
    end
    return ϕ
end

# Compute 1D shape function and derivatives in a single point
function dersbezierbasisfuns(S::BezierPatch{1}, ϕ::Vector{Float64}, dϕ::Vector{Float64})
    
    # initialize
    C = S.dec
    W = S.wpts
    P = S.cpts
    
    # NURBS are B-splines in projective three-space
    Bʷ  = broadcast(*, W, C[1]'* ϕ)
    ∇Bʷ = broadcast(*, W, C[1]'* dϕ)
    
    # compute weightfunction and derivative
    w  = sum(Bʷ)
    ∇w = sum(∇Bʷ)
    
    # compute NURBS basisfunction and derivative
    R  = Bʷ / w
    ∇R = (∇Bʷ * w - Bʷ * ∇w) / w^2
    
    # return quadtities in physical space
    return R, ∇R
end

# Compute 1D shape function and derivatives in an array of points
function dersbezierbasisfuns(S::BezierPatch{1}, ϕ::Matrix{Float64}, dϕ::Matrix{Float64})
    R, ∇R = zeros(ϕ), zeros(dϕ)
    for k in 1:size(R,2)
        R[:,k], ∇R[:,k] = dersbezierbasisfuns(S, ϕ[:,k], dϕ[:,k]) 
    end
    return R, ∇R
end

# Compute 2D shape function, derivatives and determinant in a single point
function dersbezierbasisfuns(S::BezierPatch{2}, ϕ::NTuple{2,Vector{Float64}}, dϕ::NTuple{2,Vector{Float64}})
    
    # initialize
    C = S.dec
    W = S.wpts
    P = S.cpts
    
    # compute B-spline basis functions
    ϕ  = ntuple(k-> C[k]'* ϕ[k], 2)
    dϕ = ntuple(k-> C[k]'*dϕ[k], 2)
    
    # NURBS are B-splines in projective 3-space
    Bʷ  = broadcast(*, W,  kron(ϕ[2], ϕ[1]))
    ∇Bʷ = broadcast(*, W, [kron(ϕ[2], dϕ[1]) kron(dϕ[2], ϕ[1])])
    
    # compute weightfunction
    w  = sum(Bʷ)
    ∇w = sum(∇Bʷ, 1)
    
    # compute NURBS functions and derivatives
    R  = Bʷ / w
    ∇R = (∇Bʷ * w - Bʷ * ∇w) / w^2
    
    # return quadtities in physical space
    return R, ∇R
end

# Compute 2D shape function and derivatives in an array of points
function dersbezierbasisfuns(S::BezierPatch{2}, ϕ::NTuple{2,Matrix{Float64}}, dϕ::NTuple{2,Matrix{Float64}})
    n = prod([size(ϕ[k],1) for k in 1:2]); m = ntuple(k -> size(ϕ[k],2), 2)
    R, ∇R = zeros(n,prod(m)), zeros(n,prod(m),2)
    i = 1
    for k in 1:m[2]
        for l in 1:m[1]
            R[:,i], ∇R[:,i,:] = dersbezierbasisfuns(S, (ϕ[1][:,k], ϕ[2][:,l]), (dϕ[1][:,k], dϕ[2][:,l])) 
            i+=1
        end
    end
    return R, ∇R
end



######################################################################################################
########################################## Assembly routines #########################################
######################################################################################################

typealias Forcing  Function
typealias Traction Function
typealias BoundaryDisplacement Function
typealias Material Function

function assembly(S::NURBS, D::Material, f::Forcing, g::BoundaryDisplacement, h::Traction)

    #########################  Initialize #############################
    p     = S.deg                            # polynomial degree
    nsd   = nsdim(S)                         # number of space dimensions
    ndofs = dimsplinespace(S)                # dimension of the spline space
    dims  = nbezierpatches(S)                # number of elements in each direction
    
    IEN   = S.ien                            # set IEN-array 
    ID    = S.id                             # set ID-array
    
    neq   = maximum(ID)                      # number of equations
    K     = spzeros(neq,neq)                 # allocate space for the stiffness matrix
    F     = zeros(neq)                       # allocate space for the righthandside forcing

    ########################  Precompute data #########################

    # precompute Gauss-Legendre nodes and weights
    qr = ntuple(k -> QuadRule(p[k]+1,"legendre"), nsd)

    # pre-compute univariate Bernstein polynomials and their derivatives
    ϕ = ntuple(k -> dersbezierbasisfuns(p[k], qr[k].points, 2), nsd)

    ########################  Assembly loop over elements #########################
    
    # integration over domain Ω
    for e in 1:size(IEN,2)

        # initialize
        A = IEN[:,e]                       # global node numbers
        P = ID[:,A]                          # global equation numbers
        i,j = ind2sub(dims,e)                # element subscripts

        # construct Bezier Patch
        Sb = BezierPatch(tuple(p...), S.cpts[A,:], S.wpts[A], (S.dec[1][:,:,i], S.dec[2][:,:,j]))

        # compute element contribution of the stiffness matrix and the forcing vector
        Kₑ, Fₑ = element_stiffness_and_forcing(Sb, ϕ, qr, D, f)

        # add contributions of element stiffness matrix and right hand side forcing
        I = P.!=0                  # index internal displacements
        K[P[I],P[I]] += Kₑ[I[:],I[:]]
        F[P[I]]      += Fₑ[I[:]]
        
        # subtract Dirichlet boundary displacements
        for i in 1:2
            for a in 1:size(P,2)
                if P[i,a]==0
                    F[P[I]] -= Kₑ[I,a] * g(vec(S.cpts[A[a],:]))[i]
                end
            end
        end
    end
    
    # integration over boundary ∂Ω
    for dir = 1:2
        I = 1:nsd.!=dir
        for k in 1:2
            C = ∂(S,dir,k) # boundary k, direction dir, of S
            for e in 1:size(C.ien,2)

                # initialize
                A = C.ien[:,e]                       # global node numbers
                P = ID[:,A][:]                       # global equation numbers
                J = P.!=0                            # dofs with Neumann data

                # construct Bezier Patch
                Cb = BezierPatch(tuple(C.deg...), S.cpts[A,:], S.wpts[A], (C.dec[1][:,:,e],))

                # compute element contribution of the stiffness matrix and the forcing vector
                F[P[J]] += neumann_boundary_condition(Cb, ϕ[I], qr[I], h)[J]
                
            end
        end
    end
    
    return K, F
end


# computation of element stiffness and forcing
function element_stiffness_and_forcing(Sb::BezierPatch{2}, ϕ::NTuple{2,Array{Float64,3}}, qr::NTuple{2,QuadRule}, D::Material, f::Forcing)
    
    # initialize
    nen = prod(Sb.deg+1)                 # number of element nodes
    ned = 2                              # number of element dofs per node
    Kₑ = zeros(nen*ned, nen*ned)         # allocate space for element stiffness matrix
    Fₑ = zeros(nen*ned)                  # allocate space for element forcing vector

    # loop over quadrature points
    for j in 1:dim(qr[2])
        for i in 1:dim(qr[1])

            # shape function routine
            R, ∇R = dersbezierbasisfuns(Sb,(ϕ[1][:,i,1],ϕ[2][:,j,1]), (ϕ[1][:,i,2],ϕ[2][:,j,2]))

            # compute geometrical properties
            J = ∇R' * Sb.cpts                                # Jacobian 
            x = vec(R'  * Sb.cpts)                           # quadrature point in physical space
            α = qr[1].weights[i]*qr[2].weights[j]*det(J)     # quadrature weight in physical space

            # compute strain components
            N = (inv(J) * ∇R')'
            B = Matrix{Float64}[[N[a,1]  0.0;
                                 0.0    N[a,2];
                                 N[a,2] N[a,1]] for a in 1:nen]

            # compute element stiffness contribution of quadrature point i,j
            for a in 1:nen
                p = ned*(a-1)+(1:ned)
                for b in a:nen
                    q = ned*(b-1)+(1:ned)
                    Kₑ[p,q] += (B[a]' * D(x) * B[b]) * α
                end
            end
            
            # compute element forcing contribution of quadrature point i,j
            Fₑ += ((f(x) * α) * R')[:]
        end
    end

    # Use symmetry to determine the remainder of the element stiffness matrix
    for a in 1:nen
        p = ned*(a-1)+(1:ned)
        for b in a+1:nen
            q = ned*(b-1)+(1:ned)
            Kₑ[q,p] = Kₑ[p,q]'
        end
    end
    return Kₑ, Fₑ
end

# computation of element tractions
function neumann_boundary_condition(Cb::BezierPatch{1}, ϕ::Tuple{Array{Float64,3}}, qr::Tuple{QuadRule}, h::Function)
    
    # initialize
    nen = prod(Cb.deg+1)                 # number of element nodes
    ned = 2                              # number of element dofs per node
    Fₑ = zeros(nen*ned)                  # allocate space for element forcing vector
    
    # loop over quadrature points
    for i in 1:dim(qr[1])

        # shape function routine
        R, ∇R = dersbezierbasisfuns(Cb, ϕ[1][:,i,1], ϕ[1][:,i,2])

        # compute geometrical properties
        t = vec(∇R' * Cb.cpts)                           # Tangent vector
        x = vec(R'  * Cb.cpts)                           # quadrature point in physical space
        α = qr[1].weights[i]*norm(t)                     # quadrature weight in physical space

        # compute element forcing contribution of quadrature point i,j
        Fₑ += ((h(x) * α) * R')[:]
    end
   
    return Fₑ
end


######################################################################################################
########################################## Postprocessing routines ###################################
######################################################################################################

import PyPlot.plot

# plot the NURBS geometry
function plot(S::NURBS, gran=(3,3))
    # initialize
    if isempty(S.ien)==true
        S.ien = ienarray(S)
    end
    if isempty(S.dec[1])==true
        S.dec = bezierdecompoperator(S)
    end
    dims = ntuple(i -> size(S.dec[i],3), 2)
    u    = ntuple(i -> linearspace(-1.0,1.0,gran[i]), 2)
    ϕ    = ntuple(i -> dersbezierbasisfuns(S.deg[i], u[i], 2), 2)
    
    # loop over elements
    for e in 1:prod(dims)
        i,j = ind2sub(dims,e)          # element subscripts
        A   = S.ien[:,e]                 # global node numbers
        
        # construct Bezier Patch
        Sb  = BezierPatch(tuple(S.deg...), S.cpts[A,:], S.wpts[A], (S.dec[1][:,:,i], S.dec[2][:,:,j]))
        
        # compute 2-dimensional rational basis functions
        R, ∇R = dersbezierbasisfuns(Sb, (ϕ[1][:,:,1], ϕ[2][:,:,1]), (ϕ[1][:,:,2], ϕ[2][:,:,2]))
        
        # plot surface
        X = ntuple(i -> reshape(R' * Sb.cpts[:,i], gran), size(Sb.cpts,2))
        contourf(X[1], X[2], zeros(X[1]))
        
        # plot boundary
        plot(vec(X[1][:,1]),   vec(X[2][:,1]),"k")
        plot(vec(X[1][:,end]), vec(X[2][:,end]),"k")
        plot(vec(X[1][1,:]),   vec(X[2][1,:]),"k")
        plot(vec(X[1][end,:]), vec(X[2][end,:]),"k")
    end
end 


# plot the contourlevels of stress
function plot_stress_contours(G::NURBS{2}, S::NURBS{2}, ii::Int, D::Material, levels::Vector{Float64}, gran=(4,4))

    # initialize
    p   = G.deg
    nen = prod(p+1)
    dims = ntuple(i -> size(S.dec[i],3), 2)

    # sample points
    u = ntuple(k -> linearspace(-1.0,1.0,gran[k]), 2)

    # pre-compute univariate Bernstein polynomials and their derivatives
    ϕ = ntuple(k -> dersbezierbasisfuns(p[k], u[k], 2), 2)

    # loop over elements
    for e in 1:size(S.ien,2)
        i,j = ind2sub(dims,e)          # element subscripts
        A   = G.ien[:,e]               # global node numbers

        # construct Bezier Patch
        Sb = BezierPatch(tuple(p...), G.cpts[A,:], G.wpts[A], (G.dec[1][:,:,i], G.dec[2][:,:,j]))

        # compute the stress at the sample points
        σd = zeros(gran[1],gran[2],3)
        X = zeros(gran); Y = zeros(gran)
        for j in 1:gran[2]
            for i in 1:gran[1]

                # shape function routine
                R, ∇R = dersbezierbasisfuns(Sb,(ϕ[1][:,i,1],ϕ[2][:,j,1]), (ϕ[1][:,i,2],ϕ[2][:,j,2]))

                # compute geometrical properties
                J = ∇R' * Sb.cpts                                # Jacobian 
                x = vec(R'  * Sb.cpts)                           # quadrature point in physical space

                if det(J)==0.0
                    J = eye(2)
                end

                # compute strain basis functions
                N = (inv(J) * ∇R')'
                B = Matrix{Float64}[[N[a,1]  0.0;
                                     0.0    N[a,2];
                                     N[a,2] N[a,1]] for a in 1:nen]

                # compute the stress in point x
                for a in 1:nen
                    save = D(x)*B[a]*S.cpts[A[a],:]'
                    for k in 1:3
                        σd[i,j,k] += save[k]
                    end
                end

                # save the coordinates
                X[i,j] = x[1]; Y[i,j] = x[2]
            end
        end

        cp = contourf(X, Y, σd[:,:,ii],levels)
        if e==1
            colorbar(cp)
        end
    end
end


typealias Geometry NURBS{2}
typealias Solution NURBS{2}
typealias ExactSolution Function

# Compute the L2 residial in stress	
function L2residual(G::Geometry, S::Solution, σ::ExactSolution, D::Material)


    # initialize
    p     = S.deg                            # polynomial degree
    nsd   = nsdim(S)                         # number of space dimensions
    nen   = prod(p+1)
    ndofs = dimsplinespace(S)                # dimension of the spline space
    dims  = nbezierpatches(S)                # number of elements in each direction    
    IEN   = S.ien                            # set IEN-array 

    # precompute Gauss-Legendre nodes and weights
    qr = ntuple(k -> QuadRule(p[k]+1,"legendre"), nsd)

    # pre-compute univariate Bernstein polynomials and their derivatives
    ϕ = ntuple(k -> dersbezierbasisfuns(p[k], qr[k].points, 2), nsd)

    # integration over domain Ω
    L²res = zeros(3)
    L21 = zeros(3)
    L22 = zeros(3)
    for e in 1:size(IEN,2)

        # initialize
        A = IEN[:,e]                         # global node numbers
        i,j = ind2sub(dims,e)                # element subscripts

        # construct Bezier Patch
        Sb = BezierPatch(tuple(p...), G.cpts[A,:], G.wpts[A], (G.dec[1][:,:,i], G.dec[2][:,:,j]))

        # loop over quadrature points
        for j in 1:dim(qr[2])
            for i in 1:dim(qr[1])
        
                # shape function routine
                R, ∇R = dersbezierbasisfuns(Sb,(ϕ[1][:,i,1],ϕ[2][:,j,1]), (ϕ[1][:,i,2],ϕ[2][:,j,2]))

                # compute geometrical properties
                J = ∇R' * Sb.cpts                                # Jacobian 
                x = vec(R'  * Sb.cpts)                           # quadrature point in physical space
                α = qr[1].weights[i]*qr[2].weights[j]*det(J)     # quadrature weight in physical space

                # compute strain components
                N = (inv(J) * ∇R')'

                # compute the stress in point x
                σd = zeros(3)
                for a in 1:nen
                    # strain basis function `a'
                    B = [N[a,1]    0.0;
                            0.0 N[a,2];
                         N[a,2] N[a,1]]   
                    
                    # compute the stress in quadrature point i,j
                    σd += D(x) * B * S.cpts[A[a],:]'
                end
                # compute element forcing contribution of quadrature point i,j
                L²res += (σd - σ(x)).^2 * α
            end
        end
        
    end
    return sqrt(L²res)
end


end