module Bsplines

export Degree, KnotVector

# export types
export ndgrid, allcomb, linearspace
export buildvector, uniqueknots, globalrefine, findspan
export dimsplinespace, grevillepoints, onebasisfuneval, bsplinebasiseval, dersbsplinebasisfuns, dersbsplinebasisfunsinterpolationmatrix
export knotinsertionoperator!, knotinsertionoperator, degreeelevationoperator, bezierdecompoperator
export gaussquadrule, optimalquadrule

include("functions/baseextension.jl")
include("functions/knotvector.jl")
include("functions/bsplinebasisfuns.jl")
include("functions/knotrefinement.jl")
include("functions/quadrature.jl")

end