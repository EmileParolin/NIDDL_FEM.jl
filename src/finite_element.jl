# TODO: 
# - rewrite everything without the sml2big mappings, where sml and big are two
# simplices such that dim(sml) < dim(big). to save memory on large meshes
# - can we use only the mappins big2vtx, with big = simplex of the physical
# domains (i.e. no edg2vtx if closed surface, no edg2vtx for volumic internal
# edges, no tri2vtx for volumic internal faces...) ?
# - this means that the assembly loop must be rewritten, with the main loop
# being on the largest simplex on which the forms are defined (integration
# performed)
# - and all the definitions of the finite elements as well...

# Note: the fact that the matrices are defined on quadrature nodes, leaving the
# computation of support intersection and evaluation of the quadrature to be
# done by sparse matrix product later DOES use more memory. however, this
# choice allows to define finite element matrices much more easily and also a
# lot faster. In addition, it is easier to evaluate solution on a fine grid.

##############################
# Finite Element definitions #
##############################

abstract type BFS end

""" Lagrange Finite Element P0 on edge. """
struct P0edg <: BFS
    ibfs::Array{Int64,1} # DOF indices of basis functions supported on edge (Nbf=1)
    coef::Array{Float64,1} # coef of basis functions
end

function P0edg(m::Mesh, ielt::Int64, params...)
    s = m.vtx[:,m.edg2vtx[:,ielt]] # vertices (dim=3, 2)
    l = norm(s[:,2]-s[:,1]) # length of edge
    P0edg([ielt], [1/l])
end

@inline (bfs::P0edg)(x::Array{Float64,1}) = bfs.ibfs, bfs.coef

""" Lagrange Finite Element P1 on edge. """
struct P1edg <: BFS
    ibfs::Array{Int64,1} # DOF indices of basis functions supported on edge (Nbf=3)
    s::Array{Float64,2} # vertices (dim=3, Nbf=2)
    e::Array{Float64,2}  # normalised edge (dim=3, 2)
end

function P1edg(m::Mesh, ielt::Int64, params...)
    ibfs = m.edg2nod[:,ielt] # indices of basis functions = indices of vertices
    s = m.vtx[:,m.edg2vtx[:,ielt]] # vertices (dim=3, 2)
    e = Array{Float64,2}(undef,3,2)
    for ibf in 1:2 # loop on basis functions
        # Tangent vector
        tau = s[:,[2,1][ibf]] - s[:,ibf]
        # Normalised edge
        e[:,ibf] = tau / norm(tau)^2
    end
    # creating new dof
    P1edg(ibfs, s, e)
end

@inline function (bfs::P1edg)(x::Array{Float64,1})
    val = Array{Float64,2}(undef,1,2)
    for ibf in 1:2 # loop on basis functions
        val[1,ibf] = dot(bfs.s[:,[2,1][ibf]] - x, bfs.e[:,ibf])
    end
    bfs.ibfs, val
end

struct gradP1edg <: BFS
    bfs::P1edg
end

function gradP1edg(m::Mesh, ielt::Int64, params...)
    gradP1edg(P1edg(m, ielt))
end

@inline function (bfs::gradP1edg)(x::Array{Float64,1})
    val = Array{Float64,2}(undef,3,2)
    for ibf in 1:2 # loop on basis functions
        val[:,ibf] = -bfs.bfs.e[:,ibf]
    end
    bfs.bfs.ibfs, val
end

""" Lagrange Finite Element P0 on triangle. """
struct P0tri <: BFS
    ibfs::Array{Int64,1} # DOF indices of basis functions supported on triangle (Nbf=1)
    coef::Array{Float64,1} # coef of basis functions
end

function P0tri(m::Mesh, ielt::Int64, params...)
    s = m.vtx[:,m.tri2vtx[:,ielt]] # vertices (dim=3, 3)
    cp = cross(s[:,2]-s[:,1], s[:,3]-s[:,1]) # cross product between 2 edges
    area = 0.5 * norm(cp) # triangle area
    P0tri([ielt], [1/area])
end

@inline (bfs::P0tri)(x::Array{Float64,1}) = bfs.ibfs, bfs.coef

""" Lagrange Finite Element P1 on triangle. """
struct P1tri <: BFS
    ibfs::Array{Int64,1} # DOF indices of basis functions supported on triangle (Nbf=3)
    s::Array{Float64,2} # vertices (dim=3, Nbf=3)
    nu::Array{Float64,2} # opposite edge normal vector (dim=3, Nbf=3)
    h::Array{Float64,1} # altitude (Nbf=3)
    nrm::Array{Float64,1} # normal vector (dim=3) -> curl phi
end

function P1tri(m::Mesh, ielt::Int64, params...)
    ibfs = m.tri2nod[:,ielt] # indices of basis functions = indices of vertices
    s = m.vtx[:,m.tri2vtx[:,ielt]] # vertices (dim=3, 3)
    cp = cross(s[:,2]-s[:,1], s[:,3]-s[:,1]) # cross product between 2 edges
    nrm = cp / norm(cp) # normal vector
    nu = Array{Float64,2}(undef,3,3)
    h = Array{Float64,1}(undef,3)
    for ibf in 1:3 # loop on basis functions
        # Opposite edge tangent vector
        tau = s[:,[3,1,2][ibf]] - s[:,[2,3,1][ibf]]
        tau /= norm(tau)
        # Opposite edge normal vector, pointing inside triangle
        nu[:,ibf] = cross(nrm, tau)
        # Altitude
        h[ibf] = dot(s[:,[3,1,2][ibf]] - s[:,ibf], nu[:,ibf])
    end
    # creating new dof
    P1tri(ibfs, s, nu, h, nrm)
end

@inline function (bfs::P1tri)(x::Array{Float64,1})
    val = Array{Float64,2}(undef,1,3)
    for ibf in 1:3 # loop on basis functions
        val[1,ibf] = 1 .- dot((x .- bfs.s[:,ibf]), bfs.nu[:,ibf]) ./ bfs.h[ibf]
    end
    bfs.ibfs, val
end

struct gradP1tri <: BFS
    bfs::P1tri
end

function gradP1tri(m::Mesh, ielt::Int64, params...)
    gradP1tri(P1tri(m, ielt))
end

@inline function (bfs::gradP1tri)(x::Array{Float64,1})
    val = Array{Float64,2}(undef,3,3)
    for ibf in 1:3 # loop on basis functions
        val[:,ibf] = -bfs.bfs.nu[:,ibf] / bfs.bfs.h[ibf]
    end
    bfs.bfs.ibfs, val
end

"""
    curl = grad × n
    grad = n × curl
"""
struct curlP1tri <: BFS
    bfs::gradP1tri
end

function curlP1tri(m::Mesh, ielt::Int64, params...)
    curlP1tri(gradP1tri(m, ielt))
end

@inline function (bfs::curlP1tri)(x::Array{Float64,1})
    ibfs, valgradP1 = bfs.bfs(x)
    val = Array{Float64,2}(undef,3,3)
    for ibf in 1:3 # loop on basis functions
        val[:,ibf] = cross(valgradP1[:,ibf], bfs.bfs.bfs.nrm)
    end
    ibfs, val
end

""" Lagrange Finite Element P0 on tetrahedron. """
struct P0tet <: BFS
    ibfs::Array{Int64,1} # DOF indices of basis functions supported on tet (Nbf=1)
    coef::Array{Float64,1} # coef of basis functions
end

function P0tet(m::Mesh, ielt::Int64, params...)
    s = m.vtx[:,m.tet2vtx[:,ielt]] # vertices (dim=3, 4)
    cp = cross(s[:,2]-s[:,1], s[:,3]-s[:,1]) # cross product between 2 edges
    sp = dot(s[:,4]-s[:,1], cp) # dot product with other (non colinear) edge
    vol = abs(sp) / 6 # tet volume
    P0tet([ielt], [1/vol])
end

@inline (bfs::P0tet)(x::Array{Float64,1}) = bfs.ibfs, bfs.coef

""" Lagrange Finite Element P1 on tetrahedron. """
struct P1tet <: BFS
    ibfs::Array{Int64,1} # DOF indices of basis functions supported on tet (Nbf=4)
    s::Array{Float64,2} # vertices (dim=3, Nbf=4)
    nu::Array{Float64,2} # opposite edge normal vector (dim=3, Nbf=4)
    h::Array{Float64,1} # altitude (Nbf=4)
end

function P1tet(m::Mesh, ielt::Int64, params...)
    ibfs = m.tet2nod[:,ielt] # indices of basis functions = indices of vertices
    s = m.vtx[:,m.tet2vtx[:,ielt]] # vertices (dim=3, 4)
    nu = Array{Float64,2}(undef,3,4)
    h = Array{Float64,1}(undef,4)
    for ibf in 1:4 # loop on basis functions
        # Getting 3 edges to form a basis
        e1 = s[:,[4,3,4,2][ibf]] - s[:,[2,1,1,1][ibf]]
        e2 = s[:,[3,4,2,3][ibf]] - s[:,[2,1,1,1][ibf]]
        e3 = s[:,[2,1,1,1][ibf]] - s[:,ibf]
        # Opposite face normal vector
        e12 = cross(e1,e2)
        nu[:,ibf] = e12 / norm(e12)
        # Altitude
        h[ibf] = dot(nu[:,ibf],e3)
    end
    # creating new dof
    P1tet(ibfs, s, nu, h)
end

@inline function (bfs::P1tet)(x::Array{Float64,1})
    val = Array{Float64,2}(undef,1,4)
    for ibf in 1:4 # loop on basis functions
        val[1,ibf] = 1 .- dot((bfs.s[:,ibf] - x), bfs.nu[:,ibf]) / bfs.h[ibf]
    end
    bfs.ibfs, val
end

struct gradP1tet <: BFS
    bfs::P1tet
end

function gradP1tet(m::Mesh, ielt::Int64, params...)
    gradP1tet(P1tet(m, ielt))
end

@inline function (bfs::gradP1tet)(x::Array{Float64,1})
    val = Array{Float64,2}(undef,3,4)
    for ibf in 1:4 # loop on basis functions
        val[:,ibf] = bfs.bfs.nu[:,ibf] / bfs.bfs.h[ibf]
    end
    bfs.bfs.ibfs, val
end

""" Raviart-Thomas Finite Element on triangle. """
struct RT <: BFS
    ibfs::Array{Int64,1} # DOF indices of basis functions supported on triangle (Nbf=3)
    coef::Array{Float64,1} # (signed) coefficient in RT basis function (Nbf=3)
    ov::Array{Float64,2} # opposite vertices (dim=3,Nbf=3)
    nrm::Array{Float64,1} # normal vectors -> nx phi, curl phi (dim=3)
end

function RT(m::Mesh, ielt::Int64, params...)
    ibfs = m.tri2edg[:,ielt]
    tri = m.tri2vtx[:,ielt] # triangle vertex indices (3)
    s = m.vtx[:,tri] # vertices (dim=3, 3)
    ov = s[:,[3,1,2]] # opposite vertices
    cp = cross(s[:,2]-s[:,1], s[:,3]-s[:,1]) # cross product between 2 edges
    nrm = cp / norm(cp) # normal vector
    area = 0.5*norm(cp) # triangle area
    coef = sign.(tri[[2,3,1]]-tri) * 0.5/area # /!\ consistently signed coefficient
    RT(ibfs, coef, ov, nrm)
end

@inline function (bfs::RT)(x::Array{Float64,1})
    val = Array{Float64,2}(undef,3,3)
    for ibf in 1:3 # loop on basis functions
        val[:,ibf] = bfs.coef[ibf] * (x-bfs.ov[:,ibf])
    end
    bfs.ibfs, val
end

struct nxRT <: BFS
    bfs::RT
end

function nxRT(m::Mesh, ielt::Int64, params...)
    nxRT(RT(m, ielt))
end

@inline function (bfs::nxRT)(x::Array{Float64,1})
    ibfs, valRT = bfs.bfs(x)
    val = Array{Float64,2}(undef,3,3)
    for ibf in 1:3 # loop on basis functions
        val[:,ibf] = cross(bfs.bfs.nrm, valRT[:,ibf])
    end
    ibfs, val
end

struct divRT <: BFS
    bfs::RT
end

function divRT(m::Mesh, ielt::Int64, params...)
    divRT(RT(m, ielt))
end

@inline function(bfs::divRT)(x::Array{Float64,1})
    bfs.bfs.ibfs, reshape(2 * bfs.bfs.coef, (1,3))
end

"""
If NED are along the edges in the same direction given by the node numbering,
then n×NED gives a vector pointing inside the triangle, so opposite RT (if sign
edge = +1).

            RT  = -n × NED =  NED × n
            NED =  n × RT  = -RT  × n

            curl NED = div (NED × n) = div RT
"""
struct NEDtri <: BFS
    bfs::nxRT
end

function NEDtri(m::Mesh, ielt::Int64, params...)
    NEDtri(nxRT(m, ielt))
end

@inline function (bfs::NEDtri)(x::Array{Float64,1})
    return bfs.bfs(x)
end

struct curlNEDtri <: BFS
    bfs::divRT
end

function curlNEDtri(m::Mesh, ielt::Int64, params...)
    curlNEDtri(divRT(m, ielt))
end

@inline function (bfs::curlNEDtri)(x::Array{Float64,1})
    return bfs.bfs(x)
end

""" Nedelec Finite Element on tetrahedron. """
struct NEDtet <: BFS
    ibfs::Array{Int64,1} # DOF indices of basis functions supported on tetra (Nbf=6)
    a::Array{Float64,2} # a (dim=3,Nbf=6) ϕ = a + b × x
    b::Array{Float64,2} # b (dim=3,Nbf=6)
end

function NEDtet(m::Mesh, ielt::Int64, params...)
    ibfs = m.tet2edg[:,ielt]
    tet = @view m.tet2vtx[:,ielt] # tet vertex indices (4)
    s = @view m.vtx[:,tet] # vertices (dim=3, 4)
    a = Array{Float64,2}(undef,3,6)
    b = Array{Float64,2}(undef,3,6)
    for ibf in 1:6 # loop on basis functions
        # vertex indices in local basis
        i1 = [1,1,1,2,2,3][ibf] # current edge first vertex
        i2 = [2,3,4,3,4,4][ibf] # current edge second vertex
        i3 = [3,2,2,1,1,1][ibf] # opposite edge first vertex
        i4 = [4,4,3,4,3,2][ibf] # opposite edge second vertex
        # Current edge
        cv = s[:,i1] # first vertex (we don't care about which one)
        ce = sign(tet[i2]-tet[i1])*(s[:,i2] - cv) # /!\ consistently signed tangent vector
        # Opposite edge
        ov = s[:,i3] # first vertex (we don't care about which one)
        oe = s[:,i4] - ov # tangent vector (we don't care about the sign)
        # Coefficient
        coef = 1/dot(cross(oe, (cv-ov)), ce)
        a[:,ibf] = coef * cross(oe,-ov)
        b[:,ibf] = coef * oe
    end
    NEDtet(ibfs, a, b)
end

@inline function(bfs::NEDtet)(x::Array{Float64,1})
    val = Array{Float64,2}(undef,3,6)
    for ibf in 1:6 # loop on basis functions
        val[:,ibf] = bfs.a[:,ibf] + cross(bfs.b[:,ibf], x)
    end
    bfs.ibfs, val
end

struct curlNEDtet <: BFS
    bfs::NEDtet
end

function curlNEDtet(m::Mesh, ielt::Int64, params...)
    curlNEDtet(NEDtet(m, ielt))
end

@inline function(bfs::curlNEDtet)(x::Array{Float64,1})
    val = Array{Float64,2}(undef,3,6)
    for ibf in 1:6 # loop on basis functions
        val[:,ibf] = 2 * bfs.bfs.b[:,ibf]
    end
    bfs.bfs.ibfs, val
end

"""
Scalar functions supported on an edge
"""
struct ScaEdgFunc <: BFS
    ielt::Number # Edge index
    nrm::Array{Float64,1} # normal vectors
    func::Function # function to evaluate
end

function ScaEdgFunc(m::Mesh, ielt::Int64, f::Function)
    edg = m.edg2vtx[:,ielt] # edge vertex indices (2)
    s = m.vtx[:,edg] # vertices (dim=3, 2)
    nrm = cross(s[:,2] - s[:,1], [0, 0, 1]) # normal vector
    nrm /= norm(nrm)
    ScaEdgFunc(ielt, nrm, f)
end

@inline function(bfs::ScaEdgFunc)(x::Array{Float64,1})
    1, [bfs.func(x,bfs.ielt,bfs.nrm)]
end

"""
Scalar functions supported on a triangle
"""
struct ScaTriFunc <: BFS
    ielt::Number # Triangle index
    nrm::Array{Float64,1} # normal vectors -> nx phi, curl phi (dim=3)
    func::Function # function to evaluate
end

function ScaTriFunc(m::Mesh, ielt::Int64, f::Function)
    tri = m.tri2vtx[:,ielt] # triangle vertex indices (3)
    s = m.vtx[:,tri] # vertices (dim=3, 3)
    cp = cross(s[:,2]-s[:,1], s[:,3]-s[:,1]) # cross product between 2 edges
    nrm = cp / norm(cp) # normal vector
    ScaTriFunc(ielt, nrm, f)
end

@inline function(bfs::ScaTriFunc)(x::Array{Float64,1})
    1, [bfs.func(x,bfs.ielt,bfs.nrm)]
end

"""
Vectorial functions supported on a triangle
"""
struct VecTriFunc <: BFS
    ielt::Number # Triangle index
    nrm::Array{Float64,1} # normal vectors -> nx phi, curl phi (dim=3)
    func::Function # function to evaluate
end

function VecTriFunc(m::Mesh, ielt::Int64, f::Function)
    tri = m.tri2vtx[:,ielt] # triangle vertex indices (3)
    s = m.vtx[:,tri] # vertices (dim=3, 3)
    cp = cross(s[:,2]-s[:,1], s[:,3]-s[:,1]) # cross product between 2 edges
    nrm = cp / norm(cp) # normal vector
    VecTriFunc(ielt, nrm, f)
end

@inline function(bfs::VecTriFunc)(x::Array{Float64,1})
    1, reshape(bfs.func(x,bfs.ielt,bfs.nrm),(3,1))
end

#####################
# Assembly routines #
#####################

function nbdof(fe::DataType)
    if fe in [P1edg, gradP1edg]
        return 2 # vertices in edge 
    elseif fe in [P1tri, gradP1tri, curlP1tri]
        return 3 # vertices in triangle
    elseif fe in [P1tet, gradP1tet]
        return 4 # vertices in tetrahedron
    elseif fe in [RT, nxRT, divRT, NEDtri, curlNEDtri]
        return 3 # edges in triangle
    elseif fe in [NEDtet, curlNEDtet]
        return 6 # edges in tetrahedron
    elseif fe in [P0edg, P0tri, P0tet]
        return 1 # edge, triangle or tetrahedron
    elseif fe in [ScaEdgFunc, ScaTriFunc, VecTriFunc,]
        return 1
    else
        error("Please modify nbdof function for DataType $fe")
    end
end

function dimdof(fe::DataType)
    if fe in [P1edg, gradP1edg, P1tri, gradP1tri, curlP1tri, P1tet, gradP1tet]
        return 0 # vertices in edge / triangle / tetrahedron
    elseif fe in [P0edg, RT, nxRT, divRT, NEDtri, curlNEDtri, NEDtet, curlNEDtet]
        return 1 # edges in triangle / tetrahedron
    elseif fe in [P0tri,]
        return 2 # triangle
    elseif fe in [P0tet,]
        return 3 # tetrahedron
    else
        error("Please modify dimdof function for DataType $fe")
    end
end

function fem(fe::DataType, m::Mesh)
    if fe in [P0edg, P1edg, gradP1edg, P0tri, P1tri, gradP1tri, curlP1tri, RT,
              nxRT, divRT, NEDtri, curlNEDtri, P0tet, P1tet, gradP1tet, NEDtet,
              curlNEDtet, ScaEdgFunc, ScaTriFunc, VecTriFunc]
        return m
    else
        error("Please modify fem function for DataType $fe")
    end
end

function feelts(fe::DataType)
    if fe in [P0edg, P1edg, gradP1edg, ScaEdgFunc]
        return 1 # Associated to edges
    elseif fe in [P0tri, P1tri, gradP1tri, curlP1tri, RT, nxRT, divRT, NEDtri,
              curlNEDtri, ScaTriFunc, VecTriFunc]
        return 2 # Associated to triangles
    elseif fe in [P0tet, P1tet, gradP1tet, NEDtet, curlNEDtet]
        return 3 # Associated to tetrahedrons
    else
        error("Please modify feelts function for DataType $fe")
    end
end

function femaxelts(fe::DataType, m::Mesh)
    if fe in [P0edg, P1edg, gradP1edg, ScaEdgFunc]
        return size(m.edg2vtx,2)
    elseif fe in [P0tri, P1tri, gradP1tri, curlP1tri, RT, nxRT, divRT, NEDtri,
              curlNEDtri, ScaTriFunc, VecTriFunc]
        return size(m.tri2vtx,2)
    elseif fe in [P0tet, P1tet, gradP1tet, NEDtet, curlNEDtet]
        return size(m.tet2vtx,2)
    else
        error("Please modify femaxelts function for DataType $fe")
    end
end

function femaxdofs(fe::DataType, m::Mesh)
    if fe in [P1edg, gradP1edg, P1tri, gradP1tri, curlP1tri, P1tet, gradP1tet]
        return size(m.nod2vtx,2)
    elseif fe in [P0edg, RT, nxRT, divRT, NEDtri, curlNEDtri, NEDtet, curlNEDtet]
        return size(m.edg2vtx,2)
    elseif fe in [P0tri]
        return size(m.tri2vtx,2)
    elseif fe in [P0tet]
        return size(m.tet2vtx,2)
    elseif fe in [ScaEdgFunc, ScaTriFunc, VecTriFunc,]
        return 1
    else
        error("Please modify femaxdofs function for DataType $fe")
    end
end

function dim(fe::DataType)
    if fe in [P0edg, P1edg, gradP1edg, P0tri, P1tri, P0tet, P1tet, divRT,
              curlNEDtri, ScaEdgFunc, ScaTriFunc]
        return 1 # Scalar field
    elseif fe in [gradP1tri, curlP1tri, gradP1tet, RT, nxRT, NEDtri,
                  NEDtet, curlNEDtet, VecTriFunc]
        return 3 # Vectorial field
    else
        error("Please modify dim function for DataType $fe")
    end
end

function valuetype(fe::DataType)
    if fe in [P0edg, P1edg, gradP1edg, P0tri, P1tri, gradP1tri, curlP1tri,
              P0tet, P1tet, gradP1tet, RT, nxRT, divRT, NEDtri, curlNEDtri,
              NEDtet, curlNEDtet]
        return Float64
    elseif fe in [ScaEdgFunc, ScaTriFunc, VecTriFunc,]
        return Complex{Float64}
    else
        error("Please modify valuetype function for DataType $fe")
    end
end


# TODO :
# - using metaprogramming unroll loop on fe to compute: phi, gradphi, etc tous
# les elements finis utiles ?
# - can we unroll the two inner loops? does it provided faster evaluation?
# - maybe the api should be something like: @assemble [RT, nx(RT)] m x
# with x being a list of Gauss points defined on the reference simplex
# that might default to one Gauss point at the barycenter ?

function assemble(fe::DataType, m::Mesh, Ω::Domain, q::Quadrature, params...)
    nbvalues = number_of_elements(m,Ω,dim(Ω)) * nbdof(fe) * size(q.x,1)
    iis = Array{Int64,1}(undef,nbvalues) # non-zeros element row indexes
    ijs = Array{Int64,1}(undef,nbvalues) # non-zeros element column indexes
    vas = [Array{valuetype(fe),1}(undef,nbvalues) for i in 1:dim(fe)] # non-zeros element values
    i = 1
    # Loop on tags
    for tag in tags(Ω)
        # Loop on elements on which the integration is performed
        for ielt in tag2range(m,tag,feelts(fe))
            bfs = fe(m, ielt, params...) # create basis functions
            # Gauss points real coordinates in current element
            xt = ref2real(fem(fe, m), ielt, q.x) # (dim, nbgausspts)
            for ix in 1:size(q.x,1) # loop on Gauss points in element
                rg = i:i+nbdof(fe)-1
                # evaluate each basis functions at Gauss point
                ids, vals = bfs(xt[:,ix])
                iis[rg] .= size(q.x,1)*(ielt-1)+ix # line index
                ijs[rg] .= ids # column indices
                for id in 1:dim(fe) # loop on dimension (scalar or vector)
                    vas[id][rg] = vals[id,:]
                end
                i += nbdof(fe)
            end
        end
    end
    N, M = size(q.x,1)*femaxelts(fe,m), femaxdofs(fe,m)
    mats = Vector{SparseArrays.SparseMatrixCSC{valuetype(fe),Int64}}(undef,dim(fe))
    for id in 1:dim(fe)
        mats[id] = sparse(iis,ijs,vas[id],N,M)
    end
    mats
end

"""
    femdot(vs,w,us)

...
# Arguments
- `vs`: test function
- `w`: quadrature matrix
- `us`: trial function
...

"""
function femdot(vs, w, us)
    # Sanity checks
    @assert length(vs) == length(us) "test and trial matrices dimension are
        not compatible: test=$(length(vs)) != trial=$(length(us))"
    @assert size(w,1) == size(w,2) "quadrature matrix must be square
        $(size(w,1)) != $(size(w,2))"
    # initialisation
    mats = []
    for (u,v) in zip(us, vs)
        # Sanity checks
        @assert size(u,1) == size(w,2) "trial matrix $(size(u,1)) not
            compatible with quadrature matrix $(size(w,2))"
        @assert size(v,1) == size(w,1) "test matrix $(size(v,1)) not
            compatible with quadrature matrix $(size(w,1))"
        # Compute contribution
        push!(mats, transpose(v) * w * u)
    end
    # scalar product
    sum(mats)
end


function restriction(m::Mesh, Ω::Domain, d)
    # get maximum number of DOFs
    M = number_of_elements(m, d)
    # Safety short-cut
    if length(tags(Ω,d:dim(Ω))) == 0 return spzeros(Bool,0,M) end
    # column indices with non-zero values
    j = element_indices(m,Ω,d)
    # lines indices
    i = collect(1:length(j))
    # sparse matrix
    R = sparse(i,j,ones(Bool,length(i)),length(i),M)
    return R
end


"""
    d is the dimension of the simplex on which the integration is performed
    for instance, it can be computed using feelts()
"""
function restriction(m::Mesh, Ω::Domain, q::Quadrature, d)
    # get maximum number of DOFs
    M = number_of_elements(m, d)
    # Safety short-cut
    if length(tags(Ω,d:dim(Ω))) == 0 return spzeros(Bool,0,M) end
    # Number of quadrature nodes per element
    nx = size(q.x,1)
    # column indices with non-zero values: j = nx * (elts-1) + (1:nx)
    elts = element_indices(m,Ω,d)
    j = nx * (kron(elts, ones(Int64,nx)) .- 1)
    j+= kron(ones(Int64,length(elts)), collect(1:nx))
    # lines indices
    i = collect(1:length(j))
    # sparse matrix
    R = sparse(i,j,ones(Bool,length(i)),length(i),nx*M)
    return R
end


"""
Compute the restriction matrix of the complement of the domain ω in the domain Ω.
"""
function restriction_complement(m::Mesh, ω::Domain, Ω::Domain, d)
    RΩ = restriction(m,Ω,d)
    Rω = restriction(m,ω,d)
    Iω_c = transpose(RΩ) * RΩ - transpose(Rω) * Rω
    i, j, v = findnz(Iω_c)
    Rω_c = sparse(collect(1:length(j)), j, v, length(j), size(RΩ,2))
    return Rω_c
end


function sign_correction(m::Mesh, Ω::Domain, mat, dofdim::Integer)
    N = size(mat[1],2)
    [@assert size(mat[i],2) == N for i in 1:length(mat)]
    v = ones(Int64, N) # initialisation of diagonal (identity)
    # sign changes
    v[element_indices(m,inward_boundary(Ω),dofdim)] *= -1
    # sparse matrix
    S = spdiagm(0=>v)
    # correction
    [A*S for A in mat]
end
