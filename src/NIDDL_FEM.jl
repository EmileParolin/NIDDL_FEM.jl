module NIDDL_FEM

import Base.==, Base.length, Base.in, Base.iterate
import Base.getindex, Base.lastindex, Base.isempty
import Base.union, Base.setdiff, Base.intersect, Base.issubset

using LinearAlgebra
using SparseArrays
using OrderedCollections

include("mesh.jl")
include("domain.jl")
include("quadrature.jl")
include("finite_element.jl")

export
    Mesh,
    construct_mappings!, detect_relation, tag2range,

    SingleDomain, Domain,
    ==, length, in, iterate, isempty, union, setdiff, intersect, issubset,
    remove, parts, dim, assertequaldim,
    outward_boundary, inward_boundary, embedded_boundary, boundary, skeleton,
    orientation, tag, tags, number_of_elements, element_indices,

    Quadrature,
    edgquad, triquad, tetquad,
    ref2real, weight_matrix, nodes,

    # Lagrange
    P0edg, P1edg, gradP1edg,
    P0tri, P1tri, gradP1tri, curlP1tri,
    P0tet, P1tet, gradP1tet,
    # Nedelec and Raviart-Thomas
    RT, nxRT, divRT, NEDtri, curlNEDtri,
    NEDtet, curlNEDtet,
    # Functions (to evaluate RHS)
    ScaEdgFunc, ScaTriFunc, VecTriFunc,

    nbdof, dimdof, femesh, feelts, femaxelts, femaxdofs, valuetype,
    assemble, femdot, restriction, restriction_complement, sign_correction

end # module
