"""
Data structure to represent a mesh.

In elttags, the key is of the form `(dim, tag)`, where `dim` is the dimension
of the elements of the geometrical entity associated with `tag = (dtag, flag)`
which is itself of dimension `dtag`.
"""
mutable struct Mesh
    vtx::Array{Float64,2}
    # mappings 2 vtx
    nod2vtx::Array{Int64,2}
    edg2vtx::Array{Int64,2}
    tri2vtx::Array{Int64,2}
    tet2vtx::Array{Int64,2}
    # mappings 2 nod
    edg2nod::Array{Int64,2}
    tri2nod::Array{Int64,2}
    tet2nod::Array{Int64,2}
    # mappings 2 edg
    tri2edg::Array{Int64,2}
    tet2edg::Array{Int64,2}
    # mappings 2 tri
    tet2tri::Array{Int64,2}
    # dictionary with key=(d, tag) and value=(start, number)
    elttags::OrderedDict{Tuple{Int64,Tuple{Int64,Int64}},Tuple{Int64,Int64}}
end


function Mesh(vtx, nod2vtx, edg2vtx, tri2vtx, tet2vtx, elttags; sanity_checks=true)
    # Sanity checks on inputs
    if sanity_checks
        @assert size(vtx,1) == 3 "3D coordinates should be in column."
        @assert size(nod2vtx,1) == 1 "Input nod2vtx invalid."
        @assert size(edg2vtx,1) == 2 "Input edg2vtx invalid."
        @assert size(tri2vtx,1) == 3 "Input tri2vtx invalid."
        @assert size(tet2vtx,1) == 4 "Input tet2vtx invalid."
        for t in keys(elttags)
            @assert t[1] == t[2][1] "Input elttags invalid."
        end
    end
    # Interior nodes
    nod2vtx, edg2nod, elttags = construct_mappings!(edg2vtx, nod2vtx, elttags)
    nod2vtx, tri2nod, elttags = construct_mappings!(tri2vtx, nod2vtx, elttags)
    nod2vtx, tet2nod, elttags = construct_mappings!(tet2vtx, nod2vtx, elttags)
    # Interior edges
    edg2vtx, tri2edg, elttags = construct_mappings!(tri2vtx, edg2vtx, elttags)
    edg2vtx, tet2edg, elttags = construct_mappings!(tet2vtx, edg2vtx, elttags)
    # Interior triangles
    tri2vtx, tet2tri, elttags = construct_mappings!(tet2vtx, tri2vtx, elttags)
    # Sanity checks on outputs
    if sanity_checks
        for d in 0:3
            tt = [t[2] for t in elttags if t[1][1]==d]
            elt2vtx = [nod2vtx, edg2vtx, tri2vtx, tet2vtx][d+1]
            if length(elt2vtx) > 0
                msg = "Something is wrong with the mesh tags (dim = $d)."
                @assert sum(tt[argmax(tt)]) == size(elt2vtx,2)+1 msg
                for s in sum.(tt)
                    @assert s in vcat(getindex.(tt,1), size(elt2vtx,2)+1) msg
                end
            end
        end
    end
    Mesh(vtx, nod2vtx, edg2vtx, tri2vtx, tet2vtx, edg2nod, tri2nod, tet2nod,
         tri2edg, tet2edg, tet2tri, sort(OrderedDict(elttags)))
end


"""
We assume that sml boundary elements, for which the list of the nodes gives the
orientation of the simplex, are already provided as sml2vtx.  The other sml
elements that will be constructed from big2vtx are interior sml elements. They
are added to sml2vtx with their list of nodes sorted so that they can be
uniquely identified. The numbering of their nodes is unimportant in this case
as they are shared by two big simplices with reversed numbering: no intrinsic
numbering/orientation.
"""
function construct_mappings!(big2vtx, sml2vtx, elttags)
    Nb = size(big2vtx,1) # dimension of big element: edg=2, tri=3, tet=4
    Ns = size(sml2vtx,1) # dimension of small element: nod=1, edg=2, tri=3
    # Creating all combinations of vertices of big element that can make small
    # elements, in order to respect orientation of boundary small elements
    if (Ns,Nb) == (1,2)      # Nodes in edge
        comb = [[1,], [2,]]
    elseif (Ns,Nb) == (1,3)  # Nodes in triangle
        comb = [[1,], [2,], [3,]]
    elseif (Ns,Nb) == (1,4)  # Nodes in tetrahedron
        comb = [[1,], [2,], [3,], [4,]]
    elseif (Ns,Nb) == (2,3)  # Edges in triangle
        comb = [[1,2], [2,3], [3,1]]
    elseif (Ns,Nb) == (2,4)  # Edges in tetrahedron
        comb = [[1,2], [1,3], [1,4], [2,3], [2,4], [3,4]]
    elseif (Ns,Nb) == (3,4)  # Triangles in tetrahedron
        comb = [[2,3,4], [1,4,3], [1,2,4], [1,3,2]]
    end
    # Initialisation
    is = size(sml2vtx,2) # initialisation of small index
    big2sml = Array{Int64,2}(undef,binomial(Nb, Ns),size(big2vtx,2))
    sml2vtx = hcat(sml2vtx, zeros(size(sml2vtx,1), binomial(Nb, Ns)*size(big2vtx,2))) # too large
    newelttags = copy(elttags) # Copying otherwise looping on tags is uncertain
    seen = Dict{Array{Int64,1},Int64}() # Dict using small as keys and small indicies as values
    # Adding small elements already identified in seen
    for is in 1:size(sml2vtx,2) seen[sort(sml2vtx[:,is])] = is end
    # Creating mappings
    for ((d, tag), (istart, inb)) in sort(OrderedDict(elttags)) # loop on elttags
        if d == Nb-1
            # Initialisation of tag info
            @assert !((Ns-1, tag) in keys(newelttags)) "Tag is already present."
            newelttags[(Ns-1,tag)] = (is+1, 0)
            for ib in istart:(istart+inb-1) # loop on big elements
                big = big2vtx[:,ib] # big element
                for (i, ind) in enumerate(comb) # loop on small elements contained in big one
                    sml = sort(big[ind]) # small element
                    # testing if sml elt has already been seen
                    if !(sml in keys(seen))
                        is += 1 # next small index
                        seen[sml] = is # updating small index in sml2vtx
                        sml2vtx[:,is] = sml # adding the unseen small to sml2vtx
                        newelttags[(Ns-1,tag)] = newelttags[(Ns-1,tag)].+(0, 1) # update number of elements added
                    end
                    big2sml[i,ib] = seen[sml] # adding the small index to the list of small
                end
            end
        end
    end
    elttags = newelttags
    sml2vtx[:,1:is], big2sml, elttags
end


"""
Detect relation between a surface domain and a volume domain.

Return `true` or `false` to indicate whether or not the surface domain is a
boundary of the volume domain. The second returned value is the sign of the
orientation. If the sign = 0, the surface domain is embedded in the volume
domain.
"""
function detect_relation(m::Mesh, tagΓ::Tuple{T,T}, tagΩ::Tuple{T,T}) where T <: Integer
    dΓ, flagΓ = tagΓ
    dΩ, flagΩ = tagΩ
    @assert dΓ == dΩ-1 "Incorrect input."
    # Range of elements in Γ
    rgΓ = tag2range(m, tagΓ, dΓ)
    # Data to monitor detection
    seen = zeros(Int64,length(rgΓ))
    sgns = zeros(Int64,length(rgΓ))
    # Loop on elements in Ω
    for ielt in tag2range(m, tagΩ, dΩ)
        if dΩ == 1
            elt = m.edg2nod[:,ielt]
        elseif dΩ == 2
            elt = m.tri2edg[:,ielt]
        elseif dΩ == 3
            elt = m.tet2tri[:,ielt]
        end
        # Loop on faces of element
        for (nfce, ifce) in enumerate(elt)
            if ifce in rgΓ
                # Element belongs to the boundary Γ
                seen[ifce-rgΓ[1]+1] += 1
                # Sign
                if dΩ == 1
                    elt2vtx = m.edg2vtx[:,ielt]
                    fce2vtx = m.nod2vtx[:,ifce]
                    # Sign
                    sgns[ifce-rgΓ[1]+1] += fce2vtx[1] == elt2vtx[2] ? 1 : -1
                elseif dΩ == 2
                    elt2vtx = m.tri2vtx[:,ielt]
                    fce2vtx = m.edg2vtx[:,ifce]
                    nde2vtx = setdiff(elt2vtx,fce2vtx) # last/opposite node
                    se = m.vtx[:,elt2vtx] # vertices of element
                    sf = m.vtx[:,fce2vtx] # vertices of face
                    sn = m.vtx[:,nde2vtx] # vertex of last/opposite node
                    # Normal vector (non-normalized)
                    nrm = cross(se[:,2]-se[:,1], se[:,3]-se[:,1])
                    # Tangent vector (non-normalized)
                    tau = sf[:,2]-sf[:,1]
                    # Sign
                    sgns[ifce-rgΓ[1]+1] += sign(dot(sn-sf[:,1], cross(nrm, tau)))
                elseif dΩ == 3
                    elt2vtx = m.tet2vtx[:,ielt]
                    fce2vtx = m.tri2vtx[:,ifce]
                    nde2vtx = setdiff(elt2vtx,fce2vtx) # last/opposite node
                    se = m.vtx[:,elt2vtx] # vertices of element
                    sf = m.vtx[:,fce2vtx] # vertices of face
                    sn = m.vtx[:,nde2vtx] # vertex of last/opposite node
                    # Normal vector (non-normalized)
                    nrm = cross(sf[:,2]-sf[:,1], sf[:,3]-sf[:,1])
                    # Sign
                    sgns[ifce-rgΓ[1]+1] += sign(dot(sf[:,1]-sn, nrm))
                end
            end
        end
    end
    # Interpretation of the results
    if sum(seen) == 0
        return (false, Int64(0))
    elseif sum(seen) > 0 && prod(seen) == 0
        error("Some (but not all) boundary elements of $(tagΓ) are not present
              in domain $(tagΩ).")
    else
        if sum(seen) == length(seen)
            abs(sum(sgns)) == length(seen) || "Something went wrong: Γ = $(tagΓ) and Ω = $(tagΩ)."
            abs(sum(sgns)) == length(sgns) || "Variable orientation in boundary domain $(tagΓ)."
            return (true, sgns[1])
        elseif sum(seen) == 2*length(seen)
            sum(sgns) == 0 || "Something went wrong: Γ = $(tagΓ) and Ω = $(tagΩ)."
            return (true, Int64(0))
        else
            error("Some (but not all) boundary elements of $(tagΓ) are present
                  twice in domain $(tagΩ).")
        end
    end
end


function tag2range(elttags, tag, d)
    (i, l) = elttags[(d,tag)]
    return i:(i+l-1)
end
function tag2range(m::Mesh, tag, d)
    return tag2range(m.elttags, tag, d)
end
