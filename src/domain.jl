struct SingleDomain
    tag::Tuple{Int64,Int64} # (dim, tag)
    # boundaries (domains of dimension = dim-1)
    outwardbnd::Array{SingleDomain,1} # correctly oriented boundaries
    inwardbnd::Array{SingleDomain,1}  # wrongly oriented boundaries
    # embedded interfaces (domains of dimension = dim-1)
    embeddedbnd::Array{SingleDomain,1}
    function SingleDomain(tag, outwardbnd, inwardbnd, embeddedbnd)
        @assert tag[1] in 0:3 "Incorrect input dimension."
        # Sanity checks on boundary intersections
        msg = "Non empty intersection between outward and inward boundaries."
        @assert isempty(intersect(outwardbnd, inwardbnd)) msg
        msg = "Non empty intersection between outward and embedded boundaries."
        @assert isempty(intersect(outwardbnd, embeddedbnd)) msg
        msg = "Non empty intersection between inward and embedded boundaries."
        @assert isempty(intersect(inwardbnd, embeddedbnd)) msg
        # Sanity checks on boundary elements dimension
        for Γ in union(outwardbnd, inwardbnd, embeddedbnd)
            @assert dim(Γ) == tag[1]-1 "Boundary element has incompatible dimension."
        end
        new(tag, outwardbnd, inwardbnd, embeddedbnd)
    end
end


struct Domain
    # Domain parts
    doms::Array{SingleDomain,1}
    function Domain(doms::Array{SingleDomain,1})
        # Sanity checks on input
        msg = "Cannot create a Domain with multiple occurence of a SingleDomain"
        @assert allunique(doms) msg
        for ω in doms[2:end]
            msg = "It's not possible to create a Domain from SingleDomains of
                    different dimensions."
            @assert dim(ω) == dim(doms[1]) msg
        end
        new(doms)
    end
end
Domain(ω::SingleDomain) = Domain([ω,])
Domain() = Domain(SingleDomain[])


function Base.:(==)(Ω1::SingleDomain, Ω2::SingleDomain)
    Ω1.tag == Ω2.tag || return false
    Ω1.outwardbnd == Ω2.outwardbnd || return false
    Ω1.inwardbnd == Ω2.inwardbnd || return false
    Ω1.embeddedbnd == Ω2.embeddedbnd || return false
    return true
end
function Base.:(==)(Ω1::Domain, Ω2::Domain)
    return parts(Ω1) == parts(Ω2)
end


parts(Ω::Domain) = Ω.doms


dim(Ω::SingleDomain) = Ω.tag[1]
dim(Ω::Domain) = dim(parts(Ω)[1])


"""
Check that two domains have same dimension.

If one of the domain (or both) are empty, the assertion is assumed to be true.
"""
function assertequaldim(Ω1::Domain, Ω2::Domain)
    if isempty(Ω1) || isempty(Ω2) return nothing end
    msg = "The dimension of the first domain should be equal to the dimension
    of the second domain."
    @assert dim(Ω1) == dim(Ω2) msg
end


"""
The length of a domain corresponds to the number of SingleDomain that make it.
"""
Base.length(Ω::Domain) = length(parts(Ω))


"""
Check whether a SingleDomain belongs to a Domain.
"""
Base.in(ω::SingleDomain, Ω::Domain) = Base.in(ω, parts(Ω))


Base.getindex(Ω::Domain, i) = parts(Ω)[i]
Base.lastindex(Ω::Domain) = length(Ω)


"""
Domain can be viewed as an iterable collection.
"""
Base.iterate(Ω::Domain) = Base.iterate(Ω, 1)
function Base.iterate(Ω::Domain, state)
    # Check if we are done iterating
    if state > length(Ω)
        return nothing
    end
    # Return (result, state)
    return (Ω[state], state+1)
end


"""
Check whether a Domain is empty.
"""
Base.isempty(Ω::Domain) = length(parts(Ω)) == 0


"""
Union of domains.
"""
Base.union(doms::Domain...) = Domain(Vector{SingleDomain}(vcat(parts.(doms)...)))


"""
Domain which is not in the intersection of domains Ω1 and Ω2.
"""
function Base.setdiff(Ω1::Domain, Ω2::Domain)
    #Domain(setdiff(parts(Ω1), parts(Ω2)))
    return Domain([ω for ω in parts(Ω1) if !(ω in parts(Ω2))])
end


"""
Intersection between domain Ω1 and domain Ω2.
"""
function Base.intersect(Ω1::Domain, Ω2::Domain)
    assertequaldim(Ω1, Ω2)
    Ωinter = Domain(intersect(parts(Ω1), parts(Ω2)))
    if isempty(Ωinter)
        if isempty(Ω1) || dim(Ω1) == 0
            return Domain()
        else
            return Base.intersect(boundary(Ω1), boundary(Ω2))
        end
    else
        return Ωinter
    end
end


"""
Determine whether every element of domain Ω1 is also in domain Ω2.
"""
function Base.issubset(Ω1::Domain, Ω2::Domain)
    assertequaldim(Ω1, Ω2)
    return issubset(parts(Ω1), parts(Ω2))
end


"""
Remove domain Ω1 from domain Ω2.
"""
function remove(Ω1::Domain, Ω2::Domain)
    assertequaldim(Ω1, Ω2)
    return Domain(setdiff(parts(Ω2), intersect(parts(Ω1), parts(Ω2))))
end


outward_boundary(Ω::SingleDomain) = Ω.outwardbnd
inward_boundary(Ω::SingleDomain) = Ω.inwardbnd
embedded_boundary(Ω::SingleDomain) = Ω.embeddedbnd
boundary(Ω::SingleDomain) = vcat(outward_boundary(Ω), inward_boundary(Ω))
skeleton(Ω::SingleDomain) = vcat(boundary(Ω), embedded_boundary(Ω))


function get_all_boundary_sets(Ω::Domain)
    Γout = Vector{SingleDomain}(vcat(outward_boundary.(parts(Ω))...))
    Γin = Vector{SingleDomain}(vcat(inward_boundary.(parts(Ω))...))
    Γemb = Vector{SingleDomain}(vcat(embedded_boundary.(parts(Ω))...))
    return Γout, Γin, Γemb
end


function outward_boundary(Ω::Domain)
    Γout, Γin, Γemb = get_all_boundary_sets(Ω)
    return Domain(setdiff(Γout, intersect(Γin, Γout)))
end


function inward_boundary(Ω::Domain)
    Γout, Γin, Γemb = get_all_boundary_sets(Ω)
    return Domain(setdiff(Γin, intersect(Γin, Γout)))
end


function embedded_boundary(Ω::Domain)
    Γout, Γin, Γemb = get_all_boundary_sets(Ω)
    return Domain(unique(vcat(Γemb, intersect(Γin, Γout))))
end


boundary(Ω::Domain) = union(outward_boundary(Ω), inward_boundary(Ω))


function skeleton(Ω::Domain)
    Γ = parts(boundary(Ω))
    Σ = parts(embedded_boundary(Ω))
    return Domain(unique(vcat(Γ, Σ)))
end


tag(Ω::SingleDomain) = Ω.tag
function tags(Ω::Domain, d::Int64)
    if isempty(Ω)
        return Tuple{Int64,Int64}[]
    elseif d == dim(Ω)
        return Vector{Tuple{Int64,Int64}}(vcat(tag.(parts(Ω))))
    elseif d < dim(Ω)
        return tags(skeleton(Ω),d)
    else
        error("Asking for tags with dimension > dimension of domain")
    end
end
function tags(Ω::Domain,dims::Array{Int64,1})
    tgs = Vector{Tuple{Int64,Int64}}(undef,0)
    for d in dims push!(tgs, tags(Ω,d)...) end
    return tgs
end
tags(Ω::Domain,dims::UnitRange{Int64}) = tags(Ω, collect(dims))
tags(Ω::Domain) = isempty(Ω) ? Tuple{Int64,Int64}[] : tags(Ω,dim(Ω))


function orientation(Γ::Domain, Ω::Domain)
    @assert dim(Γ) == dim(Ω)-1 "First argument is not a boundary of second argument."
    if issubset(Γ, outward_boundary(Ω))
        return 1
    elseif issubset(Γ, inward_boundary(Ω))
        return -1
    else
        error("Given boundary is not a boundary of given domain.")
    end
end


function number_of_elements(elttags::OrderedDict{Tuple{Int64,Tuple{Int64,Int64}},Tuple{Int64,Int64}},
                            Ω::Domain, d::Integer)
    if isempty(Ω)
        return 0
    else
        return sum([length(tag2range(elttags,tg,d)) for tg in tags(Ω,d:dim(Ω))])
    end
end
function number_of_elements(m::Mesh, Ω::Domain, d::Integer)
    return number_of_elements(m.elttags, Ω, d)
end
function number_of_elements(m::Mesh, d::Integer)
    if d == 0
        return size(m.nod2vtx,2)
    elseif d == 1
        return size(m.edg2vtx,2)
    elseif d == 2
        return size(m.tri2vtx,2)
    elseif d == 3
        return size(m.tet2vtx,2)
    else
        error("Only dimension <= 3 makes sense.")
    end
end


function element_indices(elttags::OrderedDict{Tuple{Int64,Tuple{Int64,Int64}},Tuple{Int64,Int64}},
                         Ω::Domain, d::Integer)
    inds = Int64[]
    if !(isempty(Ω))
        for tg in tags(Ω,d:dim(Ω))
            push!(inds, tag2range(elttags,tg,d)...)
        end
    end
    return inds
end
function element_indices(m::Mesh, Ω::Domain, d::Integer)
    return sort(element_indices(m.elttags, Ω, d))
end
