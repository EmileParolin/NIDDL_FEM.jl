struct Quadrature
    x::Array{Float64,2}
    w::Array{Float64,1}
end

# 1 Gauss points on reference edge
xedg1 = [1/2]; xedg1 = reshape(xedg1,1,1)
wedg1 = [1.0]
# 2 Gauss points on reference edge
xedg2 = [[(3+√3)/6]
         [(3-√3)/6]];  xedg2 = reshape(xedg2,2,1)
wedg2 = [1/2, 1/2]

# 1 Gauss points on reference triangle
xtri1 = [1/3 1/3]'
wtri1 = [1.0]
# 3 Gauss points on reference triangle
xtri3 = [[1/6 1/6]
         [2/3 1/6]
         [1/6 2/3]]
wtri3 = [1/3, 1/3, 1/3]

# 1 Gauss points on reference tetra
xtet1 = [1/4 1/4 1/4]'
wtet1 = [1.0]
# 4 Gauss points on reference tetra
a = (5-  √5)/20
b = (5+3*√5)/20
xtet4 = [[a a a]
         [a a b]
         [a b a]
         [b a a]]
wtet4 = [1/4, 1/4, 1/4, 1/4]

# Dictionnaries of quadrature rules for convenience
edgquad = Dict([(1, Quadrature(xedg1, wedg1)),
                (2, Quadrature(xedg2, wedg2))])
triquad = Dict([(1, Quadrature(xtri1, wtri1)),
                (3, Quadrature(xtri3, wtri3))])
tetquad = Dict([(1, Quadrature(xtet1, wtet1)),
                (4, Quadrature(xtet4, wtet4))])

@inline function ref2real(m::Mesh, s, x)
    if size(x,2) == 1 # edge
        xs = m.vtx[:,m.edg2vtx[1,s]] * (1 .- x[:,1])'
        xs+= m.vtx[:,m.edg2vtx[2,s]] * x[:,1]'
    elseif size(x,2) == 2 # triangle
        xs = m.vtx[:,m.tri2vtx[1,s]] * (1 .- x[:,1] .- x[:,2])'
        xs+= m.vtx[:,m.tri2vtx[2,s]] * x[:,1]'
        xs+= m.vtx[:,m.tri2vtx[3,s]] * x[:,2]'
    elseif size(x,2) == 3 # tetra
        xs = m.vtx[:,m.tet2vtx[1,s]] * (1 .- x[:,1] .- x[:,2] .- x[:,3])'
        xs+= m.vtx[:,m.tet2vtx[2,s]] * x[:,1]'
        xs+= m.vtx[:,m.tet2vtx[3,s]] * x[:,2]'
        xs+= m.vtx[:,m.tet2vtx[4,s]] * x[:,3]'
    end
    xs
end

function compute_vol(s)
    d = size(s,2)-1
    if d == 1 # edges
        vol = norm(s[:,2]-s[:,1]) # edge length
    elseif d == 2 # triangles
        cp = cross(s[:,2]-s[:,1], s[:,3]-s[:,1]) # cross product between 2 edges
        vol = 0.5 * norm(cp) # triangle area
    elseif d == 3 # tetrahedron
        cp = cross(s[:,2]-s[:,1], s[:,3]-s[:,1]) # cross product between 2 edges
        sp = dot(s[:,4]-s[:,1], cp) # dot product with other (non colinear) edge
        vol = abs(sp) / 6 # tet volume
    end
end

function weight_matrix(m::Mesh, Ω::Domain, q::Quadrature; coef=(args...)->1)
    d = dim(Ω)
    # Getting correct elements on which the loop is performed
    if d == 1 elt2vtx = m.edg2vtx
    elseif d == 2 elt2vtx = m.tri2vtx
    elseif d == 3 elt2vtx = m.tet2vtx
    end
    # Initialisation
    vol = Array{Complex{Float64},1}(undef,size(elt2vtx,2))
    # Loop on tags
    for tag in tags(Ω)
        # Loop on elements on which the integration is performed
        for ielt in tag2range(m,tag,d)
            s = m.vtx[:,elt2vtx[:,ielt]] # element vertices
            g = sum(s, dims=2) ./ size(s,2) # barycenter
            vol[ielt] = compute_vol(s) * coef(g)
        end
    end
    spdiagm(0 => kron(vol, q.w))
end

function nodes(m::Mesh, Ω::Domain, q::Quadrature, d)
    elts = element_indices(m,Ω,d)
    # Number of Gauss points per element
    nx = size(q.x,1)
    # Initialisation
    y = Array{Float64,2}(undef,3,nx*length(elts))
    # Loop on elements on which the integration is performed
    for (ie, ielt) in enumerate(elts)
        # Gauss points real coordinates in current element
        y[:,nx*(ie-1).+(1:nx)] = ref2real(m, ielt, q.x) # (dim, nbgausspts)
    end
    return y
end
