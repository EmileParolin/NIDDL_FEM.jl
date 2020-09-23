vtx1  = SingleDomain((0, 1), [], [], [])
vtx5  = SingleDomain((0, 5), [], [], [])
vtx9  = SingleDomain((0, 9), [], [], [])

edg2  = SingleDomain((1, 2),  [], [], [vtx1])
edg6  = SingleDomain((1, 6),  [], [], [vtx5])
edg10 = SingleDomain((1, 10), [], [], [vtx9])

tri3  = SingleDomain((2, 3),  [], [], [edg2])
tri7  = SingleDomain((2, 7),  [], [], [edg6])
tri11 = SingleDomain((2, 11), [], [], [edg10])

tet4  = SingleDomain((3, 4),      [tri3],  [],     [])
tet8  = SingleDomain((3, 8),      [tri7],  [tri3], [])
tet12 = SingleDomain((3, 12),     [tri11], [tri7], [])

utet48   = union(Domain.([tet4, tet8])...)
utet812  = union(Domain.([tet8, tet12])...)
utet4812 = union(Domain.([tet4, tet8, tet12])...)


@testset "Domains" begin
    @test tags(Domain(tet4), 0) == [(0,1)]
    @test tags(Domain(tet4), 1) == [(1,2)]
    @test tags(Domain(tet4), 2) == [(2,3)]
    @test tags(Domain(tet4), 3) == [(3,4)]
    @test tags(boundary(Domain(tet4)), 0) == [(0,1)]
    @test tags(boundary(Domain(tet4)), 1) == [(1,2)]
    @test tags(boundary(Domain(tet4)), 2) == [(2,3)]
    @test tags(Domain(tet8), 0) == [(0,5), (0,1)]
    @test tags(Domain(tet8), 1) == [(1,6), (1,2)]
    @test tags(Domain(tet8), 2) == [(2,7), (2,3)]
    @test tags(Domain(tet8), 3) == [(3,8)]
    @test tags(boundary(Domain(tet8)), 0) == [(0,5), (0,1)]
    @test tags(boundary(Domain(tet8)), 1) == [(1,6), (1,2)]
    @test tags(boundary(Domain(tet8)), 2) == [(2,7), (2,3)]
    @test tags(Domain(tet12), 0) == [(0,9), (0,5)]
    @test tags(Domain(tet12), 1) == [(1,10), (1,6)]
    @test tags(Domain(tet12), 2) == [(2,11), (2,7)]
    @test tags(Domain(tet12), 3) == [(3,12)]
    @test tags(boundary(Domain(tet12)), 0) == [(0,9), (0,5)]
    @test tags(boundary(Domain(tet12)), 1) == [(1,10), (1,6)]
    @test tags(boundary(Domain(tet12)), 2) == [(2,11), (2,7)]
    @test tags(utet48, 0) == [(0,5), (0,1)]  
    @test tags(utet48, 1) == [(1,6), (1,2)]
    @test tags(utet48, 2) == [(2,7), (2,3)]
    @test tags(utet48, 3) == [(3,4), (3,8)]
    @test tags(boundary(utet48), 0) == [(0,5)]
    @test tags(boundary(utet48), 1) == [(1,6)]
    @test tags(boundary(utet48), 2) == [(2,7)]
    @test tags(utet812, 0) == [(0,9), (0,1), (0,5)]
    @test tags(utet812, 1) == [(1,10), (1,2), (1,6)]
    @test tags(utet812, 2) == [(2,11), (2,3), (2,7)]
    @test tags(utet812, 3) == [(3,8), (3,12)]
    @test tags(boundary(utet812), 0) == [(0,9), (0,1)]
    @test tags(boundary(utet812), 1) == [(1,10), (1,2)]
    @test tags(boundary(utet812), 2) == [(2,11), (2,3)]
    @test tags(utet4812, 0) == [(0,9), (0,1), (0,5)]
    @test tags(utet4812, 1) == [(1,10), (1,2), (1,6)]
    @test tags(utet4812, 2) == [(2,11), (2,3), (2,7)]
    @test tags(utet4812, 3) == [(3,4), (3,8), (3,12)]
    @test tags(boundary(utet4812), 0) == [(0,9)]
    @test tags(boundary(utet4812), 1) == [(1,10)]
    @test tags(boundary(utet4812), 2) == [(2,11)]
end    
