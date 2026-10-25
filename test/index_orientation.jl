using Test,LinearAlgebra,Random,CovarianceSelection
using CovarianceSelection.SingleIndexGraph

@testset "Index orientation without a nonzero first coordinate" begin
    orient=CovarianceSelection.SingleIndexGraph._orient_index_direction
    @test orient([0.,-2.,1.]) ≈ [0.,2.,-1.]/sqrt(5)
    @test orient([0.,-2.,1.];anchor=2) ≈ [0.,2.,-1.]/sqrt(5)
    @test_throws DomainError orient([0.,-2.,1.];anchor=1)
    @test_throws DomainError orient(zeros(3))
    @test_throws ArgumentError orient([0.,-2.,1.];anchor=4)
    @test_throws ArgumentError orient([0.,-2.,1.];reference=zeros(3))
    @test_throws ArgumentError orient([0.,-2.,1.];anchor=2,reference=[0.,1.,0.])
    @test_throws DomainError orient([0.,1.,0.];reference=[1.,0.,0.])
    # Two nearly tied components can have opposite signs. Reference alignment
    # prevents a sign flip when their order changes between bootstrap draws.
    reference=[0.,1.,-.9]
    next=[0.,.9,-1.]
    @test dot(orient(next),reference)<0
    @test dot(orient(next;reference),reference)>0
    @test norm(orient(next;reference))≈1
    z=collect(range(-1.,1.;length=120))
    X=hcat(ones(120),z)
    score=hcat(zeros(120),z,.1z)
    idx=estimate_index(X,score;blocks=1,threshold=0,unshrink=true)
    @test idx.direction[1]==0
    @test idx.direction≈[0.,1.,.1]/sqrt(1.01)
    @test_throws DomainError estimate_index(X,score;blocks=1,threshold=0,anchor=1)
    flipped=estimate_index(X,score;blocks=1,threshold=0,reference=-idx.direction)
    @test flipped.direction≈-idx.direction
    tuning=select_index_tuning(X,score;blocks=1,block_options=[1],
        multipliers=[0.,.5],unshrink_options=[true],rng=MersenneTwister(42))
    opposite=select_index_tuning(X,score;blocks=1,block_options=[1],
        multipliers=[0.,.5],unshrink_options=[true],reference=-idx.direction,rng=MersenneTwister(42))
    @test tuning.estimate.direction[1]==0
    @test all(q->isempty(q.failure),tuning.scores)
    @test [q.alignment for q in tuning.scores]≈[q.alignment for q in opposite.scores]
    @test tuning.multiplier==opposite.multiplier
    @test dot(tuning.estimate.direction,opposite.estimate.direction)<0
end
