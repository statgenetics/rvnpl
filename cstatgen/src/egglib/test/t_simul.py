from .. import simul, data, egglib_binding



########################################################################

def test_all():

    """
    Launch all tests of this module
    """
    
    print "# Testing the simul module"

    m = simul.CoalesceFiniteAlleleMutator(5)

    print '## ... reconnect a population'
    ps = simul.CoalesceParamSet([10,10,10], M=[[None,0,0],
                                               [0,None,1],
                                               [0,1,None]])
    ps.changePairwiseMigrationRate(10, 0, 1, 0.01)
    simul.coalesce(ps, m, 100)

    print '## ... many parameters'
    ps = simul.CoalesceParamSet([10,10,10], [5, 5, 5], 0.5, 0.5, 100, 4, 0.5, 2)
    simul.coalesce(ps, m, 100)
    [(a,b),(c,d),(e,f)] = ps.samples()

    print '## ... migration matrix'
    ps = simul.CoalesceParamSet([10,10,10], [5, 5, 5], 0.5, 0.5, 100, 
                M=[[None,1,0],[0.4,None,2],[0,2,None]], N=[1,2,3])
    simul.coalesce(ps, m, 100)

    print '## ... change growth rates'
    ps = simul.CoalesceParamSet([10,10], alpha=[-2,2])
    ps.changeAllGrowthRates(5, 0)
    simul.coalesce(ps, m, 100)

    print '## ... bottleneck'
    ps = simul.CoalesceParamSet([30,30,30,30])
    ps.bottleneck(0.5, 1)
    simul.coalesce(ps, m, 100)

    print '## ... change migration rates'
    ps = simul.CoalesceParamSet([30,30,30,30])
    ps.changeAllMigrationRates(1.4, 5)
    simul.coalesce(ps, m, 100)

    print '## ... change population sizes'
    ps = simul.CoalesceParamSet([30,30,30,30])
    ps.changeAllPopulationSizes(0.5, 0.1)
    simul.coalesce(ps, m, 100)

    print '## ... change selfing rate'
    ps = simul.CoalesceParamSet(None, doubleSamples=50)
    ps.changeSelfingRate(0.8, 1)
    simul.coalesce(ps, m, 100)

    print '## ... change one growth rate'
    ps = simul.CoalesceParamSet([20, 20], alpha=[0, -8])
    ps.changeSinglePopulationGrowthRate(5, 1, 10)
    simul.coalesce(ps, m, 100)

    print '## ... change one population size'
    ps = simul.CoalesceParamSet([20, 20], alpha=[0, -8])
    ps.changeSinglePopulationSize(5, 1, 1)
    simul.coalesce(ps, m, 100)

    print '## ... single population bottleneck'
    ps = simul.CoalesceParamSet([20, 20, 20])
    ps.singlePopulationBottleneck(1, 0, 1)
    ps.singlePopulationBottleneck(1.5, 1, 10)
    simul.coalesce(ps, m, 100)

    print '## ... fusion and split'
    ps = simul.CoalesceParamSet([20, 20, 20])
    ps.populationFusion(0.2, 0, 1)
    ps.populationSplit(0.5, 2, 0.1)
    ps.changeAllMigrationRates(2, 2)
    simuls = simul.coalesce(ps, m, 100)

    print '## ... check outputs'
    assert len(simuls)==100
    for sim in simuls:
        sim.tMRCA
        sim.totLength
        sim.nMutations
        sim.nRecomb
        assert isinstance(sim, data.Align)

    simuls = simul.coalesce(ps, m, 100, forceSSR=True)
    assert len(simuls)==100
    for sim in simuls:
        sim.tMRCA
        sim.totLength
        sim.nMutations
        sim.nRecomb
        assert isinstance(sim, data.SSR)

    simuls = simul.coalesce(ps, m, 10, convert=False)
    assert len(simuls)==10
    for sim in simuls:
        sim.tMRCA
        sim.totLength
        sim.nMutations
        sim.nRecomb
        assert isinstance(sim, egglib_binding.DataMatrix)

    print '## ... mutators'

    ps = simul.CoalesceParamSet(40)

    m = simul.CoalesceFiniteAlleleMutator(2, 4, True)
    assert isinstance(simul.coalesce(ps, m, 1)[0], data.Align)

    m = simul.CoalesceFiniteAlleleMutator(1, 6, True)
    x = simul.coalesce(ps, m, 1)[0]
    assert isinstance(x, data.SSR)
    assert x.numberOfLoci()==m.numberOfMutations()

    m = simul.CoalesceFiniteAlleleMutator(0, 4, True)
    m.fixedNumberOfMutation(19)
    for i in simul.coalesce(ps, m, 100): assert i.ls()==19
    assert m.numberOfMutations()==19

    m = simul.CoalesceFiniteAlleleMutator(10, 4, True)
    m.setSites([])
    for i in range(100):
        x = simul.coalesce(ps, m, 1)[0]
        assert x.ls()==m.numberOfMutations()

    m = simul.CoalesceFiniteAlleleMutator(50, 4, True)
    m.setSites(5)
    for i in range(100):
        x = simul.coalesce(ps, m, 1)[0]
        assert x.ls()==5 <= m.numberOfMutations()

    m = simul.CoalesceFiniteAlleleMutator(50, 4, True)
    m.setSites(zip([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9, 1.0], [1,1,1,1,1,1,1,5,5]))
    for i in range(100):
        x = simul.coalesce(ps, m, 1)[0]
        assert x.ls()== 9 <= m.numberOfMutations()

    m = simul.CoalesceFiniteAlleleMutator(2, 6, True)
    m.transitionWeights([[None,2,2,2,2,5],
                         [1,None,1,1,1,5],
                         [1,1,None,3,1,1],
                         [1,1,3,None,1,1],
                         [1,1,1,1,None,1],
                         [5,5,1,1,1,None]])
    simul.coalesce(ps, m, 100)


    m = simul.CoalesceInfiniteAlleleMutator(2)
    simul.coalesce(ps, m, 100)
    
    m = simul.CoalesceInfiniteAlleleMutator(0)
    m.fixedNumberOfMutation(100)
    m.setSites(20)
    simul.coalesce(ps, m, 100)

    m = simul.CoalesceStepwiseMutator(5)
    simul.coalesce(ps, m, 100)

    m = simul.CoalesceStepwiseMutator(0)
    m.fixedNumberOfMutation(100)
    m.setSites(20)
    simul.coalesce(ps, m, 100)

    m = simul.CoalesceTwoPhaseMutator(5, 1, 0.1)
    simul.coalesce(ps, m, 100)

    m = simul.CoalesceTwoPhaseMutator(0, 0.1, 0.8)
    m.fixedNumberOfMutation(100)
    m.setSites(20)
    simul.coalesce(ps, m, 100)

