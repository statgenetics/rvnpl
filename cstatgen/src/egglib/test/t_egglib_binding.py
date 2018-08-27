import os
from .. import egglib_binding
import example_files

########################################################################

def test_ABC():

    """
    Test ABC
    """
    
    print "## Testing egglib_binding.ABC"

    f = open("ABC_input1.txt", 'w')
    f.write("""4.95194298907 0.0153325269453 # 1.15981742471 0.516166666667 10.9333333333
1.98635510461 0.491933633245 # 31.9733452219 21.6196666667 55.0666666667
3.29351499573 0.228908428875 # 20.6259423098 11.9238333333 52.2666666667
4.7437088623 0.305096623189 # 27.2713826892 13.6665 57.8666666667
2.64284692458 0.198574028268 # 15.8299405265 9.7905 48.0
1.02024068344 0.441248261522 # 17.6793791227 12.2138333333 48.6666666667
1.66583032231 0.375171886971 # 22.193263154 15.79 48.1333333333
2.47682701097 0.41098071521 # 30.7821813803 20.7763333333 56.4
2.23190710444 0.264390633197 # 19.2780463838 12.9171666667 51.3333333333
4.90752660309 0.105732060591 # 10.0621998198 5.34483333333 41.4666666667
3.08032053071 0.323599198806 # 28.8700499503 16.6708333333 57.7333333333
2.60467202517 0.0267422062183 # 2.41367410008 1.45183333333 18.0
4.79410305881 0.466240648083 # 44.1357549729 24.1845 65.0666666667
2.26126178743 0.215717101296 # 15.0149336875 10.5583333333 45.7333333333
2.85407209632 0.00567788064123 # 0.532889087031 0.4105 7.6
2.35756332194 0.488754120473 # 35.9229937493 22.2838333333 58.0
2.59873167477 0.132524627148 # 10.814513825 7.52066666667 41.8666666667
3.78029735352 0.25989812506 # 23.4471198294 13.7135 55.0666666667
3.80525557452 0.17916015372 # 16.0493654447 9.7595 50.4
3.55586791557 0.0442993587781 # 3.54214510791 2.16066666667 25.3333333333
4.91906847941 0.282067576661 # 27.7102325256 14.4543333333 59.2
3.02938569678 0.492301567372 # 42.1922771261 26.7441666667 61.4666666667
1.3042908075 0.393186747802 # 20.7199815604 13.2196666667 49.0666666667
0.231909071363 0.121422489456 # 1.00308534029 0.930666666667 10.4
1.85256018463 0.0905626391996 # 5.79908712357 4.07433333333 30.8
0.144764703453 0.0479666941369 # 0.188078501305 0.214 2.26666666667
0.320654364195 0.287756119316 # 3.94964852741 2.71366666667 27.7333333333
3.03275337793 0.498803120198 # 41.2205382027 25.9251666667 63.0666666667
4.58480058627 0.342914241632 # 30.0298673751 17.6878333333 56.5333333333
0.938088339711 0.296943758093 # 13.6983841784 10.2643333333 44.0
2.2423060314 0.401898522407 # 30.9702598816 19.4815 59.7333333333
0.921005398924 0.0577682828909 # 2.35098126631 1.88416666667 18.1333333333
1.27348469316 0.40813404689 # 18.6511180461 13.26 48.1333333333
2.13387569874 0.153659566536 # 11.5668278303 7.80933333333 44.1333333333
1.15057063398 0.371128085872 # 15.4537835239 11.9005 47.4666666667
1.46765056341 0.414746734908 # 23.54115908 16.6858333333 51.2
1.10394729051 0.339684734299 # 16.9584115343 12.3315 48.1333333333
1.3378381372 0.352703818653 # 17.0211043681 12.6275 46.4
1.19431686816 0.429554816402 # 20.5319030591 15.128 48.8
2.87289011439 0.292179557542 # 25.8294475126 16.4733333333 55.8666666667
3.49125361273 0.2136503791 # 17.0211043681 10.8833333333 49.7333333333
4.79254622735 0.217694384159 # 21.5036419826 11.2396666667 54.9333333333
4.76891534596 0.0150428341611 # 1.47328159356 0.806166666667 12.9333333333
3.52170597349 0.470239589002 # 41.7847737066 25.7671666667 62.1333333333
0.308579938981 0.278337291578 # 3.79291644299 2.729 25.8666666667
4.38703040908 0.318761162533 # 29.3715926205 16.1536666667 58.9333333333
4.49931825435 0.0233353606266 # 1.7867457624 1.02666666667 17.2
4.9014891786 0.0306204612331 # 2.5704061845 1.23883333333 18.8
3.98263248294 0.160991598118 # 14.9522408538 8.81716666667 50.1333333333
1.25708893341 0.0387669318678 # 2.19424918189 1.672 20.0
0.326615341063 0.49514820845 # 7.5231400522 6.3805 32.5333333333
4.11922609252 0.442668784601 # 41.0951525352 22.262 62.2666666667
2.6656577415 0.083868534029 # 6.30062979372 3.8895 36.1333333333
0.860622853681 0.0606174872015 # 1.9748242637 1.19616666667 18.4
3.89690292911 0.17427862282 # 15.3597442733 8.55183333333 49.6
2.92461510803 0.462439136844 # 36.9260790896 22.6243333333 60.6666666667
3.66404780628 0.234430737897 # 20.3438245578 11.8148333333 50.8
3.09959322746 0.407777141718 # 36.6126149207 20.7783333333 59.8666666667
3.1998186802 0.465569957742 # 36.9260790896 21.8086666667 62.1333333333
2.84967534561 0.112841294662 # 9.18450014707 5.998 34.6666666667
2.41147988297 0.349742554348 # 28.0236966945 16.9838333333 55.6
1.6202851128 0.494180147502 # 26.9579185204 17.835 54.2666666667
4.3025792422 0.378200913378 # 35.9543401662 21.4943333333 61.7333333333
1.98380508978 0.162098133698 # 10.3756639887 7.3515 42.4
4.62229723369 0.425949845897 # 39.6845637754 20.4628333333 63.7333333333
3.12889919064 0.276346489403 # 23.1650020774 14.1796666667 55.0666666667
0.440733340009 0.0701937662012 # 1.28520309225 0.883333333333 14.6666666667
0.212828834754 0.0280329192247 # 0.282117751958 0.177666666667 3.46666666667
4.86019594803 0.30882644587 # 27.8042717763 15.7166666667 61.8666666667
2.43751472685 0.29437756619 # 23.008269993 14.4385 52.1333333333
2.25464292428 0.0582855932469 # 3.85560927676 2.669 30.1333333333
4.45739883972 0.3354820328 # 29.4342854543 17.2266666667 58.8
2.57383701655 0.369238080742 # 30.0925602088 18.9933333333 55.6
3.77521195494 0.0471847466015 # 4.13772702871 2.14433333333 25.0666666667
4.46980910917 0.490549870633 # 44.4492191418 24.568 64.0
4.43335090973 0.329294370615 # 29.1521677023 15.635 62.0
2.19577682136 0.258246134555 # 17.0211043681 11.3206666667 48.6666666667
1.29140264848 0.497711071806 # 25.0457870905 18.0733333333 54.1333333333
0.654334969958 0.456926256532 # 14.0431947641 11.2698333333 42.2666666667
2.20217444417 0.317910889447 # 23.1650020774 15.9455 51.7333333333
4.41777837164 0.102383714864 # 8.65161106004 4.62416666667 38.0
2.46004227502 0.116104208769 # 9.34123223149 6.2765 39.7333333333
4.68816725743 0.42553313794 # 41.1891917858 22.8011666667 62.4
1.85612205818 0.0326126692392 # 2.53905976762 1.7745 21.3333333333
4.04609030198 0.0828302134053 # 7.64852571974 4.01383333333 33.8666666667
1.04244523031 0.31085060383 # 12.444527503 8.199 43.7333333333
0.533529751109 0.178973234014 # 4.79600178328 4.04433333333 28.0
4.27995103668 0.176800938253 # 16.3001367798 8.50933333333 54.5333333333
""")
    f.close()

    f = open("ABC_input2.txt", 'w')
    f.write("""0.0153325269453 # 1.15981742471 0.516166666667 10.9333333333
0.491933633245 # 31.9733452219 21.6196666667 55.0666666667
0.228908428875 # 20.6259423098 11.9238333333 52.2666666667
0.305096623189 # 27.2713826892 13.6665 57.8666666667
0.198574028268 # 15.8299405265 9.7905 48.0
0.441248261522 # 17.6793791227 12.2138333333 48.6666666667
0.375171886971 # 22.193263154 15.79 48.1333333333
0.41098071521 # 30.7821813803 20.7763333333 56.4
0.264390633197 # 19.2780463838 12.9171666667 51.3333333333
0.105732060591 # 10.0621998198 5.34483333333 41.4666666667
0.323599198806 # 28.8700499503 16.6708333333 57.7333333333
0.0267422062183 # 2.41367410008 1.45183333333 18.0
0.466240648083 # 44.1357549729 24.1845 65.0666666667
0.215717101296 # 15.0149336875 10.5583333333 45.7333333333
0.00567788064123 # 0.532889087031 0.4105 7.6
0.488754120473 # 35.9229937493 22.2838333333 58.0
0.132524627148 # 10.814513825 7.52066666667 41.8666666667
0.25989812506 # 23.4471198294 13.7135 55.0666666667
0.17916015372 # 16.0493654447 9.7595 50.4
0.0442993587781 # 3.54214510791 2.16066666667 25.3333333333
0.282067576661 # 27.7102325256 14.4543333333 59.2
0.492301567372 # 42.1922771261 26.7441666667 61.4666666667
0.393186747802 # 20.7199815604 13.2196666667 49.0666666667
0.121422489456 # 1.00308534029 0.930666666667 10.4
0.0905626391996 # 5.79908712357 4.07433333333 30.8
0.0479666941369 # 0.188078501305 0.214 2.26666666667
0.287756119316 # 3.94964852741 2.71366666667 27.7333333333
0.498803120198 # 41.2205382027 25.9251666667 63.0666666667
0.00759146872179 # 0.0 0.0 0.0
""")
    f.close()

    try:
        ABC = egglib_binding.ABC()
        ABC.number_of_statistics(3)
        ABC.add_fname('ABC_input1.txt', 2)
        ABC.add_fname('ABC_input2.txt', 1)
        ABC.get_threshold(0.5)
        ABC.sd(0)
        ABC.sd(1)
        ABC.sd(2)
        ABC.threshold()
        ABC.obs(0, 12.0)
        ABC.obs(0, 6.2)
        ABC.obs(0, 31.9)
        ABC.rejection("ABC_output1.txt", False)
        ABC.rejection("ABC_output2.txt", True)
        del ABC
        ABC = egglib_binding.ABC()
        ABC.number_of_statistics(3)
        ABC.add_fname('ABC_input1.txt', 2)
        ABC.get_threshold(0.4)
        ABC.number_of_samples()
        ABC.sd(0)
        ABC.sd(1)
        ABC.sd(2)
        ABC.threshold()
        ABC.obs(0, 0.1)
        ABC.obs(1, 6.2)
        ABC.obs(2, 31.9)
        ABC.rejection("ABC_output3.txt", False)
        ABC.rejection("ABC_output4.txt", True)
        ABC.regression("ABC_output3.txt", "ABC_output5.txt", ABC.NONE, "param1\tparam2")
        ABC.regression("ABC_output3.txt", "ABC_output6.txt", ABC.LOG, "param1\tparam2")
        ABC.regression("ABC_output3.txt", "ABC_output7.txt", ABC.TAN)
        del ABC
    finally:
        if os.path.isfile('ABC_input1.txt'): os.remove("ABC_input1.txt")
        if os.path.isfile('ABC_input2.txt'): os.remove("ABC_input2.txt")
        if os.path.isfile('ABC_output1.txt'): os.remove("ABC_output1.txt")
        if os.path.isfile('ABC_output2.txt'): os.remove("ABC_output2.txt")
        if os.path.isfile('ABC_output3.txt'): os.remove("ABC_output3.txt")
        if os.path.isfile('ABC_output4.txt'): os.remove("ABC_output4.txt")
        if os.path.isfile('ABC_output5.txt'): os.remove("ABC_output5.txt")
        if os.path.isfile('ABC_output6.txt'): os.remove("ABC_output6.txt")
        if os.path.isfile('ABC_output7.txt'): os.remove("ABC_output7.txt")


########################################################################

def test_Align():

    """
    Test Align
    """

    print "## Testing egglib_binding.Align"

    align = egglib_binding.Align()
    align = egglib_binding.Align(10, 100)
    assert align.append("name", "G"*100)==11
    assert align.removePosition(50)==99
    assert align.removePosition(0)==98
    assert align.removePosition(97)==97
    assert align.remove(10)==10
    assert align.remove(0)==9
    assert align.remove(5)==8
    align.sequence(1, "A"*97)
    assert align.ls()==97
    assert align.ls(0)==97
    assert align.ls(1)==97
    assert align.ls(7)==97
    assert align.character(7, 96)=='?'
    assert align.get(7, 96)=='?'
    align.set(0, 0, "N")
    for i in range(6): align.set(i, 7, "0")
    for i in range(6,8): align.set(i, 7, "1")
    align.binSwitch(7)
    assert ''.join([align.get(i, 7) for i in range(8)]) == '11111100'
    align2 = align.vslice(0, 4)
    align.clear()
    assert align.ns()==0
    assert align.ls()==0
    assert align2.numberOfSequences()==8
    assert align2.numberOfSites()==4
    align2.populationLabel(4)
    align2.sitePosition(3)
    align2.name(4, "flahiheup")
    assert align2.name(4) == "flahiheup"
    assert align2.find("flahiheup")==4
    assert align2.find("flahi")==-1
    assert align2.find("flahi", False)==4
    align2.group(1, 747)
    assert align2.group(1)==747
    cont = align2.hslice(2,5)
    assert cont.ns()==3
    assert cont.ls(0)==align2.ls()
    
    
########################################################################

def test_BppDiversity():

    """
    Test BppDiversity
    """
    
    print "## Testing egglib_binding.BppDiversity"

    align = egglib_binding.Align()
    align.append("name1", "AAGAAAAAAAAAAAAAAAAAA")
    align.append("name2", "AAGAAAAAAAAAAAAA-AAAA")
    align.append("name3", "AAAACAAAAAAAAAAA-AAAA")
    align.append("name4", "AAAACAAAAAGAAAAA-AAAA")
    align.append("name5", "AAAACAAAAAGAAAAAAAAAA")
    align.append("name6", "AAGACAATAAGAAAAAAAAAG", 999)

    bpp = egglib_binding.BppDiversity()
    bpp.load(align)
    assert bpp.hasOutgroup()==True
    assert bpp.S() == 3
    bpp.Sinf()
    bpp.Ssin()
    bpp.eta()
    bpp.Sext()
    bpp.He()
    bpp.He2()
    bpp.tW()
    bpp.T83()
    bpp.K()
    bpp.H()
    bpp.Ti()
    bpp.Tv()
    bpp.TiTv()
    bpp.load(align, 4)
    bpp.nstop()
    bpp.ncodon1mut()
    bpp.nsyn()
    assert bpp.tWS() != 0
    bpp.tWNS()
    bpp.PiS()
    bpp.PiNS()
    bpp.Ssites()
    bpp.NSsites()
    bpp.SS()
    bpp.SNS()
    a, b, c, d = bpp.MK()
    assert bpp.NI() != 0
    bpp.D()
    bpp.Deta()
    assert bpp.Dfl() != 0
    bpp.Dflstar()
    bpp.F()
    bpp.Fstar()
    bpp.rhoH()


########################################################################

def test_Consensus():

    """
    Test Consensus
    """
    
    print "## Testing egglib_binding.Consensus"

    align = egglib_binding.Align()
    align.append("name1_a", "AAAAAAAAAAAAAA")
    align.append("name1_b", "AAAAAAAAAAAAAA")
    align.append("name2_c", "AAACAAAAAAAAAA")
    align.append("name2_d", "AAACAAATAAAAAA")
    align.append("name3#e", "AAAAAAAAAAAAAA")
    align.append("name3#f", "AAAAAAANAAA?AA")
    consensus = egglib_binding.Consensus()
    consensus.check_sequences(align)
    align.set(5, 11, '?')
    assert align.ns()==6
    consensus.setMissing("#")
    consensus.setDisagreement("?")
    align2 = consensus.consensus(align, '_', True)
    assert align2.ns() == 4
    a, b, c, d = consensus.firstSequenceNames()
    a, b, c, d = consensus.secondSequenceNames()
    a, b, c, d = consensus.roots()
    consensus.consistentPositions()
    consensus.complementaryPositions()
    consensus.uninformativePositions()
    consensus.ambiguousPositions()
    consensus.atLeastPartiallyResolvedAmbiguities()
    consensus.inconsistentPositions()
    align3 = consensus.consensus(align, '#', True)
    assert align3.ns() == 5
    

########################################################################

def test_Container():

    """
    Test Container
    """
    
    print "## Testing egglib_binding.Container"

    cont = egglib_binding.Container()
    cont.append('name1', 'AAAAAAAAAAAA')
    cont.append('name2', 'AAAAAAAAAAAAAAAAAAAAAA')
    cont.append('name3', 'AAAAAAAAAAAACCCG')
    cont.append('name4', 'AAAAAAAAAAAAGCCCAAAAGGGGGCC')
    cont.remove(2)
    cont.name(0)
    cont.name(1, 'name2bis')
    cont.sequence(2)
    cont.sequence(2, 'AAAAAAAAAAAACCCG')
    cont.appendSequence(2, 'AAAAAAAAAAAACCCG')
    cont.set(2,31, 'U')
    assert cont.get(0,4) == 'A'
    assert cont.group(1) == 0
    cont.group(1,743)
    cont2 = cont.hslice(1,3)
    assert cont2.ns() == 2
    assert cont2.ls(0) == 22
    assert cont2.ls(1) == 32
    assert cont2.group(0) == 743
    assert cont.isEqual() == False
    cont.equalize('N')
    assert cont.isEqual() == True
    assert cont.get(0,30) == 'N'
    assert cont.find('name2') == -1
    assert cont.find('name2', False) == 1
    cont.clear()


########################################################################

def test_coalesce():

    """
    Test coalescence module classes and Random
    """
    
    print "## Testing egglib_binding coalescence module and egglib_binding.Random"

    r = egglib_binding.Random()
    r = egglib_binding.Random(1445.12, 0.1485)
    assert r.seed1()==1445.12
    assert r.seed2()==0.1485
    r.erand(4.)
    for i in range(100000):
        if r.irand(120)>=120: raise AssertionError, 'too large irand!'
    r.prand(0.5)
    r.nrand()
    r.grand(0.163)
    r.uniform()

    r = egglib_binding.Random()
    ps = egglib_binding.ParamSet()
    ps.numberOfPopulations()
    ps.addPopulation(0)
    ps.addPopulation(0)
    ps.pairwiseMigrationRate(0, 1, 0.5)
    assert ps.pairwiseMigrationRate(0, 1) == 0.5
    ps.migrationRate(1.2)
    ps.populationSize(1, 2)
    assert ps.populationSize(1)==2
    ps.growthRate(2, 1.8)
    assert ps.growthRate(2) == 1.8
    ps.numberOfSegments(1000)
    assert ps.numberOfSegments()==1000
    ps.recombinationRate(4.)
    assert ps.recombinationRate() == 4
    ps.selfingRate(0.96)
    assert ps.selfingRate()==0.96
    ps.singles(0, 20)
    ps.doubles(1, 20)
    ps.singles(2, 8)
    assert ps.singles(0)==20
    assert ps.doubles(0)==0
    assert ps.singles(1)==0
    assert ps.doubles(1)==20
    assert ps.singles(2)==8
    assert ps.doubles(2)==0
    assert ps.numberOfSamples()==68
    
    e1 = egglib_binding.PopulationFusion(0.1, 0, 1)
    e2 = egglib_binding.AllMigrationRateChange(0.3, 0.4)
    e3 = egglib_binding.PopulationSplit(0.2, 2, 0.8)
    e4 = egglib_binding.AllPopulationSizeChange(0.4, 1)
    e5 = egglib_binding.Bottleneck(0.5, 0.8)
    e6 = egglib_binding.GrowthRateChange(0.6, 0)
    e7 = egglib_binding.SelfingRateChange(0.7, 0.01)
    e8 = egglib_binding.SingleMigrationRateChange(0.8, 0, 2, 1.8)
    e9 = egglib_binding.PopulationBottleneck(0.9, 0, 0.05)
    e10 = egglib_binding.PopulationGrowthRateChange(1.0, 0, 5)
    e11 = egglib_binding.SinglePopulationSizeChange(1.1, 2, 4)
    
    ps.addChange(e1)
    ps.addChange(e2)
    ps.addChange(e3)
    ps.addChange(e4)
    ps.addChange(e5)
    ps.addChange(e6)
    ps.addChange(e7)
    ps.addChange(e8)
    ps.addChange(e9)
    ps.addChange(e10)
    ps.addChange(e11)
    
    ## we suppose that change array methods are called by Controller
    # we will also skip special methods of Controller to apply events
    # the same for Current, Edge, EdgePool, Mutation, Population
    
    controller = egglib_binding.Controller(ps, r)
    x = ps.numberOfSamples()
    while x>1:
        x = controller.step()
    
    controller.reset()
    x = ps.numberOfSamples()
    while x>1:
        x = controller.step()

    arg = controller.getArg()
    arg.time()
    arg.ageUltimateMRCA()
    for i in range(1000):
        arg.ageMRCA(i)
        arg.newick(i)
    # some methods skipped again

    m = egglib_binding.Mutator()
    m.fixedNumberOfMutations(12)
    assert m.fixedNumberOfMutations() == 12
    m.mutationModel('F')
    assert m.mutationModel() == 'F'
    m.numberOfAlleles(7)
    assert m.numberOfAlleles() == 7
    m.transitionWeight(0, 1, 2.)
    m.transitionWeight(0, 2, 2.5)
    m.transitionWeight(2, 3, 3.2)
    m.transitionWeight(2, 4, 1.2)
    assert m.transitionWeight(2, 4) == 1.2
    dm = m.mute(arg, r)
    
    m.reset()
    m.mutationRate(4.)
    m.randomAncestralAllele(True)
    m.mutationModel('I')
    dm = m.mute(arg, r)
    
    m.reset()
    m.mutationRate(12.2)
    m.randomAncestralAllele(False)
    m.mutationModel('S')
    m.numberOfSites(5)
    assert m.numberOfSites() == 5
    m.sitePosition(0, 0)
    m.sitePosition(1, 0.4)
    assert m.sitePosition(1) == 0.4
    m.sitePosition(2, 0.6)
    m.sitePosition(3, 0.8)
    m.sitePosition(3, 0.9)
    m.siteWeight(0, 4)
    m.siteWeight(1, 4)
    m.siteWeight(2, 2)
    assert m.siteWeight(2) == 2
    
    dm = m.mute(arg, r)
    m.numberOfMutations()

    m.reset()
    assert m.numberOfMutations()==0
    m.mutationRate(1.2)
    m.mutationModel('T')
    m.TPMproba(0.1)
    assert m.TPMproba()==0.1
    m.TPMparam(0.7)
    assert m.TPMparam() == 0.7
    m.numberOfSites(5)
    assert m.numberOfSites() == 5
    m.sitePosition(0, 0)
    m.sitePosition(1, 0.4)
    assert m.sitePosition(1) == 0.4
    m.sitePosition(2, 0.6)
    m.sitePosition(3, 0.8)
    m.sitePosition(3, 0.9)
    m.siteWeight(0, 4)
    m.siteWeight(1, 4)
    m.siteWeight(2, 2)
    assert m.siteWeight(2) == 2
    
    ps.setGroups(dm, True)
    ps.reset()


########################################################################

def test_Convert():

    """
    Test Convert
    """
    
    print "## Testing egglib_binding.Convert"

    r = egglib_binding.Random()
    
    ps = egglib_binding.ParamSet()
    ps.addPopulation(0)
    ps.migrationRate(1.5)
    ps.singles(0, 20)
    ps.singles(1, 20)
    
    controller = egglib_binding.Controller(ps, r)
    x = ps.numberOfSamples()
    while x>1:
        x = controller.step()
    
    m = egglib_binding.Mutator()
    m.mutationModel('F')
    m.numberOfAlleles(4)
    m.mutationRate(4.)
    arg = controller.getArg()
    dm = m.mute(arg, r)
    
    align = egglib_binding.Convert.align(dm)
    assert align.ns()==dm.numberOfSequences()
    assert align.ls()==dm.numberOfSites()
    align2 = egglib_binding.Convert.align(dm, 1000, r, True, True, True, True, "$@&%")
    assert align2.ns()==dm.numberOfSequences()
    assert align2.ls()>=dm.numberOfSites()
    align.group(0, 999)


########################################################################

def test_DataMatrix():

    """
    Test DataMatrix
    """
    
    print "## Testing egglib_binding.DataMatrix"

    dm = egglib_binding.DataMatrix()
    dm = egglib_binding.DataMatrix(100, 50)
    assert dm.numberOfSites()==50
    assert dm.numberOfSequences()==100
    dm.set(50,20, -7)
    assert dm.get(50,20) == -7
    assert dm.sitePosition(48) == 0
    dm.sitePosition(48, 14.2)
    assert dm.sitePosition(48) == 14.2
    dm.resize(105, 55)
    assert dm.sitePosition(48) == 14.2
    assert dm.get(50,20) == -7
    assert dm.numberOfSites()==55
    assert dm.numberOfSequences()==105
    dm.populationLabel(53, 100)
    assert dm.populationLabel(53) == 100
    assert dm.character(0, 0) == '0'
    dm.set(4,4, 8)
    assert dm.character(4, 4) == '8'
    dm.shift(4)
    assert dm.get(0, 20) == 11
    assert dm.get(20, 20) == 11
    assert dm.get(50, 20) == 4
    dm.clear()
    assert dm.numberOfSites()==0
    assert dm.numberOfSequences()==0
    


########################################################################

def test_Fasta():

    """
    Test Fasta
    """
    
    print "## Testing egglib_binding.Fasta"

    try:
        string = '>name1@1\nAAAAAAAAAA\n>name2@2\nCCCCCCCCCC\n>name3@2\nGGGGGGGGGG\n'
        f = open('testfile1.txt', 'w')
        f.write(string)
        f.close()
        cont = egglib_binding.Fasta.parsef('testfile1.txt')
        assert cont.ns()==3
        assert cont.name(0) == "name1@1"
        cont = egglib_binding.Fasta.parsef('testfile1.txt', True)
        assert cont.name(0) == "name1"
        assert cont.group(0) == 1
        cont = egglib_binding.Fasta.parse(string)
        assert cont.ns()==3
        assert cont.name(0) == "name1@1"
        cont = egglib_binding.Fasta.parse(string, True)
        assert cont.name(0) == "name1"
        assert cont.group(0) == 1
        egglib_binding.Fasta.formatf("testfile2.txt", cont)
        egglib_binding.Fasta.formatf("testfile2.txt", cont, True,40)
        string2 = egglib_binding.Fasta.format(cont)
        string2 = egglib_binding.Fasta.format(cont, True,100)
        assert string2==string
        string2 = egglib_binding.Fasta.format(cont, True,2)
    finally:
        if os.path.isfile('testfile1.txt'): os.remove('testfile1.txt')
        if os.path.isfile('testfile2.txt'): os.remove('testfile2.txt')


########################################################################

def test_FStatistics():

    """
    Test FStatistics and HFStatistics
    """
    
    print "## Testing egglib_binding.FStatistics/HFStatistics"

    fs = egglib_binding.FStatistics()
    fs.reserve(28)
    fs.loadIndividual(100,100,1)
    fs.loadIndividual(125,125,1)
    fs.loadIndividual(100,125,1)
    fs.loadIndividual(200,200,1)
    fs.loadIndividual(200,200,2)
    fs.loadIndividual(200,200,2)
    fs.loadIndividual(100,100,2)
    fs.loadIndividual(100,200,2)
    fs.loadIndividual(342,342,3)
    fs.loadIndividual(342,342,3)
    assert fs.populationLabel(2)==3
    assert fs.alleleValue(2)==200
    assert fs.numberOfAlleles()==4
    assert fs.numberOfPopulations()==3
    assert fs.numberOfGenotypes()==10
    assert fs.alleleFrequencyTotal(0)==6
    assert fs.alleleFrequencyPerPopulation(0,0)==3
    assert fs.alleleFrequencyPerPopulation(2,0)==0
    assert fs.genotypeFrequencyTotal(3,3)==2
    assert fs.genotypeFrequencyPerPopulation(0,3,3)==0
    assert fs.genotypeFrequencyPerPopulation(2,3,3)==2
    assert fs.populationFrequency(2) == 2
    F = fs.F()
    t = fs.theta()
    f = fs.f()
    P = fs.Vpopulation()
    I = fs.Vindividual()
    A = fs.Vallele()
    assert t - P/(P+I+A) < 0.00000001
    assert (1-f) - A/(A+I) < 0.00000001
    assert (1-F) - A/(P+I+A) < 0.00000001
    fs = egglib_binding.FStatistics()
    fs.reserve(5)
    fs.loadIndividual(100,100,1)
    fs.loadIndividual(125,125,1)
    fs.loadIndividual(100,125,1)
    fs.loadIndividual(200,200,1)
    fs.loadIndividual(200,200,2)
    fs.loadIndividual(200,200,2)
    fs.loadIndividual(100,100,2)
    fs.loadIndividual(100,200,2)
    fs.loadIndividual(342,342,3)
    fs.loadIndividual(342,342,3)
    assert fs.populationLabel(2)==3
    assert fs.alleleValue(2)==200
    assert fs.numberOfAlleles()==4
    assert fs.numberOfPopulations()==3
    assert fs.numberOfGenotypes()==10
    assert fs.alleleFrequencyTotal(0)==6
    assert fs.alleleFrequencyPerPopulation(0,0)==3
    assert fs.alleleFrequencyPerPopulation(2,0)==0
    assert fs.genotypeFrequencyTotal(3,3)==2
    assert fs.genotypeFrequencyPerPopulation(0,3,3)==0
    assert fs.genotypeFrequencyPerPopulation(2,3,3)==2
    assert fs.populationFrequency(2) == 2
    F = fs.F()
    t = fs.theta()
    f = fs.f()
    P = fs.Vpopulation()
    I = fs.Vindividual()
    A = fs.Vallele()
    assert t - P/(P+I+A) < 0.00000001
    assert (1-f) - A/(A+I) < 0.00000001
    assert (1-F) - A/(P+I+A) < 0.00000001    
    fs = egglib_binding.HFStatistics()
    fs.reserve(28)
    fs.loadIndividual(100,1)
    fs.loadIndividual(125,1)
    fs.loadIndividual(100,1)
    fs.loadIndividual(200,1)
    fs.loadIndividual(200,2)
    fs.loadIndividual(200,2)
    fs.loadIndividual(100,2)
    fs.loadIndividual(100,2)
    fs.loadIndividual(342,3)
    fs.loadIndividual(342,3)
    assert fs.populationLabel(2)==3
    assert fs.alleleValue(2)==200
    assert fs.numberOfAlleles()==4
    assert fs.numberOfPopulations()==3
    assert fs.numberOfGenotypes()==10
    assert fs.alleleFrequencyTotal(0)==4
    assert fs.alleleFrequencyPerPopulation(0,0)==2
    assert fs.alleleFrequencyPerPopulation(2,0)==0
    assert fs.populationFrequency(2) == 2
    t = fs.theta()
    T1 = fs.T1()
    T2 = fs.T2()
    assert t - T1/T2 < 0.00000001
    fs = egglib_binding.HFStatistics()
    fs.reserve(3)
    fs.loadIndividual(100,1)
    fs.loadIndividual(125,1)
    fs.loadIndividual(100,1)
    fs.loadIndividual(200,1)
    fs.loadIndividual(200,2)
    fs.loadIndividual(200,2)
    fs.loadIndividual(100,2)
    fs.loadIndividual(100,2)
    fs.loadIndividual(342,3)
    fs.loadIndividual(342,3)
    assert fs.populationLabel(2)==3
    assert fs.alleleValue(2)==200
    assert fs.numberOfAlleles()==4
    assert fs.numberOfPopulations()==3
    assert fs.numberOfGenotypes()==10
    assert fs.alleleFrequencyTotal(0)==4
    assert fs.alleleFrequencyPerPopulation(0,0)==2
    assert fs.alleleFrequencyPerPopulation(2,0)==0
    assert fs.populationFrequency(2) == 2
    t = fs.theta()
    T1 = fs.T1()
    T2 = fs.T2()
    assert t - T1/T2 < 0.00000001


########################################################################

def test_diversity():

    """
    Test diversity components
    """
    
    print "## Testing egglib_binding diversity classes"
    
    align = egglib_binding.Align()
    align.append('', 'AAAAATAAAAAAAAAAAAGAAAAAGAAAAA',0)
    align.append('', 'ATAAATATTAAAAACAAAGAAGAAGAAAAA',0)
    align.append('', 'ATAAAAATTAAAAACAAAGAAGAAGAAAAA',0)
    align.append('', 'ACAAAAATTAAACACAAAGAAGAAAAACAA',1)
    align.append('', 'AGAAAAAAAAAACACAAAAAAAAAAAACAA',1)
    align.append('', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',999)
    
    bpp = egglib_binding.BppDiversity()
    bpp.load(align, 1)
    assert bpp.hasOutgroup()
    assert bpp.S() == 10
    a = bpp.Sinf()
    a = bpp.Ssin()
    a = bpp.eta()
    a = bpp.Sext()
    a = bpp.He()
    a = bpp.He2()
    a = bpp.tW()
    a = bpp.T83()
    a = bpp.K()
    a = bpp.H()
    a = bpp.Ti()
    a = bpp.Tv()
    a = bpp.TiTv()
    a = bpp.D()
    a = bpp.Deta()
    a = bpp.Dfl()
    a = bpp.Dflstar()
    a = bpp.F()
    a = bpp.Fstar()
    a = bpp.rhoH()
    a = bpp.load(align, 4)
    assert bpp.hasOutgroup()
    a = bpp.S()
    a = bpp.Sinf()
    a = bpp.Ssin()
    a = bpp.eta()
    a = bpp.Sext()
    a = bpp.He()
    a = bpp.He2()
    a = bpp.tW()
    a = bpp.T83()
    a = bpp.K()
    a = bpp.H()
    a = bpp.Ti()
    a = bpp.Tv()
    a = bpp.TiTv()
    a = bpp.D()
    a = bpp.Deta()
    a = bpp.Dfl()
    a = bpp.Dflstar()
    a = bpp.F()
    a = bpp.Fstar()
    a = bpp.rhoH()
    a = bpp.nstop()
    a = bpp.ncodon1mut()
    a = bpp.nsyn()
    a = bpp.tWS()
    a = bpp.tWNS()
    a = bpp.PiS()
    a = bpp.PiNS()
    a = bpp.Ssites()
    a = bpp.NSsites()
    a = bpp.SS()
    a = bpp.SNS()
    a,b,c,d = bpp.MK()
    a = bpp.NI()
    
    hd = egglib_binding.HaplotypeDiversity()
    align.set(0,0, "&")
    hd.load(align, True, 1, "ACGT& MRWSYKBDHVN?-")
    a = hd.K()
    a = hd.He()
    assert hd.haplotypeIndex(0) == 0
    a = hd.haplotypeIndex(3)
    a = hd.Kst()
    a = hd.Fst()
    a = hd.Gst()
    a = hd.Hst()
    a = hd.Snn()
    
    ld = egglib_binding.LinkageDisequilibrium()
    ld.load(align, 0.4, 1, "ACGT& MRWSYKBDHVN?-")
    n = ld.numberOfPairs()
    for i in range(n):
        a = ld.d(i)
        a = ld.D(i)
        a = ld.Dp(i)
        a = ld.r(i)
        a = ld.r2(i)
        a = ld.site1(i)
        b = ld.site2(i)
        assert b>a
    correl = ld.correl()
    assert 1 == ld.Rmin(align)
    
    align2 = egglib_binding.Align()
    align2.append('', '00010000000')
    align2.append('', '00010000000')
    align2.append('', '00010000100')
    align2.append('', '00000000100')
    
    nd = egglib_binding.NucleotideDiversity()
    nd.load(align2, False, 1, 0, "01", True)
    assert nd.S()==2
    assert nd.H()!=0    
    align.set(0,0, "T")
    nd.load(align, True)
    assert nd.S()==11
    align.set(0,0, "A")
    nd.load(align, True)
    assert nd.S()==10
    align.set(0,0, "*")
    nd.load(align, True, 1, 0, "ACGT* MRWSYKBDHVN?-")
    assert nd.S()==11
    x= nd.So()
    x= nd.eta()
    assert nd.eta()>nd.S()
    x= nd.nseff()
    x= nd.lseff()
    assert nd.lseffo()>0
    assert nd.nseffo()>0
    assert nd.npop()==2
    assert nd.popLabel(1)==1
    x= nd.Pi()
    x= nd.thetaW()
    x= nd.average_Pi()
    x= nd.pop_Pi(0)
    x= nd.pop_Pi(1)
    x= nd.D()
    x= nd.thetaH()
    x= nd.thetaL()
    x= nd.H()
    x= nd.Z()
    x= nd.E()
    x= nd.FixedDifferences()
    x= nd.CommonAlleles()
    x= nd.SharedAlleles()
    x= nd.SpecificAlleles()
    x= nd.SpecificDerivedAlleles()
    x= nd.Polymorphisms(0)
    x= nd.SpecificAlleles(0)
    x= nd.SpecificDerivedAlleles(0)
    x= nd.Polymorphisms(1)
    x= nd.SpecificAlleles(1)
    x= nd.SpecificDerivedAlleles(1)
    x= nd.FixedDifferences(0,1)
    x= nd.CommonAlleles(0,1)
    x= nd.SharedAlleles(0,1)
    assert len(nd.polymorphic_positions())>0
    x= nd.singleton_positions()
    
    align.append('', 'AAAAATAAAAAAAAAAAAGAAAAAGATAAA',878)
    nd.load(align, True, 1, 0, "ACGT* MRWSYKBDHVN?-")
    for i in range(13):
        x = nd.triConfiguration(i)
    
    x = nd.get_position(4)
    site = nd.get_site(4)
    assert site.numberOfAlleles()==2
    site.allele(0)
    site.allele(1)
    site.alleleFrequency(0)
    site.alleleFrequency(1)
    site.alleleFrequency(0,0)
    site.alleleFrequency(1,0)
    site.alleleFrequency(0,1)
    site.alleleFrequency(1,1)
    site.alleleFrequency(0,2)
    site.alleleFrequency(1,2)
    site.derivedAlleleFrequency()
    site.ns()
    site.ns(0)
    site.ns(1)
    site.isOrientable()
    site.isPolymorphic(0)
    site.hasSpecificAllele(0, False)
    site.hasSpecificAllele(0, True)
    site.haveFixedDifference(0,1)
    site.haveCommonAllele(0,1)
    site.haveSharedAllele(0,1)
    
    r = egglib_binding.Random()
    ps = egglib_binding.ParamSet()
    ps.addPopulation(0)
    ps.migrationRate(1.5)
    ps.singles(0, 20)
    ps.singles(1, 20)
    controller = egglib_binding.Controller(ps, r)
    x = ps.numberOfSamples()
    while x>1:
        x = controller.step()
    m = egglib_binding.Mutator()
    m.mutationModel('F')
    m.numberOfAlleles(4)
    m.mutationRate(4.)
    arg = controller.getArg()
    dm = m.mute(arg, r)
    dm.set(4,0, 978)
    md = egglib_binding.MicrosatelliteDiversity()
    md.load(dm, 978, False)
    md.load(dm, 978, True)
    md.numberOfSites()
    for i in range(md.numberOfSites()): md.He(i)
    for i in range(md.numberOfSites()): md.numberOfAlleles(i)
    for i in range(md.numberOfSites()): md.sizeVariance(i)
    for i in range(md.numberOfSites()): md.thetaAssumingIAM(i)
    for i in range(md.numberOfSites()): md.thetaAssumingSMMfromHe(i)
    for i in range(md.numberOfSites()): md.thetaAssumingSMMfromSizeVariance(i)
   

########################################################################

def test_Ms():

    """
    Test Ms
    """
    
    print "## Testing egglib_binding.Ms"

    string = """//
segsites: 17
positions: 0.0053 0.0863 0.0985 0.2221 0.3759 0.5361 0.5504 0.5906 0.6244 0.6619 0.6762 0.7639 0.7769 0.7871 0.8197 0.9608 0.9831 
00110011001010010
00010001010010010
00010001010010010
00001100100100100
00010001010010011
00110011000011010
00010001010010010
00010001010010010
00010001010010010
10110011000010010
00001100100100100
01110011000010010
00110011001010010
00010001010010010
00010001010010010
00010001010010010
01110011000010010
00110011000010010
00010001010010010
00010001010010010

"""

    dm = egglib_binding.Ms.get(string, 20, False)
    assert dm.numberOfSites()==17

    string2 = """//
segsites: 17
positions: 0.0053 0.0863 0.0985 0.2221 0.3759 0.5361 0.5504 0.5906 0.6244 0.6619 0.6762 0.7639 0.7769 0.7871 0.8197 0.9608 0.9831 
0 0 1 1 0 0 1 1 0 0 1 0 1 0 0 1 0
0 0 0 1 0 0 0 1 0 1 0 0 1 0 0 1 0
0 0 0 0 1 1 0 0 1 0 0 1 0 0 1 0 0
0 0 0 1 0 0 0 1 0 1 0 0 1 0 0 1 1
0 0 1 1 0 0 1 1 0 0 0 0 1 1 0 1 0
0 0 0 1 0 0 0 1 0 1 0 0 1 0 0 1 0
0 0 0 1 0 0 0 1 0 1 0 0 1 0 0 1 0
0 0 0 1 0 0 0 1 0 1 0 0 1 0 0 1 0
1 0 1 1 0 0 1 1 0 0 0 0 1 0 0 1 0
0 0 0 0 1 1 0 0 1 0 0 1 0 0 1 0 0

"""

    dm = egglib_binding.Ms.get(string2, 10, True)
    assert dm.numberOfSites()==17
    assert egglib_binding.Ms.tMRCA()==-1
    
    string3 = """//
time:	0.897534	3.751144
segsites: 6
positions: 0.0424 0.1198 0.2475 0.4468 0.7608 0.8724 
100100
100100
100100
110100
100000
001000
100100
100001
100100
100100
100000
100100
100100
110110
100100
001000
100000
100100
100000
110100

"""
    dm = egglib_binding.Ms.get(string3, 20, False)
    assert dm.numberOfSites()==6
    assert egglib_binding.Ms.tMRCA()==0.897534
    dm = egglib_binding.Ms.get(string, 20, False)
    assert string.split() == egglib_binding.Ms.format(dm, False).split()
    dm = egglib_binding.Ms.get(string2, 10, True)
    assert string2.split() == egglib_binding.Ms.format(dm, True).split()


########################################################################

def test_Staden():
    
    """
    Test Staden
    """
    
    print "## Testing egglib_binding.Staden"

    align = egglib_binding.Staden.parse(example_files.staden, False)
    assert align.find("Stephane-ABP1_AF")>=0
    assert align.find("CONSENSUS")>=0
    align = egglib_binding.Staden.parse(example_files.staden, True)
    assert align.find("Stephane-ABP1_AF")>=0
    assert align.find("CONSENSUS")==-1
    
    



########################################################################

def test_all():
    
    """
    Launch all tests of this module
    """
    
    print "# Testing egglib_binding module"

    test_ABC()
    test_Container()
    test_Align()
    test_BppDiversity()
    test_Consensus()
    test_coalesce()
    test_Convert()
    test_DataMatrix()
    test_Fasta()
    test_FStatistics()
    test_diversity()
    test_Ms()
    test_Staden()

