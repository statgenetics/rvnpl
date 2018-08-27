from .. import fitmodel, simul, egglib_binding
import example_files
import os



########################################################################

def test_utilities():

    """
    Test the utilities class
    """
    
    print "## Testing utilities classes"

    ds = fitmodel.Dataset()
    ps = simul.CoalesceParamSet([20,20,20], M=1)
    m = simul.CoalesceFiniteAlleleMutator(5)
    aligns = simul.coalesce(ps, m, 100)
    for i in aligns:
        ds.add(i)
    ds.config()
    ls = []
    for align, item in ds.iterator():
        assert align.ns()==20*3
        t, (p1,p2,p3), o, lsi = item
        assert t==20*3
        assert p1==p2==p3==20
        assert o==0
        ls.append(lsi)
    config = [[60, [29,30], 1, None]] * 100
    for i,v in enumerate(config): config[i][-1] = ls[i]
    for align, item in ds.iterator(config):
        pass
    ds.sort_aligns()
    a,b,c = ds.pops()
    
    ps = fitmodel.ParamSample(4)
    assert len(ps)==4
    ps[0] = 1
    ps[1] = 2
    ps[2] = 3
    ps[3] = 4
    try: ps[4] = 1
    except IndexError: pass
    else: raise AssertionError
    ps[1]+=10
    a = ps[0]+ps[1]+ps[2]+ps[3]
    ps[:2]
    assert ps.values() == [1,12,3,4]
    
    f = open('test_posterior.txt', 'w')
    f.write("""PARAM1 PARAM2 PARAM3
4.2 0.14 2.13
8.3 0.19 1.49
1.1 0.13 1.03
18.3 1.42 4.50""")
    f.close()

    (p1,p2,p3), (data1,data2,data3) = fitmodel.import_posterior('test_posterior.txt')
    assert len(data1)==4==len(data2)==len(data3)
    
    f = open('test_posterior.txt', 'w')
    f.write("""4.2 0.14 2.13
8.3 0.19 1.49
1.1 0.13 1.03
18.3 1.42 4.50""")
    f.close()

    (p1,p2,p3), (data1,data2,data3) = fitmodel.import_posterior('test_posterior.txt')
    assert len(data1)==4==len(data2)==len(data3)
    
    os.remove('test_posterior.txt')
    

########################################################################

def test_priors():

    """
    Test the prior classes
    """
    
    print "## Testing prior classes"

    R = egglib_binding.Random()
    
    pd = fitmodel.PriorDiscrete(R)
    try: pd.parse("""[Perl] is intended to be practical (easy to
       use, efficient, complete) rather than beautiful (tiny, elegant,
       minimal).""")
    except fitmodel.PriorParseError: pass
    else: raise AssertionError

    pd.parse("""0.6   0;0.5   0;0.2   0;2.
0.2   0;0.5   0;0.2   2.;10.
0.08  0;0.5   0.2;0.6 0;2.
0.06  0;0.5   0.2;0.6 2.;10.
0.04  0.5;0.8 0;0.2   0;2.
0.02  0.5;0.8 0.2;0.6 2.;10.""")

    ps = pd.draw()
    assert len(ps)==3
    assert pd.number_of_params()==3
    pd.str()
    pd.clear()
    pd.parse("""2 -0.5;0.5
1 0.5;1
1 1;2
""")
    pd.force_positive()
    for i in range(1000): assert pd.draw()[0]>=0

    pd.clear()
    try: pd.draw()
    except RuntimeError, e: assert "cannot draw parameters" in str(e)
    else: raise AssertionError
    pd.add(1, (0.5, 1.5), (1.5, 3.))
    assert len(pd.draw())==2
    
    pd = fitmodel.PriorDumb(R)
    pd.add('U', 0, 10)
    pd.add('E', 1.2)
    pd.add('P', 1.5)
    pd.add('G', 1.3)
    pd.add('N', 8, 1)
    pd.add('F', 4.)
    assert pd.number_of_params()==6
    assert len(pd.draw())==6
    
    pd.clear()
    assert pd.number_of_params()==0
    
    pd.add('N', 0, 1)
    assert pd.number_of_params()==1

    c=0
    for i in range(100):
        [a] = pd.draw()
        if a<0:c+=1
    assert c>0

    pd.force_positive()
    c=0
    for i in range(100):
        a =     pd.draw()[0]
        if a<0:c+=1
    assert c==0
    
    assert pd.str()=="N(0;1)\n"

    pd.parse("U(0.4;0.6) N(1;0.2) F(5)")
    assert pd.number_of_params()==3
    a, b, c = pd.draw()
    
    assert fitmodel.PriorDumb in fitmodel.priors
    assert fitmodel.PriorDiscrete in fitmodel.priors
    x = fitmodel.priors[0](R)
    assert x.number_of_params()==0
    
    

########################################################################

def test_models():

    """
    Test the model classes
    """
    
    print "## Testing demographic models"

    random = egglib_binding.Random()

    print '### SNM'
    fitmodel.SNM(False).generate([[20, [20], 0, 1000]]*10, [0.001], random)
    fitmodel.SNM(True).generate([[20, [20], 0, 1000]]*10, [0.001, 0.005], random)

    print '### PEM'
    fitmodel.PEM(False).generate([[20, [20], 0, 1000]]*10, [0.001, 4.5], random)
    fitmodel.PEM(True).generate([[20, [20], 0, 1000]]*10, [0.001, 4.5, 0.005], random)

    print '### BNM'
    fitmodel.BNM(False).generate([[20, [20], 0, 1000]]*10, [0.001, 0.5, 0.2, 0.2, 1.2], random)
    fitmodel.BNM(True).generate([[20, [20], 0, 1000]]*10, [0.001, 0.5, 0.2, 0.2, 1.2, 0.005], random)

    print '### GDB'
    fitmodel.GDB(False).generate([[20, [20], 0, 1000]]*10, [0.001, 0.5, 1.5], random)
    fitmodel.GDB(True).generate([[20, [20], 0, 1000]]*10, [0.001, 0.5, 1.5, 0.005], random)

    print '### GGDB'
    fitmodel.GGDB(False).generate([[20, [20], 0, 1000]]*10, [0.001, 0.5, 1.5, 0.78], random)
    fitmodel.GGDB(True).generate([[20, [20], 0, 1000]]*10, [0.001, 0.5, 1.5, 0.78, 0.005], random)

    print '### IM'
    fitmodel.IM(False).generate([[20, [10, 0, 5, 5], 0, 1000]]*10, [0.001, 1.1], random)
    fitmodel.IM(True).generate([[20, [10, 0, 5, 5], 0, 1000]]*10, [0.001, 1.1, 0.005], random)

    print '### IMn'
    fitmodel.IMn(False).generate([[20, [10, 0, 5, 5], 0, 1000]]*10, [0.001, 1.1, 2, 3, 0.5], random)
    fitmodel.IMn(True).generate([[20, [10, 0, 5, 5], 0, 1000]]*10, [0.001, 1.1, 2, 3, 0.5, 0.005], random)

    print '### IMG'
    fitmodel.IMG(False).generate([[20, [10, 0, 5, 5], 0, 1000]]*10, [0.001, 0.9, 1.1], random)
    fitmodel.IMG(True).generate([[20, [10, 0, 5, 5], 0, 1000]]*10, [0.001, 0.9, 1.1, 0.005], random)

    print '### IMiG'
    fitmodel.IMiG(False).generate([[20, [10, 0, 5, 5], 0, 1000]]*10, [0.001, 0.4, 2, 0.1, 1, 5], random)
    fitmodel.IMiG(True).generate([[20, [10, 0, 5, 5], 0, 1000]]*10, [0.001, 0.4, 2, 0.1, 1, 5, 0.005], random)

    print '### IMiGn'
    fitmodel.IMiGn(False).generate([[20, [10, 0, 5, 5], 0, 1000]]*10, [0.001, 1.75, 2, 0.1, 1, 1, 3,5,4], random)
    fitmodel.IMiGn(True).generate([[20, [10, 0, 5, 5], 0, 1000]]*10, [0.001, 1.75, 2, 0.1, 1, 1, 3,5,4, 0.005], random)

    print '### SM'
    fitmodel.SM(False).generate([[20, [20, 20, 0], 0, 1000]]*10, [0.001, 0.2, 1.5], random)
    fitmodel.SM(True).generate([[20, [20, 20, 0], 0, 1000]]*10, [0.001, 0.2, 1.5, 0.005], random)

    print '### AM'
    fitmodel.AM(4, False).generate([[20, [20], 0, 1000]]*10, [0.001, 0.2, 1.5], random)
    fitmodel.AM(4, True).generate([[20, [20], 0, 1000]]*10, [0.001, 0.2, 1.5, 0.005], random)

    print '### MRC'
    fitmodel.MRC(False).generate([[20, [20], 0, 1000]]*10, [0.001, 0.2, 1.5, 0.5], random)
    fitmodel.MRC(True).generate([[20, [20], 0, 1000]]*10, [0.001, 0.2, 1.5, 0.5, 0.005], random)

    print '### DOM'
    fitmodel.DOM(False).generate([[20, [10,10], 0, 1000]]*10, [0.001, 10, 0.4, 0.1, 0.05, 0.2], random)
    fitmodel.DOM(True).generate([[20, [10,10], 0, 1000]]*10, [0.001, 10, 0.4, 0.1, 0.05, 0.2, 0.005], random)

    assert fitmodel.SM in fitmodel.models
    fitmodel.models[4](True).generate
    
    print '### custom model'
    
    model = example_files.CustomModel

    f = open('ErrorModel.py', 'w')
    f.write(model)
    f.close()
    fitmodel.add_model('ErrorModel')
    assert fitmodel.models[-1].name=='ErrorModel'
    fitmodel.models[-1](False).generate([[20, [10, 10], 0, 1000]]*10, [0.001, 1.75, 0.1], random)
    fitmodel.models[-1](True).generate([[20, [10, 10], 0, 1000]]*10, [0.001, 1.75, 0.1, 0.005], random)
    
    os.remove('ErrorModel.py')
    os.remove('ErrorModel.pyc')


########################################################################

def test_stats():

    """
    Test the stats classes
    """
    
    print '## Testing summary statistics sets'
    
    random = egglib_binding.Random()

    cfg = [[20, [10, 10], 0, 1000]]*10
    ds = fitmodel.SNM(False).generate(cfg, [0.001], random)
    
    print '### TPH'
    stats = fitmodel.TPH()
    stats.compute(ds, cfg)
    a, b, c = stats.stats
    
    print '### TPS'
    stats = fitmodel.TPS()
    stats.compute(ds, cfg)
    a, b1, b2, c = stats.stats

    print '### SFS'
    stats = fitmodel.SFS(2)
    stats.compute(ds, cfg)
    x, a, b = stats.stats

    stats = fitmodel.SFS(7)
    stats.compute(ds, cfg)
    assert len(stats.stats)==7+1

    print '### JFS'
    stats = fitmodel.JFS(6)
    stats.compute(ds, cfg)
    assert len(stats.stats)==2+6*6
    assert sum(stats.stats[2:])-1<0.000000000001

    fitmodel.summstats[0]().compute(ds, cfg)

    
########################################################################

def test_all():

    """
    Launch all tests of this module
    """
    
    print "# Testing the fitmodel module"

    test_utilities()
    test_priors()
    test_models()
    test_stats()
