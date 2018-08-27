from .. import wrappers, data, tools
import example_files


########################################################################

def test_ms():

    """
    Test ms wrapper
    """
    
    print '## Testing ms wrapper'
    results = wrappers.ms(20, 100, theta=5.0, segsites=10, tMRCA=True)
    assert len(results)==100
    for i in results:
        i.prob
        i.tMRCA
        
    results = wrappers.ms(20, 100, theta=5.0, tMRCA=False)
    for i in results:
        try: i.prob
        except AttributeError: pass
        else: raise AssertionError
        try: i.tMRCA
        except AttributeError: pass
        else: raise AssertionError
    
    [result] = wrappers.ms(20, 1, theta=5, T=True, tMRCA=True, segsites=10)
    assert result.prob>0
    assert result.tMRCA>0
    assert len(result.trees)==1
    assert result.trees[0].number_of_leaves()==20
    
    results = wrappers.ms(20,10, theta=5, T=True, F=2, r=(0.3, 1000),
                          c=(0.5, 50),  G=5, I=(3,10,5,5),
                          ma=(0,1,2, 2,0,3, 1,1,0),
                          eG=[(0.1, -2),(0.4, 3)],
                          eg=[(0.1, 1, 0),(0.4, 1,2)],
                          eM=[(0.1, 0),(0.4, 5),(5,10)],
                          em=[(0.1, 0, 1, 10),(0.4, 0, 1, 3)],
                          ema=[(0.1, 3, 0,1,2, 2,0,3, 1,1,0),
                               (0.4, 3, 0,1,2, 2,0,3, 1,1,0)],
                          eN=[(0.1, 0.1),(0.4, 0.5)],
                          en=[(0.1, 2, 5),(0.4, 2, 1)],
                          es=[(1, 2, 0.5),(1.2, 1, 0.1)],
                          ej=[(1.5, 1, 4),(1.6, 2, 3)],
                          n=[(0,0.5), (0,2.), (0,0.1)], g=[(1,3),(0,2)],
                          m=[(0,1,2),(2,1,0.4)]
                        )
    [i.polymorphism()['S'] for i in results]
    for i in results:
        for j in i.trees:
            assert j.number_of_leaves()==20


########################################################################

def test_blast():
    
    """
    Test BLAST wrappers
    """
    
    print "## Testing wrappers for BLAST+ programs"
    
    cds = data.Container(string=example_files.cds_medicago)
    prot = tools.translate(cds)
    db1 = wrappers.BLASTdb(cds, 'nucl')
    p1 = db1.path()
    db2 = wrappers.BLASTdb(prot, 'prot')
    p2 = db2.path()

    blast = wrappers.BLAST()
    hits = blast.blastn(cds, db1, evalue=0.1, penalty=-4, reward=5, gapopen=10, gapextend=6)
    hits[cds.name(0)][0]['qstart']
    blast.xml_results
    blast.blastp(prot, db2, evalue=0.1, gapopen=8, gapextend=2)
    blast.tblastn(prot, db1, evalue=0.1, gapopen=8, gapextend=2)
    blast.blastx(cds, db2, evalue=0.1)
    blast.tblastx(cds, db1, evalue=0.1)
    
    blast.blastn("CCGAAGGTGTTGCGCGCGAAAGGTGTGCA", db1, evalue=0.1, penalty=-4, reward=5, gapopen=10, gapextend=6)['']
    blast.blastp("ALGTKJEJDCSHSAHAHANQNEHJHRGF", db2, evalue=0.1, matrix="BLOSUM62",gapopen=8, gapextend=2)['']
    blast.tblastn("ALGTKJEJDCSHSAHAHANQNEHJHRGF", db1, evalue=0.1, matrix="BLOSUM62",gapopen=8, gapextend=2)['']
    blast.blastx("CCGAAGGTGTTGCGCGCGAAAGGTGTGCA", db2, evalue=0.1)['']
    blast.tblastx("CCGAAGGTGTTGCGCGCGAAAGGTGTGCA", db1, evalue=0.1)['']
    
    bl2seq = wrappers.BL2SEQ()
    bl2seq.blastn("CCGAAGGTGTTGCGCGCGAAAGGTGTGCA", "GTGAAAGCGCGCGAGAGAGAAAGGAGTGTTGGCGCGAGAGCGCCGCGCCGAGAAGAGAGAGT", evalue=0.1, penalty=-4, reward=5, gapopen=10, gapextend=6)[:2]
    bl2seq.blastp("ALGTKJEJDCSHSAHAHANQNEHJHRGF", "MACCCMAMTGGJHJFKSALAAKDSKFFFIDKASKAADDDKFJC", evalue=0.1, matrix="BLOSUM62",gapopen=8, gapextend=2)[:2]
    bl2seq.tblastn("ALGTKJEJDCSHSAHAHANQNEHJHRGF", "GTGAAAGCGCGCGAGAGAGAAAGGAGTGTTGGCGCGAGAGCGCCGCGCCGAGAAGAGAGAGT", evalue=0.1, matrix="BLOSUM62",gapopen=8, gapextend=2)[:2]
    bl2seq.blastx("CCGAAGGTGTTGCGCGCGAAAGGTGTGCA", "MACCCMAMTGGJHJFKSALAAKDSKFFFIDKASKAADDDKFJC", evalue=0.1)[:2]
    bl2seq.tblastx("CCGAAGGTGTTGCGCGCGAAAGGTGTGCA", "GTGAAAGCGCGCGAGAGAGAAAGGAGTGTTGGCGCGAGAGCGCCGCGCCGAGAAGAGAGAGT", evalue=0.1)[:2]


########################################################################

def test_alignment():
    
    """
    Test sequence alignment wrappers
    """
    
    print "## Testing wrappers for multiple alignment programs"

    cont = data.Container(string=example_files.aat, groups=True)

    groups = cont.groups()
    for i in groups: groups[i].sort()

    #aln = wrappers.clustal(cont, False)#too noisy and useful?
    #groups2 = aln.groups()
    #for i in groups2: groups2[i].sort()
    #assert groups==groups2

    print '### clustal'
    aln = wrappers.clustal(cont, True)
    groups2 = aln.groups()
    for i in groups2: groups2[i].sort()
    assert groups==groups2

    #aln = wrappers.muscle(cont, False)#too noisy and useful?
    #groups2 = aln.groups()
    #for i in groups2: groups2[i].sort()
    #assert groups==groups2

    print '### muscle'
    aln = wrappers.muscle(cont, True)
    groups2 = aln.groups()
    for i in groups2: groups2[i].sort()
    assert groups==groups2

    print "### with proteins"

    cds = data.Container()
    for i in cont:
        cds.append(i.name, tools.longest_orf(i.sequence))

    print "### clustal"
    aln = wrappers.clustal(cont, True)
    groups2 = aln.groups()
    for i in groups2: groups2[i].sort()
    assert groups==groups2

    print "### muscle"
    aln = wrappers.muscle(cont, True)
    groups2 = aln.groups()
    for i in groups2: groups2[i].sort()
    assert groups==groups2

    cont.name(1, cont.name(0))
    aln = wrappers.muscle(cont, nogroups=True)
    aln = wrappers.muscle(cont, nogroups=True)
    

########################################################################

def test_phylogeny():
    
    """
    Test phylogeny wrappers
    """
    
    print "## Testing wrappers for phylogeny reconstruction programs"

    cont = data.Container(string=example_files.aat, groups=True)
    aln = wrappers.muscle(cont)
    assert len(aln.groups())==2
    cds = data.Container.create([[i.name, tools.longest_orf(i.sequence), i.group] for i in cont])
    prot = tools.translate(cds)
    prot = data.Container.create([[i.name, i.sequence.rstrip('*'), i.group] for i in prot])
    prot = wrappers.muscle(prot)
    assert len(prot.groups())==2

    print '### Testing phylip'

    tree = wrappers.nj(aln)
    assert 'L0163C' in tree.all_leaves()
    tree = wrappers.nj(aln, True)
    assert 'L0163C@0' in tree.all_leaves()
    
    print "### Testing all phyml models"

    tree, lk = wrappers.phyml(aln)
    assert tree.number_of_leaves()==len(cont)

    for model in ['HKY85', 'JC69', 'K80', 'F81', 'F84', 'TN93']:
        print '#### ' + model
        wrappers.phyml(aln, model=model)

    for model in ['LG', 'WAG', 'JTT', 'MtREV', 'Dayhoff', 'DCMut',
                  'RtREV', 'CpREV', 'VT', 'Blosum62', 'MtMam', 'MtArt',
                  'HIVw', 'HIVb']:
        print '#### ' + model
        wrappers.phyml(prot, model=model)
    
    print '#### with rates'

    tree, lk = wrappers.phyml(aln, rates=8)
    tree, lk = wrappers.phyml(prot, model='Blosum62', rates=8)
    
    print "### Some bootstrap"

    tree, lk = wrappers.phyml(aln, rates=2, boot=10)
    tree, lk = wrappers.phyml(prot, model='Blosum62', rates=1, boot=10)
    
    print '### Other options'
    
    tree, lk = wrappers.phyml(aln, rates=1, boot=0, topo=tree)
    tree, lk = wrappers.phyml(aln, rates=1, boot=0, start=tree)
    tree, lk = wrappers.phyml(prot, model='Blosum62', rates=1, boot=1, topo=tree)
    tree, lk = wrappers.phyml(aln, rates=1, boot=0, search="SPR")
    tree, lk = wrappers.phyml(prot, model='Blosum62', rates=1, boot=0, search="SPR")
    tree, lk = wrappers.phyml(aln, rates=1, boot=0, search="BEST")
    tree, lk = wrappers.phyml(prot, model='Blosum62', rates=1, boot=0, search="BEST")

    print '### Codeml'
    
    seq = data.Align(string=example_files.nork)
    tre, lk = wrappers.phyml(seq)
    tre.findMonophyleticGroup(['Melilotus','truncatul','Vicia','Pisum']).set_label('$1')
    codeml = wrappers.Codeml(seq, tre)

    print '#### fixing omega'
    codeml.number_of_categories(3)
    codeml.fix_omega(1.5)
    codeml.fit('M8')

    codeml.number_of_categories(6)
    for model in ['M0', 'M1a', 'M2a', 'M8a', 'M8', 'A0', 'A', 'nW', 'b']:
        print '#### fitting %s' %model
        codeml.fit(model)
        str(codeml)

    codeml.unfix_omega(0.8)
    codeml.start_omega(1.2)

    codeml.number_of_categories(4)
    codeml = wrappers.Codeml(seq)
    for model in ['M0', 'M1a', 'M2a', 'M8a', 'M8', 'b']:
        print '#### fitting %s (star topology)' %model
        codeml.fit(model)
        str(codeml)


########################################################################

def test_primer3():
    
    """
    Test primer3 wrapper
    """
    
    print "## Testing wrapper for primer3"

    primer3 = wrappers.Primer3(example_files.Medicago_LYK_region)
    n1, n2 = primer3.find_primers()
    primer3.forward_primers()[1000:] = []
    primer3.reverse_primers()[1000:] = []
    n3 = primer3.find_pairs(200, 600)
    primer3.clean_primer_ends(3)
    primer3.pairs()[1000:] = []
    assert len(primer3.pairs())<=1000
    n4 = primer3.check_pairs()
    assert n4<=n3
    primer3.clean_pair_ends(2)
    primer3.select(100)
    primer3.sort()

    params = {'PRIMER_MAX_SIZE': 25, 'PRIMER_SALT_DIVALENT': 0.2,
              'PRIMER_SALT_MONOVALENT': 40.0,
              'PRIMER_OPT_GC_PERCENT': 52.0, 'PRIMER_MAX_GC': 60.0,
              'PRIMER_OPT_TM': 58.0, 'PRIMER_PICK_ANYWAY': 0,
              'PRIMER_MIN_GC': 35.0, 'PRIMER_GC_CLAMP': 2,
              'PRIMER_OPT_SIZE': 22, 'PRIMER_MIN_SIZE': 19,
              'PRIMER_FIRST_BASE_INDEX': 0, 'PRIMER_MAX_NS_ACCEPTED': 1,
              'PRIMER_NUM_RETURN': 8, 'PRIMER_LIBERAL_BASE': 0,
              'PRIMER_PAIR_MAX_DIFF_TM': 4, 'PRIMER_MIN_TM': 55.0,
              'PRIMER_MAX_TM': 60.0, 'PRIMER_DNTP_CONC': 1.0}

    primer3 = wrappers.Primer3(example_files.Medicago_LYK_region, **params)
    n1, n2 = primer3.find_primers()
    primer3.forward_primers()[1000:] = []
    primer3.reverse_primers()[1000:] = []
    primer3.find_pairs(200, 600)
    primer3.pairs()[1000:] = []
    primer3.clean_primer_ends(3)
    primer3.check_pairs()
    primer3.clean_pair_ends(2)
    primer3.select(100)
    primer3.sort()


########################################################################

def test_all():

    """
    Launch all tests of this module
    """
    
    print "# Testing wrappers module"
    test_primer3()
    test_ms()
    test_blast()
    test_alignment()
    test_phylogeny()

