import subprocess, sys, os, random
from .. import utils, simul, tools, data, wrappers
import example_files


########################################################################

def test_abc(stdout):

    """
    Test ABC commands
    """
    
    stdout.write('## Test ABC commands\n'); stdout.flush()

    stdout.write('### Generating data files\n'); stdout.flush()

    aligns = simul.coalesce(simul.CoalesceParamSet([20,20]), 
                            simul.CoalesceFiniteAlleleMutator(5.), 20)
   
    os.mkdir('tempdir-fas')
    for i,align in enumerate(aligns):
        align = tools.concat([align, data.Align()], 1000-align.ls(), 'A')
        align.write('tempdir-fas/%s.fas' % (str(i+1).rjust(2, '0')), True)
    
    
    try:

        stdout.write('### abc_sample\n')
        stdout.write('#### Checking help pages\n'); stdout.flush()
        
        sys.argv = ['egglib', 'abc_sample']
        utils.execute()
        
        sys.argv = ['egglib', 'abc_sample', 'model?']
        utils.execute()
        
        sys.argv = ['egglib', 'abc_sample', 'stats?']
        utils.execute()

        sys.argv = ['egglib', 'abc_sample', 'prior?']
        utils.execute()
        
    
        # small runs
        
        stdout.write('#### Performing a small run - continuous prior\n'); stdout.flush()
        sys.argv = ['egglib', 'abc_sample', 'dir=tempdir-fas',
                    'ext=fas', 'params=test.txt', 'data=test.out',
                    'prior=%N(0.001;0.1) U(-0.1,4) U(0;0.001)', 'stats=SFS:8',
                    'model=PEMR', 'post=100', 'seeds=5465,78870',
                    'max_threads=4', 'force_positive']
        utils.execute()

        stdout.write('#### Discrete prior\n'); stdout.flush()
        sys.argv = ['egglib', 'abc_sample', 'dir=tempdir-fas',
                    'ext=fas', 'params=test.txt', 'data=test.out',
                    'prior=%0.3 0;0.001 0;0.5\n0.3 0;0.001 0.5;1.\n\
0.25 0.001;0.0015 0;0.5\n0.15 0.001;0.0015 0.5;1', 'stats=SFS:4',
                    'model=PEM', 'post=100', 'max_threads=1']
        utils.execute()
        
        stdout.write('#### Continuous prior from file\n'); stdout.flush()

        f = open('prior1.txt', 'w')
        f.write('F(0.002)\nU(0;0.001)\n')
        f.close()
        sys.argv = ['egglib', 'abc_sample', 'dir=tempdir-fas',
                    'ext=fas', 'params=test-PEM.txt', 'data=test-PEM.out',
                    'prior=prior1.txt', 'stats=SFS:6',
                    'model=PEM', 'post=1000', 'max_threads=1']
        utils.execute()
        
        stdout.write('#### Discrete prior from file\n'); stdout.flush()
        f = open('prior2.txt', 'w')
        f.write('0.3 0;0.001 0;0.5\n0.3 0;0.001 0.5;1.\n\
0.25 0.001;0.0015 0;0.5\n0.15 0.001;0.0015 0.5;1\n')
        f.close()
        sys.argv = ['egglib', 'abc_sample', 'dir=tempdir-fas',
                    'ext=fas', 'params=test.txt', 'data=test.out',
                    'prior=prior2.txt', 'stats=SFS:4',
                    'model=PEM', 'post=100', 'max_threads=1']
        utils.execute()

        stdout.write('#### Model option\n'); stdout.flush()
        sys.argv = ['egglib', 'abc_sample', 'dir=tempdir-fas',
                    'ext=fas', 'params=test-AM2.txt', 'data=test-AM2.out',
                    'prior=%U(0;0.001) U(0,1) U(0.01;5)', 'stats=SFS:6',
                    'model=AM:2', 'post=1000', 'max_threads=1']
        utils.execute()

        stdout.write('#### Add model\n'); stdout.flush()
        f = open('ErrorModel.py', 'w')
        f.write(example_files.CustomModel)
        f.close()
        sys.argv = ['egglib', 'abc_sample', 'dir=tempdir-fas',
                    'ext=fas', 'params=test.txt', 'data=test.out',
                    'prior=%U(0;0.001) U(0;1) U(0;1)', 'stats=JFS:6',
                    'add_model=ErrorModel', 'model=ErrorModel',
                    'post=1000', 'max_threads=4']
        utils.execute()

        stdout.write('#### Restart option\n'); stdout.flush()
        lines = open('test.out').readlines()
        f = open('test.out', 'w')
        for i in lines[:900]: f.write(i)
        f.close()
        sys.argv = ['egglib', 'abc_sample', 'restart=test.txt']
        utils.execute()


        stdout.write('### abc_fit\n'); stdout.flush()
        sys.argv = ['egglib', 'abc_fit', 'input=test.txt', 'quiet',
                    'tolerance=0.5', 'transform=none', 'output=test.fit']
        utils.execute()
        assert os.path.isfile('test.fit')
        sys.argv = ['egglib', 'abc_fit', 'input=test-PEM.txt',
                    'tolerance=0.5', 'transform=log', 'output=test.fit']
        utils.execute()
        sys.argv = ['egglib', 'abc_fit', 'input=test.txt',
                    'tolerance=0.5', 'transform=tan', 'output=test.fit']
        utils.execute()

        stdout.write('### abc_bin\n'); stdout.flush()
        sys.argv = ['egglib', 'abc_bin', 'debug', 'input=test.fit',
                    'bins=5', 'ranges=0:0.01,0:1,0:1', 'output=test.bin']
        utils.execute()

        sys.argv = ['egglib', 'abc_bin', 'input=test.fit',
                    'bins=12', 'output=test.bin']
        utils.execute()
        assert os.path.isfile('test.bin')

        stdout.write('### abc_compare\n'); stdout.flush()
        sys.argv = ['egglib', 'abc_compare', 'tolerance=0.2',
                    'input=test-PEM.txt,test-AM2.txt']
        utils.execute()

        stdout.write('### abc_statsdisc\n'); stdout.flush()
        sys.argv = ['egglib', 'abc_statsdisc', 'input=test.bin',
                    'q=0.667']
        utils.execute()

        stdout.write('### abc_statsmarg\n'); stdout.flush()
        sys.argv = ['egglib', 'abc_statsmarg', 'input=test.fit']
        utils.execute()

        stdout.write('### abc_plot1D and abc_plot2D\n'); stdout.flush()
        sys.argv = ['egglib', 'abc_plot1D', 'input=test.bin',
                    'index=1', 'params=flap', 'root=testplots']
        utils.execute()
        sys.argv = ['egglib', 'abc_plot1D', 'input=test.bin',
                    'params=flap,flip,flop', 'root=testplots']
        utils.execute()
        sys.argv = ['egglib', 'abc_plot1D', 'input=test.bin',
                    'root=testplots']
        utils.execute()
        sys.argv = ['egglib', 'abc_plot1D', 'input=test.bin']
        utils.execute()
        sys.argv = ['egglib', 'abc_plot1D', 'input=test.bin', 'index=3']
        utils.execute()
        sys.argv = ['egglib', 'abc_plot2D', 'input=test.bin',
                    'index1=1', 'index2=2', 'param1=flap', 'param2=flop',
                    'output=hu.png']
        utils.execute()
        sys.argv = ['egglib', 'abc_plot2D', 'input=test.bin',
                    'index1=1', 'index2=2', 'CI']
        utils.execute()

        stdout.write('### abc_psimuls\n'); stdout.flush()
        sys.argv = ['egglib', 'abc_psimuls', 'model=ErrorModel',
                    'add_model=ErrorModel', 'prior=test.bin',
                    'ns=20,20', 'ls=1000', 'nrepets=100',
                    'stats=S,D,K', 'seeds=4155,5514']
        utils.execute()
        sys.argv = ['egglib', 'abc_psimuls', 'model=ErrorModel',
                    'add_model=ErrorModel', 'prior=%F(0.00001) F(2.5) F(0.1)',
                    'ns=20,20', 'ls=1000', 'nrepets=100',
                    'stats=S,D,K', 'seeds=4155,5514', 'debug']
        utils.execute()

    finally:
        for i in os.listdir('tempdir-fas'):
            os.remove(os.path.join('tempdir-fas', i))
        os.rmdir('tempdir-fas')
        
        if os.path.isfile('test.txt'): os.remove('test.txt')
        if os.path.isfile('test.out'): os.remove('test.out')
        if os.path.isfile('test-PEM.txt'): os.remove('test-PEM.txt')
        if os.path.isfile('test-PEM.out'): os.remove('test-PEM.out')
        if os.path.isfile('test-AM2.txt'): os.remove('test-AM2.txt')
        if os.path.isfile('test-AM2.out'): os.remove('test-AM2.out')
        if os.path.isfile('test.fit'): os.remove('test.fit')
        if os.path.isfile('test.bin'): os.remove('test.bin')
        if os.path.isfile('test.prior'): os.remove('test.prior')
        if os.path.isfile('prior1.txt'): os.remove('prior1.txt')
        if os.path.isfile('prior2.txt'): os.remove('prior2.txt')
        if os.path.isfile('ErrorModel.py'): os.remove('ErrorModel.py')
        if os.path.isfile('ErrorModel.pyc'): os.remove('ErrorModel.pyc')
        for i in os.listdir('.'):
            if i[-4:]=='.png': os.remove(i)

########################################################################

def test_analyzer(stdout):

    """
    Test analyzer command
    """
    
    stdout.write( '## Test analyzer command\n'); stdout.flush()

    f = open('ms_output.txt', 'w')
    f.write(example_files.ms_output)
    f.close()
    
    utils.execute('analyzer', input='ms_output.txt', config=40)
    utils.execute('analyzer', input='ms_output.txt', config='20,20')
    utils.execute('analyzer', input='ms_output.txt', config=40, stats='D,Z')
    backup = sys.stdin
    sys.stdin = open('ms_output.txt')
    utils.execute('analyzer', config='20,20', stats='D,Z,Snn')
    sys.stdin.close()
    sys.stdin = backup
    
    os.remove('ms_output.txt')

########################################################################

def test_blastgb(stdout):

    """
    Test blastgb command
    """
    
    stdout.write( '## Test blastgb command\n'); stdout.flush()

    f = open('gb_file.txt', 'w')
    f.write(example_files.gb_poplar)
    f.close()

    f = open('db_file.txt', 'w')
    f.write(example_files.db_dmi3)
    f.close()

    utils.execute('blastgb', input='gb_file.txt',
                  output='blastgb_output.txt', db='db_file.txt',
                  evalue=0.01, nresults=1)

    utils.execute('blastgb', input='gb_file.txt', nresults=10,
                  output='blastgb_output.txt', db='db_file.txt')

    db = tools.translate(data.Container(string=example_files.db_dmi3))
    db.write('db_file.txt')

    utils.execute('blastgb', 'prot', input='gb_file.txt',
                  output='blastgb_output.txt', db='db_file.txt',
                  evalue=0.01, nresults=1)
    gb = data.GenBank('blastgb_output.txt')

    db = wrappers.BLASTdb(db, 'prot')
    utils.execute('blastgb', 'prot', 'formatted', input='gb_file.txt',
                  output='blastgb_output.txt', db=db.path())

    os.remove('gb_file.txt')
    os.remove('db_file.txt')
    os.remove('blastgb_output.txt')

########################################################################

def test_c_commands(stdout):

    """
    Test commands with name in c
    """
    
    stdout.write( '## Test commands with name in c\n'); stdout.flush()

    ##########

    stdout.write( '### clean_seq\n'); stdout.flush()

    a = simul.coalesce(simul.CoalesceParamSet(40), simul.CoalesceFiniteAlleleMutator(4), 1)[0]
    a = tools.concat([a, data.Align()], 1000-a.ls(), 'A')
    for i in range(10):
        x = random.randint(0,39)
        y = random.randint(0,999)
        a.set(x, y, 'Z')
    a.write('a.fas')
    
    utils.execute('clean_seq', input='a.fas', output='b.fas')
    utils.execute('clean_seq', input='a.fas', output='c.fas', chars='CZ')

    ##########
    
    stdout.write( '### clean_tree\n'); stdout.flush()
    
    f = open('a.tre', 'w')
    f.write(example_files.tree)
    f.close()
    
    utils.execute('clean_tree', input='a.tre', output='b.tre')
    utils.execute('clean_tree', 'keep_labels', input='a.tre', output='b.tre')
    utils.execute('clean_tree', 'keep_brlens', input='a.tre', output='b.tre')
    utils.execute('clean_tree', 'keep_brlens', 'keep_labels', input='a.tre', output='b.tre')
    
    ##########

    stdout.write( '### codalign\n'); stdout.flush()

    aat = data.Container(string=example_files.aat)
    for i in aat: i.sequence = tools.longest_orf(i.sequence)
    aat.write('a.fas')

    utils.execute('codalign', 'quiet', 'debug', input='a.fas', output='b.fas')
    utils.execute('codalign', 'muscle', 'quiet', input='a.fas', output='b.fas')
    aln = data.Align('b.fas')
    aln = tools.translate(aln)
    aln.write('c.fas')

    utils.execute('codalign', input='a.fas', output='b.fas', prot='c.fas')

    ##########

    stdout.write( '### concat\n'); stdout.flush()

    nork = data.Align(string=example_files.nork)
    nork.extract(0, 200).write('a.fas')
    nork.extract(200, 600).write('b.fas')
    nork.extract(600, 1176).write('c.fas')

    utils.execute('concat', 'debug', input='a.fas,b.fas,c.fas',
                  output='d.fas', spacer='50,150', character='N')

    assert len(data.Align('d.fas'))==len(nork)

    a = nork.extract(0, 200)
    b = nork.extract(200, 600)
    c = nork.extract(600, 1176)
    for i in a: i.name += '$1'
    for i in c: i.name += '$2'
    for i in b: i.name += '$3'
    a.write('a.fas')
    b.write('b.fas')
    c.write('c.fas')

    utils.execute('concat', 'debug', input='a.fas,b.fas,c.fas',
                  output='d.fas', character='N', sep='$')

    assert len(data.Align('d.fas'))==len(nork)

    for i in range(len(nork)):
        a.name(i, str(i+1).rjust(2,'0') + ''.join(random.sample('ABCD1234567890', 6)))
        b.name(i, str(i+1).rjust(2,'0') + ''.join(random.sample('ABCD1234567890', 6)))
        c.name(i, str(i+1).rjust(2,'0') + ''.join(random.sample('ABCD1234567890', 6)))
    a.write('a.fas')
    a.write('b.fas')
    a.write('c.fas')

    utils.execute('concat', 'debug', input='a.fas,b.fas,c.fas',
                  output='d.fas', spacer=100, character='N', len=2)

    assert len(data.Align('d.fas'))==len(nork)

    for i in b: i.name = i.name[:2]
    b.write('b.fas')

    utils.execute('concat', 'partial', input='a.fas,b.fas,c.fas',
                  output='d.fas')

    assert len(data.Align('d.fas'))==len(nork)

    os.remove('d.fas')
    utils.execute('concat', input='*.fas', output='d.fas')

    ##########

    stdout.write( '### concatgb\n'); stdout.flush()

    gb = data.GenBank(string=example_files.gb_poplar)
    a = gb.extract(0,9999)
    b = gb.extract(20000,29999)
    a.write('a.gb')
    b.write('b.gb')
    
    utils.execute('concatgb', file1='a.gb', file2='b.gb',
                  output='c.gb')

    utils.execute('concatgb', file1='a.gb', file2='b.gb',
                  output='c.gb', spacer=10000)

    utils.execute('concatgb', file1='a.gb', file2='b.gb',
                  output='c.gb', spacer=10000, character='?')

    ##########

    stdout.write( '### consensus\n'); stdout.flush()

    a = simul.coalesce(simul.CoalesceParamSet(None, 20),
                       simul.CoalesceFiniteAlleleMutator(4), 1)[0]
    a = tools.concat([a, data.Align()], 1000-a.ls(), 'A')
    for i in range(20):
        a.name(i*2, str(i).rjust(2, '0') + '_1')
        a.name(i*2+1, str(i).rjust(2, '0') + '_2')
    a.write('a.fas')
    
    utils.execute('consensus', input='a.fas', output='b.fas')
    assert len(data.Align('b.fas'))==20

    for i in range(20):
        a.name(i*2, str(i).rjust(2, '0') + '%uuY')
        a.name(i*2+1, str(i).rjust(2, '0') + '%Yuu')
    a.write('a.fas')

    utils.execute('consensus', input='a.fas', output='b.fas')
    assert len(data.Align('b.fas'))==40
    
    utils.execute('consensus', input='a.fas', output='b.fas',
                  separator='%')
    assert len(data.Align('b.fas'))==20

    utils.execute('consensus', 'conservative', input='a.fas',
                  output='b.fas', separator='%', missing='N',
                  inconsistency='?')

    #########
    
    stdout.write( '### cprimers\n'); stdout.flush()
    
    a = data.Container(string=example_files.aat)
    a = wrappers.clustal(a)
    a.write('a.fas')
    
    utils.execute('cprimers', 'debug', input='a.fas')

    f = open('a.gb', 'w')
    f.write("""LOCUS       None
FEATURES             Location/Qualifiers
     CDS             join(1..203,264..443,624..1017)
BASE COUNT      226 a    128 c    127 g    296 t
ORIGIN
        1 tgaaaacgtg agctactgag cttcagaagg agggcaaaaa ggtatagtca aatttctttt
       61 ctggattata tagagaaagt tatatctgct tgacttgccc ttataggttg gttgaagtat
      121 taaatctgac ctgatcattg aaatcctttc actgaatgtg tattatattt gctatacaat
      181 ctataatgta acttggaaga ctgnnnnnnn nnnnnnnnnn nnnnnnnnnn nnnnnnnnnn
      241 nnnnnnnnnn nnnnnnnnnn nnnttattgc caactattca tagattctca tgcatttttt
      301 gtttgatcca caagtgtgac cctgttgaac aagttataga ctttagcatt gccatcatgg
      361 taatagtagt atgttacact ttgggtaagt tgtatagcta atttttacat ggaaaaatct
      421 tcatcaatgc acatatcagt tttnnnnnnn nnnnnnnnnn nnnnnnnnnn nnnnnnnnnn
      481 nnnnnnnnnn nnnnnnnnnn nnnnnnnnnn nnnnnnnnnn nnnnnnnnnn nnnnnnnnnn
      541 nnnnnnnnnn nnnnnnnnnn nnnnnnnnnn nnnnnnnnnn nnnnnnnnnn nnnnnnnnnn
      601 nnnnnnnnnn nnnnnnnnnn nnntgccata ctttccagct gcttgaaatt taacgagttt
      661 atatttatat tcttgtcgag ttgtatatcc caagatttga ttgcttaaaa ggctgaatat
      721 atgtttatgt aactataggt cttattgtgt atgctattag tttacctctt ctcatagtta
      781 agtttattca tgtttatatt cactaatgtt ggaaacccac atgctttggg acagaaacca
      841 ctgaccttcc ctcgccaggt tcggcacaac taccatattt ctacttttac taattgcccc
      901 ttgtaatttt cgtaatatga tatgctaaac atactccaca tctaattaaa ttgttagagt
      961 tcaaagtatt ttatttgttg attattctat gtggtccata aagctaacct atagatt
//
""")
    f.close()

    utils.execute('cprimers', 'debug', input='a.fas', gbin='a.gb')
    utils.execute('cprimers', input='a.fas', ndeg=6)
    utils.execute('cprimers', input='a.fas', ndeg=0)
    utils.execute('cprimers', input='a.fas', liml=50)
    utils.execute('cprimers', input='a.fas', liml=700)
    utils.execute('cprimers', input='a.fas', gbin='a.gb', liml=50, limr=700)
    utils.execute('cprimers', 'debug', 'no_check', input='a.fas', clean_ends=4, nseq=4)



    # cleaning

    os.remove('a.fas')
    os.remove('b.fas')
    os.remove('c.fas')
    os.remove('d.fas')
    os.remove('a.tre')
    os.remove('b.tre')
    os.remove('a.gb')
    os.remove('b.gb')
    os.remove('c.gb')
    os.remove('cprimers.list.txt')
    os.remove('cprimers.pairs.txt')
    os.remove('cprimers.primers.gb')

########################################################################

def test_extract_family(stdout):

    """
    Test extract, extract_gb and family commands
    """
    
    stdout.write('## Test extract, extract_gb and family commands\n'); stdout.flush()

    ##########

    stdout.write('### extract\n'); stdout.flush()

    nork = data.Align(string=example_files.nork)
    nork.write('a.fas')
    
    utils.execute('extract', input='a.fas', output='b.fas',
                                         ranges='1,10,100,1000,599,466')

    utils.execute('extract', input='a.fas', output='b.fas',
                                         ranges='1000,800-1001')

    utils.execute('extract', input='a.fas', output='b.fas',
                                         ranges='4-10,7,8,9,15,20-50,10,10')

    ##########

    stdout.write('### extract_clade\n'); stdout.flush()

    tree, lk = wrappers.phyml(nork, boot=10)
    tree.write('a.tre')
    
    utils.execute('extract_clade', sequences='a.fas', tree='a.tre',
                    output='b.fas', names='Sesbania,Lotus')
    utils.execute('extract_clade', sequences='a.fas', tree='a.tre',
                   output='b.fas', names='Sesbania,Lotus', threshold=10)
    utils.execute('extract_clade', sequences='a.fas', tree='a.tre',
                    output='b.fas', names='Sesbania,Lotus', minimum=4)
    utils.execute('extract_clade', 'monophyletic', sequences='a.fas',
                   tree='a.tre', output='b.fas', names='Sesbania,Lotus')
    utils.execute('extract_clade', 'exact', sequences='a.fas',
        tree='a.tre', output='b.fas', names='Sesbania,Lotus', minimum=5)

    ##########

    stdout.write('### family\n'); stdout.flush()

    b = data.Container(string=example_files.cds_medicago)
    a = b.slice(0, 2)
    a.write('a.fas')
    b.write('b.fas')
    
    utils.execute('family', 'debug', input='a.fas', target='b.fas',
                    output='c.fas', mode='blastn')
    utils.execute('family', 'debug', input='a.fas', target='b.fas',
                    output='c.fas', mode='blastn', evalue=10**-20)
    assert len(data.Container('c.fas'))==2
    utils.execute('family', 'debug', input='a.fas', target='b.fas',
                    output='c.fas', mode='tblastx', evalue=10**-20)
    assert len(data.Container('c.fas'))==2

    tools.translate(b).write('b.fas')
    utils.execute('family', 'debug', input='a.fas', target='b.fas',
                    output='c.fas', mode='blastx', evalue=10**-20)
    assert len(data.Container('c.fas'))==2

    tools.translate(a).write('a.fas')
    utils.execute('family', 'debug', input='a.fas', target='b.fas',
                    output='c.fas', mode='blastp', evalue=10**-20)
    assert len(data.Container('c.fas'))==2

    b.write('b.fas')
    utils.execute('family', 'debug', input='a.fas', target='b.fas',
                    output='c.fas', mode='tblastn', evalue=10**-20)
    assert len(data.Container('c.fas'))==2


    ##########
    
    os.remove('a.fas')
    os.remove('b.fas')
    os.remove('c.fas')
    os.remove('a.tre')

########################################################################

def test_convert(stdout):

    """
    Test convert commands
    """
    
    stdout.write('## Test convert commands\n'); stdout.flush()

    nork = data.Align(string=example_files.nork)
    nork.write('a.fas')

    ##########

    stdout.write('### fasta2mase\n'); stdout.flush()
    utils.execute('fasta2mase', input='a.fas', output='a.mase')

    #########

    stdout.write('### fasta2nexus\n'); stdout.flush()
    utils.execute('fasta2nexus', input='a.fas', output='a.nex', type='nucl')
    tools.translate(nork).write('b.fas')
    utils.execute('fasta2nexus', input='b.fas', output='b.nex', type='prot')

    #########

    stdout.write('### fasta2phyml\n'); stdout.flush()
    utils.execute('fasta2phyml', input='a.fas', output='a.phyml')

    #########

    open('a.txt', 'w').write('>Medicago_LYK_region\n' + 
                                    example_files.Medicago_LYK_region)
    open('b.txt', 'w').write(example_files.Medicago_LYK_region_ann)

    stdout.write('### fg2gb\n'); stdout.flush()
    utils.execute('fg2gb', 'debug',  seq='a.txt', ann='b.txt', output='c.txt')

    #########

    stdout.write('### gb2fas\n'); stdout.flush()
    utils.execute('gb2fas', input='c.txt,c.txt', output='a.fas')

    #########
    
    os.remove('a.txt')
    os.remove('b.txt')
    os.remove('c.txt')
    os.remove('a.fas')
    os.remove('b.fas')
    os.remove('a.nex')
    os.remove('b.nex')
    os.remove('a.mase')
    os.remove('a.phyml')

########################################################################

def test_i_p_commands(stdout):

    """
    Test commands with name in i through p
    """
    
    stdout.write( '## Test commands with name in i through p\n'); stdout.flush()

    ##########

    stdout.write( '### infos\n'); stdout.flush()

    nork = data.Align(string=example_files.nork)
    nork.write('a.fas')
   
    utils.execute('infos', input='a.fas')


    ##########

    stdout.write( '### interLD\n'); stdout.flush()

    nork = data.Align(string=example_files.nork)
    nork.write('a.fas')
   
    utils.execute('interLD', align1='a.fas', align2='a.fas')

    a = nork.extract(0, 500)
    b = nork.extract(500, 1176)
    a.write('a.fas')
    b.write('b.fas')

    utils.execute('interLD', align1='a.fas', align2='b.fas')
    utils.execute('interLD', align1='a.fas', align2='b.fas', permus=100, output='a.txt')
    
    ##########

    stdout.write( '### matcher\n'); stdout.flush()

    open('a.fas', 'w').write('>a\n' + example_files.Medicago_LYK_region)
    open('b.fas', 'w').write('>a\n' + example_files.Medicago_LYK_region[1000:2000])
   
    utils.execute('matcher', long='a.fas', short='b.fas', output='a.gb',
                  mode='blastn', evalue=10**-20)
    data.GenBank('a.gb')

    open('b.fas', 'w').write('>a\n' + 
        tools.translate(example_files.Medicago_LYK_region[1000:2000]))
    utils.execute('matcher', long='a.fas', short='b.fas', output='a.gb',
                  mode='blastx', evalue=10**-20)
   
    ##########

    stdout.write( '### names\n'); stdout.flush()
    data.Container(string=example_files.cds_medicago).write('a.fas')
    utils.execute('names', input='a.fas')
    utils.execute('names', 'wrap', input='a.fas')

    ##########
    
    stdout.write('### phyml\n'); stdout.flush()
    data.Align(string=example_files.nork).write('a.fas')
    utils.execute('phyml', 'quiet', input='a.fas', output='a.tre')

    stdout.write('#### search modes\n'); stdout.flush()
    utils.execute('phyml', 'quiet', input='a.fas', output='a.tre', search='SPR')
    utils.execute('phyml', 'quiet', input='a.fas', output='a.tre', search='BEST')

    stdout.write('#### rates\n'); stdout.flush()
    utils.execute('phyml', 'quiet', input='a.fas', output='a.tre', rates=6)

    stdout.write('#### bootstrap\n'); stdout.flush()
    utils.execute('phyml', 'quiet', input='a.fas', output='a.tre', boot=10, boot_trees='b.tre')

    stdout.write('#### models\n'); stdout.flush()
    utils.execute('phyml', 'quiet', input='a.fas', output='a.tre', model='JC69')
    tools.translate(data.Align(string=example_files.nork)).write('a.fas')
    utils.execute('phyml', 'quiet', input='a.fas', output='a.tre', model='LG')
    utils.execute('phyml', 'quiet', input='a.fas', output='a.tre', model='JTT')
    
    ###########
    
    os.remove('a.fas')
    os.remove('b.fas')
    os.remove('interLD.txt')
    os.remove('a.txt')
    os.remove('a.gb')
    os.remove('a.tre')
    os.remove('b.tre')

########################################################################

def test_r_s_commands(stdout):

    """
    Test commands with name in r and s
    """
    
    stdout.write( '## Test commands with name in r and s\n'); stdout.flush()

    ##########

    stdout.write( '### rename\n'); stdout.flush()

    nork = data.Align(string=example_files.nork)
    nork.write('a.fas')
    f = open('a.txt', 'w')
    f.write("""truncatul\tMedicago truncatula
Astragalu\tAstragalus sinicus
Lotus\tLotus japonicus
Melilotus\tMelilotus alba
Pisum\tPisum sativa
Sesbania\tSesbania rostrata
Vicia\tVicia hirsuta""")
    f.close()
   
    utils.execute('rename', input='a.fas', output='b.fas', list='a.txt')

    ##########

    stdout.write( '### reroot\n'); stdout.flush()
    f = open('a.tre', 'w')
    f.write(example_files.tree)
    f.close()
    utils.execute('reroot', input='a.tre', output='b.tre',
                  outgroup='MtAC146553_34,Ljchr6.CM1829.110.nc,Glyma15g10660.1,Glyma13g28390.1,Glyma07g38310.1,AtHP1')
    f = open('a.txt', 'w')
    f.write("""MtAC146553_34
Ljchr6.CM1829.110.nc
Glyma15g10660.1
Glyma13g28390.1
Glyma07g38310.1
AtHP1""")
    f.close()
    utils.execute('reroot', input='a.tre', output='b.tre',
                  outgroup='@a.txt')

    ##########

    stdout.write( '### select\n'); stdout.flush()
    utils.execute('select', input='a.fas', output='b.fas',
                  list='truncatul,Lotus,Pisum')
    f = open('a.txt', 'w')
    f.write('Lotus\nSesbania\nVicia\n')
    f.close()
    utils.execute('select', input='a.fas', output='b.fas',
                  list='@a.txt')

    ##########

    stdout.write( '### sprimers\n'); stdout.flush()

    nork = data.Align(string=example_files.nork)
    nork.group(0, 1)
    nork.group(1, 1)
    nork.group(2, 1)
    nork.write('a.fas', True)

    utils.execute('sprimers', 'debug', 'selection', 'prefer_end',
        input='a.fas',  output='a.txt',
        sizemin=200, sizemax=800, minTm=59, maxTm=65, optTm=63,
        minGc=40, optGc=55, maxGc=75, numAmb=1,
        filter1=2000, filter2=50, threshold1=0.80, threshold2=0.80,
        show=5)

    ##########

    stdout.write( '### staden2fasta\n'); stdout.flush()

    f = open('a.sta', 'w')
    f.write(example_files.staden)
    f.close()

    utils.execute('staden2fasta', 'debug', input='a.sta', output='a.fas')
    utils.execute('staden2fasta', 'debug', input='a.sta', output='a.fas', consensus='remove')
    utils.execute('staden2fasta', 'debug', input='a.sta', output='a.fas', consensus='keep')
    utils.execute('staden2fasta', 'debug', input='a.sta', output='a.fas', consensus='only')

    ##########
    
    os.remove('a.fas')
    os.remove('b.fas')
    os.remove('a.txt')
    os.remove('a.tre')
    os.remove('b.tre')
    os.remove('a.sta')


########################################################################

def test_t_w_commands(stdout):

    """
    Test commands with name in t to w
    """
    
    stdout.write( '## Test commands with name in t and w\n'); stdout.flush()

    ##########

    stdout.write( '### translate\n'); stdout.flush()

    nork = data.Align(string=example_files.nork)
    nork.write('a.fas')
   
    utils.execute('translate', 'codes')
    utils.execute('translate', input='a.fas', output='b.fas')
    utils.execute('translate', input='a.fas', output='b.fas', code=1)
    utils.execute('translate', input='a.fas', output='b.fas', code=2)
    utils.execute('translate', input='a.fas', output='b.fas', code=16)
    utils.execute('translate', input='a.fas', output='b.fas', code=21)


    ##########

    stdout.write( '### truncate\n'); stdout.flush()

    nork = data.Align(string=example_files.nork)
    nork.name(0, 'truncatul_0')
    nork.name(1, 'Astra_1')
    nork.name(2, 'Lotus')
    nork.name(3, 'Melilo_3')
    nork.name(4, 'Pisum')
    nork.name(5, 'Sesbania_4')
    nork.name(1, 'Vi_5')
    nork.write('a.fas')

    utils.execute('truncate', input='a.fas', output='b.fas', length=4)
    utils.execute('truncate', input='a.fas', output='b.fas', separator='_')
    utils.execute('truncate', input='a.fas', output='b.fas', length=4,
                                                              separator='_')
    ##########

    stdout.write( '### ungap\n'); stdout.flush()
    utils.execute('ungap', input='a.fas', output='b.fas', threshold=0)
    utils.execute('ungap', input='a.fas', output='b.fas', threshold=0.5)
    utils.execute('ungap', input='a.fas', output='b.fas', threshold=1)
    utils.execute('ungap', 'all', input='a.fas', output='b.fas')

    ##########
   
    stdout.write( '### winphyml\n'); stdout.flush()
    nork = data.Align(string=example_files.nork)
    nork.write('a.fas')
    tre1, junk = wrappers.phyml(nork)
    tre2 = wrappers.nj(nork)
    f = open('a.tre', 'w')
    f.write(str(tre1))
    f.write(str(tre2))
    f.close()
    utils.execute('winphyml', input='a.fas', trees='a.tre',
                  output='a.txt', wsize=200, wstep=150)

    ##########
    
    os.remove('a.fas')
    os.remove('b.fas')
    os.remove('a.txt')
    os.remove('a.tre')



########################################################################

def test_all():

    """
    Launch all tests of this module
    """
    
    backup1 = list(sys.argv)
    backup2 = sys.stdout
    
    
    print "# Testing utils module"

    f = open('redirection.txt', 'w')
    sys.stdout = f
    sys.argv = ['egglib']
    utils.execute()
    test_abc(backup2)

    utils.execute('ungap')
    test_analyzer(backup2)
    test_blastgb(backup2)
    test_c_commands(backup2)
    test_extract_family(backup2)
    test_convert(backup2)
    test_i_p_commands(backup2)
    test_r_s_commands(backup2)
    test_t_w_commands(backup2)
    
    sys.argv = backup1
    os.remove('redirection.txt')


