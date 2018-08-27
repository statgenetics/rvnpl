import os, random, time
from .. import tools, data
import example_files


########################################################################

def test_SpecificDataFormat():

    """
    Test functions of the Specific Data Format section
    """
    
    print '## Test tools functions for specific data formats'
    
    try:
        f = open('test_file1.txt', 'w')
        f.write(example_files.clustal_aln)
        f.close()
        aln = tools.aln2fas('test_file1.txt')
        assert len(aln)==15

        f = open('test_file2.txt', 'w')
        f.write(example_files.staden)
        f.close()
        align = tools.staden('test_file2.txt')
        ns = align.ns()
        ls = align.ls()
        align = tools.staden(string=example_files.staden)
        assert ns == align.ns()
        assert ls == align.ls()
        align = tools.staden('test_file2.txt', delete_consensus=False)
        assert ns+1 == align.ns()
        assert ls == align.ls()
        align = tools.staden(string=example_files.staden, delete_consensus=False)
        assert ns+1 == align.ns()
        assert ls == align.ls()
        
        f = open('test_file3.txt', 'w')
        f.write(example_files.fgenesh_2_5)
        f.close()
        feats = tools.get_fgenesh('test_file3.txt')

        f = open('test_file4.txt', 'w')
        f.write(example_files.genalys_snp)
        f.close()
        align = tools.genalys2fasta('test_file4.txt')

        mase = tools.Mase()
        assert mase.header==''
        assert mase.species==None
        assert len(mase.align)==0
        mase.align.append('name1', 'AAAAAAAAAAAAA', 0)
        mase.align.append('name2', 'AAAAAAAAAAACC', 0)
        assert str(mase).strip()==';;'
        mase = tools.Mase(mase.align)
        assert str(mase).strip()!=';;'
        f = open('test_file5.txt', 'w')
        f.write(example_files.mase)
        f.close()
        mase = tools.Mase('test_file5.txt')
        assert len(mase)==22
        assert mase.species=='Sorghum_bicolor'
        for header, name, sequence in mase: pass
        list(mase)


        
    finally:
        if os.path.isfile('test_file1.txt'): os.remove('test_file1.txt')
        if os.path.isfile('test_file2.txt'): os.remove('test_file2.txt')
        if os.path.isfile('test_file3.txt'): os.remove('test_file3.txt')
        if os.path.isfile('test_file4.txt'): os.remove('test_file4.txt')
        if os.path.isfile('test_file5.txt'): os.remove('test_file5.txt')


########################################################################

def test_DataAnalysis():

    """
    Test functions of the Data Analysis section
    """
    
    print '## Test tools functions for data analysis'

    align1 = data.Align.create([('name1', 'TAAA'), ('name2', 'TAAA'),
                                ('name3', 'TACC'), ('name4', 'TAAC'),
                                ('name5', 'AACC'), ('name6', 'AAAC'),
                                ('name7', 'AAAC'), ('name8', 'AAAC'),
                                ('name9', 'AAAC'), ('name0', 'AAAC'),])
    align2 = data.Align.create([('name1', 'ACCA'), ('name2', 'ACAA'),
                                ('name3', 'ACAC'), ('name4', 'AAAC'),
                                ('name5', 'AAAA'), ('name6', 'AACA'),
                                ('name7', 'AACA'), ('name8', 'AACA'),
                                ('name9', 'AACA'), ('name0', 'AACA')])
    
    n1,n2,S1,S2,K1,K2,D1,Dp1 = tools.LD(align1, align2, False)
    n1,n2,S1,S2,K1,K2,D2,Dp2 = tools.LD(align1, align2, True)


########################################################################

def test_SequenceManipulationTools():

    """
    Test functions of the Sequence Manipulation Tools section
    """
    
    print '## Test tools functions for sequence manipulation'

    nucseq = data.Container.create([('name1', 'ATGGATCATTGGGTA'),
         ('name2', 'ATCGACCATGTA'), ('name3', 'ATGGAACATTGCGTA'),
         ('name4', 'ATGCATGATCTTTGGGTA')])
    protseq = data.Align(string=""">name1
M-DHWV
>name2
I-DH-V
>name3
M-EHCV
>name4
MHDLWV""")
    align = tools.backalign(nucseq, protseq, 1)
    align = tools.backalign(nucseq, protseq, 2)

    align1 = data.Align.create([('a', 'ACGTACGT', 0)
                               ,('b', 'ACGTACGT', 0)
                               ,('c', 'ACGT----', 1)
                               ,('d', 'ACCTACCT', 1)])
    align2 = data.Align.create([('a', 'GAAAAAATATA', 0)
                               ,('b', 'GAAAATTTTTA', 0)
                               ,('d', 'GAAA---TTTA', 1)])
    align3 = data.Align.create([('a', 'TATCGGATA', 0)
                               ,('b', 'T---GAGCA', 0)
                               ,('c', 'TATATATCA', 0)])
    res = tools.concat([align1,align2,align3], spacer=0, ch='?', groupCheck=False)
    assert len(res)==4
    assert res.ls()==28
    res = tools.concat([align1,align2,align3], spacer=100, ch='?', groupCheck=False)
    assert res.ls()==228
    res = tools.concat([align1,align2,align3], spacer=50, ch='N', groupCheck=False)
    assert res.ls()==128
    align3.group(2, 1)
    res = tools.concat([align1,align2,align3], spacer=0, ch='?', groupCheck=True)
    res = tools.concat([align1,align2,align3], spacer=100, ch='?', groupCheck=True)
    res = tools.concat([align1,align2,align3], spacer=50, ch='N', groupCheck=True)

    tools.longest_orf(example_files.Medicago_LYK_region)
    tools.longest_orf(example_files.Medicago_LYK_region, full=True)
    tools.longest_orf(example_files.Medicago_LYK_region, all=True)
    tools.longest_orf(example_files.Medicago_LYK_region, mini=10)
    s = "TTATAATAACTAGATCTTTCCTCAGTAGTATAATTATAAATGCATTCAGAACTTTGATAATAATAATAATAATAATAA"
    tools.longest_orf(s)
    tools.translate(tools.longest_orf(s))
    tools.longest_orf(s, full=True)
    tools.translate(tools.longest_orf(s, full=True))
    tools.longest_orf(s, all=True, mini=0)
    tools.rc(example_files.Medicago_LYK_region)[-3:]

    s = 'ATGGATCTTGTTATTCTTTAA'
    tools.translate(s)
    c = data.Container(string=example_files.cds_medicago)
    d = tools.translate(c)
    d = tools.translate(c, strip=True)
    assert tools.translate('TTA')=='L'
    assert tools.translate('TTA', code=23)=='*'
    
    align = data.Align(string=""">
AAAAAA
>
A-A-A-
>
A-A-A-
>
AAA-A-
>@999
AAAAAA""", groups=True)
    assert tools.ungap(align, 0.1).ls() == 3
    assert tools.ungap(align, 0.5).ls() == 4
    assert tools.ungap(align, 0.9).ls() == 6
    assert tools.ungap(align, 0.1, False).ls() == 3
    assert tools.ungap(align, 0.5, False).ls() == 3
    assert tools.ungap(align, 0.9, False).ls() == 6
    
    tools.GeneticCodes.codes()
    tools.GeneticCodes.index('The Ciliate, Dasycladacean and Hexamita Nuclear Code')
    tools.GeneticCodes.is_start('ATG', 1)
    tools.GeneticCodes.is_start('ATG', 2)
    tools.GeneticCodes.is_start('ATG', 3)
    tools.GeneticCodes.is_start('ATG', 4)
    tools.GeneticCodes.is_start('ATG', 23)
    tools.GeneticCodes.is_start('TTA', 1)
    tools.GeneticCodes.is_start('TTA', 2)
    tools.GeneticCodes.is_start('TTA', 3)
    tools.GeneticCodes.is_start('TTA', 4)
    tools.GeneticCodes.is_start('TTA', 23)
    tools.GeneticCodes.translate('ATG', 1)
    tools.GeneticCodes.translate('ATG', 2)
    tools.GeneticCodes.translate('ATG', 3)
    tools.GeneticCodes.translate('ATG', 4)
    tools.GeneticCodes.translate('ATG', 23)
    tools.GeneticCodes.translate('TTA', 1)
    tools.GeneticCodes.translate('TTA', 2)
    tools.GeneticCodes.translate('TTA', 3)
    tools.GeneticCodes.translate('TTA', 4)
    tools.GeneticCodes.translate('TTA', 23)
    tools.GeneticCodes.translate('TAG', 1)
    tools.GeneticCodes.translate('TAG', 9)


########################################################################

def test_SequenceComparison():

    """
    Test functions of the Sequence Comparison section
    """
    
    print '## Test tools functions for sequence comparison'

    assert tools.compare('AGCCGTGCCAGT', 'AGCCGTGCCAGT') == True
    assert tools.compare('AGCCGTGCCAGT', 'AGCCGTGCCAGTA') == False
    assert tools.compare('AGCCG-GCCAGT', 'AGCCGTGCCAGT') == False
    assert tools.compare('AGCCG?GCCAGT', 'AGCCGTGCCAGT') == True
    assert tools.compare('AGCCGNGCCAGT', 'AGCCGTGCCAGT') == True
    assert tools.compare('MGCCGNGCCAGT', 'AGCCGTGCCAGT') == True
    assert tools.compare('MGCCGNGCCAGT', 'NGCCGTGCCAGT') == True
    assert tools.compare('MGCCGNGCCAGT', 'BGCCGTGCCAGT') == True
    assert tools.compare('MGCCGNGCCAGT', 'KGCCGTGCCAGT') == False
    assert tools.compare('A', 'H')                       == True
    assert tools.compare('AA', 'HB')                     == False
    assert tools.compare('MGCCGNGCCAGT', 'VGCCGTGCCAGT') == True
    assert len(tools.motifs('AAAAGATCTCCAAAAAAAGATATCCAAAAAAAA', 'GATCTCC'))==1
    assert len(tools.motifs('AAAAGATCTCCAAAAAAAGATATCCAAAAAAAA', 'GATCTCC', 1))==2
    assert len(tools.motifs('AAAAGAGAAAAAACTATAGAAAAAAAAAAA', 'CTATAG'))==2
    assert len(tools.motifs('AAAAGAGAAAAAACTATAGAAAAAAAAAAA', 'CTATAG', reverse=False))==1
    assert tools.locate('AAAACCAGTCAAAAA', 'CCAGT')  == 4
    assert tools.locate('AAAACCAGTCAAAAA', 'CCTGT')  == None
    assert tools.locate('AAAACCAGTCAAAAA', 'CCMGT')  == 4
    assert tools.locate('AAAACCAGTCAAAAA', 'CCAGT', 4, 10) == 4
    assert tools.locate('AAAACCAGTCAAAAA', 'CCAGT', 5, 10) == None
    assert tools.locate('AAAACCAGTCAAAAA', 'CCAGT', 2, 7) == None


########################################################################

def test_Other():

    """
    Test functions of the Tools not (directly) Related to Sequence Data section
    """
    
    print '## Test tools functions not directly related to sequence data'

    for i in range(1,101): tools.chisquare(i)
    x = [random.random() for i in range(100)]
    y = [random.random() for i in range(100)]
    a, b, c = tools.correl(x, y)
    d, e, f = tools.correl(sorted(x), sorted(y))
    g, h, i = tools.correl(x, x)
    assert g>=d>=a and h>=e>=b and i>=min(1,f)>=c
    tools.ranges([random.randint(0,1000) for i in range(100)])
    frame = tools.ReadingFrame([(100, 120, 2), (200, 220, None), (300, 330, 3)])
    frame.codons()
    assert frame.codon(200) == (120, 200, 201)
    assert frame.exon(200) == 1 and frame.codon(200) == (120, 200, 201)
    assert frame.exon(300) == 2 and frame.codon(300) == None
    assert frame.exon(150) == None and frame.codon(150) == None
    
    updater = tools.Updater(1000)
    for i in range(1000):
        time.sleep(0.001)
        updater.refresh()
    updater.close()
    updater = tools.Updater(1000)
    assert updater.closed==False
    for i in range(1000):
        time.sleep(0.001)
        updater.refresh('done: $DONE - total: $TOTAL ($REMAINING)')
    updater.format('done in $TOTAL', 0)
    updater.increment(1840)
    updater.format('done in $TOTAL - done: $DONE', 0)
    time.sleep(2)
    updater.wipe()
    assert updater.stats()['last_increment']==1840
    updater.refresh('done: $DONE - total: $TOTAL ($REMAINING)')
    updater.close()
    assert updater.closed==True
    tools.wrap('ag;lgfzvkraegjiaNJKaruegjuaarguaregjregjzdf;g', 11, 4)
    tools.wrap('ag;lgfzv kraegjiaNJKaruegjua arg uar egjre gjzdf;g', 11, 4)

    

########################################################################

def test_all():

    """
    Launch all tests of this module
    """
    
    print "# Testing tools module"
    test_SpecificDataFormat()
    test_DataAnalysis()
    test_SequenceComparison()
    test_Other()
    test_SequenceManipulationTools()

