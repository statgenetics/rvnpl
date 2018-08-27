import os
from .. import data, simul, egglib_binding
import example_files


########################################################################

def test_Container():

    """
    Test Container
    """
    
    print "## Testing data.Container"

    try:
        string = """>name1
AAAAAAAAAAAAAAAAAAAAAAAAA
>name2@0
AAAAAAAAAAAAAAAA
>name3@1
AAAAAAAAAAAAAAAAAAAAAAA
>name4@1
AAAAAAAAAAAAAA"""

        f = open('test_file1.txt', 'w')
        f.write(string)
        f.close()
        
        c = data.Container('test_file1.txt')
        assert len(c)==4
        assert len(c.groups())==1
        c = data.Container('test_file1.txt', groups=True)
        assert len(c)==4
        assert len(c.groups())==2
        c = data.Container(string=string, groups=True)
        assert c.group(0)==0
        assert c.name(2)=='name3'
        assert c[2] == ('name3', 'AAAAAAAAAAAAAAAAAAAAAAA', 1)
        for i in c:
            i.group = 4
            i.name = 'prout'
        assert len(c.groups())==1
        assert c.groups()[4] == ['prout']*4
        assert len(c)==c.ns()
        for i in range(4): c.ls(i)
        c.addSequences(c)
        assert len(c)==8
        c.addSequences([('name5', 'G'), ('name1', '', 4)])
        assert len(c)==10
        assert c.ls(9)==0
        c.append('name6', 'GGGGGGGGGGGGGGGc')
        c.append('name7', 'GGGGGGGGGGGGGGGC', 5)
        c.appendSequence(9, 'GGC')
        assert c.ls(9)==3
        assert c.contains_duplicates()
        assert c.duplicates() == ['prout']*7
        c.composition()
        names = c.names()
        code = c.encode()
        c.rename(code)
        assert c.names()==names
        assert c.find('prout')==0
        assert c.find('pro')==None
        assert c.find('pro', False)==0
        c.group(3,9999)
        assert c.group(3)==9999
        c.name(3,'hihi')
        assert c.name(3)=='hihi'
        assert c.groupByName('hihi')==9999
        assert c.isEqual()==False
        assert c.matches('huhu')==[]
        assert c.matches('hihi')==[3]
        assert c.matches('ih')==[3]
        l = len(c)
        c.remove('hihi')
        assert len(c)==l-1
        del c[0]
        assert len(c)==l-2
        c.no_duplicates()
        assert c.contains_duplicates()==False
        c.rename({'huhu':'flap'}, True)
        c.sequence(4, 'GAG')
        c.name(4, 'hook')
        assert c.sequence(4)=='GAG'
        assert c.sequenceByName('hook')=='GAG'
        c.set(4,1,'U')
        c.group(4, 999)
        c.shuffle(True)
        assert c.group(4)==999
        c.shuffle()
        c.encode()
        c2 = c.slice(1,4)
        assert len(c2)==3
        assert c2.name(0)==c.name(1)
        c.str()
        c.str(True)
        c.str(True, 22)
        c.str(lineLength=1)
        c.write('test_file2.txt')
        c.write('test_file2.txt', True)
        c.write('test_file2.txt', True, 22)
        c.write('test_file2.txt', lineLength=1)
        
        c.equalize('N')
        a = data.Align.create(c)
        assert a.get(1,a.ls()-1)=='N'
        c = data.Container.create(a)
        assert c.isEqual()
        assert [c.ls(i) for i in range(len(c))] == [a.ls()] * len(a)
        c.clear()
        assert len(c)==0
        
    finally:
        if os.path.isfile('test_file1.txt'): os.remove('test_file1.txt')
        if os.path.isfile('test_file2.txt'): os.remove('test_file2.txt')


########################################################################

def test_Align():

    """
    Test Align
    """
    
    print "## Testing data.Align"

    try:
        string = """>name1
AAAAAAAAAAAAAAAAAAAAAAAAA
>name2@0
AAAAAAAAAAAAAAAAAAAAAAAAA
>name3@1
AAAAAAAAAAAAAAAAAAAAAAAAA
>name4@1
AAAAAAAAAAAAAAAAAAAAAAAAA"""

        f = open('test_file1.txt', 'w')
        f.write(string)
        f.close()
        
        a = data.Align('test_file1.txt')
        assert len(a)==4
        assert len(a.groups())==1
        a = data.Align('test_file1.txt', groups=True)
        assert len(a)==4
        assert len(a.groups())==2
        a = data.Align(string=string, groups=True)
        assert a.group(0)==0
        assert a.name(2)=='name3'
        assert a[2] == ('name3', 'AAAAAAAAAAAAAAAAAAAAAAAAA', 1)
        for i in a:
            i.group = 4
            i.name = 'prout'
        assert len(a.groups())==1
        assert a.groups()[4] == ['prout']*4
        assert len(a)==a.ns()
        assert a.ls()==len('AAAAAAAAAAAAAAAAAAAAAAAAA')
        a.addSequences(a)
        assert len(a)==8
        a.addSequences([('name5', 'AAAAAAAAAAAAAAAAAAAAAAAAA'), ('name1', 'AAAAAAAAAAAAAAAAAAAAAAAAA', 4)])
        assert len(a)==10
        a.append('name6', 'AAAAAAAAAAAAAAAAAAAAAAAAA')
        a.append('name7', 'AAAAAAAAAAAAAAAAAAAAAAAAA', 5)
        assert a.contains_duplicates()
        assert a.duplicates() == ['prout']*7    
        a.composition()
        names = a.names()
        code = a.encode()
        a.rename(code)
        assert a.names()==names
        assert a.find('prout')==0
        assert a.find('pro')==None
        assert a.find('pro', False)==0
        a.group(3,9999)
        assert a.group(3)==9999
        a.name(3,'hihi')
        assert a.name(3)=='hihi'
        assert a.groupByName('hihi')==9999
        assert a.matches('huhu')==[]
        assert a.matches('hihi')==[3]
        assert a.matches('ih')==[3]
        l = len(a)
        a.remove('hihi')
        assert len(a)==l-1
        del a[0]
        assert len(a)==l-2
        a.no_duplicates()
        assert a.contains_duplicates()==False
        a.rename({'huhu':'flap'}, True)
        a.sequence(4, 'AAAAAAAAAAAAAAAAAAAAAAGAG')
        a.name(4, 'hook')
        assert a.sequence(4)=='AAAAAAAAAAAAAAAAAAAAAAGAG'
        assert a.sequenceByName('hook')=='AAAAAAAAAAAAAAAAAAAAAAGAG'
        a.set(4,1,'U')
        a.group(4, 999)
        a.shuffle(True)
        assert a.group(4)==999
        a.shuffle()
        a.encode()
        a2 = a.slice(1,4)
        assert len(a2)==3
        assert a2.name(0)==a.name(1)
        a.str()
        a.str(True)
        a.str(True, 22)
        a.str(lineLength=1)
        a.write('test_file2.txt')
        a.write('test_file2.txt', True)
        a.write('test_file2.txt', True, 22)
        a.write('test_file2.txt', lineLength=1)
        
        a.clear()
        assert len(a)==0
        assert a.ls()==0
        assert a.ns()==0
        
        a = data.Align.create([('', '10000'),('', '11100'),('', '00110')])
        a.binSwitch(0)
        
        a = data.Align.create([('', 'AAAAAAAAAAAAAAAAAAAAC', 1),
                               ('', 'AAAAAAAAAAAAGAAAAAAAA', 1),
                               ('', 'AAAAAAAAAAAAAAAAAACAA', 1),
                               ('', 'AACAAAAAAAAAGAAAAACAA', 2),
                               ('', 'AACAAAAATAAAGAAAAAAAA', 2),
                               ('', 'AACAAAAATAAAGAAAAAAAA', 2),
                               ('', 'AACAAAAAAAAAGAAAAACAA', 999)])
                               
        a.Rmin()
        a.Rmin(0.1, 2)
        assert a.column(1) == list('AAAAAAA')
        a.character(6,5)
        a.consensus()
        dm = a.dataMatrix('ACGT', 999)
        dm = a.dataMatrix('ACGT', 888)
        b = a.extract(5,8)
        assert b.ls()==3
        b = a.extract([1,2,3,3,2,1])
        assert b.ls()==6
        assert b.column(0)==b.column(5)
        assert b.column(1)==b.column(1)
        assert b.column(2)==b.column(2)
        x = len(a)
        
        b = data.Align.create(a)
        b.filter(1)
        assert len(b)==x
        b.filter(1, 'ACG')
        assert len(b)==x-2
        
        c = data.Align.create([('', '--AA--'),('', '--AAAA'),('', 'AAAA--')])
        c.fix_gap_ends()
         
        a.nexus()
        a.nexus(True)
        a.phylip('I')
        a.phylip('S')
        a.phyml()
        # polymorphism functions not tested thoroughly
        a.matrixLD(0.8, 1, 'ACG', 'TMRWSYKBDHVN?-')
        a.polymorphism() 
        a.polymorphism(allowMultipleMutations=True,
                minimumExploitableData=0.7, ignoreFrequency=3,
                validCharacters='ACG', missingData='MRWSYKBDHVNT?-',
                useZeroAsAncestral=False, skipDifferentiationStats=True,
                skipOutgroupBasedStats=True, skipAllHaplotypeStats=True,
                skipHaplotypeDifferentiationStats=True)
        a.polymorphismBPP(4)
        a.polymorphismBPP(5)

        
        x = a.ls()
        a.removePosition(3)
        assert a.ls()==x-1
        a.simErrors(0.2)
        for i in a.slider(20,5):
            i.polymorphism(minimumExploitableData=0.5)['D']

    finally:
        if os.path.isfile('test_file1.txt'): os.remove('test_file1.txt')
        if os.path.isfile('test_file2.txt'): os.remove('test_file2.txt')


########################################################################

def test_SSR():

    """
    Test SSR
    """
    
    print "## Testing data.SSR"
    
    ps = simul.CoalesceParamSet(None, [10,10], rho=4, nsites=4)
    m = simul.CoalesceStepwiseMutator(4*10)
    m.setSites(4)
    [ssr] = simul.coalesce(ps, m, 1)
    
    a0,b0,c0 = ssr.Fstats()
    a1,b1,c1 = ssr.Fstats(0)
    a2,b2,c2 = ssr.Fstats(1)
    a3,b3,c3 = ssr.Fstats(2)
    assert ssr.numberOfGenotypes()==20
    assert ssr.numberOfLoci()==4

    dm = egglib_binding.DataMatrix(40,4)
    for i in range(40):
        for j in range(4):
            dm.set(i, j, ssr.dataMatrix.get(i,j))
    
    ssr2 = data.SSR()
    ssr2.load(dm, [10,10])
    assert (a0,b0,c0)==ssr2.Fstats()
    assert (a1,b1,c1)==ssr2.Fstats(0)
    assert (a2,b2,c2)==ssr2.Fstats(1)
    assert (a3,b3,c3)==ssr2.Fstats(2)
    
    dm = egglib_binding.DataMatrix(40,4)
    for i in range(40):
        for j in range(4):
            dm.set(i, j, ssr.dataMatrix.get(i,j))
    
    ssr2.load(dm, [(10,0),(10,0)])
    assert (a0,b0,c0)==ssr2.Fstats()
    assert (a1,b1,c1)==ssr2.Fstats(0)
    assert (a2,b2,c2)==ssr2.Fstats(1)
    assert (a3,b3,c3)==ssr2.Fstats(2)

    dm = egglib_binding.DataMatrix()
    ssr3 = data.SSR()
    ssr3.load(dm)
    assert ssr.stats()==ssr.stats()

    ssr.str()

    string="""locus0  locus1  locus2  locus3
population1 indiv01 011/011 042/042 021/022 008/012
population1 indiv02 013/010 041/038 022/024 010/010
population1 indiv03 009/012 042/042 022/022 008/008
population1 indiv04 012/011 043/000 027/000 009/013
population1 indiv05 010/010 038/042 024/023 010/011
population1 indiv06 010/007 038/042 024/022 010/011
population1 indiv07 013/009 041/000 022/006 010/013
population1 indiv08 010/011 000/042 001/021 014/008
population1 indiv09 013/013 043/041 027/022 010/010
population1 indiv10 008/010 000/000 007/000 012/013
population2 indiv11 001/001 007/014 003/007 007/003
population2 indiv12 007/001 006/015 004/008 005/000
population2 indiv13 001/001 007/007 003/003 007/007
population2 indiv14 002/003 010/008 000/005 004/005
population2 indiv15 001/001 014/012 007/004 003/005
population2 indiv16 001/001 007/014 003/007 007/003
population2 indiv17 002/001 003/014 001/007 008/003
population2 indiv18 007/007 004/004 004/004 005/005
population2 indiv19 000/001 006/007 004/002 006/006
population2 indiv20 001/004 009/008 001/005 005/005
"""

    ssr4 = data.SSR()
    ssr4.parse(string, True, None, '/', True, '999')
    assert ssr4.numberOfGenotypes()==20
    assert ssr4.numberOfLoci()==4
    assert ssr4.dataMatrix.numberOfSequences()==40

    string="""population1@indiv01@011@042@021@008
population1@indiv02@013@041@022@010
population1@indiv03@009@042@022@008
population1@indiv04@012@043@027@009
population1@indiv05@010@038@024@010
population1@indiv06@010@038@024@010
population1@indiv07@013@041@022@010
population1@indiv08@010@000@001@014
population1@indiv09@013@043@027@010
population1@indiv10@008@000@007@012
population2@indiv11@001@007@003@007
population2@indiv12@007@006@004@005
population2@indiv13@001@007@003@007
population2@indiv14@002@010@000@004
population2@indiv15@001@014@007@003
population2@indiv16@001@007@003@007
population2@indiv17@002@003@001@008
population2@indiv18@007@004@004@005
population2@indiv19@000@006@004@006
population2@indiv20@001@009@001@005
"""

    ssr4.parse(string, False, "@", '/', False, '000')
    assert ssr4.numberOfGenotypes()==20
    assert ssr4.numberOfLoci()==4
    assert ssr4.dataMatrix.numberOfSequences()==20

    string="""                     locus0 locus1 locus2 locus3
population1 indiv01 011?011 042?042 021?022 008?012
population1 indiv02 013?010 041?038 022?024 010?010
population1 indiv03 009?012 042?042 022?022 008?008
population1 indiv04 012?011 043?000 027?000 009?013
"""

    ssr4.parse(string, True, " ", '?', True, '000')
    assert ssr4.numberOfGenotypes()==4
    assert ssr4.numberOfLoci()==4
    assert ssr4.dataMatrix.numberOfSequences()==8

    ssr.clear()
    assert ssr.numberOfGenotypes()==0
    assert ssr.numberOfLoci()==0


########################################################################

def test_TIGR():
    
    """
    Test TIGR
    """
    
    print "## Testing data.TIGR"
    
    tigr = data.TIGR(string=example_files.tigr_xml)
    gb = tigr.extract(50000,58000)
    assert len(gb) == 58000-50000+1
    
    try:
        f = open('test_file.xml', 'w')
        f.write(example_files.tigr_xml)
        f.close()
        tigr = data.TIGR('test_file.xml')
        gb = tigr.extract(42,84700)
        assert len(gb) == 84700-42+1
    finally:
        if os.path.isfile('test_file.xml'): os.remove('test_file.xml')
        

########################################################################

def test_GenBank():
    
    """
    Test GenBank
    """
    
    print '## Test data.GenBank'
    
    string = """LOCUS       AJ418369                3568 bp    mRNA    linear   PLN 15-APR-2005
DEFINITION  Medicago truncatula mRNA for nodulation receptor kinase (nork
            gene), wild type allele.
ACCESSION   AJ418369
VERSION     AJ418369.1  GI:21698782
KEYWORDS    nodulation receptor kinase; nork gene.
SOURCE      Medicago truncatula (barrel medic)
  ORGANISM  Medicago truncatula
            Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
            Spermatophyta; Magnoliophyta; eudicotyledons; core eudicotyledons;
            rosids; fabids; Fabales; Fabaceae; Papilionoideae; Trifolieae;
            Medicago.
REFERENCE   1
  AUTHORS   Endre,G., Kereszt,A., Kevei,Z., Mihacea,S., Kalo,P. and Kiss,G.B.
  TITLE     A receptor kinase gene regulating symbiotic nodule development
  JOURNAL   Nature 417 (6892), 962-966 (2002)
   PUBMED   12087406
REFERENCE   2  (bases 1 to 3568)
  AUTHORS   Kereszt,A.
  TITLE     Direct Submission
  JOURNAL   Submitted (26-OCT-2001) Institute of Genetics, Biological Research
            Center, Temesvari korut 62., Szeged 6726, HUNGARY
FEATURES             Location/Qualifiers
     source          1..3568
                     /organism="Medicago truncatula"
                     /mol_type="mRNA"
                     /db_xref="taxon:3880"
                     /chromosome="5"
                     /note="genotype A17"
     gene            1..3568
                     /gene="nork"
     CDS             376..3150
                     /gene="nork"
                     /note="Allele: wild type; NN1"
                     /codon_start=1
                     /product="nodulation receptor kinase"
                     /protein_id="CAD10808.1"
                     /db_xref="GI:21698783"
                     /db_xref="GOA:Q8L4H4"
                     /db_xref="InterPro:IPR000719"
                     /db_xref="InterPro:IPR001245"
                     /db_xref="InterPro:IPR001611"
                     /db_xref="InterPro:IPR008271"
                     /db_xref="InterPro:IPR011009"
                     /db_xref="InterPro:IPR013210"
                     /db_xref="InterPro:IPR017441"
                     /db_xref="InterPro:IPR021720"
                     /db_xref="UniProtKB/Swiss-Prot:Q8L4H4"
                     /translation="MELQVIRIFRLVVAFVLCLCIFIRSASSATKGFESIACCADSNY
                     TDPKTTLTYTTDHIWFSDKRSCRQIPEILFSHRSNKNVRKFEIYEGKRCYNLPTVKDQ
                     VYLIRGIFPFDSLNSSFYVSIGVTELGELRSSRLEDLEIEGVFRATKDYIDFCLLKED
                     VNPFISQIELRPLPEEYLHGFGTSVLKLISRNNLGDTNDDIRFPDDQNDRIWKRKETS
                     TPTSALPLSFNVSNVDLKDSVTPPLQVLQTALTHPERLEFVHDGLETDDYEYSVFLHF
                     LELNGTVRAGQRVFDIYLNNEIKKEKFDVLAGGSKNSYTALNISANGSLNITLVKASG
                     SEFGPLLNAYEILQARSWIEETNQKDLEVIQKMREELLLHNQENEALESWSGDPCMIF
                     PWKGITCDDSTGSSIITKLDLSSNNLKGAIPSIVTKMTNLQILNLSHNQFDMLFPSFP
                     PSSLLISLDLSYNDLSGWLPESIISLPHLKSLYFGCNPSMSDEDTTKLNSSLINTDYG
                     RCKAKKPKFGQVFVIGAITSGSLLITLAVGILFFCRYRHKSITLEGFGKTYPMATNII
                     FSLPSKDDFFIKSVSVKPFTLEYIEQATEQYKTLIGEGGFGSVYRGTLDDGQEVAVKV
                     RSSTSTQGTREFDNELNLLSAIQHENLVPLLGYCNEYDQQILVYPFMSNGSLLDRLYG
                     EASKRKILDWPTRLSIALGAARGLAYLHTFPGRSVIHRDVKSSNILLDQSMCAKVADF
                     GFSKYAPQEGDSYVSLEVRGTAGYLDPEYYKTQQLSEKSDVFSFGVVLLEIVSGREPL
                     NIKRPRIEWSLVEWAKPYIRASKVDEIVDPGIKGGYHAEALWRVVEVALQCLEPYSTY
                     RPCMVDIVRELEDALIIENNASEYMKSIDSLGGSNRYSIVMDKRALPSTTSTAESTIT
                     TQTLSHPQPR"
ORIGIN      
        1 taatgtttaa catcaacaat aaaactatag gaaaaaaaca taatcaacta tgcattgtac
       61 taattcaatc tctaactcgt cttccatctc tttcctagct acctcctgca gtttcctttc
      121 caggcctaaa gtcaaacacc atattttaac aatattcttt cttctacagc ttttatacaa
      181 gttcactata ttataggatt gatcagggtt cattttttct ttctttgaaa aatctctaag
      241 gggtgctgtt tccaaggcag agaatgaaat agaattcaga agaattttta tgttactata
      301 aaggaaagat gaaaagttag ttagcatgga ttcaagtttg ataaccctat ggggtaaaat
      361 ctctttcaga ttatgatgga gttacaagtt attaggatat ttagattggt tgtggcattt
      421 gttctttgtt tgtgtatatt tatcagatca gcttcttctg caactaaagg gtttgagagc
      481 atagcatgtt gtgctgattc aaattacaca gatccaaaaa ccaccctaac ttatacaaca
      541 gatcacatct ggttctctga taaaagaagt tgcagacaaa tacccgaaat tttgtttagc
      601 cacagaagca ataaaaatgt tcgaaaattt gaaatatatg aaggaaagag atgttataat
      661 ttgccaacag ttaaggatca agtatatttg ataaggggca tatttccctt tgatagttta
      721 aattcttcgt tttatgtttc gatcggggta acagaactag gcgaattaag atcgtctagg
      781 ctcgaggact tggaaattga gggagttttt agagccacca aagactacat agatttctgc
      841 ttattgaagg aggatgtcaa tcccttcatt tctcagattg aattgaggcc attacctgaa
      901 gaatacctac atggtttcgg tactagtgtt ttaaaactga taagcagaaa caatcttggt
      961 gacacaaatg atgatataag gttcccagat gaccaaaatg atagaatctg gaaacggaaa
     1021 gaaacttcaa ctccaacatc tgcccttcca ctgtctttca atgtcagcaa tgttgacctc
     1081 aaagacagtg tcacacctcc tctacaagtc ctacaaacag ctcttactca ccctgagcga
     1141 ttggagttcg tccatgatgg cctcgagacc gatgattatg aatactctgt gtttctccac
     1201 tttcttgaac taaatggcac tgtcagagca ggacaaaggg tgtttgacat ctatctaaac
     1261 aatgagatta aaaaggagaa atttgatgtt ttggctggag ggtccaagaa cagttacact
     1321 gccttgaaca tttcagcaaa tggatcactc aatataacct tagtcaaggc atctggatct
     1381 gagtttggac cccttttgaa tgcctatgaa atcctgcagg cacggtcgtg gattgaagag
     1441 accaaccaaa aagatttgga agttattcag aagatgagag aagaactgct gctgcacaac
     1501 caagaaaatg aagcattgga gagttggagt ggagaccctt gtatgatttt cccctggaaa
     1561 ggaataacat gtgatgattc aactggttca tctattatca ctaagctgga tctttcttcc
     1621 aataatctca agggagcaat tccttccatt gtcactaaga tgaccaattt acaaatactg
     1681 aacctgagcc acaaccagtt cgatatgtta ttcccctcgt ttccaccgtc ctccttgctg
     1741 atatcattgg atcttagcta caatgatctt tcaggatggc ttccagaatc cattatctca
     1801 ctgccacatt taaaatcatt atattttggc tgcaatccat ctatgagtga cgaagataca
     1861 acaaagttga acagttcact aatcaataca gattatggga gatgcaaagc aaaaaaacca
     1921 aagtttggac aagtattcgt gattggagct attacaagtg gatcactttt gattactttg
     1981 gctgttggaa ttctattttt ttgccgttat agacacaagt caattacttt ggaaggattt
     2041 ggaaagacct acccaatggc aacaaatata attttctctt tgccaagcaa agacgatttc
     2101 ttcataaagt ctgtatcagt taaaccgttc actttggagt atatagagca ggctacagaa
     2161 cagtacaaaa ctttaatagg tgaaggagga tttggctctg tttacagggg cactctagac
     2221 gatggtcaag aagtggcagt gaaagtgcgg tcatccacat caactcaggg aacccgagaa
     2281 tttgataatg agctaaacct actttcagct atacaacatg agaacctggt gcctcttctg
     2341 ggttactgta atgagtatga tcaacaaatt ctcgtgtatc ctttcatgtc caatggctct
     2401 ttgctagata gactatacgg ggaagcatca aagagaaaaa tattagactg gccaactaga
     2461 ctctctattg ctctcggtgc agctcgaggt ttggcatatc ttcacacatt tccaggacgt
     2521 tctgtaatac acagggacgt aaaatcgagc aatatactgc tggatcagag catgtgtgct
     2581 aaggttgcag attttggttt ctcaaaatat gctcctcagg aaggagacag ttatgtttcc
     2641 cttgaagtaa gaggaactgc agggtatctg gatcctgagt actacaaaac ccagcaatta
     2701 tctgaaaaaa gtgatgtttt cagctttggt gtggttcttc ttgaaattgt aagcggacgg
     2761 gaacctctca acataaagag accaaggatc gagtggagct tggttgaatg ggctaaacca
     2821 tacataagag catcaaaggt ggatgaaatt gtagatcctg gcatcaaggg aggatatcat
     2881 gcagaggcat tgtggagagt tgtggaagta gcactacaat gtctagaacc ctactcaaca
     2941 tatcggccat gcatggttga tattgtccgc gagttggagg atgctctcat tattgaaaac
     3001 aatgcatctg aatacatgaa atccatagac agccttggag gatccaaccg ctactcaatt
     3061 gttatggaca aacgggcgct gccttcaact acatctacag cagaatcaac tatcacaacc
     3121 caaaccttgt cacaccctca accgagatag taaatgggtc gatggaattc ttttgatttg
     3181 tttttgatca ttgctttagt aatatcacat tttaaatggt aaaggagaaa aatactactt
     3241 ctgattgtat ttccatccac tctatgtttc ttgaaactga atctctcttg ctcagcccca
     3301 gtttttatgg gtgaaaagaa taatttgggt caaatgcaag tgaaaccata tggtgcataa
     3361 tttaaaagcc atatcatatc atttgccaag tccaaagtaa aaatttcaca aactagttag
     3421 attgcgattt agtctataga cacttcaaca gagctatata cactatggtt gacttgcgac
     3481 taattcgctc aagcaggagg aacacatata tatgggaaac acttgtagaa ctattttgtt
     3541 tatagaaatg aaaatatttt ccgttttt
//"""

    gb0 = data.GenBank(string=string)
    assert gb0.accession == "AJ418369"
    assert gb0.definition == "Medicago truncatula mRNA for nodulation receptor kinase (nork gene), wild type allele."
    assert gb0.version == "AJ418369.1"
    assert gb0.locus == "AJ418369                3568 bp    mRNA    linear   PLN 15-APR-2005"
    assert gb0.GI == "21698782"
    assert gb0.keywords == "nodulation receptor kinase; nork gene."
    description, species, taxonomy = gb0.source
    assert description == "Medicago truncatula (barrel medic)"
    assert species == "Medicago truncatula"
    assert taxonomy == "Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliophyta; eudicotyledons; core eudicotyledons; rosids; fabids; Fabales; Fabaceae; Papilionoideae; Trifolieae; Medicago."
    assert len(gb0.references)==2
    assert len(gb0.others)==0
    
    try:
        f = open('test_file1.txt', 'w')
        f.write(string)
        f.close()
        gb = data.GenBank('test_file1.txt')
        assert gb.accession == "AJ418369"
        assert gb.definition == "Medicago truncatula mRNA for nodulation receptor kinase (nork gene), wild type allele."
        assert gb.version == "AJ418369.1"
        assert gb.locus == "AJ418369                3568 bp    mRNA    linear   PLN 15-APR-2005"
        assert gb.GI == "21698782"
        assert gb.keywords == "nodulation receptor kinase; nork gene."
        description, species, taxonomy = gb.source
        assert description == "Medicago truncatula (barrel medic)"
        assert species == "Medicago truncatula"
        assert taxonomy == "Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliophyta; eudicotyledons; core eudicotyledons; rosids; fabids; Fabales; Fabaceae; Papilionoideae; Trifolieae; Medicago."
        assert len(gb.references)==2
        assert len(gb.others)==0        
    finally:
        if os.path.isfile('test_file1.txt'): os.remove('test_file1.txt')
        
    gbf = data.GenBankFeature(gb)
    gbf.parse('''     vague_info            100..500
                     /note='for_testing'
''')
    gbf0 = gbf.copy(gb0)
    gbf.add_qualifier('comment', 'vague_info')
    assert 'comment' in gbf.qualifiers()
    assert 'note' in gbf.qualifiers()
    assert gbf.qualifiers()['comment']=='vague_info'
    assert 'comment' not in gbf0.qualifiers()
    assert len(gbf0.get_sequence())>0 and gbf0.get_sequence()==gbf.get_sequence()
    gbf.rc(len(gb))
    gbf.rc()

    loc = data.GenBankFeatureLocation("100..200")
    loc.addBaseChoice(250, 300, True, True)
    loc.addBaseRange(300,350, False, True)
    loc.addBetweenBase(410)
    loc.addSingleBase(450)
    assert loc.isRange()
    loc2 = loc.copy()
    loc.asOrder()
    assert not loc.isRange()
    assert loc2.isRange()
    loc2.rc(len(gb))
    loc2.setComplement()
    loc2.setNotComplement()
    loc2.shift(4)
    for a,b in loc: pass
    assert loc2[0][1]<loc2[1][0]
    str(loc)
    len(loc)
    
    gbf.set('none', loc2, note='this is for testing only')
    gbf.shift(-4)
    assert gbf.start()<gbf.stop()
    
    gb.add_feature(gbf)
    gb2 = gb.extract(454, 656)
    assert len(gb2)==656-454+1
    gb.number_of_features()
    gb.rc()
    gb.set_sequence('?'*len(gb))
    for f in gb:
        assert set(f.get_sequence())==set('?')
    try:
        gb.write('test_file2.txt')
        f = open('test_file3.txt', 'w')
        gb.write_stream(f)
        f.close()
    finally:
        if os.path.isfile('test_file2.txt'): os.remove('test_file2.txt')
        if os.path.isfile('test_file3.txt'): os.remove('test_file3.txt')
    str(gb)
    
    

########################################################################

def test_Tree():
    
    """
    Test Tree
    """
    
    print '## Test data.Tree'

    newick = "((a:0.1,b:0.2)100:0.4,(c:0.4,d:0.2)label,e:0.8);"
    
    try:
        f = open('test_file1.txt', 'w')
        f.write(newick)
        f.close()
        tree = data.Tree('test_file1.txt')
    finally:
        if os.path.isfile('test_file1.txt'): os.remove('test_file1.txt')
    
    tree = data.Tree(string=newick)
    tree.add_node(tree.root_node(), "f", 0.62)
    assert len(tree.all_leaves())==6
    c=0
    for i in tree: c+=1
    assert c==9
    tree2 = tree.copy()
    tree2.clean_edge_lengths()
    tree2.clean_internal_labels()
    for i in tree2:
        if i.numberOfDescendants()>0:
            assert i.get_label()==None
    assert tree.root_node().get_label()==None
    node = tree.get_node('label')
    tree.collapse(node)
    assert tree.root_node().numberOfDescendants()==5
    tree = data.Tree(string="((A,B),(C,(D,E)),((F,G),(H,I))));")
    assert tree.findGroup('ABCDE')!=None
    assert tree.findGroup('DEFG')==None
    assert tree.findMonophyleticGroup('ABCDE')==None
    assert tree.findMonophyleticGroup('HI')!=None
    
    r = egglib_binding.Random()
    ps = egglib_binding.ParamSet()
    ps.addPopulation(0)
    ps.migrationRate(0.1)
    ps.singles(0, 4)
    ps.singles(1, 4)
    controller = egglib_binding.Controller(ps, r)
    tree = data.Tree(string='(((1,2),(3,4)),(5,(6,(7,8))));')
    trees = []
    for i in range(100):
        x = ps.numberOfSamples()
        while x>1:
            x = controller.step()
        arg = controller.getArg()
        trees.append(data.Tree(string=arg.newick(0)))
        controller.reset()
    tree.frequency_nodes(trees)
    b = int(tree.root_node().descendants()[0].get_label())
    assert b>50 and b<=100
    tree.frequency_nodes(trees, True)
    b = float(tree.root_node().descendants()[0].get_label())
    assert b>0.5 and b<=1.
    tree = data.Tree(string='(((clade1_1,clade1_2)clade1_n1,(clade1_3,clade1_4)clade1_n2),(clade2_1,(clade2_2,(clade2_3,clade4_4)clade2_n1)clade2_n2)clade2_n3);')
    assert tree.get_node('clade1_1')
    assert tree.get_node_re('.+\d_n\d').get_label()=='clade1_n1'
    assert len(tree.get_nodes('clade1_n2'))==1
    assert len(tree.get_nodes_re('clade\d_n\d'))==5
    assert len(tree.get_terminal_nodes())==8
    assert tree.get_terminal_nodes()[-1].get_label()=='clade4_4'
    tree.last_node().get_label()
    tree.lateralize()

    while True:
        ps = egglib_binding.ParamSet()
        ps.singles(0, 20)
        controller = egglib_binding.Controller(ps, r)
        x = ps.numberOfSamples()
        while x>1:
            x = controller.step()
        arg = controller.getArg()
        tree = data.Tree(string=arg.newick(0))
        n = tree.root_node().descendants()[0]
        if len(n.leaves_down()) > 1: break
    L = tree.total_length()
    assert tree.total_length()==L
    tree.collapse(n)
    tree.midroot()
    assert tree.number_of_leaves()==len(tree.get_terminal_nodes())==20
    tree.number_of_nodes()
    tree.newick()
    tree.newick(True)
    tree.newick(False,True)
    tree.newick(True,True)
    tree = data.Tree(string='((c:0.2,(a:0.1,b:0.1)ab:0.1):0.1,(d:0.15,e:0.15):0.15,(f:0.2,(g:0.1,(h:0.05,i:0.05)h_i:0.05):0.1):0.1);')
    tree.remove_node(tree.get_node('f'))
    n = tree.get_node('ab')
    tree.remove_node(n)
    tree2 = tree.copy()
    n = tree2.get_node('h_i')
    tree2.root(n)
    n = tree.get_node('h_i')
    tree.root(n, 0.1)
    tree = data.Tree(string='((c:0.2,(a:0.1,b:0.1)ab:0.1):0.1,(d:0.15,e:0.15):0.15,(f:0.2,(g:0.1,(h:0.05,i:0.05)h_i:0.05):0.1):0.1);')
    assert tree.smallest_group(['a', 'b', 'c']) != None
    assert tree.smallest_group(['a', 'b', 'c', 'd', 'e']) != None
    assert tree.smallest_monophyleticGroup(['a', 'b', 'c']) != None
    assert tree.smallest_monophyleticGroup(['a', 'b', 'c', 'd', 'e']) == None
    try:
        tree.write('test_file2.txt')
        tree.write('test_file3.txt', False)
        tree.write('test_file4.txt', True, False)
        tree.write('test_file5.txt', False, False)
    finally:
        if os.path.isfile('test_file2.txt'): os.remove('test_file2.txt')
        if os.path.isfile('test_file3.txt'): os.remove('test_file3.txt')
        if os.path.isfile('test_file4.txt'): os.remove('test_file4.txt')
        if os.path.isfile('test_file5.txt'): os.remove('test_file5.txt')
    node = tree.smallest_group(['a', 'b'])
    node.add_son(label="raccoon", brlen=8.0)
    assert len(node.ascendants())==1
    assert node.numberOfAscendants()==1
    assert len(node.descendants())==3
    assert node.numberOfDescendants()==3
    assert len(tree.root_node().ascendants())==0
    assert tree.root_node().numberOfAscendants()==0
    assert len(tree.root_node().descendants())==3
    assert tree.root_node().numberOfDescendants()==3
    assert node.branch_to(node.descendants()[2])==8
    node.branch_from(node.ascendants()[0])
    node2 = data.TreeNode()
    node2.set_label("new")
    node2.add_son("X", 0.1)
    node2.add_son("Y")
    node.connect(node2, 0.4)
    node.set_branch_to(node2, 0.5)
    node2.set_branch_from(node, 0.2)
    node.get_label()
    assert node.is_descendant(node2)==True
    assert node.is_ascendant(node2)==False
    assert node2.is_descendant(node)==False
    assert node2.is_ascendant(node)==True
    assert node.numberOfRelatives()==node.numberOfAscendants()+node.numberOfDescendants()
    node.reverse(node2, True)
    node.leaves_down()
    node.leaves_up()
    node.set_label(None)
    node2.sort()
    node.remove_ascendant(node2)
    node.remove_descendant(node.descendants()[0])
    node.unlink()


########################################################################

def test_all():
    
    """
    Launch all tests of this module
    """
    
    print "# Testing data module"

    test_Container()
    test_Align()
    test_SSR()
    test_GenBank()
    test_Tree()
    test_TIGR()



