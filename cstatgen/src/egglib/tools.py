"""
This module contains tools used in some other parts of EggLib but that
might be of use for the package's users.
"""


__license__ = """
    Copyright 2008-2012 Stephane De Mita, Mathieu Siol

    This file is part of EggLib.

    EggLib is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    EggLib is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with EggLib.  If not, see <http://www.gnu.org/licenses/>.
"""

import re, os, string, StringIO, math, time, sys, random
import data, egglib_binding

RCconv = string.maketrans('ACGTMRYKBDHV','TGCAKYMRVHDB')


########################################################################

def rc(seq):
    
    """
    Reverse-complements a DNA sequence. Upper and lower-cases characters
    can be passed, the output is always upper-case. IUPAC characters
    (``ACGTMRWSYKBDHV``) are reverted. The characters ``N-?`` are
    returned as is. Other characters raise a ValueError.
    
    .. versionchanged:: 2.0.1
        Characters N, - and ? are correctly processed.

    .. versionchanged:: 2.0.2
        Reimplemented (will be faster for large sequences).
    """
    
    invalid = set(seq.upper())-set('ACTGMRWSYKBDHVN-?')
    if len(invalid):
        raise ValueError, 'cannot reverse-complement the sequence - invalid character(s): "%s"' %''.join(sorted(invalid))

    return seq[::-1].upper().translate(RCconv) 


########################################################################

def aln2fas(fname):

    """
    Imports a clustal-formatted alignment from the file name *fname* and
    returns a :class:`~egglib.Align` instance.
    """

    try:
        file = open(fname)
    except IOError:
        raise ValueError, 'cannot open this file: %s' %(fname)

    # the first line must start by CLUSTAL W or CLUSTALW
    
    line = file.readline()
    if line[:8]=='CLUSTALW': pass
    elif line[:9]=='CLUSTAL W': pass
    elif line[:8]=='CLUSTAL ': pass
    else:
        raise ValueError, 'invalid clustal format: %s' %fname

    # starts reading blocks
    container = data.Container()
    line = file.readline()

    while True:
        
        # skip empty lines
        while line!='' and line.strip()=='':
            line = file.readline()

        # detects end of file
        if line=='': break

        block_length = None

        # reads sequences
        while True:

            # conservation line
            if line[0]==' ':
                if block_length==None:
                    raise ValueError, 'invalid clustal format: %s' %fname
                line = file.readline()
                if line=='': break # end of file
                if line.strip()!='':
                    raise ValueError, 'invalid clustal format: %s' %fname
                break

            # reads sequence
            bits = line.split()

            # sequence line
            if len(bits)>3:
                raise ValueError, 'invalid clustal format: %s' %fname
                
            name = bits[0]
            seq = bits[1]
            if block_length==None:
                block_length=len(seq)
            elif block_length!=len(seq):
                print block_length, len(seq)
                print seq
                raise ValueError, 'invalid clustal format: %s' %fname
               
            # adds the sequence to the container (new sequence)
            pos = container.find(name)
            if pos==None:
                container.append(name, seq)
                pos = len(container)-1
                
            # adds the sequence (continuing old sequence)
            else:
                container.appendSequence(pos, seq)

            if len(line)==3:
                try:
                    i = int(line[2])
                except ValueError:
                    raise ValueError, 'invalid clustal format: %s' %fname
                if container.ls(pos) != i:
                    raise ValueError, 'invalid clustal format: %s' %fname

            # checks next line
            line = file.readline()
            if line=='':
                raise ValueError, 'invalid clustal format: %s' %fname
                
            # empty conservation line is caught by this line
            if line.strip()=='': break

    file.close()

    if not container.isEqual(): raise ValueError, 'invalid clustal format: %s' %fname
    return data.Align.create(container)


#######################################################################

def translate(input, code=1, strip=False):
        
    """
    Translates all sequences from nucleotide to proteins. Accepts
    sequence container instances and the return type matches the passed
    type. If *strip* is ``True``, all stop codon(s) present at the end
    of any sequence will be automatically stripped off. Setting this
    option to ``True`` will raise a :class:`ValueError` in case a
    :class:`~egglib.Align` is passed and sequences don't have all the
    same number of trailing stop codons. See the documentation of
    :meth:`GeneticCodes.translate` for documentation of the argument
    *code*. Ambiguous codons are translated if the implied possibilities
    translate all to the same codon. The IUPAC nomenclature is used.
    Note that ``N`` means ``A``, ``C``, ``G`` or ``T`` but that codons
    containing ``?`` or ``-`` will  always be translated as ``X``
    (except for ``---`` codons that are be translated as ``-``).
    """
    
    if isinstance(input, basestring):
        codons = [input[i-2:i+1] for i in range(2,len(input),3)]

        aa = [GeneticCodes.translate(i,code) for i in codons]
        prot = ''.join(aa)

        if strip:
            return prot.rstrip('*')
        else:
            return prot
    
    receptacle = input.__class__()
    
    for n,s,g in input:
        receptacle.append(n, translate(s, code, strip), g)

    return receptacle
    
#######################################################################

def longest_orf(sequence, clean=False, full=False, all=False, code=1, mini=1):
    
    """
    Finds the longest open reading frame in the sequence. By default,
    the longest sequence without stop codon (except for the trailing
    stop codon) is returned, therefore the returned ORF doesn't
    necessarily start by ATG and stops by a stop  codon. If *clean* is
    ``True``, returns the longest sequence encoding a valid protein
    sequence (without stop, without missing data, without gap). If
    *full* is ``True``, returns only genuine ORFs (starting by ATG and
    ending by a stop codon). If *all* is ``True``, returns a list of all
    ORFs (at least 3 of length), sorted by decreasing length. *code*
    specifies the genetic code; refer to the documentation of
    :meth:`GeneticCodes.translate`. *mini* specifies the minimum number
    of codons (or amino acids) of the returned ORF or ORFs (stop codons
    are not taken into account).
    
    .. versionchanged:: 2.0.1
        Added options; return the trailing stop codon when appropriate.

    .. versionchanged:: 2.1.0
        Added option *mini*. The behaviour of previous versions is
        reproduced by setting *mini* to 0.
    """
    
    # generates the sequences along the six frames
    translations = [
        sequence,
        sequence[1:],
        sequence[2:],
        rc(sequence),
        rc(sequence[1:]),
        rc(sequence[2:])
    ]

    # generates all the segments between stops
    orf = []
    for i in translations:
        
        codons = [i[j-2:j+1] for j in range(2,len(i),3)]
        
        aa = [GeneticCodes.translate(j, code) for j in codons]

        j=0
        while j<len(codons):
            
            buffer = []
            ln = 0 # orf length excluding stop codons
            
            # if needed, search for an ATG
            if full:
                while j<len(codons):
                    if codons[j].upper()=='ATG':
                        break
                    j+=1
            
            while j<len(codons):
                
                if aa[j]=='*':
                    buffer.append(codons[j])
                    break
                    
                if clean and aa[j]=='X':
                    break
                    
                # adds until a stop (included) or a break (excluded, if clean)
                buffer.append(codons[j])
                ln+=1
                j+=1
                
            # adds the orf (needs to be full if requested)
            if ln>=mini and (not full or (j<len(codons) and aa[j]=='*')):
                orf.append(''.join(buffer))
            
            j+=1
        
    orf.sort(key=len, reverse=True)
    
    if all: return orf
    elif len(orf): return orf[0]


########################################################################

_raw_genetic_codes = """1. The Standard Code (transl_table=1) Standard

    AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  Starts = ---M---------------M---------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

2. The Vertebrate Mitochondrial Code (transl_table=2) VertebrateMitochondrial

    AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG
  Starts = --------------------------------MMMM---------------M------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

3. The Yeast Mitochondrial Code (transl_table=3) YeastMitochondrial

    AAs  = FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  Starts = ----------------------------------MM----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code (transl_table=4) MoldProtozoanCoelenterateMitochondrialMycoplasmaSpiroplasma

    AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  Starts = --MM---------------M------------MMMM---------------M------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

5. The Invertebrate Mitochondrial Code (transl_table=5) InvertebrateMitochondrial

    AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG
  Starts = ---M----------------------------MMMM---------------M------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

6. The Ciliate, Dasycladacean and Hexamita Nuclear Code (transl_table=6) CiliateDasycladaceanHexamitaNuclear

    AAs  = FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  Starts = -----------------------------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

9. The Echinoderm and Flatworm Mitochondrial Code (transl_table=9) EchinodermFlatwormMitochondrial

    AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG
  Starts = -----------------------------------M---------------M------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

10. The Euplotid Nuclear Code (transl_table=10) EuplotidNuclear

    AAs  = FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  Starts = -----------------------------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

11. The Bacterial, Archaeal and Plant Plastid Code (transl_table=11) Plastid

    AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  Starts = ---M---------------M------------MMMM---------------M------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

12. The Alternative Yeast Nuclear Code (transl_table=12) AltYeastNuclear

    AAs  = FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  Starts = -------------------M---------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

13. The Ascidian Mitochondrial Code (transl_table=13) AscidianMitochondrial

    AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG
  Starts = ---M------------------------------MM---------------M------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

14. The Alternative Flatworm Mitochondrial Code (transl_table=14) AltFlatwormMitochondrial

    AAs  = FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG
  Starts = -----------------------------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

15. Blepharisma Nuclear Code (transl_table=15) BlepharismaNuclear

    AAs  = FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  Starts = -----------------------------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

16. Chlorophycean Mitochondrial Code (transl_table=16) ChlorophyceanMitochondrial

    AAs  = FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  Starts = -----------------------------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

21. Trematode Mitochondrial Code (transl_table=21) TrematodeMitochondrial

    AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG
  Starts = -----------------------------------M---------------M------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

22. Scenedesmus obliquus mitochondrial Code (transl_table=22) ScenedesmusMitochondrial

    AAs  = FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  Starts = -----------------------------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

23. Thraustochytrium Mitochondrial Code (transl_table=23) ThraustochytriumMitochondrial

    AAs  = FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  Starts = --------------------------------M--M---------------M------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""


########################################################################

_raw_genetic_codes = _raw_genetic_codes.split('\n')
_genetic_codes = {}
for name, empty1, aa, starts, base1, base2, base3 in [
    _raw_genetic_codes[i:i+7] for i in range(0, len(_raw_genetic_codes), 8) ]:
    
    # names
    match = re.match('(\d+)\.[ The]? (.+) \(transl_table=\d+\) (.+)', name)
    index = int(match.group(1))
    name = match.group(2)
    short = match.group(3)

    _genetic_codes[index] = {
        'short_name': short,     'long_name': name,
        'initiation_codons': [],
        'table': {} }

    # code
    for A,F,S,T,M in zip(*map(lambda x: x.split()[-1],[aa,base1,base2,base3,starts])):
        codon = ''.join([F,S,T])
        _genetic_codes[index]['table'][codon] = A
        if M=='M':
            _genetic_codes[index]['initiation_codons'].append(codon)
del name,empty1,aa,starts,base1,base2,base3,A,F,S,T,M


########################################################################

class GeneticCodes:
    
    """
    Holds genetic code. Instantiating this class is pointless since its
    contains only class methods.
    """

    @classmethod
    def codes(cls):
        
        """
        Gives the list of code identifiers. Each code is represented by
        three identifiers: ``(index, short, long)`` where *index* is the
        integer identifier matching NCBI nomenclature (beware that
        indices are not consecutive); *short* is a egglib-defined word
        summarizing the code which can be used as an alternative access
        means; and *long* is the full name of the genetic code.
        """
        
        return [ (index, _genetic_codes[index]['short_name'],
                _genetic_codes[index]['long_name']) for index in _genetic_codes]

    ####################################################################

    @classmethod
    def index(cls, name):
        
        """
        Tries to identify the index of the genetic code from its short
        or full name. Returns ``None`` if the string matches no model.
        The comparison is case-independent.
        """
        
        for index in _genetic_codes:
            if (_genetic_codes[index]['short_name'].lower()==name.lower() or
                _genetic_codes[index]['long_name'].lower()==name.lower()):
                    
                    return index
                    
        return None

    ####################################################################

    @classmethod
    def translate(cls, codon, code=1):
        
        """
        Translate the codon *codon* using the indicated code. *code*
        is an identifier (index, short or long name) matching NCBI
        nomenclature. Returns the one-letter amino acid code
        corresponding to codon, ``*`` for stop codon and ``X`` for any
        invalid codon (string with a length different than 3 or
        containing missing data or gaps). The codon specification is
        case-independent. Ambiguous codons might still be translated if
        the implied possibilities translate all to the same codon. The
        IUPAC nomenclature is used. Note that ``N`` means ``A``, ``C``,
        ``G`` or ``T`` but that codons containing ``?`` or ``-`` will
        always be translated as ``X``. However, ``---`` will be
        translated as ``-``.
        """
        
        if isinstance(code, basestring):
            index = GeneticCodes.index(code)
        else:
            index = code
        
        if index not in _genetic_codes:
            raise KeyError, 'unknown genetic code specification: %s' %str(code)

        codon = codon.upper()
        
        if codon == '---':
            return '-'
        
        if codon in _genetic_codes[index]['table']:
            return _genetic_codes[index]['table'][codon]

        bases = []
        for c in codon:
            if c in set('ACGT'): bases.append([c])
            elif c=='M': bases.append(['A','C'])
            elif c=='R': bases.append(['A','G'])
            elif c=='W': bases.append(['A','T'])
            elif c=='S': bases.append(['C','G'])
            elif c=='Y': bases.append(['C','T'])
            elif c=='K': bases.append(['G','T'])
            elif c=='V': bases.append(['A','C','G'])
            elif c=='H': bases.append(['A','C','T'])
            elif c=='D': bases.append(['A','G','T'])
            elif c=='B': bases.append(['C','G','T'])
            elif c=='N': bases.append(['A','C','G','T'])
            else: return 'X'
        
        aa = None
        for i in bases[0]:
            for j in bases[1]:
                for k in bases[2]:
                    codonx = ''.join((i,j,k))
                    aax = _genetic_codes[index]['table'][codonx]
                    if aa!=None:
                        if aax!=aa:
                            return 'X'
                    else:
                        aa = aax
                        
        return aa

    ####################################################################
    
    @classmethod
    def is_start(cls, codon, code=1):

        """
        Returns ``True`` if the codon is encoding one the observed
        translational start for this genetic code, ``False`` otherwise
        (including if the codon is invalid). Arguments are similar as
        for :meth:`translate`.
        """
        
        if isinstance(code, basestring):
            index = GeneticCodes.index(code)
        else:
            index = code
        
        if index not in _genetic_codes:
            raise KeyError, 'unknown genetic code specification: %s' %str(code)

        codon = codon.upper()

        return codon in _genetic_codes[index]['initiation_codons']


########################################################################

def staden(fname=None, string=None, delete_consensus=True):
    
    """
    Imports a Staden output file as an :class:`~egglib.Align` instance. The
    file should have been  generated from a contig alignment by the GAP4
    contig editor, using the command "dump contig. to file". The
    sequence named ``CONSENSUS``, if present, is automatically removed
    unless the option *delete_consensus* is ``False``.
    
    The Staden outfile file can be read from a file (using the argument
    *fname* or directly from a string (using *string*). It is required
    to pass either a file name as *fname* or a Staden string as *string*
    and it is not allowed to passb both.
    
    Staden's default convention is followed:

    * ``-`` codes for an unknown base and is replaced by ``N``.
    * ``*`` codes for an alignment gap and is replaced by ``-``.
    * ``.`` represents the same sequence than the consensus at that
      position.
    * White space represents missing data and is replaced by ``?``.
    
    .. versionadded:: 2.0.1
        Add argument *delete_consensus*.
        
    .. versionchanged:: 2.1.0
        Read from string or fname.
    """

    if fname!=None and string!=None:
        raise ValueError, 'it is not allowed to pass both fname and string to `staden()`'

    if fname==None and string==None:
        raise ValueError, '`staden()`: no data to convert'

    if fname!=None:
        f = open(fname)
        string = f.read()
        f.close()

    align = data.Align()
    align._object = egglib_binding.Staden.parse(string, delete_consensus)
    return align

########################################################################

def concat(aligns, spacer=0, ch='?', strict=True, groupCheck=True):
    
    """ Concatenates sequence alignments. A unique :class:`~egglib.Align` is
    returned. All different sequences from all passed alignments are
    represented in the final alignment. Sequences whose name match are
    matching are concatenated. In case several sequences have the same
    name in a given segment, the first one is considered and others are
    discarded. In case a sequence is missing for a particular segment,
    a stretch of non-varying characters is inserted to replace the
    unknown sequence.
    
    *aligns* must be an iterable containing :class:`~egglib.Align` instances.
    
    *spacer* specifies the length  of  unsequenced stretches
    (represented by non-varying characters)  between concatenated
    alignments. If *spacer* is a positive integer, the length of all
    stretches will be identical. If *spacer* is an iterable containing
    integers, each specifying the interval between two consecutive
    alignments (if *aligns* contains ``n`` alignments, *spacer* must be
    of length ``n-1``).
    
    *ch* gives the character to used for conserved stretches and for
    missing segments.
    
    If *strict* is ``False``, the name comparison will not extend
    further than the length of the shorter name: for example, names
    ``anaconda`` and ``anaco`` will match,  and the concatenated
    sequence will be named ``anaconda`` (regardless of which name
    appears first in the list of :class:`~egglib.Align` instances).
    
    If *groupCheck* is ``True``, an exception will be raised in case of
    a mismatch between group labels of different sequence segments
    bearing the same name. Otherwise, the group of the first segment
    found will be used as group label of the final sequence.
    
    .. versionadded:: 2.0.1
        The arguments allowing to customize function's behaviour. """

    # carefully checks arguments
    if not ('__getitem__' in dir(aligns) and '__len__' in dir(aligns)):
        raise ValueError, 'first argument of concat function must be a sequence of Align instances'
    if isinstance(spacer, int):
        spacer = [spacer] * (len(aligns) -1)
    if min(spacer)<0:
        raise ValueError, 'spacer argument of concat function must be a non-negative integer or a sequence of non-negative integers'
    if len(spacer)!=(len(aligns)-1):
        raise ValueError, 'spacer argument must have n-1 arguments, where n is the number of alignments'
    elif ('__getitem__' in dir(spacer) and '__len__' in dir(spacer)):
        pass
    else:
        raise ValueError, 'spacer argument of concat function must be a non-negative integer or a sequence of non-negative integers'
    if not (isinstance(ch, str) and len(ch)==1):
        raise ValueError, 'ch argument of concat function must be a string of length 1'
    if not isinstance(strict, bool):
        raise ValueError, 'strict argument of concat function must be a boolean'
    if not isinstance(groupCheck, bool):
        raise ValueError, 'groupCheck argument of concat function must be a boolean'
    
    # sequences will be initially stored in a Container
    conc = data.Container()
    
    c=0
    length=0
    for align in aligns:

        # checks
        if not isinstance(align, data.Align):
            raise ValueError, 'first argument of concat function must be a sequence of Align instances'

        # checks if already-present sequences are in the new segment
        lmapping = [None] * len(conc)
        rnew = []
        for i in range(len(conc)):
            j = align.find(conc.name(i), strict)
            if j!=-1:
                lmapping[i] = j
        
        # checks the reverse (to find reciprocal matches when strict is false)
        for i in range(len(align)):
            j = conc.find(align.name(i), strict)
            if j!=None:
                if len(align.name(i))>len(conc.name(j)):
                    conc.name(j, align.name(i))
                if lmapping[j]==None:
                    lmapping[j]=i
            if j==None:
                rnew.append(i)
        
        # extends the old sequences
        for i in range(len(conc)):
            if lmapping[i]==None:
                conc.appendSequence(i, ''.join([ch]*align.ls()))
            else:
                if groupCheck:
                    if conc.group(i) != align.group(lmapping[i]):
                        raise ValueError, 'group labels of sequence "%s" don\'t match between segments' %conc.name(i)
                conc.appendSequence(i, align.sequence(lmapping[i]))

        # adds new sequences
        for i in rnew:
            seq = ''.join([ch]*length) + align.sequence(i)
            conc.append( align.name(i), seq, align.group(i) )
        
        # consider spacing
        if c<len(spacer):

            # checks again
            if (not isinstance(spacer[c], int) or spacer[c]<0):
                raise ValueError, 'spacer argument of concat function must be a non-negative integer or a sequence of non-negative integers'
            
            # appends spacer
            for i in range(len(conc)):
                conc.appendSequence(i, ''.join([ch] * spacer[c]))
            
            # counts it
            length+= spacer[c]

        c+=1
        length+=align.ls()

    # convert to align
    aConc = data.Align()
    aConc.addSequences(conc)
    return aConc

########################################################################

def genalys2fasta(iname):
    
    """
    Converts Genalys-formatted sequence alignment files to fasta. The
    function imports files generated through the option *Save SNPs* of
    Genalys 2.8. *iname* if the name of the Genalys output file. Returns
    an :class:`~egglib.Align` instance.
    """
    
    file = open(iname,"r")
    res = data.Align()

    insertions = []
    flag = False
    for line in file:
        line = line.split("\t")

        if (len(line)>1 and line[0]=="Polymorphism"):
            flag = True

        if (len(line)>1 and line[0]=="IN" and flag):
            insertions.extend(line[1].split("/"))

    if len(insertions)>0:
        tp = insertions[0].split("_")
        if len(tp)==1:
            tp = tp[0].split(".")
            if len(tp)==1:
                tp.append("1")
        finsertions = [tp]
    for i in insertions:
        i= i.split("_")
        if len(i)==1:
            tp = tp[0].split(".")
            if len(tp)==1:
                i.append("1")
        if i[0]!=finsertions[-1][0]:
            finsertions.append(i)
        finsertions[-1][1] = i[1]
    
    if len(insertions)>0:
        insertions = finsertions
    
    file.close()
    file = open(iname)
    
    noms = []
    sequences = []
    maxlen = 0
    
    for line in file:
        line = line.split("\t")
            
        if len(line)>1:
            bidon= re.match(".+\.ab1$",line[1])
            if bidon!=None:
                noms.append(line[1])
                sequences.append("")
                index = 6
                for i in range(10):
                    if line[i]=="F" or line[i]=="R":
                        index = i+1
                        break
                if line[index]!="":
                    debut = int(line[index])-1
                    for i in insertions:
                        if int(i[0])<=debut:
                            debut= debut+ int(i[1])
                        else:
                            break
                    for i in range(debut):
                        sequences[-1]= sequences[-1]+ "?"
                sequences[-1]= sequences[-1]+line[-1].rstrip("\n")
                if len(sequences[-1])>maxlen:
                    maxlen = len(sequences[-1])
    
    for i in range(len(sequences)):
        sequences[i] = sequences[i].replace("_","-")
        for j in range(len(sequences[i]),maxlen):
            sequences[i]= sequences[i]+"?"
        res.append(noms[i], sequences[i])

    return res

########################################################################

def get_fgenesh(fname):
    
    """
    Imports fgenesh output. *fname* must be the name of a file
    containing fgenesh ouput. The feature definition are parsed an
    returned as a list of ``gene`` and ``CDS`` features represented by
    dictionaries. Note that 5' partial features might not be in the
    appropriate frame and that it can be necessary to add a ``codon_start``
    qualifier.
    """
    
    # checking
    if not os.path.isfile(fname):
        raise Exception, 'error, file not found: %s' %fname
        return []

    # define the locus name
    locus = os.path.basename(fname).split('.')[0]

    # import the raw data
    f = open(fname)
    data = f.read()
    f.close()

    # supports for mac/windows files
    data= data.replace('\r\n', '\n')
    data= data.replace('\r', '\n')

    # gets the feature table
    try:
        data_sub= data.split('   G Str   Feature   Start        End    Score           ORF           Len\n')[1].split('Predicted protein(s):\n')[0]
    except IndexError:
        raise ValueError, 'invalid fgenesh format: %s' %fname
    data_sub= data_sub.split('\n\n')

    # edit
    del data_sub[-1]
    data_sub[0]= '\n'.join(data_sub[0].split('\n')[1:])

    # iteratively grabs the features
    features = {}
    for i in data_sub:
        pos = []
        start = 1
        rank = '---'
        strand = '---'
        for j in i.split('\n'):
            a =re.search(' ?[0-9]+ ([+|-])      (TSS|PolA) +([0-9]+)', j)
            b =re.search(' ?([0-9]+) ([+|-]) + ([0-9])+ CDS(o|f|i|l) +([0-9]+) - +([0-9]+) +[-\.0-9]+ + ([0-9]+)', j)
            if b:
                if (b.group(3)=="1"):
                    if (int(b.group(5))==int(b.group(7))): start= 1
                    elif (int(b.group(5))==(int(b.group(7))-1)): start= 2
                    elif (int(b.group(5))==(int(b.group(7))-2)): start= 3
                    else:  sys.exit('error in file '+fname)
                pos.append( [int(b.group(5))-1, int(b.group(6))-1 ] )
                rank = b.group(1)
                if (b.group(2)=='+'): strand = 'plus'
                else: strand = 'minus'

        features['cds'+rank] ={
            'gene': locus+'_'+rank,
            'strand': strand,
            'pos': pos,
            'type': 'CDS',
            'note': 'fgenesh prediction'
        }

        features['gene'+rank] ={
            'gene': locus+'_'+rank,
            'strand': strand,
            'pos': [[ pos[0][0], pos[-1][1] ]],
            'type': 'gene',
            'note': 'fgenesh prediction'
        }
        
    # gets the sequence section
    try:
        data_sub= data.split('   G Str   Feature   Start        End    Score           ORF           Len\n')[1].split('Predicted protein(s):\n')[1].split('>')
    except IndexError:
        raise ValueError, 'invalid fgenesh format: %s' %fname
    del data_sub[0]

    if ( (2*len(data_sub)!=len(features)) and
           (len(data_sub)!=len(features)) ) : raise IOError, 'cannot import %s' %fname

    # returns the sequences as a table
    return [ features[i] for i in features ]
    
########################################################################

class Mase(list):

    """
    Minimal implementation of the mase format (allowing input/output
    operations). This class emulates a list of dictionaries, each
    dictionary representing a sequence and describing the keys *header*,
    *name* and *sequence*. However, the string formatter (``str(mase)``
    or ``print mase``, where *mase* is a :class:`~egglib.tools.Mase` instance)
    generates a mase-formatted string. Object attributes are *header*
    (a string with file-level information), *species* (the species of
    the ingroup), *align* (an :class:`~egglib.Align` instance corresponding
    to the data contained in the instance, and created upon construction).
    Modifying this instance has no effect.
    """

    def __init__(self, input=None):
        
        """
        The constructor takes an optional argument that can be a string
        giving the path to a mase-formatted file, or a :class:`~egglib.Align`
        instance. The constructor is currently unable to import
        population labels, and only sequences marked as *ingroup* are
        imported.
        
        .. versionchanged:: 2.0.1
            An :class:`IOError` is raised upon file formatting error.
        """
        
        if isinstance(input, basestring):
            self._parse(input)

        elif isinstance(input, data.Align):
            self._import(input)

        elif input==None:
            self.header= ''
            self.species=None
            self.align= data.Align()

        else:
            raise ValueError, 'mase initialized with invalid type'

    ####################################################################

    def _parse(self, fname):

        self.header= ''
        self.align= data.Align()

        if not os.path.isfile(fname):
            raise IOError, 'cannot open '+fname
        f=open(fname)
        stuff= f.read()
        f.close()

        # foreign file support
        stuff= stuff.replace('\r\n', '\n')
        stuff= stuff.replace('\r', '\n')
        stuff= stuff.split('\n')
        stuff= [ i+'\n' for i in stuff ]

        while (stuff[0][:2]==';;'):
            self.header+= stuff[0][2:]
            del stuff[0]

        stuff = ''.join(stuff)
        stuff= stuff.split(';')
        if not len(stuff):
            raise IOError, 'error in mase file '+fname
        del stuff[0]
        
        buff = []

        for i in stuff:
            i= i.split('\n')
            if (len(i)<3):
                self.clear()
                self.header= ''
                raise IOError, 'error in mase file '+fname
            buff.append({'header': i[0],
                          'name': i[1],
                          'sequence': ''.join(i[2:]).upper()})                

        obj= re.search('@ of species = (\d+)  INGROUP_(.+)\n([0-9,]+)',self.header)
        indices = map(int, obj.group(3).split(','))
        if len(indices) != int(obj.group(1)):
            raise IOError, 'error in mase file '+fname
        if len(indices)>len(buff):
            raise IOError, 'error in mase file '+fname
        self.species = obj.group(2).strip()
        for i in indices:
            if i-1>=len(buff):
                raise IOError, 'error in mase file '+fname
            self.align.append( buff[i-1]['name'], buff[i-1]['sequence'] )
            self.append(buff[i-1])

    ####################################################################

    def _import(self, alignment):
        
        self.header= ''
        self.align= data.Align()
        for i in alignment:            
            self.append({'header': '',
                         'name': i[0],
                         'sequence': i[1]})
            self.align.append(*i)

    ####################################################################

    def __str__(self):

        string= StringIO.StringIO()
        string.write(';;%s\n' %self.header.strip())
        for i in self:
            string.write(';%s\n%s\n' %(i['header'].strip(), i['name']))
            c=0
            for j in i['sequence']:
                string.write(j)
                c+=1
                if (c==60):
                    string.write('\n')
                    c=0
            string.write('\n')
        return string.getvalue()


########################################################################

class ReadingFrame:
    
    """
    Handles reading frame positions.
    """
    
    def __init__(self, frame):
        
        """
        *frame* must be a sequence of ``(start, stop, codon_start)``
        sequences where *start* and *stop* gives the first and last
        position of an exon and  *codon_start* is 1 if the first position
        of the exon is the first position of a codon (e.g. ``ATG ATG``),
        2 if the first position of the segment is the second position of
        a codon (e.g. ``TG ATG``), 3 if the first position of the
        segment is the third position a of codon (e.g. ``G ATG``), or
        ``None`` if the reading frame is continuing the previous exon.
        If *codon_start* of the first segment is ``None``, 1 will be
        assumed. It is not possible to modify the codon positions held
        by the instance after construction.
        """
        
        self._frame = []
        try:
            for start,stop,codon_start in frame:
              self._frame.append((start,stop,codon_start)) # deep copy  
        except ValueError,TypeError:
            raise ValueError, 'invalid reading frame specification'
        self._cached_codons = self.codons()

    ####################################################################

    def exon(self, x):
        
        """
        Returns the exon index of a position. Returns ``-1`` if the
        position falls outside specified segments (out of ranges or
        in introns).
        """
        
        for i,v in enumerate(self._frame):
            if x>=v[0] and x<=v[1]: return i
        return None
        
    ####################################################################

    def codon(self, x):
        
        """
        If the position ``x`` falls in a complete codon, returns the
        three positions of that codon. If ``x`` fall outside of defined
        segments, or in a codon that appears not to be completely
        available, returns ``None``.
        
        .. note::
            the codon positions are cached at build time. As a result,
            the result of this method will be incorrect if frame
            positions are changed after the creation of the instance.
        """
        
        for triplet in self._cached_codons:
            if x in triplet: return triplet
        return None
        
    ####################################################################

    def codons(self):

        """
        Returns the list of complete codons (as triplets of absolute
        positions).
        """
        
        ret = []
        cache = []
        for a,b,c in self._frame:
            if c==2:
                i = a+2
            elif c==3:
                i = a+1
            else:
                i = a
            if c!=None:
                cache = []
            
            while i<=b:
                cache.append(i)
                i+=1
                if len(cache)==3:
                    ret.append(tuple(cache))
                    cache = []
            
        return ret

########################################################################

def _mycmpnuc(nuc1, nuc2):
    nuc1= nuc1.upper()
    nuc2= nuc2.upper()

    if (nuc1==nuc2):
        return True

    if (nuc1=='?' or nuc2=='?'):
        return True

    if (nuc1=='A'):
        if nuc2 in set('AMRWDHVN'):
            return True

    if (nuc1=='C'):
        if nuc2 in set('CMSYBHVN'):
            return True

    if (nuc1=='G'):
        if nuc2 in set('GRSKBDVN'):
            return True

    if (nuc1=='T'):
        if nuc2 in set('TWYKBDVN'):
            return True

    if (nuc1=='M'):
        if nuc2 in set('ACMRWSYBDHVN'):
            return True

    if (nuc1=='R'):
        if nuc2 in set('AGMRWSKBDHVN'):
            return True

    if (nuc1=='W'):
        if nuc2 in set('ATMRWYKBDHVN'):
            return True

    if (nuc1=='S'):
        if nuc2 in set('CGMRSYKBDHVN'):
            return True

    if (nuc1=='Y'):
        if nuc2 in set('CTMWSYKBDHVN'):
            return True

    if (nuc1=='K'):
        if nuc2 in set('GTRWSYKBDHVN'):
            return True

    if (nuc1=='B'):
        if nuc2 in set('CGTMRWSYKBDHVN'):
            return True

    if (nuc1=='D'):
        if nuc2 in set('AGTMRWSYKBDHVN'):
            return True

    if (nuc1=='H'):
        if nuc2 in set('ACTMRWSYKBDHVN'):
            return True

    if (nuc1=='V'):
        if nuc2 in set('ACGMRWSYKBDHVN'):
            return True

    if (nuc1=='N'):
        if nuc2 in set('ACGTRYSWKMBDHVN'):
            return True

    return False

########################################################################

def compare(seq1, seq2):
    
    """
    Compares two sequences. Sequences are different if they have
    different lengths or if they differ by at least one position. The
    comparison supports IUPAC ambiguity characters (for example, A and
    M are not considered to be different). Furthermore, partially
    overlapping ambiguity characters (for example, M and R) are
    not taken as different. Returns ``True`` is sequences are
    identical (or differ only by overlapping IUPAC characters), ``False``
    otherwise.
    """
    
    if (len(seq1)!=len(seq2)): return False
    for i in range(len(seq1)):
        if not _mycmpnuc(seq1[i],seq2[i]): return False
    return True

########################################################################

def locate(sequence, motif, start=0, stop=-1):
    
    """
    Locates the position of the *motif* in *sequence*. *motif* and
    *sequence* should be DNA sequences Ambiguity characters (``M``,
    ``R``, ``W``, ``S``, ``Y``, ``K``, ``B``, ``D``, ``H``, ``V`` and
    ``N``) are recognized and match the appropriate characters. ``?``
    matches any character. Note that the meaning of ``N`` (``A``, ``C``,
    ``G`` or ``T``) is very different to ``?`` (any character). 
    *start* and *stop* allow to restrict search to a given subset of
    *sequence* (the returned index is still given with respect to the
    full sequence). The function returns the position of the first
    exact match or, if there is not exact match, the position of the
    first matching position allowing ambiguity charactor, or, if there
    is no match at all, ``None``.
    
    if *hasAmb* is ``True``, ambiguities will be supported in the target
    sequence (*sequence*). With that mode on, ambiguities of the motif
    sequence (*motif*) will only be considered as a match if the target
    sequence account for all 
    
    .. versionchanged:: 2.1.0
       Supports ambiguity characters in *sequence*. Returns exact
       matches first.
    """

    sequence = sequence.upper()
    motif = motif.upper()

    if (stop<0): stop = len(sequence)+stop
    sequence= sequence[start : stop+1]

    # first search without considering  ambiguity characters
    r= re.search(motif, sequence)
    if r!=None: return r.start()+start


    motif = motif.replace('A', '##01##')
    motif = motif.replace('C', '##02##')
    motif = motif.replace('B', '##03##')
    motif = motif.replace('D', '##04##')
    motif = motif.replace('G', '##05##')
    motif = motif.replace('H', '##06##')
    motif = motif.replace('K', '##07##')
    motif = motif.replace('M', '##08##')
    motif = motif.replace('N', '##09##')
    motif = motif.replace('S', '##10##')
    motif = motif.replace('R', '##11##')
    motif = motif.replace('T', '##12##')
    motif = motif.replace('W', '##13##')
    motif = motif.replace('V', '##14##')
    motif = motif.replace('Y', '##15##')

    motif = motif.replace('##01##', '[ADHMNRWV?]')
    motif = motif.replace('##02##', '[CBHMNSVY?]')
    motif = motif.replace('##03##', '[CBGKNSTY?]')
    motif = motif.replace('##04##', '[ADGKNRTW?]')
    motif = motif.replace('##05##', '[BDGKNSRV?]')
    motif = motif.replace('##06##', '[ACHMNTWY?]')
    motif = motif.replace('##07##', '[BDGKNT?]')
    motif = motif.replace('##08##', '[ACHMNV?]')
    motif = motif.replace('##09##', '[ACBDGHKMNSRTWVY?]')
    motif = motif.replace('##10##', '[CBGNSV?]')
    motif = motif.replace('##11##', '[ADGNRV?]')
    motif = motif.replace('##12##', '[BDHKNTWY?]')
    motif = motif.replace('##13##', '[ADHNTW?]')
    motif = motif.replace('##14##', '[ACGMNSRV?]')
    motif = motif.replace('##15##', '[CBHNTY?]')

    motif = motif.replace('?', '.')
    motif = motif.replace('-', '[-?]')

    r= re.search(motif, sequence)

    if not r: return None
    return r.start()+start
    
########################################################################

def motifs(sequence, motif, mismatches=0, reverse=True):

    """
    Locates motifs in a nucleotide sequence. Standard ambiguity
    characters are supported (as explained in :func:`compare`
    documentation). *sequence* and *motif* are nucleotide sequence
    strings. *mismatches* gives the number of nucleotide differences
    allowed for motif match. If *reverse* is ``True``, both strands are
    examined (otherwise, only the forward strand is considered. Returns
    a list of hits. Each hit is represented by a dictionary containing
    keys: *start*: starting position of the hit, *sequence*: sequence of
    the matching region, *mismatches*: number of mismatches in the hit,
    *reverse*: ``True`` if the hit is on the reverse strand. The hit
    position and the found motif are always given with respect to the
    passed sequence, even when the motif was found on the reverse hit.
    """

    hits= []
    i=0
    while ((i+len(motif))<=len(sequence)):
        subject= sequence[i:i+len(motif)]
        c=0
        for j in range(len(motif)):
            c+=(not _mycmpnuc(motif[j], subject[j]))
        
        if (c<=mismatches):
            hits+=[{
                'position': i,
                'sequence': subject,
                'mismatches': c,
                'reverse': False
            }]
        i+=1
    if reverse:
        rcseq = rc(sequence)
        rchits= motifs(rcseq, motif, mismatches, False)
        for i in rchits:
            i['reverse'] = True
            i['position'] = len(sequence)-i['position']-len(motif)
            hits.append(i)
    return hits          

########################################################################

def backalign(nucseq, protseq, code=1):
    
    """
    Alignement of coding sequences based on aligned predicted products.
    Conceptual translations of DNA sequences must match exactly passed
    protein sequences (except for gaps). Stop codons are not supported.
    *nucseq* is a :class:`~egglib.Container` instance containing raw
    coding sequence. *protseq* is a :class:`~egglib.Align` instance
    containing align amino acid sequences. *code* specifies the genetic
    code; refer to the documentation of :meth:`~GeneticCodes.translate`.
    Returns a :class:`~egglib.Align` instance containing aligned coding
    sequences.
    """
        
    # checks that all nucleotide sequences are in the protein alignment and that the aa sequences are the same
    for n,s,g in nucseq:
        if n not in protseq:
            raise ValueError, 'sequence '+n+' not found in protein sequences'
        s1= s
        s2=''
        j=3
        while(j<=len(s1)):
            s2=s2 + GeneticCodes.translate(s1[j-3:j], code)
            j=j+3
            
        # tries to make a very useful error message when needed
        if (protseq.sequenceByName(n).replace('-','') != s2):
            x= protseq.sequenceByName(n).replace('-','')
            a=''
            a+= 'the protein encoded by '+n+' does not match\n'
            a+= 'predicted: "'+s2+'"\n'
            a+= 'protein:   "'+x+'"\n'

            if len(s2)!=len(x): a+= 'not the same length'
            else:
                diff=''
                for j in range(len(x)):
                    if s2[j]==x[j]: diff+='-'
                    else: diff+='X'
                a+=diff
            raise ValueError, a

    # inserts --- codons for gaps in aa sequences, and outputs
    niou = data.Container()
    for n,s,g in nucseq:
        oldnuc= s
        prot= protseq.sequenceByName(n)
        newnuc = ''
        j= 0
        for i in prot:
            if (i!='-'):
                newnuc +=  oldnuc[j]+oldnuc[j+1]+oldnuc[j+2]
                j=j+3
            else:
                newnuc += '---'
        niou.append(n, newnuc)
    aln= data.Align()
    aln.addSequences(niou)

    return aln

########################################################################

def chisquare(ddl):

    """
    Returns the 5% critical value of the chi-square distribution
    with *ddl* degrees of liberty (maximum: 100).
    """

    table=(   3.841,   5.991,   7.815,   9.488,  11.07,   12.592,  14.067,  15.507,  16.919,  18.307,
             19.675,  21.026,  22.362,  23.685,  24.996,  26.296,  27.587,  28.869,  30.144,  31.41,
             32.671,  33.924,  35.172,  36.415,  37.652,  38.885,  40.113,  41.337,  42.557,  43.773,
             44.985,  46.194,  47.4,    48.602,  49.802,  50.998,  52.192,  53.384,  54.572,  55.758,
             56.942,  58.124,  59.304,  60.481,  61.656,  62.83,   64.001,  65.171,  66.339,  67.505,
             68.669,  69.832,  70.993,  72.153,  73.311,  74.468,  75.624,  76.778,  77.931,  79.082,
             80.232,  81.381,  82.529,  83.675,  84.821,  85.965,  87.108,  88.25,   89.391,  90.531,
             91.67,   92.808,  93.945,  95.081,  96.217,  97.351,  98.484,  99.617, 100.749, 101.879,
            103.01,  104.139, 105.267, 106.395, 107.522, 108.648, 109.773, 110.898, 112.022, 113.145,
            114.268, 115.39,  116.511, 117.632, 118.752, 119.871, 120.99,  122.108, 123.225, 129.561)

    return table[ddl-1]

########################################################################

def ranges(values):
    
    """
    Identifies continuous ranges among the iterable *values*. *values*
    must be iterable but needs not to be sorted and can contain
    duplicates. The function returns a list of ``(start,stop)`` tuples,
    where *start* and *stop* defines a continuous range.
    """
    
    if not len(values): return []
    stack= [i for i in values]
    stack.sort()
    res=[[stack[0], stack[0]]]
    del stack[0]
    while len(stack):
          if (stack[0]==res[-1][1]+1): res[-1][1]+=1
          else: res.append([stack[0], stack[0]])
          del stack[0]
    return res

########################################################################

def stats(data):

    """
    Computes basic distribution statistics from the list *data*. Returns
    a tuple ``(m, V, sd, se)`` where *m* is the mean, *V* the unbiased
    variance, *sd* the unbiased standard deviation and *se* the
    unbiased standard error.
    
    .. note::
        Returns a tuple of zeros if data is empty.
    """

    n= len(data)
    if not n: return (0., 0., 0., 0.)
    SUM=0.
    SUM2=0.
    for i in data:
        SUM+=i
        SUM2+=i**2
    m= SUM/n
    v= SUM2/n - m**2
    sd= math.sqrt(math.fabs(v))
    se= sd/math.sqrt(n)
    return (m,v,sd,se)

########################################################################

def correl(x, y):

    """
    Computes correlation coefficients. *x* is a sequence giving the
    values of the explanatory variable. *y* is a sequence giving the
    values of the response variable. Returns a tuple ``(r, r**2, a)``
    where *r* is the correlation coefficient and *a* the regression
    coefficient.
    """

    if len(x)!=len(y):
        raise ValueError, 'correl expects two lists of the same size!'
    n=len(x)

    X=Y=SSDx=SSDy=SJD=0.

    X = 1. * sum(x) / n
    Y = 1. * sum(y) / n

    for i in range(n):
        SSDx+= (x[i]-X)**2
        SSDy+= (y[i]-Y)**2
        SJD+=  (y[i]-Y)*(x[i]-X)

    a = SJD/SSDx                                 # regression coefficient (regression of y by x)
    r = SJD/(math.sqrt(SSDx)*math.sqrt(SSDy))    # correlation coefficient
    
    return (r,r**2,a)

########################################################################

class Updater():
    
    """
    Helper designed to monitor progress of long-running tasks. In
    principle, :class:`~egglib.tools.Updater` should be coupled to a repetitive
    process with a fixed and known number of steps to perform (*target*)
    each requiring the same amount of time. :class:`~egglib.tools.Updater` should be
    updated regularly at reasonnable intervals (not to short to keep it
    from being itself a resource load).
    
    The class can be used as in the following examples::

        >>> import egglib
        >>> updater = egglib.tools.Updater(1000)
        >>> for i in range(1000):
        >>>     # ... time-consuming task here ...
        >>>     updater.refresh()
        >>> updater.close()
        >>> 
        >>> # The number of iterations is not known - custom display
        >>> import random
        >>> updater = egglib.tools.Updater()
        >>> maxi = 0
        >>> while True:
        >>>     X = random.random()
        >>>     if X>maxi:
        >>>         maxi = X
        >>>     updater.refresh('$ELAPSED, max. value: %f' %maxi)
        >>>     if X>0.99999:
        >>>         break
        >>> updater.refresh('$DONE tries, time $ELAPSED, got %f' %maxi)
        >>> updater.close()
    """
    
    ####################################################################

    def get_closed(self):
        return self._closed

    closed = property(get_closed, doc='``True`` if the instance has been \
    closed using the :meth:`close` method.')

    ####################################################################
    
    def __init__(self, target=None):
        
        """
        Constructor's argument *target* gives the number of iterations
        to perform. If ``None``, this information is not available.
        """
        
        self._closed = False
        self._target = target
        self._done = 0
        self._time_point = time.time()
        self._elapsed = 0
        self._last_elapsed = 0
        self._last_increment = 0
        
        self.length_max = 80
        self._printed = 0
        
        if target==None:
            self.template = '$DONE - $ELAPSED'
        else:
            self.template = '$DONE/$TARGET (remaining: $REMAINING)'

        # invariable formatting variable

        if self._target != None:
            self._TARGET = str(self._target)
            self._L = len(self._TARGET)
        else:
            self._TARGET = '?'
            self._L = 1

        self._refresh_cache = self._time_point
        self._refresh_cache_template = self.template
    
    ####################################################################
    
    def stats(self):
        
        """
        Returns a dictionary with the current values of counters.
        """
        
        return {
            'target': self._target,
            'done': self._done,
            'elapsed': self._elapsed,
            'last_elapsed': self._last_elapsed,
            'last_increment': self._last_increment
        }
    
    ####################################################################
    
    def increment(self, number):
        
        """
        Adds *number* steps and update elapsed time and estimated
        running (if *target* was given). *number* might be negative.
        """
        
        now = time.time()
        interval = now - self._time_point
        self._time_point = now
        self._elapsed += interval
        self._last_elapsed = interval
        self._last_increment = number
        
        self._done += number
        
    ####################################################################
    
    def format(self, template=None, increment=1):
        
        """
        Returns a string providing feedback about the run's progress. 
        *template* gives the template of the string to return. Actual
        values will substitute to the following special strings:

            - ``$DONE``: Number of steps done.
            - ``$TARGET``: Total number of steps to do.
            - ``$TODO``: Number of steps left to do.
            - ``$ELAPSED``: Time used since object creation.
            - ``$REMAINING``: Estimated time to complete the task.
            - ``$LREMAINING``: Like ``$REMAINING`` but computed from the last time point.
            - ``$TOTAL``: Estimated total time (computed as ``$REMAINING`` + ``$ELAPSED``)
            - ``$PERCENT``: Percentage done (including ``%`` symbol).

        If *template* is ``None``, a template defined at construction
        time is used. This template is ``$DONE|$ELAPSED`` when *target*
        is ``None`` and ``$DONE/$TARGET (remaining: $REMAINING)`` if *target*
        is specified. It is stored at the object attribute
        :attr:`template` and can be modified dynamically.
        
        If *increment* is different than zero, :meth:`increment` is
        called and this number is passed before formatting the string.
        """
        
        if increment!=0:
            self.increment(increment)

        # defines variables

        values = {}
        values['DONE'] = str(self._done).rjust(self._L, '0')
        values['TARGET'] = self._TARGET
        values['ELAPSED'] = self._tformat(self._elapsed)
        
        if self._target != None and self._target>0:
            values['TODO'] = str(self._target - self._done).rjust(self._L, '0')
            values['PERCENT'] = '%.2f%%' %(100.*self._done/self._target)
            
            if self._done==0:
                values['REMAINING'] = 'inf.'
                values['LREMAINING'] = 'inf.'
                values['TOTAL'] = 'inf.'
            elif self._done>self._target:
                values['REMAINING'] = self._tformat(0)
                values['LREMAINING'] = self._tformat(0)
                values['TOTAL'] = self._tformat(self._elapsed)
            else:
                est = (1.*self._elapsed/self._done) * (self._target-self._done) 
                est_L = (1.*self._last_elapsed/self._last_increment) * (self._target-self._done)

                values['REMAINING'] = self._tformat(est)
                values['LREMAINING'] = self._tformat(est_L)
                values['TOTAL'] = self._tformat(est + self._elapsed)

        else:
            values['REMAINING'] = '?'
            values['LREMAINING'] = '?'
            values['PERCENT'] = '?%'
            values['TODO'] = '?'
            values['TOTAL'] = '?'
            
        # substitutes variables
        if template == None:
            template = self.template
        for i in values:
            template = template.replace('$%s' %i, values[i])
        
        return template
        
    ####################################################################

    def _tformat(self, elapsed):
        
        sec= int(elapsed)
        mn= sec//60
        sec -= mn*60
        h = mn//60
        mn-= h*60
        d = h//24
        h-= d*24
        
        string=''
        if d==1: string+='1 day'
        elif d>1: string+='%d days' %d
        if d>0 or h>0: string+= ' %s:' %str(h).rjust(2, '0')
        string+= '%s:%s' %(str(mn).rjust(2, '0'),str(sec).rjust(2, '0'))

        return string.strip()

    ####################################################################

    def wipe(self):
        
        """
        Writes an empty line of the maximal possible length, therefore
        clearing the line from characters printed by another process
        (provided that these characters are not too many).
        """
        
        sys.stdout.write('\10'*self.length_max)
        self._printed = self.length_max

    ####################################################################

    def refresh(self, template=None, increment=1, grain=1.0):
        
        """
        This method generates a string exactly as :meth:`format` does,
        but writes the string to :class:`sys.stdout` instead. If the
        same object has already wrote anything, an equivalent number
        of backspaces are written, in principle allowing to overwrite
        the previous string and making the string appear to update
        itself. The result might no be so nice if something else is
        written to the console in the mean time or if the console
        doesn't support backspaces. This method doesn't write a newline,
        but the object will upon destruction or call to :meth:`close`.
        The string will be stripped is it is longer than the object
        attribute :attr:`length_max`, which can be changed dynamically.
        If less than the number given by *grain* (in seconds) has
        occurred since the last refresh, nothing is printed.
        """
        
        if self._closed:
            raise ValueError, 'refresh operation cannont be applied to \
            a closed `Updater` instance'

        string = self.format(template, increment)
        self._refresh_cache_template = template
        
        if self._time_point > ( self._refresh_cache + grain):

            string = string[:self.length_max]
            
            sys.stdout.write('\10'*self._printed)

            if len(string) > self._printed:
                sys.stdout.write(string)
                self._printed = len(string)
            else:
                sys.stdout.write(string.ljust(self._printed))

            sys.stdout.flush()
            self._refresh_cache = self._time_point
    
    ####################################################################

    def close(self):
        
        """
        If anything was written using :meth:`refresh`, writes any cached
        refresh data and writes a new line. Otherwise, does nothing.
        This method is automatically called upon object destruction.
        """
        
        if self._printed>0:
            self.refresh(self._refresh_cache_template, 0, 0.)
            self._printed = 0
            sys.stdout.write('\n')

        self._closed = True

    ####################################################################
    
    def __del__(self):
        
        self.close()

########################################################################

def wrap(string, length, indent=0):
    
    """
    Formats the string *string*  to ensures the line lengths are not
    larger than *length*. The optional argument *indent* specifies the
    number of spaces to insert at the beginning of all lines except the
    first. The line breaks are inserted at spaces.
    
    An example is given below::
        >>> import egglib
        >>> string = "Lekrrjf djdhs eeir djs ehehf bnreh eurvz rhffdvfu dksgta."
        >>> print egglib.tools.wrap(string, 20, 4)
        Lekrrjf djdhs eeir
            djs ehehf bnreh
            eurvz rhffdvfu
            dksgta.

    
    """
        
    res = ''
        
    cache = []
    mark = 0
    for i in string:
        if i==' ' and set(cache)==' ': continue
        cache.append(i)
        if len(cache)<length and i=='\n':
            res += ''.join(cache)
            cache = []
            
        if len(cache)>=length:
            cache = ''.join(cache)
            pos = cache[mark:].rfind(' ')
            if pos<indent:
                write = cache
                cache = list(' '*indent)
            else:
                pos += mark
                write = cache[:pos]
                cache = list(' '*indent + cache[pos+1:].lstrip())
                mark = indent
                
            res += write + '\n'
    res += ''.join(cache)
    return res

########################################################################

class Bin:
    
    """
    Discretizes a multivariate distributions. ``len(Bin)`` returns
    the number of observations.
    """
    
    ####################################################################
    
    def __init__(self, data, ranges=None):
        
        """
        *data* must be a list of lists of parameters. Each
        first-dimension list (row) defines a parameter. To avoid
        confusion, note that datasets should contain a small number of
        rows and a large number of columns, and therefore that 
        ``len(data)`` should be a small number and that
        ``len(data[i]))`` should return a large number. All lists must
        have the same length. Only a reference to the  data is stored,
        therefore the data should not be changed while :class:`~egglib.tools.Bin` is
        using it. *ranges* allows to define the limits of the space
        region where the points were sampled. If specified, *ranges*
        must be a list of ``(bot,top)`` tuples giving the bottom and top
        bounds of the distribution. *bot* and *top* must be
        (respectively) lower than or equal to the smallest point and
        higher than or equal to the largest point in the distribution.
        Whenever this condition is not respected, binning the
        distribution will result in data loss. If *ranges* is not
        specified, the extreme points will be used to define an
        empirical region.
        """
        
        if len(data) == 0:
            raise ValueError, 'cannot instantiate Bin: empty data'

        lens = set(map(len, data))
        if len(lens) != 1:
            raise ValueError, 'Bin initialized with invalid data (lists of different lengths)'
                
        self._n = lens.pop()
        self._data = data
                
        self._ranges = None

        if ranges!=None:
            
            if len(ranges) != len(self._data):
                raise ValueError, '`ranges` must match `data` in length'
            
            self._ranges = list(ranges)
            
        elif self._n != 0:
            
            self._ranges = []
            for i in self._data:
                self._ranges.append( (min(i), max(i)) )
                
    ####################################################################

    def __len__(self):

        return self._n
        
    ####################################################################
    
    def ranges(self):
        
        """
        Returns a reference to the list of ``(min,max)`` tuples (one
        for each parameter) defining the region in which the bin was
        defined.
        """
        
        return self._ranges
        
    ####################################################################
            
    def extract(self, dim, bot, top):
        
        """
        Returns a new :class:`~egglib.tools.Bin` instance with only values falling
        between *bot* and *top* (excluding *top*) of dimension *dim*.
        """
        
        if dim>=len(self._data): raise ValueError, 'Bin: invalid index: %d (maximum %d)' %(dim, len(self._data)-1)
        
        data = [[] for i in self._data]

        for i in range(self._n):
            
            if self._data[dim][i] >= bot and self._data[dim][i] < top:

                for j in range(len(self._data)):
                    data[j].append(self._data[j][i])

        ranges = list(self._ranges)
        ranges[dim] = (bot, top)

        return Bin(data, ranges)

    ####################################################################

    def slice(self, dim, ncat=8, bot=None, top=None):
        
        """
        Performs automatic (and constant-interval) binarization of
        dimension *dim*. *ncat* gives the number of intervals. *bot* and
        *top* give the area to discretize (if ``None``, use the
        empirical bounds of the distribution). Returns a list of 
        :class:`~egglib.tools.Bin` instances.
        """
        
        if dim>=len(self._data): raise ValueError, 'Bin: invalid index: %d (maximum %d)' %(dim, len(self._data)-1)

        smalldiff= 0.0000000000001

        # escape route
        #if not len(self):
        #    return []

        # auto-defines the bounds

        if bot==None:
            bot = self._ranges[dim][0]

        if top==None:
            top = self._ranges[dim][1]
        
        # defines limits

        def mklim(x):
            return bot + x*(1.*top-bot)/ncat

        lim= map(mklim, range(1, ncat))
        lim = [bot] + lim + [top+smalldiff]

        # splits itself
        res = [self.extract(dim, lim[i], lim[i+1]) for i in range(ncat)]

        # done
        return res

########################################################################

def LD(align1, align2, shuffle):
        
    """
    Computes linkage disequilibrium statistics between two :class:`~data.Align`
    instances *align1* and *align2*. If *shuffle* is ``True``, randomly
    shuffles sequences of *align2* (without altering the original
    instance), emulating the hypothesis of linkage equilibrium. Returns
    a ``(n1,n2, S1,S2,K1,K2,D,Dp)`` tuple, where *n1* is the number of
    used sequences of the first alignment, *S1* is the number of
    polymorphic sites of the first alignment, *K1* is the number of
    unique haplotypes of the first alignment, *D* is the standard
    estimator of linkage disequilibrium and *Dp* is Lewontin's estimator
    (bound by 0 and 1).
    """

    # trick to sort sequences

    copy1 = [i for i in align1]
    copy2 = [i for i in align2]
    copy1.sort(lambda x,y: cmp(x[0], y[0]))
    copy2.sort(lambda x,y: cmp(x[0], y[0]))

    # makes internal alignments (removes all ougroups)
    
    _align1 = data.Align()
    _align2 = data.Align()

    for n,s,g in copy1:
        if g!=999:
            _align1.append(n, s, g)

    for n,s,g in copy2:
        if g != 999:
            _align2.append(n, s, g)


    if sorted(_align1.names())!=sorted(_align2.names()):
        raise ValueError, 'cannot compute linkage disequilibrium statistics: alignment names must match'
    if _align1.contains_duplicates():
        raise ValueError, 'cannot compute linkage disequilibrium statistics: alignments must not contain duplicates'

    # shuffle, if requested

    if shuffle:
        sequences= [s for n,s,g in _align2]
        random.shuffle(sequences)
        for i,v in enumerate(sequences):
            _align2.sequence(i, v)
    
    # computes polymorphism    

    pol1 = _align1.polymorphism()
    pol2 = _align2.polymorphism()
    
    S1 = pol1['S']
    S2 = pol2['S']
    K1 = pol1['K']
    K2 = pol2['K']
    
    N = len(align1)
    
    # finds (multi-locus) haplotypes

    haplotypes = [ '%d%d' %(pol1['alleles'][i], pol2['alleles'][i]) 
                                                    for i in range(N) ]
    
    # computes frequencies (genotype frequencies are p, haplotype frequencies are P)

    p1 = [ 1.*pol1['alleles'].count(i)/N for i in range(K1) ]
    p2 = [ 1.*pol2['alleles'].count(i)/N for i in range(K2) ]
    P = [ [ 1.*haplotypes.count('%d%d' %(i,j))/N for j in range(K2) ] for i in range(K1) ]

    # computes LD

    totD = 0.
    totDp = 0.
    for i in range(K1):
        for j in range(K2):
            D = P[i][j] - p1[i]*p2[j]
            Dmax = min(p1[i]*(1-p2[j]), (1-p1[i])*p2[j])
            Dmin = min(p1[i]*p2[j], (1-p1[i])*(1-p2[j]))
            if D>=0.: Dp = 1.*D/Dmax
            else: Dp = 1.*D/Dmin
            totD += abs(D) * p1[i]*p2[j]
            totDp += abs(Dp) * p1[i]*p2[j]
    
    return len(_align1), len(_align2), S1, S2, K1, K2, totD, totDp

########################################################################

def ungap(align, freq, includeOutgroup=True):
    
    """
    Builds a new :class:`~egglib.Align` instance containing all
    sequences of  *align* and only the columns for which the frequency
    of gaps (``-`` symbols) is less than the value given by *freq*.
    
    If *includeOutgroup* is ``True``, the sequences with group label 999
    (if any) are considered for computing the frequency of gaps. These
    sequences are however always exported to the returned alignment).
    
    .. versionchanged:: 2.1.0
       Added option *includeOutgroup*.
    """
    
    indivs = [(i, v.group) for (i,v) in enumerate(align)
                            if (includeOutgroup or v.group!=999)]
                            
    if len(indivs)<1:
        raise ValueError, 'cannot ungap: not enough sequences in alignment'
    
    good = []
    
    for pos in range(align.ls()):
        
        col = align.column(pos)
        col = [col[i] for (i,j) in indivs]
        
        gaps = 1. * col.count('-') / len(col)
        
        if gaps < freq:
            good.append(pos)
        
    return align.extract(good)
    
########################################################################

def ungap_all(align):
    
    """
    Removes all gaps (``-`` symbols) from all sequences of the :class:`~data.Align`
    instance *align* and returns the resulting ~data.Container instance.
    """
    
    result = data.Container()
    
    for n,s,g in align:
        result.append(n, s.translate(None, '-'), g)

    return result

########################################################################

def ungap_triplets(align, freq):
    
    """
    Like :meth:`~egglib.tools.ungap` but remove triplet positions
    instead of single sites. A triplet of containing at least one gap
    (``-`` symbol) is held as missing. Triplet positions for which the
    proportion of missing triplets is larger than or equal to *freq* are
    completely removed. The :class:`~data.Align` instance passed as
    *align* must have a length multiple of 3. Returns the resulting
    :class:`~data.Align` instance.

    .. versionadded:: 2.1.0
    """
    
    if align.ls()%3 != 0:
        raise ValueError, 'tools.ungap_triplet: invalid Align instance (length should be a multiple of 3)'
    
    good = []
    
    for i in range(0, align.ls(), 3):
        
        col1 = align.column(i)
        col2 = align.column(i+1)
        col3 = align.column(i+2)
        
        gaps = 1. * sum(map(lambda x: '-' in x, zip(col1, col2, col3))) / len(align)
        
        if gaps < freq:
            good += [i, i+1, i+2]
        
    return align.extract(good)




########################################################################

class Fstats(object):

    """
    This is a prototype for generic computation of F-statistics according
    to the Weir and Cockerham (1984) method -- multi-locus version
    
    The two-level statistics have yet to be validated, use with care.
    
    Instances of this class have a length: the number of of analyzed
    loci.
    """

    none = 0
    first = 1
    second = 2
    
    ####################################################################

    def __init__(self, hierarchical, level=none):

        """
        Constructor arguments (will be applied to all analyzed loci):

            - *hierarchical*: a boolean indicating whether data are 
              hierarchical (two nested levels of structure) and should
              be analyzed as such.

            - *level*: if data have actually two levels and
              *hierarchical* is ``False``, this option allows to specify
              which level of structure should be considered. The only
              accepted values are :member:`Fstats.none` (or 0),
              :member:`Fstats.first` (or the integer 1),
              :member:`Fstats.second` (or the integer 2).
              :member:`Fstats.none` must be used when data are not
              hierarchical and therefore collapsing is not required.
              :member:`Fstats.first` should be used to collapse/ignore
              the inner level (clusters, or regions, are considered) and
              :member:`Fstats.second` shoudl be used to collapse/ignore
              the outer level (local demes are considered. This option
              is ignored if *hierarchical* is ``True``.
        """

        self.loci = []
        self.hierarchical = hierarchical
        self.level = level

    ####################################################################

    def compute(self, data):
        
        """
        Doc tba
        """
        
        if self.hierarchical==False:
            if self.level==self.none:
                pass
            elif self.level==self.first:
                data = [reduce(list.__add__, i) for i in data]
            elif self.level==self.second:
                data = reduce(list.__add__, data)
            else:
                raise ValueError, ('invalid level option value (%s)'
                                                       %str(self.level))
        
        self.loci.append(Fstats1(data, self.hierarchical))
        
        if len(self.loci)>1:
            self.check_structure(self.loci[0], self.loci[-1])

    ####################################################################

    def __len__(self):
        return len(self.loci)

    ####################################################################

    def fstats(self):
        
        if self.hierarchical:

            a = sum([sum(i.a) for i in self.loci])
            b1 = sum([sum(i.b1) for i in self.loci])
            b2 = sum([sum(i.b2) for i in self.loci])
            c = sum([sum(i.c) for i in self.loci])
      
            Fis, Fsc, Fst, Fct, Fit = Fstats1.get_fstats(
                                        pop=a, dem=b2, ind=b1, gam=c)
                
            return {
                'Fis': Fis,
                'Fsc': Fsc,
                'Fst': Fst,
                'Fct': Fct,
                'Fit': Fit
            }

        else:

            a = sum([sum(i.a) for i in self.loci])
            b = sum([sum(i.b) for i in self.loci])
            c = sum([sum(i.c) for i in self.loci])
      
            Fis, Fst, Fit = Fstats1.get_fstats(pop=a, ind=b, gam=c)
            
            return {
                'Fis': Fis,
                'Fst': Fst,
                'Fit': Fit
            }


    ####################################################################

    def check_structure(self, locus1, locus2):
        
        flag = True
        flag &= locus1.r == locus2.r
        flag &= locus1.mi == locus2.mi
        flag &= locus1.m == locus2.m

        # (sample size per locus is allowed to vary)
        


########################################################################

class Fstats1(object):
    
    """
    This is a prototype for generic computation of F-statistics according
    to the Weir and Cockerham (1984) method
    """

    
    ####################################################################

    def __init__(self, data, hierarchical):
        
        """
        To compute F-statistics, you must pass a list of genotypes.
        Structure of input data tba.
        
        Resulting object structure tba
        
        Single locus version
        """
    
        self.hierarchical = hierarchical
        self.get_structure(data)

        self.alleles = sorted(self.alleles)
        
        assert len(self.alleles)>1
        if len(self.alleles)==2:
            self.get_variance(data, self.alleles[0])
        else:
            for allele in self.alleles:
                self.get_variance(data, allele)

        if self.hierarchical==True:

            a = sum(self.a)
            b1 = sum(self.b1)
            b2 = sum(self.b2)
            c = sum(self.c)
      
            Fis, Fsc, Fst, Fct, Fit = self.get_fstats(
                                           pop=a, dem=b2, ind=b1, gam=c)
       
            self.Fis = Fis, self.Fis
            self.Fsc = Fsc, self.Fsc
            self.Fst = Fst, self.Fst
            self.Fct = Fct, self.Fct

        else:

            a = sum(self.a)
            b = sum(self.b)
            c = sum(self.c)
      
            Fis, Fst, Fit = self.get_fstats(pop=a, ind=b, gam=c)
       
            self.Fis = Fis, self.Fis
            self.Fst = Fst, self.Fst
            self.Fit = Fit, self.Fit

    
    ####################################################################
    
    def clear(self):
        
        """
        Resets the instance
        """
        
        self.r = None
        self.mi = None
        self.m = None
        self.nij = None
        self.ni = None
        self.n = None
        self.nbar = None
        self.nc = None
        self.pij = []
        self.pi = []
        self.p = []
        self.pbar = []
        self.s2 = []
        self.hi = []
        self.hij = []
        self.hbar = []
        self.n1 = []
        self.n2 = []
        self.n3 = []
        self.MSP = []
        self.MSD = []
        self.MSI = []
        self.MSG = []
        self.a = []
        self.b = []
        self.b1 = []
        self.b2 = []
        self.c = []
        self.Fis = []
        self.Fsc = []
        self.Fst = []
        self.Fct = []
        self.Fit = []
        self.alleles = set()

        
    ####################################################################


    def get_structure(self, data):
        
        """
        Resets the instance and get population structure
        """
        
        self.clear()

        if self.hierarchical: self.get_structure_hierarchical(data)
        else:  self.get_structure_standard(data)
            
            
            
            
    ####################################################################

    def get_structure_hierarchical(self, data):

        # number of populations
        self.r = len(data)

        # number of demes per population
        self.mi = map(len, data)

        # number of samples per deme per population (nij)
        # list of all alleles (alleles)
            
        self.nij = []
        for pop in data:    
            self.nij.append([])
            for deme in pop:
                self.nij[-1].append(len(deme))
                for i in deme:
                    self.get_indiv(i)

        if len(self.alleles)<2:
            raise ValueError, 'cannot compute F-statistics - the locus is fixed'

        # sums

        self.m = sum(self.mi)
        self.ni = map(sum, self.nij)
        self.n = sum(self.ni)

        # summary statistics

        self.n1 = 0.
        for i in range(self.r):
            for j in range(self.mi[i]):
                self.n1 += ((self.n - self.ni[i]) * self.nij[i][j]**2) / (self.ni[i]*self.n)
        self.n1 /= (self.r - 1)
        
        self.n2 = 0.
        for i in range(self.r):
            self.n2 += self.ni[i] ** 2
        self.n2 = (self.n - self.n2/self.n) / (self.r-1)
        
        self.n3 = 0.
        for i in range(self.r):
            for j in range(self.mi[i]):
                self.n3 += self.nij[i][j]**2 / self.ni[i]
        self.n3 = (self.n - self.n3) / (self.m - self.r)
        
        
        

    ####################################################################

    def get_structure_standard(self, data):

        # number of populations
        self.r = len(data)

        # number of individuals per population and total number
        self.ni = map(len, data)
        self.n = sum(self.ni)

        # list of all alleles (alleles)
            
        for pop in data:
            for i in pop:
                self.get_indiv(i)

        if len(self.alleles)<2:
            raise ValueError, 'cannot compute F-statistics - the locus is fixed'

        # summary statistics

        self.nbar = 1. * self.n / self.r
        sumsq = sum([i**2 for i in self.ni])
        self.nc = ( (self.r * self.nbar - 1.*sumsq/(self.r*self.nbar)) /
                     (self.r - 1) )



    ####################################################################
    
    def get_indiv(self, i):

        try: a, b = i
        except TypeError: raise ValueError, 'invalid population structure (possibly, not hierarchical)'
        except ValueError: raise ValueError, 'invalid population structure (possibly, not diploid)'
        if not isinstance(a, int): raise ValueError, 'invalid allele type: %s (must be int)' %a
        if not isinstance(b, int): raise ValueError, 'invalid allele type: %s (must be int)' %b

        self.alleles.add(a)
        self.alleles.add(b)        
        

    ####################################################################

    def get_variance(self, data, allele):
        
        """
        DOESN'T reset the instance and compute allele frequency and
        F-statistics for a given allele
        """

        if self.hierarchical: self.get_variance_hierarchical(data, allele)
        else: self.get_variance_standard(data, allele)



    ####################################################################

    def get_variance_standard(self, data, allele):

        # allele frequency
        
        pi = []
        for pop in data:
            pop = reduce(list.__add__, pop)
            x = 1. * pop.count(allele) / len(pop)
            pi.append(x)

        # average frequency

        pbar = 0
        for p, n in zip(pi, self.ni):
            pbar += p*n
        pbar /= self.r * self.nbar
        
        # sample variance
        
        s2 = 0
        for p, n in zip(pi, self.ni):
            s2 += n * (p-pbar) ** 2
        s2 /= (self.r - 1) * self.nbar

        # heterozygote frequency
        
        hi = []
        for pop in data:
            x = 1. * [i.count(1) for i in pop].count(1) / len(pop)
            hi.append(x)
            
        # average heretozygote frequency
        
        hbar = 0.
        for h, n in zip(hi, self.ni):
            hbar += n * h
        hbar /= self.r * self.nbar

        # variance components
        
        a = (self.nbar / self.nc) * (s2 - (1/(self.nbar - 1)) * (
             pbar * (1-pbar) - ((self.r - 1.)/self.r) * s2 - hbar / 4))
             
        b = (self.nbar / (self.nbar - 1)) * (pbar*(1-pbar) -
             ((self.r - 1.) / self.r) * s2 - 
             ((2*self.nbar - 1)/(4*self.nbar)) * hbar)
             
        c = hbar / 2

        # f-statistics

        Fis, Fst, Fit = self.get_fstats(pop=a, ind=b, gam=c)

        # export data
        
        self.a.append    ( a    )
        self.b.append    ( b    )
        self.c.append    ( c    )
        self.pi.append   ( pi   )
        self.pbar.append ( pbar )
        self.s2.append   ( s2   )
        self.hi.append   ( hi   )
        self.hbar.append ( hbar )
        self.Fis.append  ( Fis  )
        self.Fst.append  ( Fst  )
        self.Fit.append  ( Fit  )
        




    ####################################################################

    def get_variance_hierarchical(self, data, allele):
        
        # frequency of the allele per deme per population

        pij = []
        for pop in data:
            pij.append([])
            for deme in pop:
                deme = reduce(list.__add__, deme)
                x = 1. * deme.count(allele) / len(deme)
                pij[-1].append(x)
        
        # frequency per population

        pi = []
        for i in range(self.r):
            pi.append(0.)
            for j in range(self.mi[i]):
                pi[-1] += self.nij[i][j] * pij[i][j]
            pi[-1] /= self.ni[i]
            ### ERROR HERE IF EMPTY CLUSTER (empirical)
        
        # frequency total

        p = 0.
        for i in range(self.r):
            for j in range(self.mi[i]):
                p += self.nij[i][j] * pij[i][j]
        p /= self.n
        
        if p==0 or p==1: raise ValueError, 'cannot compute: this allele is fixed'

        # frequency of heterozygotes per deme per population

        hij = []
        for pop in data:
            hij.append([])
            for deme in pop:
                x = 1. * [i.count(allele) for i in deme].count(1) / len(deme)
                hij[-1].append(x)
                
        # mean square deviations

        MSP = 0.
        for i in range(self.r):
            MSP += self.ni[i] * (pi[i] - p)**2
        MSP /= self.r - 1
        MSP *= 2.
        
        MSD = 0.
        for i in range(self.r):
            for j in range(self.mi[i]):
                MSD += self.nij[i][j] * (pij[i][j] - pi[i])**2
        MSD /= self.m - self.r
        MSD *= 2.
        
        MSI_1 = 0.
        MSI_2 = 0.
        for i in range(self.r):
            for j in range(self.mi[i]):
                MSI_1 += self.nij[i][j] * pij[i][j] * (1-pij[i][j])
                MSI_2 += self.nij[i][j] * hij[i][j]
        MSI = (2.*MSI_1 - MSI_2/2.) / (self.n - self.m)

        MSG = 0.
        for i in range(self.r):
            for j in range(self.mi[i]):
                MSG += self.nij[i][j] * hij[i][j]
        MSG /= 2 * self.n

        # variance components

        a = ((self.n3 * MSP - self.n1 * MSD - (self.n3 - self.n1)*MSI) /
             (2. * self.n2 * self.n3))
        b2 = (MSD - MSI) / (2. * self.n3)
        b1 = (MSI - MSG) / 2.
        c = MSG
        
        # fixation indices
        
        Fis, Fsc, Fst, Fct, Fit = self.get_fstats(
                                        pop=a, dem=b2, ind=b1, gam=c)

        self.pij.append (  pij )
        self.pi.append  (  pi  )
        self.p.append   (  p   )
        self.hij.append (  hij )
        self.MSP.append (  MSP )
        self.MSD.append (  MSD )
        self.MSI.append (  MSI )
        self.MSG.append (  MSG )
        self.a.append   (  a   )
        self.b1.append  (  b1  )
        self.b2.append  (  b2  )
        self.c.append   (  c   )
        self.Fis.append (  Fis )
        self.Fsc.append (  Fsc )
        self.Fst.append (  Fst )
        self.Fct.append (  Fct )
        self.Fit.append (  Fit )

    ####################################################################
    
    @classmethod
    def get_fstats(cls, pop=None, dem=None, ind=None, gam=None):
        
        # hierarchical
        
        if pop!=None and dem!=None and ind!=None and gam!=None:
            a=pop
            b2=dem
            b1=ind
            c=gam
        
            if c+b1 == 0:
                Fis = None
            else:
                Fis = 1 - c / (c + b1)

            if b1 + b2 + c == 0:
                Fsc = None
            else:
                Fsc = b2 / (b2 + b1 + c)
            
            if a + b2 + b1 + c == 0:
                Fct = None
                Fst = None
                Fit = None
            else:
                Fct = a / (a + b2 + b1 + c)
                Fst = (a + b2) / (a + b2 + b1 + c)
                Fit = 1 - c / (a + b2 + b1 + c)
            
            return Fis, Fsc, Fst, Fct, Fit

        # standard

        if pop!=None and dem==None and ind!=None and gam!=None:
            a=pop
            b=ind
            c=gam

            if c+b == 0:
                Fis = None
            else:
                Fis = 1 - c / (c + b)

            if a + b + c == 0:
                Fst = None
                Fit = None
            else:
                Fst = a / (a + b + c)
                Fit = 1 - c / (a + b + c)
            
            return Fis, Fst, Fit


        raise ValueError, 'invalid combination of variance components'



########################################################################

