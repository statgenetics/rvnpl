"""
This module contains wrappers to external applications. They must be
available on the user's system and detected properly at installation
time. To detect a new application *a posteriori*, one needs to relaunch
the detection procedure by typing ``python setup.py build_apps`` from
the egglib-py directory, and then re-install the configuration file
by typing ``python setup.py install``.

"""


__license__ = """
    Copyright 2009-2011 Stephane De Mita, Mathieu Siol

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


import os, xml.dom.minidom, subprocess, sys, tempfile, math, re
import data, tools, egglib_binding

########################################################################

apps = {}
"""
Mapping of the external applications commands.

This dictionnary is set as module initialization time from the data
specified during EggLib installation. It can be modified dynamically if
needed.

Each entry of the dictionary is one of the external commands possibly
used by :mod:`egglib`. If the value is ``None``, :mod:`egglib` will
consider that the command is not available. Otherwise, the value will be
used as command to launch the program. The mapping is initialized at
load time from the file ``apps.conf`` normally generated at build time.
Each line of the file represent a single application and the format of
the line is ``application@command`` where *application* is a name
identifying the application usually the name of the command and *command*
the exact command (including absolute path whenever necessary) to run
the program. An empty string (``application@``) or the command name *None*
(``application@None``) are replaced by ``None`` , what is interpreted as
the command being not available. 
"""

f=open(os.path.join(os.path.dirname( sys.modules[__name__].__file__ ), 'apps.conf'))
for line in f:
    line= line.strip().split('@')
    if len(line)==2:
        app = line[0].strip()
        apps[app] = line[1].strip()
        if not len(apps[app]): apps[app] = None
        if apps[app] == 'None': apps[app] = None
f.close()
del f

########################################################################

def ms(nsam, howmany, theta=None, segsites=None, T=False, F=False,
      r=False, c=False, G=False, I=False, n=False, g=False, m=False, 
      ma=False, eG=False, eg=False, eM=False,
      em=False, ema=False, eN=False, en=False, es=False, ej=False,tMRCA=False):

    """
    Runs the program ``ms`` to generate random datasets by coalescence.
    ``ms`` must be installed in the system. Arguments are for the
    command ``ms`` (refer to the program's documentation for details).
    Note that all options starting by ``e`` (past demographic changes)
    as well a ``m``, ``n`` and ``g``
    expect a list of tuples (at least one), each tuple containing the
    appropriate number of arguments. Note also that the options are
    processed in the same order as in the function's signature. Note
    that if the *tMRCA* argument is set to ``True``, the returned
    alignments will contain a ``tMRCA`` member. If both ``theta`` and
    ``segsites`` are specified to positive values, the returned
    alignments will contain a ``prob`` member. If the *T* flag is sets,
    the returned alignments will contain a ``tree`` member (that will be
    a :class:`~egglib.Tree` instance).

    .. versionadded:: 0.1
       Created to provide a closer wrapper of ``ms``.

    .. versionchanged:: 2.1.0
       Exported alignment might contain a ``prob`` and/or a ``trees`` member.
    """

    # parses argments and makes the command line
    cmd = []

    if apps['ms']==None: raise RuntimeError, 'application ms is required for this feature'
    cmd.append( apps['ms'] )

    # nsam
    if not isinstance(nsam, int) or nsam<2:
        raise ValueError, 'ms: invalid nsam'
    cmd.append( str(nsam) )

    # howmany
    if not isinstance(howmany, int) or howmany<1:
        raise ValueError, 'ms: invalid howmany'
    cmd.append( str(howmany) )

    # Theta and segsites

    if theta==None and segsites==None:
        raise ValueError, 'ms: at least one of theta and segsites is required'

    if theta!=None:
        if not isinstance(theta, (int,float)):
            raise ValueError, 'ms: invalid theta'
        if theta<0.00000000000000000001:
            raise ValueError, 'ms: theta must not be zero'
        cmd.append( '-t' )
        cmd.append( str(theta) )

    if segsites!=None:
        if not isinstance(segsites,int):
            raise ValueError, 'ms: invalid segsites'
        if segsites==0:
            raise ValueError, 'ms: S must not be zero'
        cmd.append('-s')
        cmd.append(str(segsites))

    # T
    if T: cmd.append('-T')

    # F
    if F:
        if not isinstance(F, int): raise ValueError, 'ms: invalid F'
        else:
            cmd.append('-F')
            cmd.append(str(F))

    # r
    if r:
        try: a, b  = r
        except ValueError: raise ValueError, 'ms: invalid r'
        cmd.append('-r')
        cmd.append(str(r[0]))
        cmd.append(str(r[1]))

    # c
    if not c: pass
    elif not r: raise ValueError,'ms: -c option requires -r'
    elif not isinstance(c, (tuple, list)): raise ValueError, 'ms: invalid c'
    else:
        cmd.append('-c')
        cmd.append(str(c[0]))
        cmd.append(str(c[1]))

    # G
    if not G: pass
    elif not isinstance(G, (float, int)): raise ValueError, 'ms: invalid G'
    else:
        cmd.append('-G')
        cmd.append(str(G))

    # I
    if not I: pass
    elif not ('__iter__' in dir(I)):
        raise ValueError, 'ms: invalid I (must be a list or a tuple)'
    else:
        cmd.append('-I')
        for i in I: cmd.append(str(i))

    # m
    if m!=False:
        for i in m:
            try:
                a, b, c = i
            except ValueError: raise ValueError, 'ms: invalid m argument'
            a+=1
            b+=1
            cmd += ['-m', str(a), str(b), str(c)]

    # g
    if g!=False:
        for i in g:
            cmd.append('-g')
            cmd += map(str, i)
        
    # n
    if n!=False:
        for i in n:
            cmd.append('-n')
            cmd += map(str, i)

    # ma
    if ma!=False:
        if not I: raise ValueError, 'ms: ma option requires -I'
        else:
            ma = list(ma)
            try:
                ma[1] += 1
                ma[2] += 1
            except ValueError: raise ValueError, 'ms: invalid ma argument'
            cmd.append('-ma')
            for i in ma: cmd.append(str(i))

    # es
    if es!=False:
        for i in es:
            i = list(i)
            try:  i[1]+=1
            except ValueError: raise ValueError, 'ms: invalid es argument'
            cmd.append('-es')
            cmd += map(str, i)

    # ej
    if ej!=False:
        for i in ej:
            i = list(i)
            try:
                i[1]+=1
                i[2]+=1
            except ValueError: raise ValueError, 'ms: invalid ej argument'
            cmd.append('-ej')
            cmd += map(str, i)


    # eG
    if eG!=False:
        for i in eG:
            cmd.append('-eG')
            cmd += map(str, i)

    # eg
    if eg!=False:
        for i in eg:
            i = list(i)
            try:  i[1]+=1
            except ValueError: raise ValueError, 'ms: invalid eg argument'
            cmd.append('-eg')
            cmd += map(str, i)

    # eM
    if eM!=False:
        for i in eM:
            cmd.append('-eM')
            cmd += map(str, i)
    
    # em
    if em!=False:
        for i in em:
            i = list(i)
            try:
                i[1]+=1
                i[2]+=1
            except ValueError: raise ValueError, 'ms: invalid em argument'
            cmd.append('-em')
            cmd += map(str, i)

    # ema
    if ema!=False:
        for i in ema:
            cmd.append('-ema')
            cmd += map(str, i)
            
    # eN
    if eN!=False:
        for i in eN:
            cmd.append('-eN')
            cmd += map(str, i)
            
    # en
    if en!=False:
        for i in en:
            i = list(i)
            try:
                i[1]+=1
            except ValueError: raise ValueError, 'ms: invalid en argument'
            cmd.append('-en')
            cmd += map(str, i)
            
    # tMRCA
    if tMRCA:
        cmd.append('-L')
    
    # creates a temporary directory to run ms
    tempdir = tempfile.mkdtemp()

    # goes to the temp dir
    path= os.getcwd()
    os.chdir(tempdir)

    # runs ms and retrieves the output
    try:
        pipe= subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        result = pipe.communicate()
    except OSError, e:
        os.chdir(path)
        for root, dirs, files in os.walk(tempdir, topdown=False):
            for name in files:
                    os.remove(os.path.join(root, name))
            for name in dirs:
                    os.rmdir(os.path.join(root, name))
        os.rmdir(tempdir)
        raise RuntimeError, 'system error running ms, saying: '+ str(e)
    finally:
        if os.getcwd()==path:pass
        else:os.chdir(path)
        for root, dirs, files in os.walk(tempdir, topdown=False):
            for name in files:
                    os.remove(os.path.join(root, name))
            for name in dirs:
                    os.rmdir(os.path.join(root, name))
        os.rmdir(tempdir)
    
    # parse the results
    if len(result[1]):
        raise RuntimeError, result[1]
    if len(result[0])==0: raise RuntimeError, 'ms: no output to parse, arguments might be invalid'
    RawStrings=result[0].split('//')
    if len(RawStrings)-1!=howmany: raise RuntimeError, 'ms: number of detected alignment inconsistent with argument howmany'
    del RawStrings[0]
    RawStrings = ['//'+i for i in RawStrings]
    
    # DataMatrix conversion
    DataMat = []
    if tMRCA: tMRCAlist = []
    if theta>0 and segsites>0: problist = []
    if T==True: treelist = []
    for i in RawStrings:
        try:
            DataMat.append(egglib_binding.Ms.get(i,nsam))
        except IOError,e:
            raise RuntimeError, 'format error Ms, saying: '+ str(e)
        if tMRCA:
            x = egglib_binding.Ms.tMRCA()
            if x!=-1:
                tMRCAlist.append(x)
            else:
                raise RuntimeError, 'ms did not generate tMRCA'
        if theta>0 and segsites>0:
            x = egglib_binding.Ms.prob()
            if x!=-1:
                problist.append(x)
            else:
                raise RuntimeError, 'ms did not generate prob value'

        if T==True:
            x = egglib_binding.Ms.trees()
            if x!="":
                trees = x.split(';')
                if len(trees)<2: raise RuntimeError, 'ms did generate invalid trees'
                if trees[-1] != '': raise RuntimeError, 'ms did generate invalid trees'
                trees = [data.Tree(string=j+';') for j in trees[:-1]]
                treelist.append(trees)
            else:
                raise RuntimeError, 'ms did not generate trees'


    del RawStrings

    # creation of a Random object
    Random = egglib_binding.Random()

    # conversion into Align objects
    sAligns = [egglib_binding.Convert.align(i,0,Random) for i in DataMat]
    for i in sAligns:
        for j in range(nsam):
            i.name(j,'seq%d'%(j+1))
            # appends group information if several populations
            if I:
                counter=0
                for u,v in enumerate(I[1:]):
                    if u==I[0]:break
                    if j in range(counter,counter+(v)):
                        i.group(j,u+1)
                        break
                    else:counter+=v

    # conversion into Align objects
    Aligns = [data.Align() for i in range(len(sAligns))]
    
    for i in range(len(Aligns)):
        Aligns[i]._object = sAligns[i]
        if tMRCA:
            Aligns[i].tMRCA = tMRCAlist[i]
        if theta>0 and segsites>0:
            Aligns[i].prob = problist[i]
        if T==True:
            Aligns[i].trees = treelist[i]
    return Aligns


########################################################################

class BLASTdb(object):
    
    """
    Handles a local BLAST database saved as temporary files. This class
    is most useful when a database is needed only temporarily. The
    database will be available for BLAST applications as long as
    the instance lives.
    """
    
    _dbname = 'db'
    dirname = None
    
    ####################################################################
    
    def __init__(self, container, type):
        
        """
        The constructor expects two (mandatory) arguments: *container*
        must be a :class:`~egglib.Container` or :class:`~egglib.Align` instance and
        *type* must be either ``'nucl'`` (for nucleotides) or ``'prot'``
        (for proteins) and specify the appropriate data type. For
        protein sequences, trailing stop codons are automatically
        trimmed.
        """

        # safety checkings
        if not isinstance(container, data.BaseContainer):
            raise TypeError, 'cannot make a BLAST database; invalid sequence container'
        if not apps['makeblastdb']:
            raise RuntimeError, 'the external application "makeblastdb" is required to create a local database'
        if type not in ('nucl', 'prot'):
            raise ValueError, 'invalid type argument for BLASTdb'

        # claims a directory
        self.dirname = tempfile.mkdtemp()
        rememberdir = os.getcwd()

        # writes the database on the disk
        handle, fname = tempfile.mkstemp()
        os.close(handle)
        container2 = data.Container()
        for n,s,g in container:
            container2.append(n, s.rstrip('*'))
        container2.write(fname)

        # builds the command line
        args = (apps['makeblastdb'], '-in', fname, '-dbtype', type, '-out', self._dbname)
        
        # runs makeblastdb
        os.chdir(self.dirname)
        try:
            pipe = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = pipe.communicate()
        finally:
            os.chdir(rememberdir)
            os.remove(fname)
        
        # intercepts error
        if len(stderr):
            raise RuntimeError, 'error running `makeblastdb`: %s' %stderr

    ####################################################################

    def path(self):

        """
        Returns the full path name of the local database, as required
        by BLAST programs. It is usually not required to use this method
        directly.
        """

        return os.path.join(self.dirname, self._dbname)

    ####################################################################

    def __del__(self):

        if self.dirname!=None and os.path.isdir(self.dirname):
            for i in os.listdir(self.dirname):
                if os.path.isfile(os.path.join(self.dirname, i)):
                    os.remove(os.path.join(self.dirname, i))
            os.rmdir(self.dirname)


########################################################################

class BLAST:
    
    """
    Provides NCBI Basic Local Alignment Search Tools for finding
    homologues of query sequences against a local database. All proposed
    methods return a dictionary of processed BLAST results. This
    dictionary stores the hits for each sequences, indexed by its name
    string. If only one sequence is passed as a string, the output
    dictionary will always one item indexed by an empty string. For a
    given sequences, the results are presented as a list of HSPs (there
    can be several HSPs on a single hit sequence), each hit being a
    dictionary storing the following information: *subject* (the name
    of the hit sequence), *bitScore* (the bit score value), *score*
    (the raw score value), *eValue* (the expectation value), *qstart*
    and *qend* (the start and end positions on the query sequence),
    *hstart* and *hend* (the start and end positions on the hit sequence,
    *qframe* and *hframe* (the frame in which is locate the hit in
    respectively the query and the hit sequence), *identity* (the number
    of matching positions), *gaps* (the number of gapped positions),
    *length (the length of the hit), *qseq* (the sequence of the query
    sequence at the HSP), *hseq* (the sequence of the hit sequence at
    the HSP) and *midline* (the string summarizing the quality of the
    local alignment, indicating matching positions). The full XML
    document is nonetheless accessible as instance's member *xml_results*
    after each call to any of the method. *xml_results* is ``None`` by
    default.
    """
    
    def __init__(self):
        
        self.xml_results = None

    ####################################################################

    def __del__(self):
        
        if self.xml_results:
            self.xml_results.unlink()
            
    ####################################################################

    def _search(self, program, query, target, evalue=None, **params):

        if not apps[program]:
            raise RuntimeError, 'the external application %s is not available' %program

        if isinstance(target, BLASTdb):
            target = target.path()
        
        # process the query
        if isinstance(query, basestring):
            queryContainer = data.Container()
            queryContainer.append('', query,  0)
        elif isinstance(query, data.BaseContainer):
            queryContainer = query
        else:
            raise TypeError, 'the query argument must be a Container or a string'

        # check that no duplicated name in query container
        if queryContainer.contains_duplicates():
            raise ValueError, 'cannot perform BLAST search - ensure no names are duplicated in sequences'

        # defines parameters
        if evalue==None:
            evalue = math.exp(-6)
        parameters = {'evalue': evalue}
        
        if program=='blastn':
            if 'penalty' not in params or params['penalty']==None:
                params['penalty'] = -2
            parameters['reward'] = 1
            parameters['penalty'] = params['penalty']
            
        else:
            L = 0
            for n,s,g in queryContainer:
                L+=len(s.replace('-', ''))
            L = 1. * L / len(queryContainer)
            if L<35:
                matrix = 'PAM30'
                gapopen = 9
                gapextend = 1
            elif L<50:
                matrix = 'PAM70'
                gapopen = 10
                gapextend = 1
            elif L<85:
                matrix = 'BLOSUM80'
                gapopen = 10
                gapextend = 1
            else:
                matrix = 'BLOSUM62'
                gapopen = 10
                gapextend = 1

            parameters['matrix'] = matrix
            
            if program!='tblastx':
                parameters['gapopen'] = gapopen
                parameters['gapextend'] = gapextend
        
        # sets automatic parameters
        parameters.update(params)
            
        # sets parameters
        args = [ apps[program], '-query', '-', '-out', '-',
                 '-outfmt', '5', '-db', target  ]
        for i in parameters:
            args += [ '-%s' %i, str(parameters[i]) ]
                                      
        if program in ('blastn', 'blastp'):
            args+= ['-task', program]
            
        # launch the blast search
        
        pipe = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = pipe.communicate(str(queryContainer))

        if len(stderr):
            raise RuntimeError, 'error running blast program: %s' %stderr

        # parses output
        if not len(stdout):
            self.xml_results = None
            return dict.fromkeys(queryContainer.names(), None)

        self.xml_results = xml.dom.minidom.parseString(stdout)

        def getL(node, tag):
            return node.getElementsByTagName(tag)
        def get1(node,tag):
            try: return node.getElementsByTagName(tag)[0]
            except IndexError: raise RuntimeError, 'invalid blast output; cannot get %s' %tag
            
        results = {}
        iterations = get1(self.xml_results,'BlastOutput')
        iterations = get1(iterations,'BlastOutput_iterations')
        iterations = getL(iterations,'Iteration')
        if len(iterations)!=len(queryContainer):
            raise RuntimeError, 'invalid blast output: cannot find blast results for all queries'
        c=0
        for iteration in iterations:
            name=queryContainer.name(c)
            c+=1
            if name!='' and name!=get1(iteration,'Iteration_query-def').firstChild.data:
                raise RuntimeError, 'invalid blast output: problem with query "%s".\nCannot find exact name in blast output - check sequence names!' %name
            results[name] = []
            message = iteration.getElementsByTagName('Iteration_message')
            if len(message) and message[0].firstChild.data=='No hits found':
                hits = []
            else:
                hits = get1(iteration,'Iteration_hits')
                hits = getL(hits, 'Hit')
            for hit in hits:
                subject = get1(hit,'Hit_def').firstChild.data
                hsps = get1(hit,'Hit_hsps')
                for hsp in getL(hsps, 'Hsp'):
                    try:
                        results[name].append({
                            'subject': subject,
                            'bitScore':float(get1(hsp,'Hsp_bit-score').firstChild.data),
                            'score':   int(get1(hsp,'Hsp_score').firstChild.data),
                            'eValue':  float(get1(hsp,'Hsp_evalue').firstChild.data),
                            'qstart':  int(get1(hsp,'Hsp_query-from').firstChild.data)-1,
                            'qend':    int(get1(hsp,'Hsp_query-to').firstChild.data)-1,
                            'hstart':  int(get1(hsp,'Hsp_hit-from').firstChild.data)-1,
                            'hend':    int(get1(hsp,'Hsp_hit-to').firstChild.data)-1,
                            'qframe':  int(get1(hsp,'Hsp_query-frame').firstChild.data),
                            'hframe':  int(get1(hsp,'Hsp_hit-frame').firstChild.data),
                            'identity':int(get1(hsp,'Hsp_identity').firstChild.data),
                            'gaps':    int(get1(hsp,'Hsp_gaps').firstChild.data),
                            'length':  int(get1(hsp,'Hsp_align-len').firstChild.data),
                            'qseq':    get1(hsp,'Hsp_qseq').firstChild.data,
                            'hseq':    get1(hsp,'Hsp_hseq').firstChild.data,
                            'midline': get1(hsp,'Hsp_midline').firstChild.data
                        })
                    except ValueError:
                        raise RuntimeError, 'invalid blast output; invalid value'
        
        return results

    ####################################################################
    
    def blastn(self, query, target, evalue=None, penalty=None, **params):
        
        """
        Searches a nucleotide database using nucleotide queries. *query*
        can be a string, a :class:`~egglib.Container` or :class:`~egglib.Align` instance.
        In the latter cases, all sequences will be processed. *target*
        must refer to a valid database of the correct data type, either
        represented by its file system path or by a :class:`~egglib.wrappers.BLASTdb`
        instance. *evalue* is the expectaction value (expected number of
        random hits by chance alone, depending on the database size).
        The default value is e\ :sup:`-6` (therefore much less, and more
        stringent, than ``blastn``'s default value which is 10). *penalty*
        is the penalty to apply for nucleotide mismatch (the default
        reward for nucleotide match is +1). The default value is -2. The
        value must be negative, and should be increased to account for
        most distant homologies. "A ratio of 0.33 (1/-3) is appropriate
        for sequences that are about 99% conserved; a ratio of 0.5 (1/-2)
        is best for sequences that are 95% conserved; a ratio of about
        one (1/-1) is best for sequences that are 75% conserved" (from
        BLAST online documentation). All other BLAST parameters can be
        set as keywords arguments. Keyword arguments are passed as is to
        the ``blastn`` program and can overwrite arguments default values
        of *evalue* and *penalty*. For example it is possible to set
        *reward* as a keyword argument as in ``reward=5 penalty=-4``.
        """
        
        return self._search('blastn', query, target, evalue, penalty=penalty, **params)

    ####################################################################

    def blastp(self, query, target, evalue=None, **params):

        """
        Searches a protein database using protein queries. Arguments are
        as for :meth:`blastn` with the exception that *reward* and
        *penalty* are not applicable. Parameters *matrix*, *gapopen* and
        *gapextend* are defined automatically based on the average
        length of query sequences. These automatic settings can be
        overriden by keyword arguments.
        """

        return self._search('blastp', query, target, evalue, **params)

    def tblastn(self, query, target, evalue=None, **params):

        """
        Searches a translated nucleotide database using protein queries.
        Arguments are as for :meth:`blastp`.
        """

        return self._search('tblastn', query, target, evalue, **params)

    def blastx(self, query, target, evalue=None, **params):

        """
        Searches a protein database using translated nucleotide queries.
        Arguments are as for :meth:`blastp`.
        """

        return self._search('blastx', query, target, evalue, **params)

    def tblastx(self, query, target, evalue=None, **params):

        """
        Searches a translate nucleotide database using translated
        nucleotide queries. Arguments are as for :meth:`blastp`.
        """

        return self._search('tblastx', query, target, evalue, **params)


########################################################################

class BL2SEQ:

    """
    Provides NCBI Basic Local Alignment Search Tools for aligning two
    sequences by local alignment. The proposed methods return all a list
    of dictionaries representing all HSPs. The items of the dictionaries
    corresponding to the following variables: *qstart*, *qend*, *sstart*,
    *send*, *evalue*, *bitscore*, *score*, *length*, *nident*, *qframe*,
    *sframe*, *gaps*, *qseq*, *sseq* and *midline*. The latest is given
    only as identity/mismatch marks. Parameters are defined as for
    :class:`~egglib.wrappers.BLAST` except that *query* and *subject* must both be
    sequence strings.
    """

    def _search(self, program, query, subject, evalue=None, penalty=None, **params):

        if not apps[program]:
            raise RuntimeError, 'the external application %s is not available' %program

        # defines parameters
        if evalue==None:
            evalue = math.exp(-6)
        parameters = {'evalue': evalue}
        
        if program=='blastn':
            if penalty==None:
                penalty = -2
            parameters['reward'] = 1
            parameters['penalty'] = penalty
            
        else:
            L = len(query)
            if L<35:
                matrix = 'PAM30'
                gapopen = 9
                gapextend = 1
            elif L<50:
                matrix = 'PAM70'
                gapopen = 10
                gapextend = 1
            elif L<85:
                matrix = 'BLOSUM80'
                gapopen = 10
                gapextend = 1
            else:
                matrix = 'BLOSUM62'
                gapopen = 10
                gapextend = 1

            if program!='tblastx':
                parameters['gapopen'] = gapopen
                parameters['gapextend'] = gapextend
        
        # sets automatic parameters
        parameters.update(params)
            
        # write the subject file
        handle, fname = tempfile.mkstemp()
        os.write(handle, subject)
        os.close(handle)

        # sets parameters
        format = 'qstart qend sstart send evalue bitscore score length nident qframe sframe gaps qseq sseq'
        args = [ apps[program], '-query', '-', '-out', '-',
                 '-outfmt', '6 %s' %format, '-subject', fname    ]
 
        for i in parameters:
            args += [ '-%s' %i, str(parameters[i]) ]
            
        if program in ('blastn', 'blastp'):
            args+= ['-task', program]
            
        # launch the blast search
        try:
            pipe = subprocess.Popen(tuple(args), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = pipe.communicate(query)
        finally:
            # cleans the subject file if needed
            if os.path.isfile(fname):
                os.remove(fname)

        if len(stderr):
            raise RuntimeError, 'error running `blastn`: %s' %stderr

        results = []
        stdout = stdout.strip().split()
        format = format.split()
        while len(stdout):
            if len(stdout)<len(format):
                raise RuntimeError, 'invalid blast output'
            results.append({})
            for i in range(len(format)):
                item= stdout[0]
                del stdout[0]
                try: item = int(item)
                except ValueError:
                    try: item = float(item)
                    except ValueError:
                        pass
                results[-1][format[i]] = item
            midline = ''
            for i in range(results[-1]['length']):
                try:
                    if results[-1]['qseq'][i]==results[-1]['sseq'][i]:
                        midline+='|'
                    else:
                        midline+=' '
                except IndexError:
                    raise RuntimeError, 'invalid blast output'
            results[-1]['midline'] = midline
            results[-1]['qstart'] -= 1
            results[-1]['qend'] -= 1
            results[-1]['sstart'] -= 1
            results[-1]['send'] -= 1
        return results

    ####################################################################
    
    def blastn(self, query, subject, evalue=None, penalty=None, **params):
        
        """
        Align a nucleotide query to a nucleotide subject. Arguments are
        as for :meth:`BLAST.blastn` except that *query* and *subject* 
        must both be sequence strings.
        """
        
        return self._search('blastn', query, subject, evalue, penalty, **params)

    ####################################################################
    
    def blastp(self, query, subject, evalue=None, **params):

        """
        Align a protein query to a protein subject. Arguments are
        as for :meth:`BLAST.blastp` except that *query* and *subject*
        must both be sequence strings.
        """

        return self._search('blastp', query, subject, evalue, None, **params)

    ####################################################################

    def tblastn(self, query, subject, evalue=None, **params):

        """
        Align a translated nucleotide query to a protein subject.
        Arguments are as for :meth:`BLAST.tblastn` except that *query*
        and *subject* must both be sequence strings.
        """

        return self._search('tblastn', query, subject, evalue, None, **params)

    ####################################################################
    
    def blastx(self, query, subject, evalue=None, **params):

        """
        Align a protein query to a translated nucleotide subject.
        Arguments are as for :meth:`BLAST.blastx` except that *query*
        and *subject* must both be sequence strings.
        """

        return self._search('blastx', query, subject, evalue, None, **params)

    ####################################################################

    def tblastx(self, query, subject, evalue=None, **params):

        """
        Align a translated nucleotide query to a translated nucleotide
        subject. Arguments are as for :meth:`BLAST.tblastx` except that
        *query* and *subject* must both be sequence strings.
        """

        return self._search('tblastx', query, subject, evalue, None, **params)


########################################################################

def clustal(container, quiet=True, nogroups=False):
    
    """
    Performs multiple alignment using CLUSTALW. *container* might be
    a :class:`~egglib.Container` or :class:`~egglib.Align` instance. If
    *quiet* is ``True``, the standard output (but not standard error)
    of the wrapped program will be intercepted and discarded. Returns
    a :class:`~egglib.Align` instance. By default, the function
    preserves group labels. However, if the container contains
    duplicates (even if they belong to the same group), this operation
    will fail (with a :class:`ValueError`). To process containers
    containing duplicates and for which group label information is not
    important, set the flag *nogroups* to ``True``.
    """

    if not apps['clustalw']:
        raise RuntimeError, 'the external application "clustalw" is required to perform multiple alignment'

    # stores group information
    if not nogroups:
        if container.contains_duplicates():
            raise ValueError, 'clustal: alignment contains duplicates - you might want to use the flag `nogroups`'
        groups = container.groups()
            
    # gets tempfile names
    handle, fname = tempfile.mkstemp(suffix='.fas')
    os.close(handle)
        
    # writes the sequences
    container.write(fname)

    # runs clustal
    cmd= (apps['clustalw'], '-align', '-infile=%s' %fname)
    if quiet:
        Pipe= subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr= Pipe.communicate()
    else:
        Pipe= subprocess.Popen(cmd, stderr=subprocess.PIPE)
        stdout, stderr= Pipe.communicate()
    
    # guesses output fnames
    dndf= fname[:-4]+'.dnd'
    alnf= fname[:-4]+'.aln'
    
    # processes them accordingly
    if os.path.isfile(fname): os.remove(fname)
    if os.path.isfile(dndf): os.remove(dndf)
    
    if len(stderr):
        if os.path.isfile(alnf):
            os.remove(alnf)
        raise RuntimeError, 'unable to align using CLUSTALW, error message:\n %s' %stderr
            
    if not os.path.isfile(alnf):
        if quiet:
            message = 'unable to align using CLUSTALW - last message was:\n%s' %stdout.strip().split('\n')[-1]
        else:
            message = 'unable to align using CLUSTALW'
        raise RuntimeError, message

    try:
        aln = tools.aln2fas(alnf)
    finally:
        os.remove(alnf)
    
    # restores group information
    if not nogroups:
        for group in groups:
            for i in groups[group]:
                j = aln.find(i)
                aln.group(j, group)
     
    return aln


########################################################################

def muscle(container, quiet=True, nogroups=False):

    """
    Performs multiple alignment using MUSCLE. *container* might be
    a :class:`~egglib.Container` or :class:`~egglib.Align` instance. If
    *quiet* is ``True``, progress information will not be shown If
    *quiet* is ``False``, progress information will be shown but the
    function might not be able to detect automatically errors reported
    by the wrapped program. Returns a :class:`~egglib.Align` instance.
    By default, the function preserves group labels. However, if the
    container contains duplicates (even if they belong to the same
    group), this operation will fail (with a :class:`ValueError`). To
    process containers containing duplicates and for which group label
    information is not important, set the flag *nogroups* to ``True``.
    """

    if not apps['muscle']:
        raise RuntimeError, 'the external application "muscle" is required to perform multiple alignment'

    # stores group information
    if not nogroups:
        if container.contains_duplicates():
            raise ValueError, 'muscle: alignment contains duplicates - you might want to use the flag `nogroups`'
        groups = container.groups()

    # since muscle is (currently) unable to write to standard output, creates a temporary file
    handle, name = tempfile.mkstemp()
    os.close(handle)

    try:
        # runs muscle
        args = [apps['muscle'], '-out', name ]
        if quiet:
            args.append('-quiet')
            pipe =subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = pipe.communicate(str(container))
        else:
            pipe =subprocess.Popen(args, stdin=subprocess.PIPE)
            pipe.communicate(str(container))
            stderr = []

        # checks error
        if len(stderr) and quiet:
            raise RuntimeError, 'unable to align using muscle, error message:\n %s' %stderr

        if not os.path.isfile(name):
            raise RuntimeError, 'muscle failed to produce any alignment (reason unknown)'

        # imports data
        align = data.Align(name)

        if len(container) and not len(align):
            raise RuntimeError, 'an unexpected error occurred while running muscle'

    finally:
        os.unlink(name)

    # restores group information
    if not nogroups:
        for group in groups:
            for i in groups[group]:
                j = align.find(i)
                align.group(j, group)

    return align


########################################################################

def phyml(input, model='GTR', rates=1, boot=0, topo=None, start=None,
    search='NNI', quiet=True):

    """
    Reconstructs phylogeny using maximum likelihood through the PhyML
    software. *input* should be a :class:`~egglib.Align` instance. 
    *model* indicates the model to use. Accepted values are ``HKY85``,
    ``JC69``, ``K80``, ``F81``, ``F84``, ``TN93`` and ``GTR`` for
    nucleotides and ``LG``, ``WAG``, ``JTT``, ``MtREV``, ``Dayhoff``,
    ``DCMut``, ``RtREV``, ``CpREV``, ``VT``, ``Blosum62``, ``MtMam``,
    ``MtArt``, ``HIVw`` and ``HIVb`` for protein sequences. *rates*
    gives the number of discrete categories of evolutionary rate. *boot*
    sets the number of bootstrap repetitions. Values of -1, -2 and -3
    activates one the test-based branch support evaluation methods that
    provide faster alternatives to bootstrap repetitions. A value of 0
    will provide no branch support at all. *topo* allows to fix the tree
    topology. *start* allows to set the starting topology (it is
    illegal to set both *topo* and *start* to non-``None`` values.
    *search* can be ``NNI`` (fastest), ``SPR`` or ``BEST``(best of both
    methods). For *topo* or *start*, a :class:`~egglib.Tree` instance
    must be passed. If present, branch lengths and branch labels will
    be ignored. If *quiet* is ``True``, the standard output of the
    wrapped program will be intercepted and discarded. The function
    returns a tuple ``(tree, loglk)`` where *tree* is a
    :class:`~egglib.Tree` instance and *loglk* the log-likelihood
    reported by PhyML.
    """

    if not apps['phyml']:
        raise RuntimeError, 'the external application "PhyML" is required to perform phylogenetic reconstruction'

    # checks arguments
    dna= 'HKY85 JC69 K80 F81 F84 TN93 GTR'.split()
    aa= 'LG WAG JTT MtREV Dayhoff DCMut RtREV CpREV VT Blosum62 MtMam MtArt HIVw HIVb'.split()
    if model not in  dna+aa: raise ValueError, 'phyml3: unsupported model'
    if not isinstance(rates, (float,int)) or rates<1:
        raise ValueError, '`rates` argument of function phyml() is invalid'
    if not isinstance(boot, int) or boot<-3:
        raise ValueError, '`boot` argument of function phyml() is invalid'

    try:

        # writes the phyml-formatted data in a temporary file
        fd, seqfile = tempfile.mkstemp()
        try:
            os.write(fd, input.phyml() )
        except AttributeError:
            raise ValueError, 'first argument of function phyml() has incorrect type'
        finally:
            os.close(fd)

        # writes the tree if needed
        if topo and start: raise ValueError, 'phyml: it is illegal to use both `topo` and `start` options'
        if topo or start:
            try:
                [towrite] = filter(None, [topo, start])
            except ValueError:
                raise ValueError, 'an unexpected error occurred in phyml(), please report it to the author'
            fd, treefile = tempfile.mkstemp()
            try:
                os.write(fd, towrite.newick(False,False))
            except AttributeError:
                raise ValueError, 'argument `topo` of function phyml() has incorrect type'
            finally:
                os.close(fd)

    
        # defines the argument list0
        args=[apps['phyml']]
        args+=['-i', seqfile]
        if model in aa: args+=['-d', 'aa']
        args+=['-m', model]
        args+=['-c', str(rates)]
        if rates > 1:
            args+=['-a', 'e']
        if topo or start:
            args+=['-u', treefile]
            if topo:
                args+=['-o', 'lr']
        args+=['-b', str(boot)]

        if search not in ['NNI' , 'SPR', 'BEST']:
            raise ValueError, 'invalid `search` argument'
        args += ['-s', search]

        # runs phyml
        if quiet:
            pipe= subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout,stderr = pipe.communicate('Y')
        else:
            pipe= subprocess.Popen(args, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout,stderr = pipe.communicate('Y')

        # checks error
        stderr=stderr.strip()
        if len(stderr):
            raise RuntimeError, 'the following error occurred while running PhyML:\n%s' %stderr

        # collects results
        if not os.path.isfile(seqfile+'_phyml_tree.txt'):
            raise RuntimeError, 'error while running PhyML: cannot open output tree file'
        tree = data.Tree(seqfile+'_phyml_tree.txt')

        # intercepts no-output error
        if tree.number_of_leaves()==0:
            string = 'error while running PhyML: cannot read tree'
            if stdout!=None:
                lines = stdout.strip().split('\n')
                if len(lines)>3:
                    string += ' - last line written in output was was:\n'
                    string += lines[-3]
            raise RuntimeError, string

        # gets log-likelihood
        if not os.path.isfile(seqfile+'_phyml_stats.txt'):
            raise RuntimeError, 'error while running PhyML: cannot open statistics report file'
        f = open(seqfile+'_phyml_stats.txt')
        loglk = None
        for line in f:
            if line[:18]=='. Log-likelihood: ':
                try:
                    loglk = float(line.split()[-1])
                except ValueError:
                    raise RuntimeError, 'error while running PhyML: cannot read log-likelihood in this line: %s' %line.strip()
        f.close()

        # intercepts error
        if loglk==None:
            string = 'error while running PhyML: cannot read likelihood'
            if stdout!=None:
                lines = stdout.strip().split('\n')
                if len(lines)>3:
                    string += ' - last line written in output was was:\n'
                    string += lines[-3]
            raise RuntimeError, string

    finally:
        if os.path.isfile(seqfile):
            os.remove(seqfile)
        if topo and os.path.isfile(treefile):
            os.remove(treefile)

        # cleaning
        kill = [seqfile, seqfile+'_phyml_tree.txt', seqfile+'_phyml_stats.txt']
        if boot>0: kill += [seqfile+'_phyml_boot_trees.txt', seqfile+'_phyml_boot_stats.txt']
        if topo: kill.append(treefile)
        for i in kill:
            if os.path.isfile(i): os.unlink(i)

    return tree, loglk

########################################################################

def nj(input, groups=False, quiet=True):

    """
    Constructs a neighbor-joining tree using programs from the PHYLIP
    package (``dnadist`` and ``neighbor``). *input* should be a
    :class:`~egglib.Align` instance containing DNA sequences only.
    If *group* is ``True``, the group labels will be appended to
    sequence names and therefore will appear in the final tree. If
    *quiet* is ``True``, the standard output of the wrapped program will
    be intercepted and discarded. The function returns a
    :class:`~egglib.Tree` instance.
    """

    if not apps['dnadist']:
        raise RuntimeError, 'the external application "dnadist" is required to perform neighbor-joining'

    if not apps['neighbor']:
        raise RuntimeError, 'the external application "neighbor" is required to perform neighbor-joining'

    try:
        mapping = input.encode()
    except AttributeError:
        raise ValueError, 'first argument of function nj() has incorrect type'

    try:
        # create a temporary directory
        temppath = tempfile.mkdtemp()
        curpath = os.getcwd()
        os.chdir(temppath)

        # creates the input file
        f = open('infile', 'w')
        f.write( input.phylip('S') )
        f.close()

        # runs dnadist
        if os.path.isfile('outfile'):
            os.remove('outfile')
        sp= subprocess.Popen(apps['dnadist'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        stdout, stderr = sp.communicate('Y\n\n\n')
        
        if 'ERROR: ' in stdout:
            message = stdout[stdout.index('ERROR: '):].strip()
            raise RuntimeError, 'dnadist has reported the following error message:\n%s' %message

        if os.path.isfile('infile'):
            os.remove('infile')
        if os.path.isfile('outfile'):
            os.rename('outfile', 'infile')
        else:
            raise RuntimeError, 'dnadist has failed to generate a distance matrix'

        # runs neighbor
        if os.path.isfile('outtree'):
            os.remove('outtree')

        sp= subprocess.Popen(apps['neighbor'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        stdout, stderr = sp.communicate('Y\n\n\n')

        if 'ERROR: ' in stdout:
            message = stdout[stdout.index('ERROR: '):].strip()
            raise RuntimeError, 'dnadist has reported the following error message:\n%s' %message

        # gets the tree
        if os.path.isfile('outtree'):
            njtree= data.Tree('outtree')
        else:
            raise RuntimeError, 'neighbor has failed to a generate neighbor-joining tree'
    
    finally:

        # restoring
        input.rename(mapping)

        # cleaning
        if os.path.isfile('infile'): os.remove('infile')
        if os.path.isfile('outfile'): os.remove('outfile')
        if os.path.isfile('outtree'): os.remove('outtree')

        os.chdir(curpath)


    for leaf in njtree.get_terminal_nodes():
        label = leaf.get_label()
        if label not in mapping:
            raise RuntimeError, 'unknown temporary code %s in neighbor-joining tree' %label
        name = mapping[label]

        if groups:
            index = input.find(name)
            if index==-1:
                raise RuntimeError, 'unknown leaf name %s in neighbor-joining tree' %name
            group = input.group(index)
            leaf.set_label('%s@%d' %(name, group))
        else:
            leaf.set_label(name)

    return njtree


########################################################################

class Codeml:
    
    """
    Interface to non-synonymous/synonymous substitution rate analyses
    available in the ``codeml`` program of the PAML package. The
    sequences and tree are loaded at construction time. The results can
    be accessed through the return value of :meth:`fit` or as a
    pre-formatted string by calling ``str(codeml)`` (where *codeml* is a
    :class:`~egglib.wrappers.Codeml` instance). Default options of influencal parameters
    are: start omega value of 0.4, omega not fixed and 10 discrete omega
    categories. They can be changed using the appropriate accessors.
    After running :meth:`fit`, the instance caches the control file as
    *controlfile*, the ``codeml`` main output file as *outputfile* and
    ``codeml`` standard output (where you might be able to read error
    messages) as *standardoutput*. This is done to allow manual
    inspection in case of errors.
    """

    def __init__(self, aln, tree=None):
        
        """
        Constructor arguments: *aln*, a :class:`~egglib.Align` instance,
        *tree*, a :class:`~egglib.Tree` instance. The names of *aln* and
        *tree* must match (except that ``$x`` and ``$x`` labels -where
        *x* is an integer- are ignored at the end of tree leaf names).
        If *tree* is ``None``, a star topology will be used. Branch
        lengths from the tree are discarded.
        """
        
        self._tempdir = tempfile.mkdtemp()
        self.controlfile = ''
        self.outputfile = ''
        self.rstfile = ''
        self.standardoutput = ''
        self._trees = None
        self._align = None
        self._results={}
        self._set_align(aln)
        self._set_tree(tree)
        self._startw= 0.4
        self._fixw= False
        self._ncat= 10
        self._results= {}
        self._init_results()

        if not apps['codeml']:
            raise RuntimeError, 'the external application `codeml` is required to perform this task'

    ####################################################################

    def fix_omega(self, value):
        
        """
        Fixes omega to *value*. It is not required to call this method for
        fitting models that require a fixed value of omega.
        """
        
        self._fixw = True
        self._startw = value
        
    ####################################################################
    
    def unfix_omega(self, value):
        
        """
        Releases omega from a previous call to :meth:`fix_omega` and set
        the start value to *value*.
        """
        
        self._fixw = False
        self._startw = value
        
    ####################################################################
    
    def start_omega(self, value):
        
        """
        Sets the start value of omega to *value*. It is not legal to call
        this method when omega is fixed.
        """
        
        if self._fixw:
            raise ValueError, 'it is illegal to call `start_omega` after `fix_omega` on a `Codeml` instance'
        self._startw = value

    ####################################################################
    
    def number_of_categories(self, value):
        
        """
        Sets the number of discrete omega categories.
        """
        
        self._ncat = value
        
    ####################################################################

    def _init_results(self):
        
        self._results['model']= None
        self._results['lnL']= None
        self._results['np']= None
        self._results['kappa']= None
        self._results['omega']= None
        self._results['freq']= None
        self._results['beta']= None
        self._results['trees']= None
        self._results['site_method']= None
        self._results['site_proba']= None
        self._results['site_class']= None
        self._results['site_omega']= None
        self._results['site_error']= None        

    ####################################################################

    def _set_tree(self, tree):
        
        if tree:
            leaves = []
            for i in tree.all_leaves():
                match= re.search('(.+)[\$\#]\d+$', i)
                if match:
                    leaves.append(match.group(1))
                else:
                    leaves.append(i)

            if sorted(leaves)!=sorted(self._align.names()):
                raise ValueError, 'codeml (class): the names in the alignment and the tree don\'t match'
            self._tree = tree

        else:
            self._tree = None

    ####################################################################

    def _set_align(self, aln):

        self._align= data.Align()
        self._align.addSequences(aln)
        self._chck_seq()

    ####################################################################

    def _chck_seq(self):
        if not self._align.ls() or self._align.ls()%3:
            raise ValueError, 'to run `codeml`, sequence length must be a multiple of 3'
        for a,b,c in self._align:
            aa= tools.translate(b)
            if '*' in aa:
                raise ValueError, 'stop codon found at protein position %d of sequence %s' %(aa.find('*')+1, a)

    ####################################################################

    _models = {
        'M0':  (0, 0, 0),
        'M1a': (1, 0, 0),
        'M2a': (2, 0, 0),
        'M8a': (8, 0, 1),
        'M8':  (8, 0, 0),
        'A0':  (2, 2, 1),
        'A':   (2, 2, 0),
        'nW':  (0, 2, 0),
        'b':   (0, 1, 0)
    }

    ####################################################################

    def fit(self, MODEL, quiet=True):

        """
        Fits a given model and collects the result.
        
        *MODEL* must be only of the following model specifications:
            - ``M0``: fixed-omega model.
            - ``M1a``: nearly-neutral model.
            - ``M2a``: positive selection model.
            - ``M8a``: nearly-neutral beta model.
            - ``M8``: positive selection beta model.
            - ``A0``: branch-site null model.
            - ``A``: branch-site model.
            - ``nW``: independent omegas in subtrees.
            - ``b``: branch-independent model.

        The results are returned as a dictionary containing these keys
        (note that keys irrelevant to me fitted model will not be
        exported):
            - model: the model fitted.
            - lnL: the log-likelihood.
            - np: number of parameters.
            - kappa: the transition/transversion ratio.
            - omega: omega estimate, as a single value M0, a list of two values for M1a, three values for M2a, eleven values for M8a and M8, k values for nW (where k is the number of clades in the tree), alist of four tuples of two values for A0 and A.
            - freq: estimates of the frequency of the different categories, None for M0 and nW, a list of two values for M1a, three values for M2a, eleven values for M8a and M8 and four values for A0 and A.
            - beta: a tuple for p and q (beta distribution parameters), None for all models but M1a and M2a.
            - site_method: method used to estimate posterior site.
            - trees: the trees found in the results (in order) as a list of tree instances.
            - site_proba: list with one list per site, each list contains the posterior probability of the site under each omega category, for models M1, M1a, M2a, M8a, M8, A0, A.
            - site_class: the list of highest probability class for each site, for models M1a, M2a, M8a, M8, A0, A.
            - site_omega: the list of average posterior omega for each site, for models M1, M1a, M2a, M8a, M8.
            - site_error: the standard deviation of posterior omega for each site, for models M2a and M8.
                          
        If *quiet* is ``True``, the standard output is intercepted and
        discarded.
        """
        
        ## INITIALIZES ##
        self._results= {}
        if MODEL not in self._models:
            raise ValueError, '`Codeml doesn\'t know this model: %s' %str(MODEL)

        self._init_results()
        self._results['model']= MODEL

        if MODEL in ['A0', 'A', 'nW']:
            labels = set()
            if self._tree==None:
                raise ValueError, 'this model requires a tree!'
            for node in self._tree:
                if node.get_label()!=None:
                    try:
                        labels.add(
                            int(re.search('[\#|\$](\d+)', node.get_label()).group(1)))
                    except (ValueError, AttributeError):
                        pass
            if len(labels)==0:
                raise ValueError, 'this model requires that the tree contains labels'

        ## PREPARES CONTROLE FILE ##
        NSsites, model, fixw = self._models[MODEL]
        if fixw==0: setw= self._startw
        if fixw==1: setw= 1.

        self.controlfile="""
      seqfile = seqfile
     treefile = treefile
      outfile = outfile
        noisy = 3
      verbose = 0
      runmode = 0
      seqtype = 1
    CodonFreq = 2
        clock = 0
       aaDist = 0
        model = %s
      NSsites = %s
        icode = 0
        Mgene = 0
    fix_kappa = 0
        kappa = 2
    fix_omega = %s
        omega = %s
    fix_alpha = 1
        alpha = 0.
       Malpha = 0
        ncatG = %s
        getSE = 0
 RateAncestor = 0
   Small_Diff = .5e-6
    cleandata = 0
       method = 0
""" %(str(model), str(NSsites), str(fixw), str(setw), str(self._ncat))
        self.controlfile = self.controlfile.strip()

        ## GOES TO THE TEMP DIR ##
        path= os.getcwd()
        os.chdir(self._tempdir)
        self.outputfile = ''
        self.standardoutput = ''
        try:

            ## PERFORMS THE RUN ##
            if self._tree:
                self._tree.write('treefile', brlens=False)
            else:
                f=open('treefile', 'w')
                f.write( '(%s);\n' %(','.join([i.replace(' ', '_')
                                        for i in self._align.names()])) )
                f.close()
            f=open('seqfile', 'w')
            f.write( self._align.phyml() )
            f.close()
            f=open('controlfile', 'w')
            f.write( self.controlfile )
            f.close()

            commandline = (apps['codeml'], 'controlfile')
            
            if quiet:
                pipe= subprocess.Popen(commandline, stdout=subprocess.PIPE)
                self.standardoutput, stderr = pipe.communicate()
            else:
                pipe= subprocess.Popen(commandline)
                stdout, stderr = pipe.communicate()
                self.standardoutput = ''
            
            if quiet:
                if self.standardoutput.count('is stop'):
                    raise RuntimeError, 'codeml error: stop codon in sequence'
                if self.standardoutput.split()[0]!='CODONML':
                    raise IOError, 'codeml error: unknown output'

            f = open('outfile')
            self.outputfile = f.read()
            f.close()
            f = open('rst')
            rst = f.readlines()
            self.rstfile = ''.join(rst)
            f.close()

        except OSError, e:
            raise RuntimeError, 'system error running codeml, saying: '+ str(e)
            
        finally:
            os.chdir(path)

        ## PARSES THE MAIN OUTPUT FILE ##
        for line in self.outputfile.split('\n'):
            match= re.match('lnL\(ntime: ?.+  np: ?(.+)\): *(.+) +\+.+', line)
            if match:
                self._results['np'] = int(match.group(1))
                self._results['lnL'] = float(match.group(2))
            
            match= re.match('kappa \(ts/tv\) = +(.+)', line)
            if match: self._results['kappa'] = float(match.group(1))

            match = re.match('tree length for dN: +(.+?)$', line)
            if match: self._results['treeLength_dN'] = float(match.group(1))

            match = re.match('tree length for dS: +(.+?)$', line)
            if match: self._results['treeLength_dS'] = float(match.group(1))

            match= re.match('omega \(dN/dS\) = +(.+)', line)
            if match: self._results['omega'] = float(match.group(1))

            match= re.match('p: +(.+)', line)
            if match:
                try: self._results['freq']= [float(i) for i in match.group(1).split()]
                except ValueError: pass

            match= re.match('w: +(.+)', line)
            if match:
                # needs a specific patch
                items= match.group(1).split()
                if not self._results['freq']:
                    raise RuntimeError, 'in codeml output, expects omega values after their frequencies'
                if len(items)==(len(self._results['freq'])-1):
                    items.append(items[-1][-9:])
                    items[-2]= items[-2][:9]
                # let's hope it works
                try:
                    self._results['omega']= [float(i) for i in items]
                except ValueError:
                    pass
            
            match= re.match('  p0 ?= +.+  p ?= +(.+) q ?= +(.+)', line)
            if match:
                self._results['beta'] = float(match.group(1)), float(match.group(2))
            
            match= re.match('w \(dN/dS\) for branches: (.+)', line)
            if match:
                try:
                    self._results['omega']= [float(i) for i in match.group(1).split()]
                except ValueError:
                    pass

            match= re.match('proportion       +(.+)', line)
            if match:
                try:
                    self._results['freq']= [float(i) for i in match.group(1).split()]
                except ValueError:
                    pass

            match= re.match('background w     +(.+)', line)
            if match:
                try:
                    self._results['omega']= [float(i) for i in match.group(1).split()]
                except ValueError:
                    pass

            match= re.match('foreground w     +(.+)', line)
            if match:
                try:
                    self._results['omega']= zip(self._results['omega'], [float(i) for i in match.group(1).split()])
                except ValueError:
                    pass

            match= re.match('^(\(.+\);)$', line.strip())
            if match:
                if not self._results['trees']:
                    self._results['trees'] = []
                self._results['trees'].append( data.Tree(string=match.group(1)) )


        ## PARSES THE RST OUTPUT FILE (FOR MODELS M1a, M2a, M8a, M8, A0 and A) ##
        self._results['site_method'] = None
        self._results['site_proba'] = None
        self._results['site_class'] = None
        if MODEL in ('M1a', 'M2a', 'M8a', 'M8', 'A0', 'A'):
            if MODEL=='M1a':
                seekline= 'Naive Empirical Bayes (NEB) probabilities'
            if MODEL in ('M2a', 'M8a', 'M8', 'A0', 'A'):
                seekline= 'Bayes Empirical Bayes (BEB) probabilities'

            pos= None
            for i, line in enumerate(rst):
                if MODEL in ('M1a', 'M8a', 'A0') and re.match('^Naive Empirical Bayes \(NEB\) probabilities', line):
                    pos= i
                    self._results['site_method'] = 'NEB'
                if MODEL in ('M2a', 'M8', 'A') and re.match('^Bayes Empirical Bayes \(BEB\) probabilities', line):
                    pos= i
                    self._results['site_method'] = 'BEB'
            
            if not pos:
                raise RuntimeError, 'cannot find empirical Bayes probabilities in the rst output file of codeml'
            
            self._results['site_proba'] = []
            self._results['site_class'] = []
            self._results['site_proba'] = []
            if MODEL in ('M1a', 'M2a', 'M8a', 'M8'):
                self._results['site_omega'] = []
                if MODEL in ('M2a', 'M8'):
                    self._results['site_error'] = []

            for line in rst[pos+3:]:
                if MODEL in ['A0', 'A']:
                    match= re.match('^\d+ [-A-Z*] +(.+) (.+) (.+) (.+) \( ?(\d+)\)$', line.strip())
                    if not match: break
                    self._results['site_proba'].append([float(match.group(i+1)) for i in range(4)])
                    self._results['site_class'].append(int(match.group(5)))
                if MODEL in ['M1a', 'M8a']:
                    match= re.match('^\d+ [-A-Z*] +(.+) \( ?(\d+)\) +(.+)$', line.strip())
                    if not match: break
                    self._results['site_proba'].append([float(i) for i in match.group(1).split()])
                    self._results['site_class'].append(int(match.group(2)))
                    self._results['site_omega'].append(float(match.group(3)))
                if MODEL in ['M2a', 'M8']:
                    match= re.match('^\d+ [-A-Z*] +(.+) \( ?(\d+)\) +(.+) \+- +(.+)$', line.strip())
                    if not match: break
                    self._results['site_proba'].append([float(i) for i in match.group(1).split()])
                    self._results['site_class'].append(int(match.group(2)))
                    self._results['site_omega'].append(float(match.group(3)))
                    self._results['site_error'].append(float(match.group(4)))

        if not self._results['lnL']:
            raise RuntimeError, 'could not import codeml\'s output'

        # checkings
        if MODEL=='M1a':
            if not self._results['freq'] or len(self._results['freq'])!=2:
                raise RuntimeError, 'codeml: inconsistent results'
            if not self._results['omega'] or len(self._results['omega'])!=2:
                raise RuntimeError, 'codeml: inconsistent results'
        if MODEL=='M2a':
            if not self._results['freq'] or len(self._results['freq'])!=3:
                raise RuntimeError, 'codeml: inconsistent results'
            if not self._results['omega'] or len(self._results['omega'])!=3:
                raise RuntimeError, 'codeml: inconsistent results'
        if MODEL=='A0' or MODEL=='A':
            if not self._results['freq'] or len(self._results['freq'])!=4:
                raise RuntimeError, 'codeml: inconsistent results'
            if not self._results['omega'] or len(self._results['omega'])!=4:
                raise RuntimeError, 'codeml: inconsistent results'
            for i in range(4):
                if len(self._results['omega'][i])!=2:
                    raise RuntimeError, 'codeml: inconsistent results'
        if MODEL=='M8a' or MODEL=='M8':
            if not self._results['freq'] or len(self._results['freq'])!=self._ncat+1:
                raise RuntimeError, 'codeml: inconsistent results (invalid number of freq)'
            if not self._results['omega'] or len(self._results['omega'])!=self._ncat+1:
                raise RuntimeError, 'codeml: inconsistent results (invalid number of omega)'
        
        if (MODEL=='M8a' or MODEL=='M8') and not self._results['beta']:
                raise RuntimeError, 'codeml: inconsistent results (no beta parameters)'

        if self._results['site_method']:
            if not self._results['site_proba']:
                raise RuntimeError, 'codeml: inconsistent results (no site_proba)'
            NSITES= len(self._results['site_proba'])
            if not self._results['site_class'] or len(self._results['site_class'])!=NSITES:
                raise RuntimeError, 'codeml: inconsistent results (no site_class or incorrect number of items)'
            if MODEL in ['M1a', 'M2a', 'M8a', 'M8'] and (not self._results['site_omega'] or len(self._results['site_omega']))!=NSITES:
                raise RuntimeError, 'codeml: inconsistent results (error with number of sites)'
            if MODEL in ['M2a', 'M8'] and (not self._results['site_error'] or len(self._results['site_omega']))!=NSITES:
                raise RuntimeError, 'codeml: inconsistent results (error with site_error)'

        return dict(self._results)

    ####################################################################

    def __str__(self):

        string= ''
        if self._results['model']: string+= 'model:     %s\n' %(self._results['model'])
        if self._results['lnL']:   string+= 'lnL:      %f\n' %(self._results['lnL'])
        if self._results['np']:    string+= 'np:       %f\n' %(self._results['np'])
        if self._results['kappa']: string+= 'kappa:     %f\n' %(self._results['kappa'])
        if self._results['omega']:
            if isinstance(self._results['omega'], float): omega= str(self._results['omega'])
            if isinstance(self._results['omega'], list):
                omega = ''
                for i in self._results['omega']:
                    if isinstance(i, float): omega+=' %f' %i
                    else: omega+= ' |%f %f|' %i
                    omega= omega.strip()
            string+= 'omega:     %s\n' %(omega)
        if self._results['freq']:  string+= 'frequency: %s\n' %(' '.join([str(i) for i in self._results['freq']]))
        if self._results['beta']:  string+= 'beta par.: %f %f\n' %self._results['beta']
        if self._results['trees']:
            for i in self._results['trees']: string+= 'tree: %s\n' %str(i)
        if self._results['site_method']:
            string+= 'Site probabilities (method %s):\n' %self._results['site_method']
            for i in range(len(self._results['site_proba'])):
                string+= '      %s %s [%d]' %(str(i+1).rjust(4, '0'), ' '.join([str(j).ljust(7, '0') for j in self._results['site_proba'][i]]), self._results['site_class'][i])
                if self._results['site_omega']: string+= ' %f' %self._results['site_omega'][i]
                if self._results['site_error']: string+= ' +- %f' %self._results['site_error'][i]
                string+= '\n'

        return string

    ####################################################################

    def __del__(self):

        if os.path.isdir(self._tempdir):
            if tempfile.gettempdir()!=self._tempdir[:len(tempfile.gettempdir())]:
                sys.stderr.write('WARNING: refusing to clean directory %s. Reason: doesn\'t look like a temporary directory' %self.tempdir)
            else:
                for fname in os.listdir(self._tempdir):
                    os.remove(os.path.join(self._tempdir, fname))
                os.rmdir(self._tempdir)

########################################################################     

class Primer3:

    """
    Primer design using the program ``PRIMER3``. The constructor takes a
    sequence and optional parameters. The list of parameters and default
    values can be accessed through the class-level attribute dictionary
    :attr:`.default_parameters`.
    """
    
    _tempname = 'egglib.primer3.temp'

    default_parameters = {
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_MAX_SIZE': 27,
        'PRIMER_MIN_TM': 57.0,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_SALT_DIVALENT': 0.0,
        'PRIMER_DNTP_CONC': 0.0,
        'PRIMER_MIN_GC': 30.0,
        'PRIMER_OPT_GC_PERCENT': 50.0,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MAX_NS_ACCEPTED': 2,
        'PRIMER_PICK_ANYWAY': 1,
        'PRIMER_GC_CLAMP': 1,
        'PRIMER_PAIR_MAX_DIFF_TM': 3,
        'PRIMER_LIBERAL_BASE': 1,
        'PRIMER_NUM_RETURN': 5,
        'PRIMER_FIRST_BASE_INDEX': 1
    }
    
    """
    Class-level dictionary holding default values for all run parameters.
    """
    
    ####################################################################

    def __init__(self, sequence, **params):
        
        """
        *sequence* must be a nucleotide sequence. Parameter values can
        be passed as keyword arguments. Parameter default values are
        taken from *Primer.default_parameters*. Parameters are restricted
        to the default list, such as spelling errors might result in a
        crash later, during primer search.
        """
        
        if not apps['primer3']:
            raise RuntimeError, 'the external application `primer3` is required to build an instance of the class Primer3'
        
        self._sequence = sequence.upper()
        self._Fprimers = []
        self._Rprimer = []
        self._pairs = []
        
        self._params = dict(Primer3.default_parameters)
        self._params.update(params)

    ####################################################################

    def find_primers(self):
        
        """
        Finds primers. Returns a tuple ``(nf, nr)`` where *nf* is the
        number of forward primers found and *nr* the number of reverse
        primers found.
        """
        
        self._Fprimers = []
        self._Rprimers = []
        
        # will work in a temporary directory
        origin = os.getcwd()
        dest = tempfile.mkdtemp()
        os.chdir(dest)
        
        try:
            string = ''
            # creates the parameter strin
            string += 'SEQUENCE_ID=%s\n' %self._tempname
            string += 'SEQUENCE_TEMPLATE='+self._sequence + '\n'
            string += 'PRIMER_PICK_LEFT_PRIMER=1\n'
            string += 'PRIMER_PICK_LEFT_PRIMER=2\n'
            string += 'P3_FILE_FLAG=1\n'
            for key,value in self._params.items():
                string += '%s=%s\n' %(key,str(value))
            string +='=\n'
            
            # runs the program
            pipe= subprocess.Popen([apps['primer3'], '-strict_tags'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = pipe.communicate(string)
            
            if len(stderr):
                raise RuntimeError, 'error running primer3: %s' %stderr

            mess = re.search('ERROR=(.+)', stdout)
            if mess:
                raise RuntimeError, 'primer3 error: message: %s'%mess.group(1)

            # import the data
            if (os.path.isfile(self._tempname+'.for') and os.path.isfile(self._tempname+'.rev')):

                pfile= open(self._tempname+'.for')
                pfile.readline()
                pfile.readline()
                pfile.readline()
                for j in pfile:
                    bits = j.strip().split()
                    if len(bits)!=10:
                        raise IOError, 'invalid primer3 output - this line:\n%s' %j
                    index, seq, loc, length, namb, gc, tm, any, end, q = bits
                    
                    p = int(loc) - self._params['PRIMER_FIRST_BASE_INDEX']
                    if tools.compare(seq, self._sequence[p:p+len(seq)]) == False:
                        raise IOError, 'cannot locate primer %s generated by PRIMER3' %seq
                    
                    self._Fprimers.append({'seq': seq, 'pos': p,
                                'GC%': float(gc), 'Tm': float(tm),
                                'Q': float(q), 'ANY': float(any),
                                'END': float(end)})
                pfile.close()

                pfile= open(self._tempname+'.rev')
                pfile.readline()
                pfile.readline()
                pfile.readline()
                for j in pfile:
                    bits = j.strip().split()
                    if len(bits)!=10:
                        raise IOError, 'invalid primer3 output - this line:\n%s' %j
                    index, seq, loc, length, namb, gc, tm, any, end, q = bits
                    rcseq = tools.rc(seq)
                    
                    p = int(loc) - (len(seq)-1) - self._params['PRIMER_FIRST_BASE_INDEX']
                    if tools.compare(rcseq, self._sequence[p:p+len(seq)]) == False:
                        raise IOError, 'cannot locate primer %s generated by PRIMER3' %seq
                    
                    self._Rprimers.append({'seq': seq, 'pos': p,
                                'GC%': float(gc), 'Tm': float(tm),
                                'Q': float(q), 'ANY': float(any),
                                'END': float(end)})
                pfile.close()

            else:
                raise RuntimeError, 'primer3 failed to generate ouput files for an unknown reason'

        finally:
            if os.path.isfile(self._tempname+'.for'): os.remove(self._tempname+'.for')
            if os.path.isfile(self._tempname+'.rev'): os.remove(self._tempname+'.rev')
            os.chdir(origin)
            os.rmdir(dest)

        self._Fprimers.sort(lambda x,y: cmp(x['pos'],y['pos']))
        self._Rprimers.sort(lambda x,y: cmp(x['pos'],y['pos']))

        return [ len(self._Fprimers), len(self._Rprimers) ]

    ###################################################################

    def forward_primers(self):
        
        """
        Returns a reference to the list of forward primers (that must
        have been previously detected using :meth:`find_primers`). Each item of
        the list represents a primer as a dictionary containing the
        following keys: *seq* (the primer sequence, given in the 5' to
        3' orientation), *pos* (the position of nucleotide at the 5' end),
        *GC%*, *Tm*, *Q* (the quality value), *END* (not defined in
        ``PRIMER3`` documentation, might be the maximum secondary structure
        stability) and *ANY* (also not documented in ``PRIMER3``, might be
        maximum misannealing stability with respect to the provided
        sequence). The last two parameters should be minimized, and their
        definition will be confirmed as soon as possible.
        """
        
        return self._Fprimers
        
    ###################################################################

    def reverse_primers(self):
        
        """
        Equivalent to :meth:`forward_primers`, except that the *pos*
        value is the position of the nucleotide at the 3' end of the
        primer, therefore the first nucleotide when reading in the
        original orientation of the provided sequence.
        """

        return self._Rprimers

    ###################################################################
    
    def find_pairs(self, mini=70, maxi=9999):
        
        """
        Finds primer pairs. Primers must have been previously designed.
        *mini* and *maxi* gives the range of accepted products. This
        method doesn't involve any call to ``PRIMER3``.
        """

        self._pairs = []

        if not len(self._Fprimers) or not len(self._Rprimers):
           return 0

        for i in self._Fprimers:
            for j in self._Rprimers:
                stop = j['pos']+len(j['seq'])-1
                size = stop-i['pos']+1
                if (size<mini or size>maxi): continue
                elem={}
                elem['F'] = i
                elem['R'] = j
                elem['start'] = i['pos']
                elem['end'] = stop
                elem['size'] = size
                self._pairs.append(elem)

        return len(self._pairs)

    ####################################################################
    
    def pairs(self):
        
        """
        Returns the list of primer pairs found by :meth:`find_pairs`.
        Each item is a directory with the values: *F* (the forward
        primer), *R* (the reverse primer), *start*, *end* and *size*.
        The primers are the same as given by :meth:`forward_primers` and
        :meth:`reverse_primers`.
        """

        return self._pairs

    ####################################################################

    def check_pairs(self):
        
        """
        Checks primer pairs defined using :meth:`find_pairs` and discards
        the pairs that fail to pass the test. This method includes a
        second call to the ``PRIMER3`` application.
        """

        # will work in a temporary directory
        origin = os.getcwd()
        dest = tempfile.mkdtemp()
        os.chdir(dest)
        
        try:
        
            # creates the parameter string
            string = ''
            for key,value in self._params.items():
                string += '%s=%s\n' %(key,str(value))
            string += 'PRIMER_EXPLAIN_FLAG=1\n'
            for n,i in enumerate(self._pairs):
                string += 'SEQUENCE_ID=test_pair_%d\n' %n
                string += 'SEQUENCE_TEMPLATE='+self._sequence+'\n'
                string += 'SEQUENCE_PRIMER='+i['F']['seq']+'\n'
                string += 'SEQUENCE_PRIMER_REVCOMP='+i['R']['seq']+'\n'
                string += 'PRIMER_MAX_NS_ACCEPTED='+str(self._params['PRIMER_MAX_NS_ACCEPTED'])+'\n'
                string += '=\n'
            string +='=\n'

            # runs the program
            pipe= subprocess.Popen([apps['primer3'], '-strict_tags'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = pipe.communicate(string)
            
            # detect errors
            if len(stderr):
                raise RuntimeError, 'error running primer3: %s' %stderr

            # import the data
            temppairs= []
            c=0
            for i in stdout.split('\n'):
                a= re.search('PRIMER_PAIR_EXPLAIN=considered 1, ok (0|1)',i)
                if a:
                    if a.group(1)=='1':
                        temppairs.append(self._pairs[c])
                    c+=1
            self._pairs= temppairs
            
        finally:
            os.chdir(origin)
            os.rmdir(dest)
            
        return len(self._pairs)

    ####################################################################

    def select(self, number):
        
        """
        Sorts the list of primer pairs (based on the sum of primer
        qualities) and select best primers. *number* gives the number
        of primer pairs to retain. If there is less pairs, they will
        all be retained, but still be sorted. Returns the number of
        pairs retained. The lists of forward and reverse primers are not
        affected.
        """
        
        self._pairs.sort(lambda x,y: cmp((x['F']['Q']+x['R']['Q']), (y['F']['Q']+y['R']['Q'])))
        self._pairs = self._pairs[:number]
        return len(self._pairs)

    ####################################################################

    def _lastCleanQ(self, sequence, number):

        """
        Returns ``True`` if the last *number* positions of the sequence
        *sequence* is made of A, C, G and T (case-independent) only.
        """

        return set(sequence[-number:].upper()) <= set('ACGT')

    ####################################################################

    def clean_pair_ends(self, number):

        """
        Deletes all primer pairs that contain at least one invalid
        character (fully resolved, not missing). All primer pairs with
        a least one primer containing a character other than A, C, G and
        T (case-independent) close to the 3' end are removed. *number*
        gives the number of characters to consider. If the number of
        larger than the length of the primer, the complete sequence is
        considered. Returns the number of pairs.
        """

        tp = []
        for i in self._pairs:
            if (self._lastCleanQ(i['F']['seq'],number) and
                self._lastCleanQ(i['R']['seq'],number)):
                    tp.append(i)
        self._pairs = tp
        return len(self._pairs)

    ####################################################################

    def clean_primer_ends(self, number):
        
        """
        Similar to :meth:`clean_pair_ends` except that the lists of
        forward and reverse primers are considered. The pairs, if they
        were generated, are not affected. Returns a tuple ``(nf, nr)``
        with *nf* and *nr* the numbers of forward and reverse primers,
        respectively.
        """
        
        tp = []
        for i in self._Fprimers:
            if (self._lastCleanQ(i['seq'],number)):
                tp.append(i)
        self._Fprimers = tp
        tp = []
        for i in self._Rprimers:
            if (self._lastCleanQ(i['seq'],number)): tp.append(i)
        self._Rprimers = tp
        return [len(self._Fprimers),len(self._Rprimers)]

    ####################################################################

    def sort(self):
        
        """
        Sorts all class attributes (forward and reverse primers and
        primer pairs) based on quality (sum of both primer qualities for
        pairs). 
        """
        
        self._Fprimers.sort(cmp=lambda x,y: cmp(x['Q'], y['Q']))
        self._Rprimers.sort(cmp=lambda x,y: cmp(x['Q'], y['Q']))
        self._pairs.sort(cmp=lambda x,y: cmp(x['R']['Q']+x['F']['Q'], y['R']['Q']+y['F']['Q']))

########################################################################
