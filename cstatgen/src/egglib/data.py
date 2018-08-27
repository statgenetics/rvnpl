__license__ = """
    Copyright 2008-2011 Stephane De Mita, Mathieu Siol

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


import os, sys, random, xml.dom.minidom, re, cStringIO as StringIO
import egglib_binding
import tools



#######################################################################

class SequenceItem(object):
    
    """
    Item managing the ``name``, ``sequence`` and ``group`` values of a
    given index of a :class:`~egglib.Container` or :class:`~egglib.Align` instance. Any
    change applied to the :class:`~egglib.SequenceItem` instance are immediately
    propagated to the :class:`~egglib.Container` or :class:`~egglib.Align` instance
    (generating the corresponding exception in case of misuse). It is
    important to note that some of errors might be generated when
    attempting to access data and not upon object creation. The
    ``print item`` statement (where ``item`` is a :class:`~egglib.SequenceItem`
    instance) returns a specially formatted string
    ``"name",sequence,group`` where ``name``, ``sequence`` and ``group``
    are name string, sequence string and group index for the index
    corresponding to the instance (note the double quotes around the
    name and the commas separating the three items). The instance also
    supports iteration and index-based accessing, but note that
    :class:`~egglib.SequenceItem` instance contain always three items: the 
    name, the sequence and the group (in that order).
    :class:`~egglib.SequenceItem` also supports indexing (item 0 is the name,
    item 1 the sequence and item 2 the group index).
    
    .. versionadded:: 2.0.1
    """
    
    ####################################################################

    def __init__(self, parent, index):
        
        """
        *parent* must be a :class:`~egglib.Container` or :class:`~egglib.Align` instance
        and *index* an index lying in the range of ``len(parent)``.
        """
        
        self._parent = parent
        self._index = index
        
        if not isinstance(index, int):
            raise TypeError, '`index` must be an integer'
        if index <0 or index>len(parent):
            raise IndexError, 'invalid index: %d (range 0-%d)' %(index, len(parent))

    ####################################################################

    def __iter__(self):
        
        """
        Iterates over the elements of the ``(name, sequence, group)``
        tuple
        """
            
        for i in (self.name, self.sequence, self.group):
            yield i

    ####################################################################

    def __getitem__(self, index):

        if index==0: return self.name
        if index==1: return self.sequence
        if index==2: return self.group
        raise IndexError, 'SequenceItem instances contain only three item'

    ####################################################################
    
    @property
    def name(self):
        """ Name string """
        return self._parent.name(self._index)

    @name.setter
    def name(self, string):
        self._parent.name(self._index, string)

    ####################################################################

    @property
    def sequence(self):
        """ Sequence string """
        return self._parent.sequence(self._index)

    @sequence.setter
    def sequence(self, string):
        self._parent.sequence(self._index, string)

    ####################################################################

    @property
    def group(self):
        """ Group label """
        return self._parent.group(self._index)

    @group.setter
    def group(self, integer):
        self._parent.group(self._index, integer)

    ####################################################################

    def __str__(self):
        
        """
        Returns a string representation of the values corresponding to
        this item.
        """
        
        return '"%s",%s,%d' %(self.name, self.sequence, self.group)


#######################################################################

class BaseContainer(object):

    """ Provides Python-orientated methods arounds classes Container and
    Align, as wrapped using SWIG. Internally, the build instance holds a
    Container or an Align instance as attribute *_object*.
    
    .. versionchanged:: 2.0.1
        The ``[]`` operators accept only indices.
        :meth:`.sequenceByName()` fulfils the dictionary-like
        behaviour. :meth:`append()`, :meth:`extend()` and
        :meth:`__iadd__()` (operator ``+=``) are removed. """
    
    ####################################################################

    @classmethod
    def create(cls, obj):
        
        """
        Creates an instance from the object `obj`. The created instance
        will match the type from which the method is called (
        ``Container.create(obj)`` will return a :class:`~egglib.Container`, and
        ``Align.create(obj)`` will return a :class:`~egglib.Align`, and the
        same goes if the method is called on an object). In the case of
        :class:`~egglib.Align`, the restriction of coherent sequence lengths
        applies (there is not automatic correction). *obj* is a priori a
        :class:`~egglib.Container` or a :class:`~egglib.Align`, but the method supports
        any iterable returning ``(name,sequence,group)`` or
        ``(name,sequence)`` tuples (in the latter case, groups will be
        initialized to ``0``. For example, the following is valid::
        
            import egglib
            data = []
            data.append( ("sequence1", "AAAAAAAAA",     0) )
            data.append( ("sequence2", "GGGG",          2) )
            data.append( ("sequence3", "AAAAAAAAAAAAAA") )
            container = egglib.Container.create(data)
            
        .. versionadded:: 2.0.1
        """

        result = cls()
        
        if isinstance(obj, egglib_binding.Container):
            for i in xrange(obj.ns()):
                result.append(obj.name(i),obj.sequence(i),obj.group(i))

        else: # assume a (name,sequence[,group) iterable            
            result.addSequences(obj)

        return result

    ####################################################################

    def __init__(self):

        """ This class cannot be instantiated. """

        raise RuntimeError, 'base class of Container and Align cannot be instantiated'
        
    ####################################################################

    def str(self, exportGroupLabels=False, lineLength=50):

        """ Formats the instance as a fasta string. *exportGroupLabels*:
        if ``True``, exports group/population membership as ``@x`` tags
        placed at the end of sequence names (where ``x`` is any positive
        integer). *lineLength* gives the number of characters to place
        on a single line in the fasta output. If ``0``, no newlines are
        inserted within sequences."""

        return self.__str__(exportGroupLabels, lineLength)

    ####################################################################

    def __str__(self, exportGroupLabels=False, lineLength=50):

        """ Identical to str() """

        return egglib_binding.Fasta.format(self._object, exportGroupLabels, lineLength)
    
    ####################################################################

    def write(self, fname, exportGroupLabels=False, lineLength=50):
        
        """ Writes the sequences to a fasta-formatted file. *fname* is
        the name of the file to create. Other arguments are as for
        :meth:`.str`.
        
        .. versionadded:: 2.0.1 """
        
        egglib_binding.Fasta.formatf(fname, self._object, exportGroupLabels, lineLength)

    ####################################################################

    def __len__(self):

        """ Returns the number of sequences of the instance. Can be
        used as in:
        
        >>> align = egglib.Align("file.fas")
        ... len(align)
        12 """

        return self.ns()

    ####################################################################

    def __getitem__(self, index):
        
        """ Operator []: returns the tuple (name, sequence, group)
        corresponding to the passed *index* """
        
        if (index<0): index= len(self)+index
        if (index<0 or index>=self.ns()): raise IndexError, 'invalid index'
        return (self.name(index), self.sequence(index), self.group(index))

    ####################################################################

    def sequenceByName(self, name, strict=True):
        
        """ Returns the sequence string corresponding to the first
        match of `name`. If the name is not found, raises a ``KeyError``.
        If *strict* is ``True``, seeks an exact match. If ``False``,
        compares only until the end of the requested name (for example:
        ``'ATCFF'`` will match ``'ATCFF_01'`` if *strict* is ``False``).

        .. versionadded:: 2.0.1 """

        p= self.find(str(name), strict)
        if (p==None): raise KeyError, 'no such sequence: '+str(name)
        return self.sequence(p)

    ####################################################################

    def groupByName(self, name, strict=True):
        
        """ Returns the group label corresponding to the first match of
        `name`. If the name is not found, raises a :class:`KeyError`. If *strict*
        is ``True``, seeks an exact match. If ``False``, compares only
        until the end of the requested name (for example: ``'ATCFF'``
        will match ``'ATCFF_01'`` if *strict* is ``false``).

        .. versionadded:: 2.0.1 """

        p= self.find(str(name), strict)
        if (p==None): raise KeyError, 'no such sequence: '+str(name)
        return self.group(p)

    ####################################################################

    def __setitem__(self, index, value):
        
        """ Operator []: sets an element of the instance. The assignment
        is by index. If *value* is a :class:`str`, it must be the
        sequence. If the object is of type :class:`~egglib.Align`, the sequence
        must match the alignment length. It must be a sequence with two
        or three items: *name*, *sequence* and optionally *group*.
        sequence, [group]). The index must fit in the range of the
        instance (only overwriting). If omitted, the group will be
        set to 0. """
            
        if (index<0): index= self.ns()+index
        if (index<0 or index>=self.ns()): raise IndexError, 'invalid index'
        if (isinstance(value, basestring)):
            self.sequence(index, value)
        elif (len(value)==2 and
                isinstance(value[0], basestring) and
                isinstance(value[1], basestring)):
                self.name(index, value[0])
                self.sequence(index, value[1])
                self.group(index,0)
        elif (len(value)==3 and
                isinstance(value[0], basestring) and
                isinstance(value[1], basestring) and
                isinstance(value[2], (int, long))): 
                self.name(index, value[0])
                self.sequence(index, value[1])
                self.group(index,value[2])
        else: raise TypeError, 'invalid value'

    ####################################################################

    def __delitem__(self, index):
        
        """ Operator [] - deletes an element of the instance. """
            
        if (index<0): index= self.ns()+index
        if (index<0 or index>=self.ns()): raise IndexError
        return self._object.remove(index)

    ####################################################################

    def __contains__(self, key):
        
        """ Tests whether a sequence name occurs in the instance.
        Corresponds to expressions such as "name in fasta" and
        returns a boolean. """
            
        if (self.find(str(key))==None): return False
        else: return True

    ####################################################################

    def __iter__(self):
        
        """
        Iterates over :class:`~egglib.SequenceItem` instances corresponding to
        all indices of the instance
        """
            
        for index in xrange(len(self)): yield SequenceItem(self,index)

    ####################################################################

    def addSequences(self, seqs):
        
        """ Appends repetitively (name, sequence, group) tuples to the
        end of the object (passed the last sequence. *seqs* must be
        an iterable returning (name, sequence, group) tuples (such
        as a :class:`~egglib.Container` or :class:`~egglib.Align` instance). (the group
        item is optional and tuples can be of length 2.) Returns the
        number of sequences after the operation.

        .. versionadded:: 2.0.1 """
            
        for item in seqs: self.append(*item)
        return self.ns()

    ####################################################################

    def composition(self):
        
        """ Gets the composition in characters of each sequence. Returns
        a dictionary with the sequence names as key. Each entry is
        itself a dictionary giving the absolute frequency of each
        character found in the corresponding sequences. """
            
        ret= {}
        for i in self:
            ret[i[0]]= {}
            for j in i[1]:
                j= j.upper()
                if ret[i[0]].has_key(j): ret[i[0]][j]+=1
                else: ret[i[0]][j]= 1
        return ret

    ####################################################################

    def names(self):

        """ Returns the list of sequence names """

        ret = [i[0] for i in self]
        return ret

    ####################################################################

    def duplicates(self):
        
        """ Returns the list of sequence names found more than once in
        the instance. """
        
        at_least_one= []
        more_than_one= []
        for i in self:
            if i[0] in at_least_one: more_than_one.append(i[0])
            at_least_one.append(i[0])
        return more_than_one

    ####################################################################

    def contains_duplicates(self):
        
        """ ``True`` if the instance contains at least one duplicate. """
        
        at_least_one= []
        for i in self:
            if i[0] in at_least_one: return True
            at_least_one.append(i[0])
        return False

    ####################################################################

    def no_duplicates(self):

        """ Discards all duplicates: for all sequences with the same
            name, the one with the largest index is removed. """

        i=0
        while(i<len(self)):
            j=i+1
            while(j<len(self)):
                if (self.name(i)==self.name(j)): del self[j]
                else: j+=1
            i+=1

    ####################################################################

    def groups(self):

        """ Gets the group structure. Returns a dictionary with the
        group labels (as ``int``) as keys. Values are the lists of
        sequence names corresponding to each group. """

        lst= {}
        for i in self:
            if i[2] in lst: lst[i[2]].append(i[0])
            else: lst[i[2]] = [i[0]]
        return lst

    ####################################################################

    def shuffle(self, maintain_outgroup= True):
        
        """ Randomly reassigns group labels. Modifies the current object
        and returns nothing. If *maintain_outgroup* is ``True``, doesn't
        reassign the outgroup (group label ``999``). """

        # take labels
        labels= [i[2] for i in self if (not maintain_outgroup or i[2]!=999)]

        # shuffle labels
        random.shuffle(labels)

        # reaffect labels
        for i in range(len(self)):
            if maintain_outgroup and self.group(i)==999: continue
            self.group(i, labels[0])
            if not len(labels): raise Exception, 'sorry, bug in shuffle method'
            labels= labels[1:]
        if len(labels): raise Exception, 'sorry, bug in shuffle method'

    ####################################################################

    def append(self, name, sequence, group=0):

        """ Adds a sequence to the object. *name* is the sequence name,
        *sequence* the sequence string and *group* is the population
        label. Note that the length of *sequence* must match the length
        of the alignment, if *self* is of type :class:`~egglib.Align`. Returns
        the number of sequences after the operation. """

        return self._object.append(str(name), str(sequence), int(group))
    
    ####################################################################
    
    def clear(self):

        """ Deletes all content of the current instance. """

        self._object.clear()

    ####################################################################
        
    def remove(self, name):
        
        """
        Removes the first sequence having name *name*. If no sequence
        has this name, a :class:`KeyError` is raised. A workaround is
        easy to implement::

            >>> index = align.find(name)
            >>> if index!=None:
            >>>     del align[index]
        
        .. versionchanged:: 2.0.1
            New meaning.
        """
        
        index = self.find(str(name))
        if index==None:
            raise KeyError, 'cannot delete Align sequence with name %s' %str(name)
        self._object.remove(index)

    ####################################################################

    def name(self, pos, name=None):
        
        """ Sets/gets the name of the sequence at index *pos*. If *name*
        is ``None``, returns the current name. Otherwise changes the
        name and returns nothing. """
        
        if pos<0: raise IndexError, 'invalid sequence index'
        if name==None: return self._object.name(pos)
        self._object.name(pos, name)

    ####################################################################

    def sequence(self, pos, sequence=None):
        
        """ Sets/gets the sequence string at index *pos*. If *sequence*
        is ``None``, returns the current sequence. Otherwise changes the
        sequence and returns nothing. If the object is an :class:`~egglib.Align`,
        the sequence length must match the alignment length. """
        
        if pos<0: raise IndexError, 'invalid sequence index'
        if sequence==None: return self._object.sequence(pos)
        self._object.sequence(pos, sequence)

    ####################################################################

    def group(self, pos, group=None):
        
        """ Sets/gets the group label of the sequence at index *pos*.
        If *group* is ``None``, returns the current group label.
        Otherwise changes the group label and returns nothing. If not
        ``None``, *group* must be a positive integer. """
        
        if pos<0: raise IndexError, 'invalid sequence index'
        if group==None: return self._object.group(pos)
        self._object.group(pos, group)

    ####################################################################
    
    def appendSequence(self, pos, sequence):
        
        """ Appends the sequence string *sequence* to the sequence at
        position *pos*. """

        self._object.appendSequence(pos, sequence)
        
    ####################################################################
    
    def set(self, sequence, position, ch):
        
        """ Sets the character value at string position *position* of
        the sequence at index *sequence* to value *ch*. """
        
        self._object.set(sequence, position, ch)
    
    ####################################################################

    def get(self, s, p):
        
        """ Gets the character value of the sequence *s* at position
        *p*. """
        
        return self._object.get(s, p)

    ####################################################################

    def ns(self):
        
        """ Returns the number of sequences contained in the instance. """

        return self._object.ns()

    ####################################################################

    def find(self, string, strict=True):
        
        """
        Returns the index of the first sequence with the name specified
        by *string*. If *strict* is ``False``, then the comparison
        ignores names that are longest than *string*. In other words,
        the name ``Alphacaga_tada1`` will be recognized if :meth:`.find`
        is called with *string* ``Alphacaga`` and *strict* = ``False``.
        If the name is not found, returns ``None``.
        
        .. versionchanged:: 2.1.0
            Returns None instead of -1 if the name is not found.
        """
        
        POS = self._object.find(str(string), strict)
        if POS==-1: return None
        else: return POS

    ####################################################################

    def matches(self, format):
        
        """ Returns the :class:`list` of indices matching the passed format.
        The format is passed as-is the the :mod:`re` module using the
        function :meth:`.search` (which doesn't necessarily match the
        beginning of the string). If no sequence name matches the passed
        format an empty list is returned. """
        
        return [i for i in range(len(self)) if re.search(format, self.name(i))]

    ###################################################################
    
    def encode(self, nbits=10):
        
        """
        Renames all sequences using a random mapping of names to unique
        keys of lenght *nbits*. *nbits* cannot be lower than 4 (which
        should allow renaming several millions of sequence) or larger
        than 63 (which is the number of different characters available).
        Returns a dictionary mapping all the generated keys to the
        actual sequence names. The keys are case-dependent and
        guaranteed not to start with a number. The returned mapping can
        be used to restored the original names using :meth:`.rename`

        .. versionadded:: 2.0.1
        """

        code = 'ABCDEDGHIJKLMNOPQRSTUVWXYZ'
        code += code.lower()
        code += '0123456789_'
        
        if nbits<4 or nbits>len(code):
            raise ValueError, 'invalid value of `nbits` argument of method `Align.encode`'

        mapping = {}
        for i in range(len(self)):
            name = self.name(i)
            while True:
                key = ''.join(random.sample(code, nbits))
                if key[0] not in '0123456789_' and key not in mapping:
                    break
            mapping[key] = name
            self.name(i, key)
        return mapping

    ###################################################################

    def rename(self, mapping, liberal=False):
        
        """
        Rename all sequences of the instance using the passed *mapping*.
        If *liberal* is `False` and a name does not appear in *mapping*,
        a :class:`ValueError` is raised. If *liberal* is ``True``, names
        that don't appear in *mapping* are left unchanged.
        
        .. versionadded:: 2.0.1
        """
        
        for i in range(len(self)):
            name = self.name(i)
            if name in mapping:
                self.name(i, mapping[name])
            else:
                if not liberal:
                    raise ValueError, 'cannot rename sequence name: %s' %self.name(i)

########################################################################

class Container(BaseContainer) :

    """ Holds sequences, without requiring that they have the same
    length. This class is a C++-implemented class providing performant
    storage and access utilies, wrapped within at Python layer that
    interfaces several operations. In particular it allows direct
    instanciation from a fasta-formatted file or from a string stored
    in a Python :class:`str` instance (see constructor's signature
    below).
    
    :class:`~egglib.Container` also allow subscript indexing (as in
    ``container[0]``) and iteration (as in ``for i in container``).
    Returned items are :class:`~egglib.SequenceItem` instances that can be
    either converted in ``(name, sequence, group)`` tuples or modified
    to modify the underlying instance. For example, the following code
    resets all group indices of the :class:`~egglib.Container` instance
    ``container``:
    
    >>> for i in container:
    ...     i.group = 0
    
    :class:`~egglib.Container` supports call ``str()`` and expressions such as
    ``print container``. In both cases, the result of the :meth:`.str`
    method (with default arguments) is returned. The result is a
    fasta-formatted string. Consider using the :meth:`.str` method to
    customize and :meth:`.write` to export the instance to a file on the
    disk.
    
    :class:`~egglib.Container` instances have a ``len()`` (the result of
    :meth:`.ns` is returned) supports expressions such as ``name in
    container`` which return ``True`` if ``name`` is the name of one of
    the sequences contained in the instance.
    
    .. versionchanged:: 2.0.1
       The ``[]`` operators accept only indices.
       :meth:`.sequenceByName()` fulfils the dictionary-like behaviour.
       :meth:`.append()`, :meth:`.extend()` and  :meth:`.__iadd__()`
       (operator ``+=``) are removed. """

    ####################################################################
    
    def __init__(self, fname=None, string=None, groups=False):

        """ .. rubric:: Constructor arguments

        :param fname: the path of a fasta-formatted file or ``None``.
        :param string: a string containing fasta sequences or ``None``.
        :param groups: whether to import group labels. The labels should appear as strings ``@0``, ``@1``, etc. in the input file.

        If *fname* and *string* are ``None``, the value of *groups*
        is ignored and an empty instance is built. If both *fname*
        and *string* are specified, an error is thrown.
        
        .. versionchanged:: 2.0.1
            Doesn't accept simultaneous values for *fname*  and
            *string*.
        
        .. rubric:: Methods """

        # checks
        if fname!=None and string!=None:
            raise ValueError, 'the constructor of Container accepts only one value for fname and string arguments'

        self._object = egglib_binding.Container()
        
        # imports from file
        if fname:
            if not os.path.isfile(fname): raise ValueError, 'cannot open %s' %fname
            egglib_binding.Fasta.parsef(fname, self._object, groups)
        
        # imports from string
        elif string:
            egglib_binding.Fasta.parse(string, self._object, groups)
            

    ####################################################################

    def slice(self, a, b):

        """ Extracts a selection of sequences. Sequences with indices
        *a* to *b*-1 are extracted and returned as a
        new instance. If *a* is smaller than 0,
        0 is used instead. If *b* is larger than the number of
        sequences, the latter is used instead. If *b* is not larger
        than *a*, the returned instance is empty. """

        r= Container()
        i = max(0, a)
        while (i<b and i<len(self)):
            r.append(self.name(i), self.sequence(i),self.group(i))
            i+=1
        return r

    ####################################################################

    def appendSequence(self, pos, sequence):

        """ Appends the *sequence* string to the end of the sequence at
        position *pos* of the instance. """

        self._object.appendSequence(pos, sequence)

    ####################################################################

    def ls(self, pos):
        
        """ Returns the length of the sequence stringat position *pos*.
        
        .. versionadded:: 2.0.1 """

        return self._object.ls(pos)

    ####################################################################

    def equalize(self, ch='?'):
        
        """ Appends character *ch* to the end of sequences such as all
        sequences have the same length. The length of all sequences will
        be the length of the longest sequence before call. This value is
        returned by the method. """

        if len(ch)!=1: raise ValueError, '%s: invalid equalization character' %ch
        return self._object.equalize(ch)
        
    ####################################################################
    
    def isEqual(self):
        
        """ Returns ``True`` if all sequences have the same length,
        ``False`` otherwise. """

        if not len(self): return True
        for i in range(1, len(self)):
            if self.ls(i)!=self.ls(0): return False
        return True

########################################################################

class Align(BaseContainer):

    """ Holds sequences and ensures that they have the same length.
    This class is a C++-implemented class providing performant storage
    and access utilies, wrapped within at Python layer that interfaces
    several operations. In particular it allows direct instanciation
    from a fasta-formatted file or from a string stored in a Python
    :class:`str` instance (see constructor's signature below).
    
    :class:`~egglib.Align` also allow subscript indexing (as in
    ``align[0]``) and iteration (as in ``for i in align``).
    Returned items are :class:`~egglib.SequenceItem` instances that can be
    either converted in ``(name, sequence, group)`` tuples or modified
    to modify the underlying instance. For example, the following code
    resets all group indices of the :class:`~egglib.Align` instance
    ``align``:
    
    >>> for i in align:
    ...     i.group = 0
    
    :class:`~egglib.Align` supports calls to both ``str()`` (and, as a result,
    expressions such as ``print align``). In both cases, the result of
    the :meth:`.str` method (with default arguments) is returned. The
    result is a fasta-formatted string. Consider using the :meth:`.str`
    method to customize and :meth:`.write` to export the instance to a
    file on the disk.
    
    :class:`~egglib.Align` instances have a ``len()`` (the result of
    :meth:`.ns` is returned) supports expressions such as ``name in
    align`` which return ``True`` if ``name`` is the name of one of
    the sequences contained in the instance.
    
    .. versionchanged:: 2.0.1
       The ``[]`` operators accept only indices.
       :meth:`.sequenceByName()` fulfils the dictionary-like behaviour.
       :meth:`.append()`, :meth:`.extend()` and  :meth:`.__iadd__()`
       (operator ``+=``) are removed.
    """
    
    ####################################################################

    def __init__(self, fname=None, string=None, groups=False):
        
        """ .. rubric:: Constructor arguments

        :param fname: the path of a fasta-formatted file or ``None``.
        :param string: a string containing fasta sequences or ``None``.
        :param groups: whether to import group labels. The labels should appear as strings ``@0``, ``@1``, etc. in the input file.

        If *fname* and *string* are ``None``, the value of *groups*
        is ignored and an empty instance is built. If both *fname*
        and *string* are specified, an error is thrown.
        
        .. versionchanged:: 2.0.1
            Doesn't accept simultaneous values for *fname*  and
            *string*.
        
        .. rubric:: Methods """

        # check
        if fname!=None and string!=None:
            raise ValueError, 'the constructor of Align accepts only one value for fname and string arguments'

        self._object = egglib_binding.Align()

        # imports from a fasta file
        if fname!=None:
            if not os.path.isfile(fname):
                raise ValueError, 'cannot open %s' %fname
            egglib_binding.Fasta.parsef(fname, self._object, groups)

        # imports from a string
        elif string!=None:
            egglib_binding.Fasta.parse(string, self._object, groups)

    ####################################################################

    def slice(self, a, b):

        """ Extracts a selection of sequences. Sequences with indices
        ``a`` to ``b-1`` are extracted and returned as an
        :class:`~egglib.Align` instance. If ``a`` is smaller than ``0``,
        ``0`` is used instead. If ``b`` is larger than the number of
        sequences, the latter is used instead. If ``b`` is not larger
        than ``a``, the returned instance is emptry. """

        r= Align()
        for i in range(max(0,a),min(b,len(self))):
            r.append(self.name(i), self.sequence(i), self.group(i))
        return r

    ####################################################################

    def nexus(self, prot=False):
        
        """ Generates a simple nexus-formatted string. If *prot* is
        ``True``, adds ``datatype=protein`` in the file, allowing it to
        be imported as proteins (but doesn't perform further checking).
        Returns a nexus-formatted string. Note: any spaces and tabs in
        sequence names are replaced by underscores. This nexus
        implementation is minimal but will normally suffice to export
        sequences to programs expecting nexus. """
            
        string=  '#NEXUS\n'
        string+= 'begin data;\n'
        string+= 'dimensions ntax=%d nchar=%d;\n' %(self.ns(), self.ls())
        if prot: type='prot'
        else: type='dna'
        string+= 'format datatype=%s interleave=no gap=-;\n'%type
        string+= 'matrix\n'
        L=0
        for i in self:
            if (len(i[0])>L): L=len(i[0])
        for i in self: string+= '%s  %s\n' %(i[0].ljust(L).replace(' ', '_').replace('\t', '_'), i[1])
        string+= ';\nend;\n'
        return string

    ####################################################################

    def filter(self, ratio, valid='ACGT'):
        
        """ Removes the sequences with too few valid sites. *ratio* is
        the limit threshold (relative to the sequence with the largest
        number of valid characters). The user can specify the list of
        valid states through the argument *valid*. The comparison is
        case-independent. This method modifies the current instance and
        returns nothing. """

        validu=valid.upper()
        def rlen(seq):
            return len([i for i in seq.upper() if i in validu])

        max= 0
        kill= []
        for i in self:
            if rlen(i[1])>max: max= rlen(i[1])
        for i in range(len(self)):
            if 1.*rlen(self.sequence(i))/max <ratio: kill.append(i)
        kill.reverse()

        for i in kill: del self[i]

    ####################################################################

    def column(self, pos):
        
        """ Extracts the alignment column at position *pos* as a list of
        characters. """
        
        if (pos<0 or pos>=self.ls()): raise Exception, 'invalid index'
        return [i[1][pos] for i in self]

    ####################################################################

    def consensus(self):
        
        """ Generates a consensus of the object, assuming nucleotide
        sequences. The consensus is generated based on standard
        ambiguity (IUPAC) codes. A ``-`` character is inserted if any
        sequence has a ``-``. A ``?`` character is inserted if any
        sequence has a ``?``. Returns the consensus string. """
        
        # consensus loop
        sequence = ''
        for i in range(self._object.ls()):
            data=[]
            for j in self: data.append(j[1][i].lower())
            sequence+=self.consensus_base(data)
        return sequence.upper()

    ####################################################################

    def consensus_base(self, data):

        # Returns the consensus of a sequence of characters
        #
        # data: a sequence of nucleotide (IUPAC) characters.
        # Returns: The most resolved character accounting for all
        # characters passed.

        # (this might be recoded more elegantly with sets) (this is very slow!)

        if len(set(data) - set('acgt-nmrwsykbdhv?')):
            raise ValueError, 'consensus can only generated with IUPAC-compliant sequences'

        if set(data)==set('-'): return '-'
        if '-' in data: return '?'
        if '?' in data: return '?'

        if 'n' in data:
            data.remove('n')
            data+= ['a','c', 'g', 't']
        if 'm' in data:
            data.remove('m')
            data+= ['a','c']
        if 'r' in data:
            data.remove('r')
            data+= ['a','g']
        if 'w' in data:
            data.remove('w')
            data+= ['a','t']
        if 's' in data:
            data.remove('s')
            data+= ['c','g']
        if 'y' in data:
            data.remove('y')
            data+= ['c','t']
        if 'k' in data:
            data.remove('k')
            data+= ['g','t']
        if 'b' in data:
            data.remove('b')
            data+= ['c','g','t']
        if 'd' in data:
            data.remove('d')
            data+= ['a','g','t']
        if 'h' in data:
            data.remove('h')
            data+= ['a','c','t']
        if 'v' in data:
            data.remove('v')
            data+= ['a','c','g']

        if 'a' in data:
            if 'c' in data:
                if 'g' in data:
                    if 't' in data: return 'n'
                    else: return 'v'
                elif 't' in data: return 'h'
                else: return 'm'
            elif 'g' in data:
                if 't' in data: return 'd'
                else: return 'r'
            elif 't' in data: return 'w'
            else: return 'a'
        elif 'c' in data:
            if 'g' in data:
                if 't' in data: return 'b'
                else: return 's'
            elif 't' in data: return 'y'
            else: return 'c'
        elif 'g' in data:
            if 't' in data: return 'k'
            else: return 'g'
        elif 't' in data: return 't'
        else: return '?'

    ####################################################################

    def ls(self):
        
        """ Returns the length of the alignment. """
        
        return self._object.ls()
        
    ####################################################################
        
    def removePosition(self, pos):
        
        """ Removes character at position *pos* of all sequences
        (effectively removing a column of the alignment. Returns the
        new length of the alignment. """
        
        return self._object.removePosition(pos)

    ####################################################################

    def character(self, s, p):
        
        """ Fast accessor to a character. Returns character at position
        *p* of sequence *s*. This accessor is faster than :meth:`.get`
        because it doesn't perform out-of-bound check. """
        
        return self._object.character(s, p)

    ####################################################################
        
    def binSwitch(self, pos):
        
        """ Takes all characters at position *pos*; replaces ``0`` by
        ``1`` and all the way around, raises an exception if another
        character is found. This method doesn't have a return value. """
        
        self._object.binSwitch(pos)
        
    ####################################################################
        
    def extract(self, *args):
        
        """ Extract given positions (or columns) of the alignment and
        returns a new alignment. There are two ways of using this
        method. The first is by passing a range specification as in
        ``align.extract(100, 200)``. The bounds will be passed as it to
        the slice operator on all sequences. The above example will
        extract columns ``100`` to ``199``. As a result, out of bound
        values will be silently supported. The second use of the method
        is as in ``align.extract([80, 143, 189, 842, 967])``. The single
        argument must be an iterable containing positions indices, that
        might contain repetitions and needs not to be sorted. The
        positions will be extracted in the specified order.
        
        .. versionadded:: 2.0.1 """
        
        extracted = Align()

        if len(args)==1:
            positions = map(int, args[0])
            for n,s,g in self:
                extracted.append(n, ''.join([s[i] for i in positions]), g)

        elif len(args)==2:
            for n,s,g in self:
                extracted.append(n, s[int(args[0]):int(args[1])], g)

        else:
            raise ValueError, 'invalid argument for method extracting alignment positions'
        
        return extracted

    ####################################################################

    def polymorphism(self, allowMultipleMutations=False,
                     minimumExploitableData=1., ignoreFrequency=0,
                     validCharacters="ACGT", missingData="MRWSYKBDHVN?-",
                     useZeroAsAncestral=False,
                     skipDifferentiationStats=False,
                     skipOutgroupBasedStats=False,
                     skipAllHaplotypeStats=False,
                     skipHaplotypeDifferentiationStats=False):
                         
        """
        Computes nucleotide and haplotype diversity statistics.

        Arguments:

          * *minimumExploitableData* sites where the non-missing data
            (as defined by mapping strings, see below) are at a
            frequency smaller than this value will be removed from the
            analysis. Use 1. to take only 'complete' sites into account
            and 0. to use all sites (the option is not considered for
            haplotype-based statistics).
            
          * *allowMultipleMutations*: if ``False``, only sites with 1
            or 2 alleles are considered, and sites with more alleles
            are considered as missing data. The sum of the frequencies
            of all alleles not matching the outgroup will treated as the
            derived allele frequency (for orientable sites).
            
          * *ignoreFrequency*: removes sites that are polymorph because
            of an allele at absolute frequency (as an integer: number of copies) smaller than or equal to
            this value. If ignoreFrequency=0, no sites are removed,
            if ignoreFrequency=1, singleton sites are ignored. Such
            sites are completely removed from the analysis (not counted
            in lseff). Note that if more than one mutation is allowed,
            the site is removed only if all the alleles but one are
            smaller than or equal to this value. For example, an
            alignment column AAAAAAGAAT is ignored with an
            ignoreFrequency of 1, but AAAAAAGGAT is conserved
            (including the third allele T which is a singleton).

          * *validCharacters*: a string giving the list of characters
            that should be considered as valid data.
            
          * *missingData*: characters indicating missing data, that is
            tolerated but ignored. All characters that are neither in
            *validCharacters* nor *missingData*, but found in the data,
            will cause an error.
            
          * *useZeroAsAncestral*: if true, all outgroups (if present)
            will be ignored and the character "0" will be considered as
            ancestral for all sites, whatever the character mapping. 
            
          * *skipDifferentiationStats*, *skipOutgroupBasedStats*,
            *skipAllHaplotypeStats*, *skipHaplotypeDifferentiationStats*
            allow the user to skip part of the analysis (in order to
            save time).
            
        The method returns a dictionary containing the diversity
        statistics. Some of the statistics will be computed only in
        presence of more than one group in the alignment, or in the
        presence of an outgroup, or depending on the value of other
        statistics and or if skip flags were activated (otherwise, they
        will have a ``None`` value).
        
        These statistics are always computed:
        
            * ``nseff``: Average number of analyzed sequences per
              analyzed site. It equals to ``ns()`` minus the number
              of outgroup sequences unless *minimumExploitableData* is
              set to a value < 1. In the latter case it can be a fraction.
            * ``lseff``: Number of analyzed sites.
            * ``npop``: Number of populations detected in the alignment.
            * ``S``: Number of polymorphic sites.
            * ``eta``: Minimal number of mutations (ignores
              *allowMultipleMutations*).
            * ``sites``: List of :class:`SitePolymorphism` instances
              (one for each polymorphic site).
            * ``siteIndices``: List of alignment position of each
              polymorphic site.
            * ``singletons``: List of positions of singletons.

        These statistics are computed only is ``lseff`` is > 0:

            * ``thetaW``: Theta estimator of Watterson (*Theor. Popul. Biol.* **7**:256-276, 1975).
            * ``Pi``: Nucleotide diversity.

        This statistic is computed only is ``S`` is > 0:
        
            * ``D``: Tajima statistic (*Genetics* **123**:585-595, 1989)
            
        These statistics are computed only if *skipAllHaplotypeStats* is
        ``False``:
        
            * ``He``: Haplotypic diversity.
            * ``K``: Number of distinct haplotypes.
            * ``alleles``: A list giving the haplotype index for each
              sequence of the alignment, or -1 when not applicable, such
              as for the outgroup. 

        This statistic is computed only if *skipAllHaplotypeStats*
        and *skipHaplotypeDifferentiationStats* are ``False``, ``npop``
        is > 1:

            * ``Snn``: Nearest neighbor statistics (Hudson *Genetics*
              **155**:2011-2014, 2000).

        These statistics are computed only if *skipAllHaplotypeStats*
        and *skipHaplotypeDifferentiationStats* are ``False``, ``npop``
        is > 1 and ``S`` is > 0:

            * ``Fst``: Population differentiation, based on nucleotides
              (Hudson *et al.* *Genetics* **132**:583-589, 1992).
            * ``Gst``: Population differentiation, based on haplotypes
              (Nei version, Hudson *et al.* *Mol. Biol. Evol.*
              **9**:138-151, 1992).
            * ``Hst``: Population differentiation, based on haplotypes
              (Hudson *et al.* *Mol. Biol. Evol.* **9**:138-151, 1992).
            * ``Kst``: Population differentiation, based on nucleotides
              (Hudson *et al.* *Mol. Biol. Evol.* **9**:138-151, 1992).

        These statistics are computed only if *skipDifferentiationStats*
        is ``False`` and ``npop`` is > 1:
        
            * ``pair_CommonAlleles``: For each pair of populations,
              number of sites with at least one allele shared by the
              two populations. Alleles that are fixed in one or both
              populations are taken into account, provided that they are
              polymorphic over the whole sample.
            * ``pair_FixedDifferences``: For each pair of populations,
              number of sites with a fixed differences between the two
              populations.
            * ``pair_SharedAlleles``: For each pair of populations,
              number of sites with at least one allee shared by the two
              populations. Only alleles that are segregating in both
              populations are taken into account.
            * ``pop_Polymorphisms``: For each population, number of
              polymorphic sites in this population.
            * ``pop_SpecificAlleles``: For each population, number of
              sites with at least one allele specific to this
              population.
            * ``pop_SpecificDerivedAlleles``: For each population,
              number sites with at least one derived allele specific to
              this population.
            * ``CommonAlleles``: Number of sites with at least one
              allele shared among at least two populations.
            * ``FixedDifferences``: Number of sites with at least one
              difference fixed between two populations.
            * ``SharedAlleles``: Number of sites with at least one
              allele shared by at least two populations.
            * ``SpecificAlleles``: Number of sites with at least one
              allele specific to one population.
            * ``SpecificDerivedAlleles``: Number of sites with at least
              one derived allele specific to one population.

        These statistics are computed only if *skipDifferentiationStats*
        is ``False``, ``npop`` is > 1 and ``lseff`` > 0:

            * ``average_Pi``: Average of nucleotide diversity per
              population.
            * ``pop_Pi``: Vector of nucleotide diversity per population.

        This statistic is computed only if *skipDifferentiationStats* is
        ``False`` and ``npop`` is 3:

            * ``triConfigurations``: A list of 13 numbers counting the
              number of sites falling in the possible configurations in
              three populations. Only diallelic loci are considered and
              rooting is not considered. The possible configurations are
              explained using an example below. The order is as in the
              returned list (remember that indices start from 0 in
              Python). Each line gives the allele(s) present in each
              population: A means one allele, G the other allele and A/G
              both alleles (polymorphism within this population). A and
              G can be substituted (A/G, A, A is the same as A/G, G, G).
              
                    * 0: A/G, A, A
                    * 1: A/G, A, G
                    * 2: A, A/G, A
                    * 3: A, A/G, G
                    * 4: A, A, A/G
                    * 5: A, G, A/G
                    * 6: A/G, A/G, A
                    * 7: A/G, A, A/G
                    * 8: A, A/G, A/G
                    * 9: A/G, A/G, A/G
                    * 10: A G G
                    * 11: A G A
                    * 12: A A G
               
        These statistic are computed only if *skipOutgroupBasedStats* is
        ``False``:
        
            * ``lseffo``: Number of oriented sites that were analyzed.
            * ``So``: Number of polymorphic sites among oriented sites.
            
        These statistics are computed only if *skipOutgroupBasedStats*
        is ``False`` and *So* is > 0:
        
            * ``thetaH``: Theta estimator of Fay and Wu (*Genetics*
              **155**:1405-1413, 2000).
            * ``thetaL``: Theta estimator of Zeng *et al.* (*Genetics*
              **174**:1431-1439, 2006).
            * ``H``: Fay and Wu statistic (*Genetics* **155**:1405-1413, 2000)
            * ``Z``: Fay and Wu statistic standardized by Zeng *et al.*
              *Genetics* **174**:1431-1439, 2006).
            * ``E``: Zeng *et al.* statistic *Genetics* **174**:1431-1439, 2006).
              
        The returned dictionary also contains a nest dictionary ``options``
        which feedbacks the values used at function call.
              
        .. versionchanged:: 2.0.2
           ``Polymorphisms`` is renamed ``pop_Polymorphisms``. The
           following statistics are added: ``pair_CommonAlleles``,
           ``pair_FixedDifferences``, ``pair_SharedAlleles``,
           ``pop_SpecificAlleles``, ``pop_SpecificDerivedAlleles``. The
           following statistic are now computed only if ``So`` > 0:
           ``thetaH``, ``thetaL``, ``E``, ``H`` and ``Z``. The following
           statistics are now  computed only if ``lseff`` > 0: ``thetaW``,
           ``Pi``, ``pop_Pi`` and ``average_Pi``. The following statistic
           are computed only if ``S`` > 0: ``D``, ``Fst``, ``Gst``,
           ``Hst`` and ``Kst``. ``npop`` is always returned. For
           consistency, outgroup-based  statistics are computed even if
           ``lseffo`` is 0 (except those who require that ``So`` > 0).

        .. versionchanged:: 2.1.0
            The statistics not computed are now exported and set to
            ``None``.
        """
        
        if not isinstance(ignoreFrequency, int):
            raise TypeError, 'argument ignoreFrequency must be an integer'
                
        results = {}
        
        results['options'] = {
            'allowMultipleMutations': allowMultipleMutations,
            'minimumExploitableData': minimumExploitableData,
            'ignoreFrequency': ignoreFrequency,
            'validCharacters': validCharacters,
            'missingData': missingData,
            'useZeroAsAncestral': useZeroAsAncestral
        }

        for key in ['nseff', 'lseff', 'S', 'eta', 'singletons', 'sites',
                'siteIndices', 'thetaW', 'Pi', 'D', 'npop',
                'average_Pi', 'pair_CommonAlleles',
                'pair_FixedDifferences', 'pair_SharedAlleles', 'pop_Pi',
                'pop_Polymorphisms', 'pop_SpecificAlleles',
                'pop_SpecificDerivedAlleles', 'CommonAlleles',
                'FixedDifferences', 'SharedAlleles', 'SpecificAlleles',
                'SpecificDerivedAlleles', 'triConfigurations', 'lseffo',
                'So', 'thetaH', 'thetaL', 'H', 'E', 'Z', 'He', 'K',
                'alleles', 'Fst', 'Gst', 'Hst', 'Kst', 'Snn']:
            results[key] = None

        nucleotideDiversity = egglib_binding.NucleotideDiversity()
        nucleotideDiversity.load(self._object, allowMultipleMutations,
                                minimumExploitableData, ignoreFrequency,
                                '%s %s' %(validCharacters, missingData),
                                useZeroAsAncestral)

        results['nseff'] = nucleotideDiversity.nseff()
        results['lseff'] = nucleotideDiversity.lseff()
        results['S'] = nucleotideDiversity.S()
        results['eta'] = nucleotideDiversity.eta()
        results['singletons'] = nucleotideDiversity.singleton_positions()
        results['sites'] = [egglib_binding.SitePolymorphism(nucleotideDiversity.get_site(i))
                                for i in range(nucleotideDiversity.S())]
        results['siteIndices'] = nucleotideDiversity.polymorphic_positions()
        if results['lseff']>0:
            results['thetaW'] = nucleotideDiversity.thetaW()
            results['Pi'] = nucleotideDiversity.Pi()
        if results['S']>0:
            results['D'] = nucleotideDiversity.D()

        results['npop'] = nucleotideDiversity.npop()

        if not skipDifferentiationStats:
            if results['npop']>1:

                if results['lseff']>0:
                    results['average_Pi'] = nucleotideDiversity.average_Pi()

                CommonAlleles = []
                FixedDifferences = []
                SharedAlleles = []

                for i in range(results['npop']):

                    for j in range(i+1, results['npop']):
                        
                        CommonAlleles.append(nucleotideDiversity.CommonAlleles(i,j))
                        FixedDifferences.append(nucleotideDiversity.FixedDifferences(i,j))
                        SharedAlleles.append(nucleotideDiversity.SharedAlleles(i,j))

                results['pair_CommonAlleles'] = CommonAlleles
                results['pair_FixedDifferences'] = FixedDifferences
                results['pair_SharedAlleles'] = SharedAlleles
                if results['lseff']>0:
                    results['pop_Pi'] = [nucleotideDiversity.pop_Pi(i) for i in range(results['npop'])]
                results['pop_Polymorphisms'] = [nucleotideDiversity.Polymorphisms(i) for i in range(results['npop'])]
                results['pop_SpecificAlleles'] = [nucleotideDiversity.SpecificAlleles(i) for i in range(results['npop'])]
                results['pop_SpecificDerivedAlleles'] = [nucleotideDiversity.SpecificDerivedAlleles(i) for i in range(results['npop'])]
                results['CommonAlleles'] = nucleotideDiversity.CommonAlleles()
                results['FixedDifferences'] = nucleotideDiversity.FixedDifferences()
                results['SharedAlleles'] = nucleotideDiversity.SharedAlleles()
                results['SpecificAlleles'] = nucleotideDiversity.SpecificAlleles()
                results['SpecificDerivedAlleles'] = nucleotideDiversity.SpecificDerivedAlleles()
            if results['npop']==3:
                results['triConfigurations'] = [nucleotideDiversity.triConfiguration(i) for i in range(13)]

        if not skipOutgroupBasedStats:
            results['lseffo'] = nucleotideDiversity.lseffo()
            results['So'] = nucleotideDiversity.So()
            if results['So']>0:
                results['thetaH'] = nucleotideDiversity.thetaH()
                results['thetaL'] = nucleotideDiversity.thetaL()
                results['H'] = nucleotideDiversity.H()
                results['E'] = nucleotideDiversity.E()
                results['Z'] = nucleotideDiversity.Z()
        
        del nucleotideDiversity
        
        if not skipAllHaplotypeStats:
            haplotypeDiversity = egglib_binding.HaplotypeDiversity()
            haplotypeDiversity.reserve(results['S'])
            haplotypeDiversity.load(self._object, allowMultipleMutations,
                    ignoreFrequency,'%s %s' %(validCharacters, missingData))
        
            results['He'] = haplotypeDiversity.He()
            results['K'] = haplotypeDiversity.K()
            
            results['alleles'] = []
            c = 0
            for i in range(len(self)):

                if self.group(i)==999:
                    allele = -1
                else:
                    allele = haplotypeDiversity.haplotypeIndex(c)
                    c+=1

                results['alleles'].append( allele )
            
            if not skipHaplotypeDifferentiationStats:
                if results['npop']>1:
                    if results['S']>0:
                        results['Fst'] = haplotypeDiversity.Fst()
                        results['Gst'] = haplotypeDiversity.Gst()
                        results['Hst'] = haplotypeDiversity.Hst()
                        results['Kst'] = haplotypeDiversity.Kst()
                    results['Snn'] = haplotypeDiversity.Snn()
        
        return results

    ####################################################################

    def polymorphismBPP(self, dataType=1):
        
        """
        Computes diversity statistics using tools provided through the
        Bio++ libraries. Note that attempting to call this method from
        an EggLib module compile without Bio++ support will result in a
        RuntimeError.
        
        Arguments:
        
            * *dataType*: 1 for DNA, 2 for RNA, 3 for protein sequences,
              4 for standard codons, 5 for vertebrate mitochondrial
              codons, 6 for invertebrate mitochondrial codons and 7 for
              echinoderm mitochondrial codons.

        The method returns a dictionary containing the diversity
        statistics. Some keys will be computed only in the presence of
        an outgroup, or if sequences were specified as coding or
        depending on the value of other statistics (otherwise, they will
        be ``None``).
        
        The following statistics are always computed:
        
            * ``S``: Number of polymorphic sites.
            * ``Sinf``: Number of parsimony informative sites.
            * ``Ssin``: Number of singleton sites.
            * ``eta``: Minimal number of mutations.
            * ``thetaW``: Theta estimator (Watterson *Theor. Popul.
              Biol.* **7**:256-276, 1975).
            * ``T83``: Theta estimator (Tajima *Genetics* **105**:437-460, 1983)
            * ``He``: Heterozygosity.
            * ``Ti``: Number of transitions.
            * ``Tv``: Number of transversions.
            * ``K``: Number of haplotypes.
            * ``H``: Haplotypic diversity.
            * ``rhoH``: Hudson's estimator of rho (*Genet. Res.*
              **50**:245-250, 1987).

        The following statistic is computed only if ``Tv`` > 0:

            * ``TiTv``: Transition/transversion ratio.
            
        The following statistic is computed only if ``S`` > 0:

            * ``D``: Tajima statistic (*Genetics* **123**:585-595, 1989).

        The following statistics are computed only if ``eta`` > 0:

            * ``Deta``: Tajima's D computed with ``eta`` instead of ``S``.
            * ``Dflstar``: Fu and Li's D* (without outgroup; *Genetics* **133**:693-709).
            * ``Fstar``: Fu and Li's F* (without ougroup; *Genetics* **133**:693-709).

        The following statistic is computed only if an outgroup is found:
        
            * ``Sext``: Mutations on external branches.
            
        The following statistics are computed only if an outgroup is
        found and ``eta`` > 0:
        
            * ``Dfl``: Fu and Li's D (*Genetics* **133**:693-709).
            * ``F``: Fu and Li's F (*Genetics* **133**:693-709).

        The following statistics are computed only if sequences are
        coding *dataType* = 4-7:

            * ``ncodon1mut``: Number of codon sites with exactly one mutation.
            * ``NSsites``: Average number of non-synonymous sites.
            * ``nstop``: Number of codon sites with a stop codon.
            * ``nsyn``: Number of codon sites with a synonymous change.
            * ``PiNS``: Nucleotide diversity computed on non-synonymous sites.
            * ``PiS``: Nucleotide diversity computed on synonymous sites.
            * ``SNS``: Number of non-synonymous polymorphic sites.
            * ``SS``: Number of synonymous polymorphic sites.
            * ``Ssites``: Number of synonymous sites.
            * ``tWNS``: Watterson's theta computed on non-synonymous sites.
            * ``tWS``: Watterson's theta computed on synonymous sites.

        The following statistics are computed only if sequences are
        coding *dataType* = 4-7 and an outgroup is found:
        
            * ``MK``: McDonald-Kreitman test table (*Nature*
              **351**:652-654, 1991).
            * ``NI``: Neutrality index (Rand and Kann *Mol. Biol. Evol.*
              **13**:735-748).

        The returned dictionary also contains a nest dictionary ``options``
        which feedbacks the values used at function call.

        .. versionchanged:: 2.0.2
            The following statistics are now computed only if ``S`` > 0:
            ``D``, ``Deta``, ``Dflstar``, ``Fstar``, ``Dfl``, ``F``.

        .. versionchanged:: 2.1.0
            The statistics not computed are now exported and set to
            ``None``.

        """

        if dataType not in (1,2,3,4,5,6,7):
            raise ValueError, 'invalid argument dataType'
        
        if 'BppDiversity' not in dir(egglib_binding):
            raise RuntimeError, 'Bio++ tools are not available'
        
        results = {}
        
        results['options'] = {
            'dataType': dataType
        }
        
        for key in ['S' 'Sinf', 'Ssin', 'eta', 'thetaW', 'He', 'T83',
        'Ti', 'Tv', 'H', 'K', 'rhoH', 'Tv', 'TiTv', 'D', 'Deta',
        'Dflstar', 'Fstar', 'Sext', 'Dfl', 'F', 'ncodon1mut', 'NSsites',
        'nstop', 'nsyn', 'PiNS', 'PiS', 'SNS', 'SS', 'Ssites', 'tWNS',
        'tWS', 'MK', 'NI']:
            results[key] = None
            
        bppDiversity = egglib_binding.BppDiversity()
        bppDiversity.load(self._object, dataType)

        results['S'] = bppDiversity.S()
        results['Sinf'] = bppDiversity.Sinf()
        results['Ssin'] = bppDiversity.Ssin()
        results['eta'] = bppDiversity.eta()
        results['thetaW'] = bppDiversity.tW()
        results['He'] = bppDiversity.He()
        results['T83'] = bppDiversity.T83()
        results['Ti'] = bppDiversity.Ti()
        results['Tv'] = bppDiversity.Tv()
        results['H'] = bppDiversity.H()
        results['K'] = bppDiversity.K()
        results['rhoH'] = bppDiversity.rhoH()

        if results['Tv']:
            results['TiTv'] = bppDiversity.TiTv()

        if results['S']:
            results['D'] = bppDiversity.D()

        if results['eta']:
            results['Deta'] = bppDiversity.Deta()
            results['Dflstar'] = bppDiversity.Dflstar()
            results['Fstar'] = bppDiversity.Fstar()
                
        if bppDiversity.hasOutgroup():
            results['Sext'] = bppDiversity.Sext()
            if results['eta']:
                results['Dfl'] = bppDiversity.Dfl()
                results['F'] = bppDiversity.F()
            
        if dataType>3:
            results['ncodon1mut'] = bppDiversity.ncodon1mut()
            results['NSsites'] = bppDiversity.NSsites()
            results['nstop'] = bppDiversity.nstop()
            results['nsyn'] = bppDiversity.nsyn()
            results['PiNS'] = bppDiversity.PiNS()
            results['PiS'] = bppDiversity.PiS()
            results['SNS'] = bppDiversity.SNS()
            results['SS'] = bppDiversity.SS()
            results['Ssites'] = bppDiversity.Ssites()
            results['tWNS'] = bppDiversity.tWNS()
            results['tWS'] = bppDiversity.tWS()

        if dataType>3 and bppDiversity.hasOutgroup():
            results['MK'] = bppDiversity.MK()
            results['NI'] = bppDiversity.NI()
            

        return results

    ####################################################################

    def matrixLD(self, minimumExploitableData=1., ignoreFrequency=0, 
                    validCharacters="ACGT", missingData="MRWSYKBDHVN?-"):
                        
        """
        Generates the matrix of linkage disequilibrium between all pairs
        of polymorphic sites. The options have the same meaning as for
        :meth:`.polymorphism`.
        
        Returns a dictionary containing the following keys:
            - ``minimumExploitableData`` (value of input parameter),
            - ``ignoreFrequency`` (value of input parameter),
            - ``n`` (number of pairs of sequences),
            - ``d`` (alignment distance between polymorphic sites),
            - ``D`` (*D* linkage disequilibrium statistic),
            - ``Dp`` (*D'* linkage disequilibrium statistic),
            - ``r`` (Pearson's correlation coefficient),
            - ``r2`` (square Pearson's correlation coefficient).

        ``D``, ``Dp``, ``r`` and ``r2`` are four possible measures of
        linkage disequilibrium. ``minimumExploitableData``,
        ``ignoreFrequency`` and ``n`` are provided as integer values.
        ``d``, ``D``, ``Dp``, ``r`` and ``r2`` are provided as nested
        dictionaries containing the matrix of values of the
        corresponding statistic between all pairs of polymorphic sites.
        The individual values can be accessed as this (example given for
        ``D``): ``matrixLD()['D'][i][j]`` where *i* is the alignment
        position of the first site and *j* is the alignment position of
        the second site such as *i* < *j*. The polymorphic sites as the
        same as those returned by ``polymorphism()['siteIndices']``
        (called on the same object with the same configuration options).
        """
        
        if not isinstance(ignoreFrequency, int):
            raise TypeError, 'argument ignoreFrequency must be an integer'

        linkageDisequilibrium = egglib_binding.LinkageDisequilibrium()
        linkageDisequilibrium.load(self._object, minimumExploitableData,
            ignoreFrequency, '%s %s' %(validCharacters,missingData))

        results = {
            'minimumExploitableData': minimumExploitableData,
            'ignoreFrequency': ignoreFrequency,
            'D': {}, 'd': {}, 'Dp': {}, 'r': {}, 'r2': {},
            'n': linkageDisequilibrium.numberOfPairs()
        }
                
        for i in range(results['n']):
            s1= linkageDisequilibrium.site1(i)
            s2= linkageDisequilibrium.site2(i)
            d = linkageDisequilibrium.d(i)
            D = linkageDisequilibrium.D(i)
            Dp= linkageDisequilibrium.Dp(i)
            r = linkageDisequilibrium.r(i)
            r2= linkageDisequilibrium.r2(i)
            
            if s1 not in results['d']:
                results['d'][s1] = {}
                results['D'][s1] = {}
                results['Dp'][s1] = {}
                results['r'][s1] = {}
                results['r2'][s1] = {}
            results['d'][s1][s2] = d
            results['D'][s1][s2] = D
            results['Dp'][s1][s2] = Dp
            results['r'][s1][s2] = r
            results['r2'][s1][s2] = r2
            
        return results

    ####################################################################

    def Rmin(self, minimumExploitableData=1., ignoreFrequency=0, 
                    validCharacters="ACGT", missingData="MRWSYKBDHVN?-"):
                        
        """
        Computes the minimal number of recombination events
         
        The computation is performed as described in Hudson, RR and
        NL Kaplan. 1985. Statistical properties of the number of
        recombination events in the history of a sample of DNA
        sequences. Genetics 111: 147-164. The returned parameter is
        the minimal number of recombination events, given by the
        number of non-overlapping pairs of segregating sites violating
        the rule of the four gamete. Only sites with two alleles are
        considered. Note that homoplasy (multiple mutations) mimicks
        recombination. The result of this function is not stored
        in this instance, and re-computed at each call.
        """

        if not isinstance(ignoreFrequency, int):
            raise TypeError, 'argument ignoreFrequency must be an integer'

        linkageDisequilibrium = egglib_binding.LinkageDisequilibrium()
        linkageDisequilibrium.load(self._object, minimumExploitableData,
            ignoreFrequency, '%s %s' %(validCharacters,missingData))

        return linkageDisequilibrium.Rmin(self._object)

    ####################################################################

    def fix_gap_ends(self):
        
        """
        Replaces all leading or trailing alignment gaps (``-``) by
        missing data symbols (``?``). Internal alignment gaps (those
        having at least one character other than ``-`` and ``?`` at each
        side) are left unchanged.
        """
        
        for i in range(len(self)):
            
            # determines end points
            L=0
            while L<self.ls() and (self.get(i,L)=='-' or self.get(i,L)=='?'):
                L+=1
            R=self.ls()-1
            while R>=0 and R>L and (self.get(i,R)=='-' or self.get(i,R)=='?'):
                R-=1  # no overlap
            
            # modifies characters
            for j in range(L):
                self.set(i,j,'?')
            for j in range(R+1, self.ls()):
                self.set(i,j,'?')

    ####################################################################

    def phyml(self):
        
        """
        Returns a phyml-formatted string representing the content of the
        instance. The phyml format is suitable as input data for PhyML
        and PAML software. Raises a :class:`ValueError` is any name of
        the instance contains at least one character of the following
        list: "``()[]{},;``" as well as spaces, tabs, newlines and
        linefeeds. Group information is never exported.
        """
        
        for n,s,g in self:
            if len( set('(){}[],; \n\r\t').intersection(n) ):
                raise ValueError, 'for phyml format conversion, invalid character in name: %s' %n

        return '%d %d\n' %(self.ns(), self.ls()) + '\n'.join([ '%s  %s' %(n, s) for n,s,g in self ])

    ####################################################################
    
    def phylip(self, format='I'):
        
        """
        Returns a phyml-formatted string representing the content of the
        instance. The phyml format is suitable as input data for PhyML
        and PAML software. Raises a :class:`ValueError` is any name of
        the instance contains at least one character of the following
        list: "``()[]{},;``" as well as spaces, tabs, newlines and
        linefeeds. Group labels are never exported. Sequence names
        cannot be longer than 10 characters. A :class:`ValueError` will
        be raised if a longer name is met. *format* must be 'I' or 'S'
        (case-independent), indicating whether the data should be
        formatted in the sequential (S) or interleaved (I) format (see
        PHYLIP's documentation for definitions).
        """
        
        BLOCK = 10
        NBLOCKS = 6
        
        for n,s,g in self:
            if len( set('(){}[],; \n\r\t').intersection(n) ):
                raise ValueError, 'phylip format conversion error, invalid character in name: %s' %n
            if len(n)>10:
                raise ValueError, 'phylip format conversion error, this name is too long: %s' %n

        if format.upper() not in set('IS'): 
            raise ValueError, 'unknown value for option `format`: %s' %str(format)

        if format.upper()=='I':
            string = '  {0:d} {1:d} I\n'.format(self.ns(), self.ls())
        else:
            string = '  {0:d} {1:d}\n'.format(self.ns(), self.ls())
            
        c = 0

        for n,s,g in self:
                    
            ci = c
            string += '{0}{1}'.format(n.ljust(10), s[:BLOCK-10])
            ci += BLOCK-10
            n = 0
            while n < (NBLOCKS-1) and ci<self.ls():
                string += ' {0}'.format(s[ci:ci+BLOCK])
                ci += BLOCK
                
                # if interleaved, only 1 line
                if format=='I' or format=='i': n += 1

            string += '\n'
        c = ci

        # if sequential, c should be full
        while c<self.ls():
            
            string += '\n'
            for n,s,g in self:
                    
                ci = c
                n = 0
                while n < NBLOCKS and ci<self.ls():
                    if n!=0: string += ' '
                    string += '{0}'.format(s[ci:ci+BLOCK])
                    ci += BLOCK
                    n += 1
                string += '\n'
            c = ci
            
        return string

    ####################################################################
    
    def slider(self, wwidth, wstep):
        
        """
        Provides a means to perform sliding-windows analysis over the
        alignment. This method returns a generator that can be used as
        in ``for window in align.slider(wwidth,wstep)``, where each
        step *window* of the iteration will be a :class:`~egglib.Align` instance
        of length *wwidth* (or less if not enough sequence is available
        near the end of the alignment). Each step moves forward
        following the value of *wstep*.
        """
        
        cache = 0
        for i in range(0, self.ls(), wstep):
            
            # avoids redundant windows
            if min(self.ls(), i+wwidth) == cache:
                break
                
            yield self.extract(i, i+wwidth)
            cache = min(self.ls(),i+wwidth)

    ####################################################################
    
    def dataMatrix(self, mapping='ACGT', others=999):
        
        """
        Returns a copy of the current instance as a :class:`~egglib.egglib_binding.DataMatrix`
        Mapping must be a string of type 'ACGT' indicating valid characters
        that will be encoded by their position in the string (ie 0,1,2,3).
        *others* gives the index to affect to characters not found in
        the mapping string.
        """
        
        data = egglib_binding.DataMatrix(self.ns(), self.ls())
        for i,v in enumerate(self):
            data.populationLabel(i, self.group(i))
            for j in range(self.ls()):
                c = self.get(i,j)
                if c in mapping: index = mapping.index(c)
                else: index = others
                data.set(i,j, index)
        return data

    ####################################################################

    def simErrors(self, rate):

        """ Randomly introduces missing data. *rate* is the desired
        proportion of missing data. Replaces random valid positions by
        ``N``. There should be not missing data in the original object.
        Note that the module :mod:`numpy` is required, and that this
        method might be inefficient for large error rates. Changes the
        current object and returns nothing.
        
        .. versionchanged:: 2.1.0
           Restricted to :class:`~egglib.Align instances.
        """

        # finds the number of errors
        try: import numpy
        except ImportError: raise Exception, 'the package numpy is required for random errors in alignments'
        B= numpy.random.binomial(self.ns()*self.ls(), rate)
        
        # introduces the errors at random positions
        while B:
            X= numpy.random.randint(self.ns())
            Y= numpy.random.randint(self.ls())
            if self.sequence(X)[Y].upper() not in 'ACGT': continue
            self.set(X,Y,'N')
            B-=1


########################################################################

class SSR(object):
    
    """
    SSR data container. This class is essentially a wrapper for the
    lower-level class :class:`~egglib.egglib_binding.DataMatrix` with
    parser/formatter methods.
    
    The user doesn'nt normally need to manipulate directly the class
    attributes described below.
    
    The class attribute :attr:`.dataMatrix` holds the :class:`~egglib.egglib_binding.DataMatrix`
    instance contained by this instance. :attr:`loci` contains the list
    (in the correct order) of locus names. :attr:`individuals` stores
    the names of individuals (or automatically-generated names if they
    were not available) and maps them to 1 or 2 (for diploid data)
    indices of dataMatrix. :attr:`individuals` maps individuals to
    populations names. :attr:`populations` stores the population
    names. The labels used in the :class:`~egglib.egglib_binding.DataMatrix`
    instance are the indices of this   list. :attr:`.missing` stores the
    value used to identify missing data (``None`` if missing data are
    not allowed).
    
    .. versionadded:: 2.0.1

    .. versionchanged:: 2.1.0
        Individual to population mapping and string formatting added.
    """
    
    
    def __init__(self):
        
        """ The constructor initializes an empty instance. """
        
        self.dataMatrix = egglib_binding.DataMatrix()
        self.loci = []
        self.individuals = {}
        self.populations = []
        self.indiv2pop = {}
        self.missing = None

    ####################################################################

    def clear(self):
        
        """ Removes all data from the instance. """
        
        self.dataMatrix.clear()
        self.loci = []
        self.populations = []
        self.individuals.clear()
        self.indiv2pop.clear()


    ####################################################################

    def load(self, dataMatrix, sampleConfiguration=None):
        
        """
        Imports the data present in the
        :class:`~egglib.egglib_binding.DataMatrix` instance passed as
        *dataMatrix*. Note that the :class:`SSR` instance is supposed to
        take ownership of the :class:`~egglib.egglib_binding.DataMatrix`
        instance that should not be modified outside the class and that
        will be clear if the :class:`SSR` is cleared. The appropriate
        behaviour is to delete any outside reference to *dataMatrix*
        after passing it to this method. The argument
        *sampleConfiguration* is an iterable indicating the number of
        samples per population. There must be one item per population.
        If items are integers, they give the number of diploid samples
        (one random chromosome per individual), otherwise they must be a
        sequence of two integers giving the number of diploid and
        haploid samples (in this order). The total sum of given integers
        must match the number of genotypes of *dataMatrix*. If
        *sampleConfiguration* is ``None``, it it assumed that all
        samples are haploid samples. The population labels from the
        passed :class:`~egglib.egglib_binding.DataMatrix` are discarded
        unless *sampleConfiguration* is ``None``. Below are some
        examples of accepted values for the second argument. 20
         individuals from 4 populations with boths chromosomes sampled: 
        ``[10, 10, 10, 10]``. 10 individuals from 4 populations with
        one chromosome sampled: ``[(0,20), (0,20), (0,20), (0,20)]``.
        A mixture of samples: ``[10, (0,20), (5,10), (1,18)]``. In all 
        these three examples, 20 chromosomes are sampled from each
        population, suming up to a total of 80 samples.
        
        .. versionchanged:: 2.1.0
            Population names were previously integers, they now are
            converted to strings.
        """

        self.clear()
        self.dataMatrix = dataMatrix
        
        # creates automatically locus names
        if dataMatrix.numberOfSites(): L = len(str(dataMatrix.numberOfSites()))
        else: L = 0
        self.loci = ['locus%s' %str(i).rjust(L, '0') for i in range(dataMatrix.numberOfSites())]
        
        # imports the instance's labels
        for i in range(dataMatrix.numberOfSequences()):
            label = dataMatrix.populationLabel(i)
            if label not in self.populations:
                self.populations.append(label)
        if len(self.populations): L = len(str(max(self.populations)))
        else: L = 0
        self.populations = ['population%s' %str(i).rjust(L, '0') for i in self.populations]
        
        # sets the individual assignations and population labels from the sample configuration (is passed)
        if sampleConfiguration!=None:
            L = len(str(len(sampleConfiguration)))
            chr=0
            ind=0
            pop=0
            poplabel = 'population%s' %"1".rjust(L, '0')
            try:
                for i in sampleConfiguration:
                    if isinstance(i, int):
                        double = i
                        single = 0
                    else:
                        try:
                            double, single = i
                        except ValueError:
                            raise ValueError, '`SSR.load` received invalid value for argument `sampleConfiguration`'
                    
                    for j in range(double):
                        if chr>=self.dataMatrix.numberOfSequences():
                            raise ValueError, 'argument `sampleConfiguration` of `SSR.load` defines too many samples'
                        self.individuals[ind] = (chr, chr+1)
                        self.indiv2pop[ind] = poplabel
                        self.dataMatrix.populationLabel(chr,   pop)
                        self.dataMatrix.populationLabel(chr+1, pop)
                        chr+=2
                        ind+=1
                    for j in range(single):
                        if chr>=self.dataMatrix.numberOfSequences():
                            raise ValueError, 'argument `sampleConfiguration` of `SSR.load` defines too many samples'
                        self.individuals[ind] = chr
                        self.indiv2pop[ind] = poplabel
                        self.dataMatrix.populationLabel(chr, pop)
                        chr+=1
                        ind+=1
                    pop+=1
                    poplabel = 'population%s' %str(pop+1).rjust(L, '0')
            except TypeError:
                raise TypeError, '`SSR.load` received invalid value for argument `sampleConfiguration`'

            if chr!=self.dataMatrix.numberOfSequences():
                raise ValueError, 'argument `sampleConfiguration` of `SSR.load` doesn\'t define the right number of samples'

        # assign by default each genotype to a different individual
        else:
            for i in range(self.dataMatrix.numberOfSequences()):
                self.individuals[i] = i
                self.indiv2pop[i] = self.dataMatrix.populationLabel(i)

        # changes individuals' names to strings
        L = len(str(len(self.individuals)))
        for i in self.individuals.keys():
            name = 'indiv%s' %str(i+1).rjust(L, '0')
            self.individuals[name] = self.individuals[i]
            self.indiv2pop[name] = self.indiv2pop[i]
            del self.individuals[i]
            del self.indiv2pop[i]


    ####################################################################
    
    def parse(self, string, diploid=True, genotypeSeparator=None,
                          alleleSeparator='/', header=True, missing='999'):
        
        """ Imports data from the string *string*. The data should
        follow the following format: one line containing locus names
        (if *header* is ``True``) and then one line per individual. The
        header needs not to be aligned with the data matrix. It is only
        required that the number of items on the first line matches the
        number of genotypes given for each individual. Each line is made
        of a population name (if *groups* is ``True``), the individual
        name followed by the appropriate number of genotype values. A
        given genotype is coded by two (if *diploids* is ``True``) or
        one (otherwise) integer. If the latter case, the two values must
        be separated by *alleleSeparator*. *genotypeSeparator* gives the
        separator between values for one individual. The default value
        of *genotypeSeparator* matches all white spaces (including space
        and tabulation). The same separator is used between population
        and individual labels. The user must specify another value for
        *genotypeSeparator* to be able to import data with that have
        spaces in names. Unless using the default separator, ensure that
        the separator is not duplicated in the string (e.g. two spaces
        in a row between header items): this is supported only for the
        default (refer to the standard library :meth:`str.split` method
        for more details. Missing data are coded by the string or the
        integer given by *missing*.
        
        Example (to be read with default argument values)::

                    locus1  locus2   locus3
          pop1 ind1 001/001 001/002  001/003
          pop1 ind2 002/002 001/001  001/001
          pop1 ind3 002/002 001/001  002/002
          pop1 ind4 001/002 001/002  003/002
          pop2 ind5 001/002 003/003  000/000
          pop2 ind6 003/004 004/004  002/003
        """

        # transforms the missing symbol to a string
        missing=str(missing)
        self.missing=999

        # empties the instance
        self.clear()

        # enforces foreign file support
        string=string.replace('\r', '\n')
        string=string.replace('\n\n', '\n')

        # gets the number and names of loci / splits the strings to lines
        string=string.strip()
        if header:
            # cleans any leading separators
            if genotypeSeparator:
                string = string.lstrip(genotypeSeparator)
            else:
                string = string.lstrip()
                
            # imports the name of loci
            string=string.split('\n')
            if not len(string[0]):
                raise IOError, 'SSR: no data'
            else: self.loci=string[0].split(genotypeSeparator)
            del string[0]
        else:
            string=string.split('\n')
            if not len(string[0]): raise IOError, 'SSR: no data'
            nloci=len(string[0].split(genotypeSeparator))-2
            if nloci<1: raise Exception, 'SSR: file format error: not enough loci'
            nbits= len(str(nloci))
            self.loci= ['locus%s'%str(i).rjust(nbits, '0') for i in range(1, nloci+1)]
            
        # resizes the data matrix
        if diploid: 
            self.dataMatrix.resize(len(string)*2, len(self.loci))
        else:
            self.dataMatrix.resize(len(string), len(self.loci))
            
        # loops over the different individuals
        c=0
        for i in string:
            
            # gets the items
            j=i.split(genotypeSeparator)
            if len(j)!=(len(self.loci)+2):
                raise IOError, 'SSR: inconsistent data matrix'

            # sets the group
            if j[0] not in self.populations:
                self.populations.append(j[0])
            self.dataMatrix.populationLabel(c, int(self.populations.index(j[0])))
            self.indiv2pop[j[1]] = j[0]

            # keeps name
            if diploid:
                self.individuals[j[1]] = (c, c+1)
            else:
                self.individuals[j[1]] = c

            # reads and stores the genotype for each locus
            d=0
            for k in j[2:]:
                
                if diploid:
                    k = k.split(alleleSeparator)
                    if len(k)!=2:
                        raise IOError, 'SSR: if `diploid` is `True`, expects two alleles for each genotype'
                    for i in range(2):
                        if k[i]=='999' and missing!='999':
                            raise 'SSR: 999 is not accepted as allele value, because 999 is used internally as symbol for missing data'
                        if k[i]==missing:
                            k[i] = 999
                        try:
                            k[i] = int(k[i])
                        except ValueError:
                            raise IOError, 'SSR: all allele values (except missing data) must be integers (got %s)' %k[i]
                        self.dataMatrix.set(c+i, d, k[i])
            
                else:
                    if k=='999' and missing!='999':
                        raise 'SSR: 999 is not accepted as allele value, because 999 is used internally as symbol for missing data'
                    if k==missing:
                        k = 999
                    try:
                        k = int(k)
                    except ValueError:
                        raise IOError, 'SSR: all allele values (except missing data) must be integers (got %s)' %k
                    self.dataMatrix.set(c, d, k)
                d+=1
            
            if diploid: c+=2
            else: c+=1

    ####################################################################

    def stats(self):
        
        """ Computes diversity statistics from the currently loaded
        data. Returns a dictionary containing the following values:
        
        * *k*, the number of alleles.
        * *V*, the variance of allele size.
        * *He*, the expected heterozygosity.
        * *thetaI*, theta assuming IAM.
        * *thetaHe*, theta from *He* (assuming SMM).
        * *thetaV*, theta from *V* (assuming SMM).
        
        Each entry is a list of the corresponding statistics computed
        for the corresponding locus. """
        
        MD = egglib_binding.MicrosatelliteDiversity()
        if self.missing:
            MD.load(self.dataMatrix, self.missing, False)
        else:
            MD.load(self.dataMatrix, 0, True)
        
        return {
            'k': [MD.numberOfAlleles(i) for i in range(self.dataMatrix.numberOfSites())],
            'V': [MD.sizeVariance(i) for i in range(self.dataMatrix.numberOfSites())],
            'He': [MD.He(i) for i in range(self.dataMatrix.numberOfSites())],
            'thetaI': [MD.thetaAssumingIAM(i) for i in range(self.dataMatrix.numberOfSites())],
            'thetaHe': [MD.thetaAssumingSMMfromHe(i) for i in range(self.dataMatrix.numberOfSites())],
            'thetaV': [MD.thetaAssumingSMMfromSizeVariance(i) for i in range(self.dataMatrix.numberOfSites())]
            }


    ####################################################################

    def Fstats(self, locus=None):
        
        """ Computes F-statistics from the currently loaded data, using
        Weir and Cockheram (1984) method. If *locus* is an integer, the
        statistics for that locus are returned. If *locus* is ``None``,
        the multi-locus version of the statistics are returned. This
        method returns a ``(Fis, Fst, Fit)`` tuple. If one or several
        values cannot be computed (due to lack of one or more components
        of the variance), the corresponding value is replaced by
        ``None``. This method requires that all genotypes have two
        alleles. In case of missing data, complete genotypes (i.e. data
        for one individual at a given locus) are removed."""

        # function to compute components of variance
        def site_Fstats(locus):
            FS = egglib_binding.FStatistics()
            FS.reserve(len(self.individuals))
            for I in self.individuals:
                try:
                    G1, G2 = self.individuals[I]
                except TypeError:
                    raise ValueError, 'method `SSR.Fstats` requires that all data are diploid'

                A1 = self.dataMatrix.get(G1, locus)
                A2 = self.dataMatrix.get(G2, locus)
                P = self.dataMatrix.populationLabel(G1)

                # skips g enotypes with at least one missing allele
                if A1==self.missing or A2==self.missing:
                    continue

                FS.loadIndividual(A1, A2, P)

            return (FS.Vallele(), FS.Vindividual(), FS.Vpopulation())


        if locus==None:
            FisTop = 0.
            FisBot = 0.
            FstTop = 0.
            FstBot = 0.
            FitTop = 0.
            FitBot = 0.
            for i in range(len(self.loci)):
                VA, VI, VP = site_Fstats(i)
                if VA+VI+VP:
                    FitTop+= VA
                    FitBot+= (VA+VI+VP)
                    FstTop+= VP
                    FstBot+= (VA+VI+VP)
                if VA+VI:
                    FisTop+= VA
                    FisBot+= (VA+VI)
            if FisBot!=0.:
                Fis = 1 - FisTop/FisBot
            else:
                Fis = None
            if FstBot!=0.:
                Fst = FstTop/FstBot
            else:
                Fst = None
            if FitBot!=0.:
                Fit = 1 - FitTop/FitBot
            else:
                Fit = None
                
            return (Fis, Fst, Fit)

        else:
            if locus>=len(self.loci):
                raise IndexError, 'SSR: invalid locus index'
            VA, VI, VP = site_Fstats(locus)
            if VA+VI+VP:
                Fit = 1. - 1. * VA / (VA+VI+VP)
                Fst =      1. * VP / (VA+VI+VP)
            else:
                Fit = None
                Fst = None
            if VA+VI:
                Fis = 1. - 1. * VA / (VA+VI)
            else:
                Fis = None
            
            return (Fis, Fst, Fit)


    ####################################################################

    def numberOfLoci(self):
        
        """ Returns the number of loci of the data currently loaded. """
        
        return len(self.loci)
        
    
    ####################################################################

    def numberOfGenotypes(self):
        
        """ Returns the number of genotypes of the data currently
        loaded. """
        
        return len(self.individuals)

    ####################################################################

    def __str__(self):

        return self.str()

    ####################################################################

    def str(self):
        
        """
        Returns a string representation of the object. The string can
        be parsed by the method :meth:`.parse` using default options.
        
        .. versionadded: 2.1.0
        """
        
        Lpop = max(map(len, self.populations)) if len(self.populations) else 0
        Lind = max(map(len, self.individuals)) if len(self.individuals) else 0
        Lloc = max(map(len, self.loci)) if len(self.loci) else 0
        maxallele = 0

        genotypes = []
        for indiv in sorted(self.individuals):
            geno = self.individuals[indiv]
            genotypes.append([])
            if isinstance(geno, int):
                for j in range(len(self.loci)):
                    A = self.dataMatrix.get(geno,j)
                    if A>maxallele: maxallele = A
                    genotypes[-1].append([A])
            elif len(geno)==2:
                for j in range(len(self.loci)):
                    A1,A2 = self.dataMatrix.get(geno[0],j), self.dataMatrix.get(geno[1],j)
                    genotypes[-1].append([A1, A2])
                    if A1>maxallele: maxallele = A1
                    if A2>maxallele: maxallele = A2
            else:
                raise RuntimeError, 'invalid number of genotypes in SSR instance'

        maxallele = max([maxallele, 999])
        Lall = len(str(maxallele))
        Lloc = max([Lloc, Lall*2+1])
        
        header = ' '.join([' '*(Lpop+Lind+1)] + [i.center(Lloc) for i in self.loci])
        lines = [header]
        i=0

        for indiv in sorted(self.individuals):
            line = [self.indiv2pop[indiv], indiv]
            for j in range(len(self.loci)):
                item = '/'.join([str(a).rjust(Lall, '0') for a in genotypes[i][j]])
                item = item.center(Lloc)
                line.append(item)
            lines.append(' '.join(line))
            i+=1
            
        return '\n'.join(lines)
        


########################################################################

class TIGR(object):
    
    """ Automatic wrapper around the TIGR XML format for storing genome
    annotation data. """
    
    ####################################################################
    
    def __init__(self, fname=None, string=None):
        
        """ Object initialization: ``TIGR(fname=None, string=None)``
        
        *fname* must be the name of a file containing TIGR-formatted
        XML data. *string* should be directly a string containing
        TIGR-formatted XML data. It is not allowed to specify both
        *fname*  and *string* to non-``None`` values, but at least one
        of the two arguments must be specified (it is currently
        impossible to create an empty instance). """
        
        self.document = None
        
        if fname==None and string==None:
            
            # default initialization
            #DOMimplementation = xml.dom.minidom.getDOMImplementation()
            #self.document = DOMimplementation.createDocument(None, 'TIGR', None)
            
            raise ValueError, 'TIGR class\'s constructor expects one argument exactly'
            
        else:
            
            # initializes the stream
            stream = None
            if fname!=None and string!=None:
                raise ValueError, 'TIGR class\'s constructor expects one argument exactly'
            if fname:
                self.document = xml.dom.minidom.parse(fname)
            if string:
                self.document = xml.dom.minidom.parseString(string)

            # check of document structure
            if (len(self.document.getElementsByTagName('TIGR'))!=1 or
                len(self.document.getElementsByTagName('TIGR')[0].getElementsByTagName('ASSEMBLY'))!=1):
                    raise ValueError, 'TIGR class expects TIGR-formatted XML data'

            # makes a link to the assembly level
            self.assembly = self.document.getElementsByTagName('TIGR')[0].getElementsByTagName('ASSEMBLY')[0]

            # gets information
            self.start, self.stop = self._getCoordset(self.assembly.getElementsByTagName('COORDSET')[0], 0)
            self.length = self.stop-self.start+1


    ####################################################################

    def _getCoordset(self, coordset, shift):

        END5 = coordset.getElementsByTagName('END5')[0].firstChild.data
        END3 = coordset.getElementsByTagName('END3')[0].firstChild.data
        return int(END5)-1+shift, int(END3)-1+shift
        
    ####################################################################

    def __del__(self):
        
        if self.document: self.document.unlink()

    ####################################################################
        
    def extract(self, start, stop):
        
        """ Returns a GenBank instance containing the sequence of the
        range [*start*, *stop*] and the features that are completely
        included in that range. Note that positions must be expressed
        in the TIGR system own coordinate system. """
        
        
        if start>=stop or start<self.start or stop>self.stop:
            raise ValueError, 'invalid range for extracting TIGR data'

    
        shift = self.start - start # it is therefore negative

        gb = GenBank()

        # gets the sequence
        
        sequence = self.assembly.getElementsByTagName('ASSEMBLY_SEQUENCE')[0].firstChild.data[start-self.start:stop+1-self.start]

        gb.set_sequence(sequence)
        
        # gets the features
        
        for TU in self.assembly.getElementsByTagName('GENE_LIST')[0].getElementsByTagName('PROTEIN_CODING')[0].getElementsByTagName('TU'):
            a,b = self._getCoordset(TU.getElementsByTagName('COORDSET')[0], 0)
                       
            if a>=start and b<=stop:
                            
                # collects the gene-level information
                
                name = TU.getElementsByTagName('FEAT_NAME')[0].firstChild.data
                com_names = [i for i in TU.getElementsByTagName('COM_NAME')
                                if i.getAttribute('IS_PRIMARY') == '1']
                if len(com_names)==0: com_name = name
                else: com_name = com_names[0].firstChild.data
                
                models= TU.getElementsByTagName('MODEL')
                if len(models)!=1:
                    raise ValueError, '%s has %d models' %(name, len(models))
                imgag_comment = models[0].getAttribute('COMMENT')
                
                # collect the exon and cds positions

                mRNA = []
                CDS = []
                
                for exon in models[0].getElementsByTagName('EXON'):
                    mRNA.append(self._getCoordset(exon.getElementsByTagName('COORDSET')[0], shift))
                    cds = exon.getElementsByTagName('CDS')
                    if len(cds)>1:
                        raise ValueError, 'an exon of %s has more than one CDS fragment' %name
                    if len(cds)==1:
                        CDS.append(self._getCoordset(cds[0].getElementsByTagName('COORDSET')[0], shift))

                forward= set(map(lambda (i,j): i<j, mRNA+CDS+[(a,b)]))
                if len(forward)!=1:
                    raise IOError, 'inconsistent orientation for gene %s at positions %d-%d' %(name, a,b)
                forward=forward.pop()
                if not forward:
                    c=b
                    b=a
                    a=c
                    mRNA = map(sorted, mRNA)
                    mRNA.sort(cmp=lambda x,y: cmp(x[0], y[0]))
                    CDS = map(sorted, CDS)
                    CDS.sort(cmp=lambda x,y: cmp(x[0], y[0]))
                
                # adds the features    

                geneF = GenBankFeature(gb)
                geneLoc = GenBankFeatureLocation()
                geneLoc.addBaseRange(a+shift, b+shift)
                if not forward: geneLoc.setComplement()
                geneF.set('gene', geneLoc, gene=name, product=com_name, note=imgag_comment)
                gb.add_feature(geneF)
                
                mRNAF = GenBankFeature(gb)
                mRNALoc = GenBankFeatureLocation()
                for i,j in mRNA: mRNALoc.addBaseRange(i, j)
                if not forward: mRNALoc.setComplement()
                mRNAF.set('mRNA', mRNALoc, gene=name, product=com_name, note=imgag_comment)
                gb.add_feature(mRNAF)
                
                CDSF = GenBankFeature(gb)
                CDSLoc = GenBankFeatureLocation()
                for i,j in CDS: CDSLoc.addBaseRange(i, j)
                if not forward: CDSLoc.setComplement()
                CDSF.set('CDS', CDSLoc, gene=name, product=com_name, note=imgag_comment)
                gb.add_feature(CDSF)

        return gb


########################################################################

class GenBankFeatureLocation(object):
    
    """ Holds the location of a GenBank feature. Supports various forms
    of location as defined in the GenBank format specification. The
    constructor contains a parser working from a GenBank-formatted
    string. By default, features are on the forward strand and segmented
    features are ranges (not orders). :class:`~egglib.GenBankFeatureLocation`
    supports iteration and allows to iterate over ``(first,last)``
    segments regardless of their types (for a single-base segment a
    position *position*, the tuple ``(position,position)`` is returned;
    similar 2-item tuples are returned for other types of segment as
    well). :class:`~egglib.GenBankFeatureLocation` also supports access (but
    not assignation nor deletion) thought the ``[]`` operator. A
    ``(first,last)`` tuple is returned as for the iterator. Finally, the
    instance can be GenBank-formatted using :func:`.str`. The length of
    the instance is the number of segments. """
    
    ####################################################################

    def __iter__(self):
        """ Iterator. """
        for pos in self._pos: yield pos
        
    def __len__(self):
        """ len(x) """
        return len(self._pos)

    def __getitem__(self, index):
        """ [] accessor. """
        return self._pos[index]

    ####################################################################

    def __str__(self):
        
        if not len(self._pos):
            raise RuntimeError, 'cannot format Genbank feature\'s location: no positions were loaded'
        
        string = []
        
        for i in range(len(self._pos)):
            if self._type[i]=='S':
                string.append( str(self._pos[i][0]  +1  ) )
            elif self._type[i]=='R' or self._type[i]=='C':
                item = ''
                if self._ends[i][0]:
                    item += '<'
                item += str(self._pos[i][0]  +1  )
                item += '.'
                if self._type[i]=='R':
                    item +='.'
                if self._ends[i][1]:
                    item += '>'
                item += str(self._pos[i][1]  +1  )
                string.append( item )
            elif self._type[i]=='B':
                string.append(  '%d^%d' %(self._pos[i][0] +1 , self._pos[i][1] +1) )
            else: raise RuntimeError, 'internal error in GenBankFeatureLocation.__str__'

        string = ','.join(string)

        if len(self._pos)>1:
            if self._range:
                string = 'join(%s)' %string
            else:
                string = 'order(%s)' %string
                
        if self._complement:
            string = 'complement(%s)' %string
            
        return string

    ####################################################################

    def copy(self):
        
        """ Returns a deep copy of the current instance. """
        
        gbfl = GenBankFeatureLocation()
        gbfl._complement = self._complement
        gbfl._range = self._range
        gbfl._type = str(self._type)
        gbfl._pos = [(i,j) for i,j in self._pos]
        gbfl._ends = [(i,j) for i,j in self._ends]

        return gbfl

    ####################################################################

    def __init__(self, string=None):
        self._complement = False # True or False
        self._range = True # True or False
        self._type = '' # string of S|R|B|C for Single, Range, Between and Choice
        self._pos = [] # list of start/stop positions, always in increasing number
        self._ends = [] # list of boolean tuples (True if the segment is 5'/3' partial)
        
        if string!=None:
            self._parse(string)
            
    ####################################################################

    def _parse(self, string):
        
        # process top-level location groups
            
        # complement not incompatible with the rest
        if string[:11]=='complement(':
            if string[-1]!=')':
                raise IOError, 'invalid GenBank feature location: %s' %string
            self._complement = True
            string = string[11:-1]

        # range or order (only remove it)
        if string[:5] == 'join(':
            if string[-1]!=')':
                raise IOError, 'invalid GenBank feature location: %s' %string
            string = string[5:-1]
        elif string[:6] == 'order(':
            if string[-1]!=')':
                raise IOError, 'invalid GenBank feature location: %s' %string
            self._range = True
            string = string[6:-1]
            
        # extracts the segment(s)
        
        segments = string.split(',')

        for segment in segments:
                
            # position formats
            single = re.match('^(\d+)$', segment)
            rangeChoice   = re.match('^(\<?)(\d+)(\<?)(\.\.?)(\>?)(\d+)(\>?)$', segment)
            between = re.match('^(\d+)\^(\d+)$', segment)
            
            if single:
                try:
                    pos = int(single.group(1))  -1
                except ValueError:
                    raise IOError, 'invalid position in GenBankFeature location: %s' %segment
                self.addSingleBase(pos)

            elif rangeChoice:

                # gets positions
                try:
                    start = int(rangeChoice.group(2))  -1
                    stop = int(rangeChoice.group(6))   -1
                except ValueError:
                    raise IOError, 'invalid position in GenBankFeature location: %s' %segment
                    
                # gets partial marks
                Ls = rangeChoice.group(1)+rangeChoice.group(3)
                if Ls=='':
                    Lp = False
                elif Ls=='<':
                    Lp = True
                else:
                    raise IOError, 'invalid specification in GenBankFeature location: %s' %segment
                Rs = rangeChoice.group(5)+rangeChoice.group(7)
                if Rs=='':
                    Rp = False
                elif Rs=='>':
                    Rp = True
                else:
                    raise IOError, 'invalid specification in GenBankFeature location: %s' %segment
    
                # loads appropriate feature
                if rangeChoice.group(4)=='..':
                    self.addBaseRange(start, stop, Lp, Rp)
                elif rangeChoice.group(4)=='.':
                    self.addBaseChoice(start, stop, Lp, Rp)
                else:
                    raise RuntimeError, 'error in GenBankFeatureLocation constructor'

            elif between:
                try:
                    start = int(between.group(1))  -1
                    stop = int(between.group(2))   -1
                except ValueError:
                    raise IOError, 'invalid position in GenBankFeature location: %s' %segment
                if stop!=start+1:
                    raise IOError, 'invalid between-base feature specification: %s' %segment
                self.addBetweenBase(start)
            
            else:
                raise IOError, 'invalid GenBank feature segment positions: %s' %segment
        
    ####################################################################

    def setComplement(self):
        """ Places the feature on the complement strand. """
        self._complement = True

    def setNotComplement(self):
        """ Places the feature on the forward (not complement) strand,
        which is the default. """
        self._complement = False
    
    def isComplement(self):
        """ True if the feature is on the complement strand. """
        return self._complement
        
    ####################################################################
    
    def asOrder(self):
        """ Defines the feature as an order instead of a range. """
        self._range = False

    def asRange(self):
        """ Defines the features as a range, with is the default. """
        self._range = True
    
    def isRange(self):
        """ True if the feature is a range (the default), False if it is
        an order. """
        return self._range

    ####################################################################
    
    def shift(self, shift):
        """ Shift all positions according to the (positive of negative)
        argument. """
        self._pos = map(lambda x: [x[0]+shift, x[1]+shift], self._pos)

    ####################################################################

    def addSingleBase(self, position):
        
        """ Adds a single-base segment to the feature. If no segments
        were entered previously, set the unique segment location.
        *position* must be an integer. All entered positions must be
        larger than any positions entered previously. """
        
        if not isinstance(position, int):
            raise TypeError, 'GenBankFeatureLocation positions must be of type int'
        if len(self._pos) and position<self._pos[-1][1]:
            raise ValueError, 'GenBankFeatureLocation positions must be entered in increasing order'
        
        self._pos.append((position,position))
        self._type += 'S'
        self._ends.append((False,False))
        
    ####################################################################

    def addBetweenBase(self, position):
        
        """ Adds a segment lying between two consecutive  bases. If no
        segments were entered previously, set the unique segment
        location. *position* must be an integer. The feature will be set
        between *position* and *position* + 1. If the feature is intended
        to be placed on the complement strand between positions, say,
        1127 and 1128, one must use ``addBetweenBase(1127)`` in
        combination with ``setComplement()``. All entered positions must
        be larger than any positions entered previously. """
        
        if not isinstance(position, int):
            raise TypeError, 'GenBankFeatureLocation positions must be of type int'
        if len(self._pos) and position<self._pos[-1][1]:
            raise ValueError, 'GenBankFeatureLocation positions must be entered in increasing order'
        
        self._pos.append((position,position+1))
        self._type += 'B'
        self._ends.append((False,False))
        
    ####################################################################

    def addBaseRange(self, first, last, left_partial=False, right_partial=False):
        
        """ Adds a base range the feature. If no segments were
        previously enter, set the unique segment location. *first* and
        *last* must be integers. The feature will be set between *first*
        and *last* positions, including both limits. If the feature
        is intended to be placed on the complement strand between
        positions, say, 1127 and 1482, one must use ``addBaseRange(1127,1482)``
        in combination with ``setComplement()``. All entered positions 
        must be larger than any positions entered previously and *last*
        must be larger than *first* (but can be equal). *left_partial*
        and/or *right_partial* must be set to ``True`` if, respectively,
        the real start of the segment lies 5' of *first* and/or the real
        end of the segment lies beyond *last* (relatively to the forward
        strand and consistently with the numbering system). """
        
        if not isinstance(first, int) or not isinstance(last, int):
            raise TypeError, 'GenBankFeatureLocation positions must be of type int'
        if len(self._pos) and first<self._pos[-1][1]:
            raise ValueError, 'GenBankFeatureLocation positions must be entered in increasing order'
        if (last<first):
            raise ValueError, 'GenBankFeatureLocation positions must be entered in increasing order'
        
        self._pos.append((first,last))
        self._type += 'R'
        self._ends.append(( bool(left_partial), bool(right_partial) ))
        
    ####################################################################

    def addBaseChoice(self, first, last, left_partial=False, right_partial=False):
        
        """ Adds a segment corresponding to a single base chosen within
        a base range. If no segments were previously enter, set the
        unique segment location. *first* and *last* must be integers.
        The feature will be set between *first* and *last* positions,
        including both limits. If the feature is intended to be placed
        on the complement strand between positions, say, 1127 and 1482,
        one must use ``addBaseChoice(1127,1482)`` in combination with
        ``setComplement()``. All entered positions must be larger than
        any positions entered previously and *last* must be strictly
        larger than *first*. *left_partial* and/or *right_partial* must
        be set to ``True`` if, respectively, the real start of the
        segment lies 5' of *first* and/or the real end of the segment
        lies beyond *last* (relatively to the forward strand and
        consistently with the numbering system). """
        
        if not isinstance(first, int) or not isinstance(last, int):
            raise TypeError, 'GenBankFeatureLocation positions must be of type int'
        if len(self._pos) and first<self._pos[-1][1]:
            raise ValueError, 'GenBankFeatureLocation positions must be entered in increasing order'
        if (last<first):
            raise ValueError, 'GenBankFeatureLocation positions must be entered in increasing order'
        if (first==last):
            raise ValueError, 'invalid use of GenBankFeatureLocation.addBaseChoice: start=stop (use addSingleBase instead)'
        
        self._pos.append((first,last))
        self._type += 'C'
        self._ends.append(( bool(left_partial), bool(right_partial) ))
        
    ####################################################################

    def rc(self, length):
        
        """
        Reverse the feature positions: positions are modified to be
        counted from the end. The *length* of the complete sequence must
        be passed.
        """

        self._complement = not self._complement
        self._pos = self._pos[::-1]
        self._ends = self._ends[::-1]
        self._type = self._type[::-1]
        for i in range(len(self._pos)):
            a,b = self._pos[i]
            self._pos[i] = (length-b-1, length-a-1)
            self._ends[i] = ( self._ends[i][1], self._ends[i][0] )

########################################################################

class GenBankFeature(object):
        
    """
    GenBankFeature contains a feature associated to a GenBank instance.
    Instances of this class should not be instantiated or used
    separatedly from a GenBank instance. The constructor creates an
    empty instance (athough a :class:`~egglib.GenBank` instance must be passed
    as *parent*) and either :meth:`.set` or :meth:`.parse` must be used
    subsequently.
    """
        
    ####################################################################
        
    def __init__(self, parent):
        self._parent = parent
        self._type = ''
        self._location = GenBankFeatureLocation()
        self._qualifiers = []

    ####################################################################

    def type(self):
        
        """
        Returns the type string of the instance.
        """
        
        return self._type
        
    ####################################################################

    def qualifiers(self):
        
        """
        Returns a dictionary with all qualifier values. This method
        cannot be used to change data within the instance.
        
        .. versionchanged:: 2.1.0
           Meaning changed.        
        """
        
        return dict(self._qualifiers)
        
    ####################################################################

    def add_qualifier(self, key, value):
        
        """
        Adds a qualifier to the instance's qualifiers.
        """
        
        self._qualifiers.append((key, value))
        
    ####################################################################

    def set(self, type, location, **qualifiers):
            
        """ Updates feature information: *type* is a string identifying
        the feature type (such as gene, CDS, misc_feature, etc.); 
        location must be a :class:`~egglib.GenBankFeatureLocation` instance
        giving the feature's location. Other qualifiers must be passed
        as keyword arguments. Note that *type* can be any string and
        that it is not allowed to use "type" as a qualifier keyword. """
            
        if 'type' in qualifiers:
            raise ValueError, 'cannot use "type" as custom qualifier in GenBankFeature\'s constructor'
            
        self._qualifiers = [(i,qualifiers[i]) for i in qualifiers]
        self._location = location
        self._type = type

    ####################################################################

    def parse(self, string):
        
        """ Updates feature information from information read in a
        GenBank-formatted string. """
        
        try:
            self._type = string.split()[0]
            locstring = string.split()[1]
            string = string.split()[2:]
        except IndexError:
            raise IOError, 'invalid GenBank feature string'
        while len(string) and string[0][0]!='/':
            locstring += string.pop(0)
        qualifiers=string
        
        # now we have the feature's position in a genuine string (locstring)
        self._location = GenBankFeatureLocation(locstring)
        
        # and the rest of qualifiers in a list (qualifiers)
        self._qualifiers = []
        if not len(qualifiers):
            return

        # fuses qualifiers where needed
        i=1
        while i<len(qualifiers):
            if qualifiers[i][0]!='/' or '=' not in qualifiers[i]:
                qualifiers[i-1]+=' '+qualifiers[i]
                del qualifiers[i]
            else:
                i+=1
        
        # processes the qualifiers
        for qualifier in qualifiers:
            if qualifier[0]!='/':
                raise IOError, 'invalid qualifier string: %s' %qualifier
            try:
                pos = qualifier.index('=')
            except ValueError:
                raise IOError, 'invalid qualifier string: %s' %qualifier
            key = qualifier[1:pos]
            value = qualifier[pos+1:]
            if key=='translation':
                value = ''.join(value.split())
            
            self._qualifiers.append((key,value))

    ####################################################################

    def get_sequence(self):
        
        """ Returns the string corresponding to this feature. If the
        positions pass beyond the end of the parent' sequence, a
        RuntimeError (instead of IndexError) is raised. """
        
        seq = ''
        for i,j in self._location:
            if j>len(self._parent):
                raise RuntimeError, 'GenBank feature exceeds sequence length'
            seq += self._parent._sequence[i:j+1]
        if self._location.isComplement(): return tools.rc(seq)
        else: return seq

    ####################################################################

    def start(self):
        
        """ Returns the first position of the (first) segment, such
        as ``start()`` is always smaller than ``stop()``. """
        
        return self._location[0][0]
        
    ####################################################################

    def stop(self):
        
        """ Returns the first position of the (first) segment, such
        as ``start()`` is always smaller than ``stop()``. """
        
        return self._location[-1][1]

    ####################################################################

    def copy(self, genbank):
        
        """
        Returns a copy of the current instance, connected to the
        :class:`~egglib.GenBank` instance *genbank*.
        """
        
        feature = GenBankFeature(genbank)
        feature.set(
            self._type, self._location.copy(), **dict(self._qualifiers))

        return feature
            
    ####################################################################

    def shift(self, shift):
            
        """ Shift all positions according to the (positive of
        negative) argument. """
        
        self._location.shift(shift)
        
    ####################################################################

    def __str__(self):
        
        """ GenBank formatting. """
        
        string = '     %s %s\n'%(self._type.ljust(15), str(self._location))
        string = self._parent._wrap(string, 21) + '\n'
        
        for key, value in self._qualifiers:
            if '/' in value:
                if value[0]+value[-1]!='""' and value[0]+value[-1]!='\'\'':
                    value = '"%s"' %value
            string+= self._parent._wrap('                     /%s=%s' %(key,value), 21)
            string+= '\n'

        return string
            
    ####################################################################
    
    def rc(self, length=None):
        
        """
        Reverse-complement the feature: apply it to the complement
        strand and reverse positions counting from the end. The *length*
        argument specifies the length of the complete sequence and is
        usually not required.
        """

        if length==None:
            length = len(self._parent)
        self._location.rc(length)


########################################################################

class GenBank(object):

    """
    :class:`~egglib.GenBank` represents a GenBank-formatted DNA sequence
    record.
    """
        
    ####################################################################
    
    def __init__(self, fname=None, string=None):
        
        """ Constructor signature: ``GenBank(fname=None, string=None)``.
        Only one of the two arguments *fname* and *string* can be
        non-``None``. If both are ``None``, the constructor generates
        an empty instance with sequence of length 0. If *fname* is
        non-``None``, a GenBank record is read from the file with this
        name. If ``string`` is non-``None``, a GenBank record is read
        directly from this string. The following variables are read from
        the parsed input if present: *accession*, *definition*, *title*,
        *version*, *GI*, *keywords*, *source*, *references* (which is a
        list), *locus* and *others*. Their default value is ``None``
        except for *references* and *others* for which it is an empty
        list. *source* is a (*description*, *species*, *taxonomy*)
        tuple. Each of *references* is a (*header*, *raw reference*)
        tuple and each of *others* is a (*key*, *raw*) tuple."""
        
        self._sequence = ''
        self._features = []
        self.accession = None
        self.definition = None
        self.title = None
        self.locus = None
        self.version = None
        self.GI = None
        self.source = (None,None,None)
        self.references = []
        self.keywords = None
        self.others = []
        
        if fname!=None and string!=None:
            raise ValueError, 'GenBank constructor expects at most one argument'
        
        stream = None
        
        if fname:
            stream = open(fname)
        if string:
            stream = StringIO.StringIO(string)
            
        if stream:
            self._parse(stream)

    ####################################################################

    def add_feature(self, feature):
        
        """
        Pushes a feature to the instance. The argument *feature* must
        be a well-formed :class:`~egglib.GenBankFeature` instance.
        """
        
        self._features.append(feature)
    
    ####################################################################

    def number_of_features(self):
        
        """
        Gives the number of features contained in the instance.
        """
        
        return len(self._features)

    ####################################################################
    
    def set_sequence(self, string):
        
        """ Sets the sequence string. Note that changing the record's
        string might obsolete the features. """
    
        self._sequence = string
        
    ####################################################################
    
    def get_sequence(self):
        
        """ Access to the sequence string. """
        
        return self._sequence
    
    ####################################################################

    def __iter__(self):
        
        """
        Iterates over the GenBankFeature instances contained in this
        GenBank instance.
        """
            
        for feature in self._features:
            yield feature

    ####################################################################

    def extract(self, from_pos, to_pos):
        
        """
        Returns a new GenBank instance representing a subset of the
        current instance, from position *from_pos*  to *to_pos*. All
        features that are completely included in the specified range are
        exported.
        """
        
        if from_pos < 0 or to_pos>=len(self):
            raise ValueError, 'invalid positions for GenBank extraction'
            
        gb = GenBank()
        gb.set_sequence( self._sequence[from_pos:to_pos+1] )
        for feature in self:
            if feature.start() >= from_pos and feature.stop() <= to_pos:
                    clone = feature.copy(gb)
                    clone.shift(-from_pos)
                    gb._features.append( clone )
        return gb

    ####################################################################
    
    def _parse(self, stream):
        
        # identifies blocks marked by a capitalized word at the very
        # beginning of the line, and send them to the dynamic parser
        
        key = None  # current block key
        buffer = [] # block being read
        
        while True:

            line = stream.readline()
            if not len(line): break

            if line.strip()=='ORIGIN':
                match = re.match('(.+)()', line.strip())
            else:
                match = re.match('^([A-Z]+) (.+)', line)
            if match:
                if key:
                    self._parse_block(key, buffer)
                key = match.group(1)
                buffer = [ match.group(2) ]
                if key=='ORIGIN':
                    break
            else:
                buffer.append(line)
                
        if key!='ORIGIN':
            raise IOError, 'GenBank records lacks sequence'
        
        self._parse_sequence(stream)
        
                
    ####################################################################
    
    def _parse_block(self, key, content):
        
        # alias "the dynamic parser". Processes the different GenBank
        # blocks appropriately
        
        def merge(items):
            return ' '.join([' '.join(i.split()) for i in items])
        
        if key=='LOCUS':
            if len(content)>1:
                raise IOError, 'GenBank record exhibits an invalid LOCUS line'
            self.locus = content[0].strip()
        
        elif key=='DEFINITION':
            self.definition = merge(content)

        elif key=='ACCESSION':
            self.accession = merge(content)

        elif key=='VERSION':
            content = merge(content)
            match = re.match('([^ ]+) +GI\:(\d+)', content)
            if not match:
                raise IOError, 'Incorrect VERSION/GI line in GenBank record'
            self.version = match.group(1)
            self.GI = match.group(2)

        elif key=='KEYWORDS':
            self.keywords = merge(content)
            
        elif key=='SOURCE':
            if len(content)<3:
                raise IOError, 'Incorrect SOURCE section in GenBank record'
            else:
                source= content[0].strip()
                organism= content[1].strip()
                if organism[:8] != 'ORGANISM':
                    raise IOError, 'Incorrect SOURCE section in GenBank record'
                organism= organism[10:]
                taxonomy= merge(content[2:])
            self.source = source, organism, taxonomy

        elif key=='REFERENCE':
            if len(content)<2:
                raise IOError, 'Incorrect REFERENCE section in GenBank record'
            self.references.append((content[0].strip(), ''.join(content[1:])))
            
        elif key=='FEATURES':
            self._parse_features(content[1:])
                        
        else:
            self.others.append(( key, content[0]+'\n'+''.join(content[1:]).rstrip() ))

    ####################################################################
    
    def _parse_features(self, lines):

        buffer = ''
        for line in lines:
            # new feature
            if not re.match('^ {8}', line):
                if len(buffer):
                    feature = GenBankFeature(self)
                    feature.parse(buffer)
                    self._features.append(feature)
                    buffer=''
            buffer+=line
                
        if len(buffer):
            feature = GenBankFeature(self)
            feature.parse(buffer)
            self._features.append(feature)

    ####################################################################

    def _parse_sequence(self, stream):
        
        self._sequence = ''
        
        while True:
            line= stream.readline()
            if not len(line):
                raise IOError, 'GenBank sequence doesn\'t terminated appropriately'
            line=line.strip()
            match = re.match(' *[0-9]+ ([A-Za-z ]+)', line)
            if match:
                fragment = ''.join( match.group(1).split())
                self._sequence += fragment.upper()
            elif line=='//':
                break
            else:
                raise IOError, 'invalid sequence line in GenBank record (reproduced below)\n%s' %line

    ####################################################################

    def __len__(self):
        return len(self._sequence)

    ####################################################################

    def write(self, fname):
        
        """
        Create a file named *fname* and writes the formatted record in.
        """
        
        f = open(fname, 'w')
        try:
            self.write_stream(f)
        finally:
            f.close()

    ####################################################################

    def __str__(self):

        stream = StringIO.StringIO()
        self.write_stream(stream)
        return stream.getvalue()

    ####################################################################

    def write_stream(self, stream):
        
        """ Writes the content of the instance as a Genbank-formatted
        string within the passed file (or file-compatible) stream. """
        
        if 'write' not in dir(stream):
            raise TypeError, 'GenBank write_stream formatter expects an argument with a "write" attribute'

        stream.write('LOCUS       %s\n' %self.locus)
        if self.definition!=None:
            stream.write(self._wrap('DEFINITION  %s' %self.definition, 12) + '\n')
        if self.accession!=None:
            stream.write('ACCESSION   %s\n' %self.accession)
        if self.version!=None or self.GI!=None:
            stream.write('VERSION     %s  GI:%s\n' %(self.version, self.GI))
        if self.keywords!=None:
            stream.write(self._wrap('KEYWORDS    %s' %self.keywords, 12)+ '\n')
        if (self.source[0]!=None or self.source[1]!=None or self.source[2]!=None):
            stream.write(self._wrap('SOURCE      %s' %str(self.source[0]), 12)+ '\n')
        if self.accession!=None:
            stream.write(self._wrap('  ORGANISM  %s' %str(self.source[1]), 12)+ '\n')
            stream.write(self._wrap('            %s' %str(self.source[2]), 12)+ '\n')
        
        for reference in self.references:
            stream.write('REFERENCE   %s\n%s' %reference)
        
        for other in self.others:
            # ignores "BASE COUNT"
            if other[0]=='BASE' and other[1].split()[0]=='COUNT':
                continue
            # otherwise writes
            stream.write('%s %s\n' %other)    

        stream.write('FEATURES             Location/Qualifiers\n')
        for feature in self:
            stream.write(str( feature ))

        stream.write('ORIGIN')
        
        if len(self) > 999999999:
            raise IOError, 'cannot export a GenBank instance with a sequence longest than 999999999 base pairs'

        c=0
        while c<len(self._sequence):
            stream.write('\n' +str(c+1).rjust(9))
            for i in range(6):
                stream.write(' '+self._sequence[c:c+10].lower())
                c+=10
                if c>=len(self._sequence):
                    break
        stream.write('\n//\n')

    ####################################################################

    LINEWIDTH = 80
    MAXBREAK = 40

    def _wrap(self, string, spacing):
        buffer = ''.join([i for i in string])
        space=0
        comma=0
        c=0
        res = ''
        while c<len(buffer):
            if buffer[c]==' ':
                space=c
            if buffer[c]==',':
                comma=c
            c+=1
            if c==(self.LINEWIDTH-1):
                if c<len(buffer) and buffer[c]==' ':
                    space = c
                if space>=self.MAXBREAK:
                    res+=buffer[:space]
                    res+= '\n'
                    buffer = ''.join([' ']*spacing) + buffer[space+1:]
                    c=0
                    space=0
                    comma=0
                elif comma>=self.MAXBREAK:
                    res+=buffer[:comma+1]
                    res+= '\n'
                    buffer = ''.join([' ']*spacing) + buffer[comma+1:]
                    c=0
                    space=0
                    comma=0
                else:
                    res+=buffer[:c]
                    res+= '\n'
                    buffer = ''.join([' ']*spacing) + buffer[c:]
                    c=0
                    comma=0
                    space=0
        res += buffer
        return res.rstrip()

    ####################################################################

    def rc(self):
        
        """
        Reverse-complement the instance (in place). All features
        positions and the sequence will be reverted and applied to the
        complementary strand. The features will be sorted in increasing
        start position (after reverting). This method should be applied
        only on genuine nucleotide sequences.
        """
        
        self._sequence= tools.rc(self._sequence)
        for feature in self: feature.rc(len(self))
        self._features.sort(lambda x,y: cmp(x.start(),y.start()))
        

########################################################################

class TreeNode(object):
    
    """
    This class provides an interface to a :class:`~egglib.Tree` instance's nodes
    and allows access and modification of data attached to a given node
    as well as the tree descending from that node. A node must be
    understood as the point *below* a branch. Edges (connections between
    nodes) have a direction: they go *from* a node *to* another node.
    Nodes have therefore *descendants*  and *ascendants*. Connecting a
    node to itself or making a two-way edge (to edges connecting the
    same two nodes in opposite directions) is not explicitly forbidden.
    Duplicate edges (between the same two nodes and in the same
    direction) are however illegal.
    """
    
    _cache = set()
    
    ####################################################################

    def __init__(self):
        
        """
        The constructor instantiates a tree node with default value.
        """
        
        self._label= None
        self._descendants = []
        self._ascendants = []

    ####################################################################

    def str(self, labels=True, brlens=True):
        
        """
        Formats the node and the subtree descending from is as a
        newick string. If *labels* is ``False``, omits the internal
        branch labels. If *brlens* is ``False``, omits the branch
        lengths. Doesn't support closed path.
        """
        
        if not len(TreeNode._cache):
            # this instance is top level
            TreeNode._cache.add(self)
            TopLevel = True
        else:
            # this instance is sub level
            if self in TreeNode._cache:
                raise ValueError, 'this tree has a closed path and cannot be formatted as newick'
            TreeNode._cache.add(self)
            TopLevel = False
        
        try:
            string = ''
            if len(self._descendants):
                string+='('
                for node,L in self._descendants:
                    string+=node.str(labels, brlens)
                    if (L!=None and brlens):
                        string+=':'+str(L)
                    string+=','
                string=string[:-1]+')'

            if (self._label!=None and (labels or not len(self._descendants))):
                string+=str(self._label)

        finally:
            if TopLevel:
                TreeNode._cache.clear()
        
        return string

    ####################################################################

    def leaves_down(self):
        
        """
        Recursively gets all leaf labels descending from that node. This
        method supports closed paths (networks) and nodes are never
        processed more than once.
        """
        
        if not len(TreeNode._cache):
            # this instance is top level
            TreeNode._cache.add(self)
            TopLevel = True
        else:
            # this instance is sub level
            if self in TreeNode._cache:
                return [] # don't process twice the same subtree
            TreeNode._cache.add(self)
            TopLevel = False

        try:

            leaves = []
            if len(self._descendants):
                for node,L in self._descendants:
                    leaves += node.leaves_down()
            else:
                leaves.append(self._label)

        finally:
            if TopLevel:
                TreeNode._cache.clear()

        return leaves

    ####################################################################
    
    def leaves_up(self):
        
        """
        Recursively gets all leaf labels contained on the other side of
        the tree. In case of a network, the results of this method and
        :meth:`leaves_down` might overlap. However, there can be no
        redundancy within the results returned by one of the methods
        """
        
        if not len(TreeNode._cache):
            # this instance is top level
            TreeNode._cache.add(self)
            TopLevel = True
        else:
            # this instance is sub level
            if self in TreeNode._cache:
                return [] # don't process twice the same subtree
            TreeNode._cache.add(self)
            TopLevel = False

        try:

            leaves = []
            for father,bidon in self._ascendants:
                for brother,bidon in father._descendants:
                    if brother != self:
                        leaves += brother.leaves_down()
                leaves += father.leaves_up()

        finally:
            if TopLevel:
                TreeNode._cache.clear()

        return leaves

    ####################################################################

    def unlink(self):
        
        """
        Clears all references to other :class:`~egglib.TreeNode` instances
        contained in this instance.
        """
        
        self._ascendants = []
        self._descendants = []

    ####################################################################

    def add_son(self, label=None, brlen=None):
        
        """
        Generates a new :class:`~egglib.TreeNode` instance descending
        from the current instance. *label* is to be applied to the new
        node. *brlen* is the length of the edge connecting the two nodes.
        Note that each node will refer to the other, generating a
        circular reference loop and preventing garbage collection of the
        node instances. It is therefore required to disconnect all nodes
        using the method :meth:`~data.TreeNode.disconnect`. Return the
        newly created node.
        """
        
        son = TreeNode()
        son._label = label
        son._ascendants.append((self, brlen))
        self._descendants.append((son, brlen))
            # note: no checks are required since the son is brand new

        return son

    ####################################################################

    def connect(self, node, brlen=None):
        
        """
        Connect this node to an other, existing, node. The orientation
        of the link is *from* the current instance *to* the passed
        instance. *brlen* is the length of the newly created edge. Note
        that each node will refer to the other, generating a circular
        reference loop and preventing garbage collection of the node
        instances. It is therefore required to disconnect all nodes
        using the method :meth:`~data.TreeNode.disconnect`.
        """

        if self in node._ascendants:
            raise ValueError, 'trying to duplicate an already existing edge'

        self._descendants.append((node, brlen))
        node._ascendants.append((self, brlen))
        
    ####################################################################

    def is_ascendant(self, node):
        
        """
        ``True`` if the :class:`~egglib.TreeNode` instance *node* is one of this
        node's ascendants.
        """
    
        for treeNode,L in self._ascendants:
            if node==treeNode: return True
        return False
    
    ####################################################################

    def is_descendant(self, node):
        
        """
        ``True`` if the :class:`~egglib.TreeNode` instance *node* is one of this
        node's descendants.
        """
    
        for treeNode,L in self._descendants:
            if node==treeNode: return True
        return False
    
    ####################################################################

    def ascendants(self):
        
        """
        Returns the list of all ascendants as :class:`~egglib.TreeNode`
        instances.
        """
        
        return [node for (node,L) in self._ascendants]

    ####################################################################

    def descendants(self):
        
        """
        Returns the list of all descendants as :class:`~egglib.TreeNode`
        instances.
        """
        
        return [node for (node,L) in self._descendants]

    ####################################################################
    
    def _get_index(self, list, node):
        
        index=-1
        for (i,(treeNode,L)) in enumerate(list):
            if treeNode==node:
                index=i
                break
        if index==-1: return None
        return index
    
    ####################################################################

    def remove_ascendant(self, node):
        
        """
        Removes the edges between this node and the node represented by
        the :class:`~egglib.TreeNode` instance *node*. *node* must be one of
        this node's ascendants. Note that this method removes also this
        node from *node*'s descendants.
        """

        index = self._get_index(self._ascendants, node)
        if index==None:
            raise ValueError, 'trying to remove non-existing edge'
        del self._ascendants[index]

        index = self._get_index(node._descendants, self)
        if index==None:
            raise RuntimeError, 'inconsistent edge'
        del node._descendants[index]
        
    ####################################################################
    
    def remove_descendant(self, node):
        
        """
        Removes the edges between this node and the node represented by
        the :class:`~egglib.TreeNode` instance *node*. *node* must be one of
        this node's descendants. Note that this method removes also this
        node from *node*'s ascendants.
        """

        index = self._get_index(self._descendants, node)
        if index==None:
            raise ValueError, 'trying to remove non-existing edge'
        del self._descendants[index]

        index = self._get_index(node._ascendants, self)
        if index==None:
            raise RuntimeError, 'inconsistent edge'

        del node._ascendants[index]
   
    ####################################################################

    def set_branch_to(self, node, value):
        
        """
        Sets the length of the branch connecting this node to *node*.
        *node* but be a :class:`~egglib.TreeNode` instance present amongst this
        node's descendants. *value* might be a float or ``None``. Note 
        that this methods affects both nodes.
        """

        index = self._get_index(self._descendants, node)
        if index==None:
            raise ValueError, 'trying to set the length of a non-existing edge'
        self._descendants[index] = (node, value)

        index = node._get_index(node._ascendants, self)
        if index==None:
            raise RuntimeError, 'inconsistent edge'
        node._ascendants[index] = (self, value)

    ####################################################################

    def set_branch_from(self, node, value):
        
        """
        Sets the length of the branch connecting *node* to this node.
        *node* but be a :class:`~egglib.TreeNode` instance present amongst this
        node's ascendants. *value* might be a float or ``None``.
        """

        index = self._get_index(self._ascendants, node)
        if index==None:
            raise ValueError, 'trying to set the length of a non-existing edge'
        self._ascendants[index] = (node, value)

        index = node._get_index(node._descendants, self)
        if index==None:
            raise RuntimeError, 'inconsistent edge'
        node._descendants[index] = (self, value)

    ####################################################################

    def branch_to(self, node):
        
        """
        Returns the length of the branch connecting this node to *node*.
        *node* but be a :class:`~egglib.TreeNode` instance present amongst this
        node's descendants. This method returns ``None`` if the value is
        not defined.
        """

        index = self._get_index(self._descendants, node)
        if index==None:
            raise ValueError, 'trying to get the length of a non-existing edge'
        return self._descendants[index][1]

    ####################################################################

    def branch_from(self, node):
        
        """
        Returns the length of the branch connecting *node* to this node.
        *node* but be a :class:`~egglib.TreeNode` instance present amongst this
        node's ascendants. This method returns ``None`` if the value is
        not defined.
        """

        index = self._get_index(self._ascendants, node)
        if index==None:
            raise ValueError, 'trying to get the length of a non-existing edge'
        return self._ascendants[index][1]

    ####################################################################

    def get_label(self):
        """ Returns the node's label. """        
        return self._label
        
    def set_label(self, value):
        """ Change the label *value*. """
        self._label = value

    def numberOfDescendants(self):
        """ Gets the number of nodes descending from this one. """
        return len(self._descendants)
    
    def numberOfAscendants(self):
        """ Gets the number of nodes descending to this one. """
        return len(self._ascendants)
    
    def numberOfRelatives(self):
        """ Gets the number of nodes connected to this one. """
        return len(self._descendants)+len(self._ascendants)
        
    ####################################################################
    
    def reverse(self, node, exchange_labels):
        
        """
        Reverse an edge's orientation between this node and the node
        given by *node*, as a :class:`~egglib.TreeNode` instance. The two nodes
        must be connected by exactly one edge. If *exchange_labels* is
        ``True``, the node labels are exchanged.
        """
        
        if self.is_descendant(node) and self.is_ascendant(node):
            raise ValueError, 'cannot reverse relationship between two TreeNode instances connected by more than one edge'
            
        if self.is_descendant(node):
            father = self
            son = node
        elif self.is_ascendant(node):
            father = node
            son = self
        else:
            raise ValueError, 'cannot reverse relationship between two TreeNode instances that are not connected by any edge'

        L = father.branch_to(son)
        father.remove_descendant(son)
        son.connect(father, L)
        
        if exchange_labels:
            label = father.get_label()
            father.set_label( son.get_label() )
            son.set_label( label )

    ####################################################################
    
    def sort(self):
        
        """
        Sorts the descendants based on their number of leaves.
        """

        def cmpnodes(node1, node2):
            a = len(node1[0].leaves_down())
            b = len(node2[0].leaves_down())
            return cmp(a,b)

        self._descendants.sort(cmpnodes)


########################################################################

class Tree(object):

    """
    Handles phylogenetic trees. A tree is a linked collection of nodes
    which all have at least one ascendant and any number of descendants.
    Nodes are implemented as :class:`~egglib.TreeNode` instances. A node without
    descendants is a *leaf*. A node with exactly one ascendant and one
    descendant is generally meaningless, but is allowed. All nodes
    (internal nodes as well as leaves) have a *label* which in the case
    of leaves can be used as leaf name. It is not possible to apply a
    name *and* a label to leaf node, accordingly to the newick format.
    All connections are oriented and lengthed (although the lengths can
    be omitted) but note that labels are applied to nodes, not edges
    (aka branches). All :class:`~egglib.Tree` instances have at least one root
    node which is the only one allowed not to have an ascendant. This
    class allows network-like structures, but note that some operations
    are available only for genuine trees (ie without closed paths).
    Import and export to/from strings and files are in the bracket-based
    newick format, and is concerned by this limitation. :class:`~egglib.Tree`
    instances can be exported using the built-in :func:`str` function,
    and the methods :meth:`.newick` and :meth:`.write`. :class:`~egglib.Tree`
    instances are iterable. Each step yields a :class:`~egglib.TreeNode`
    instance, starting with the root node but without a defined order.
    """

    ####################################################################
    
    def __del__(self):
        
        """
        To prevent bad surprises with circular references (networks)
        this destructor disconnects explicitely all TreeNode instances.
        """
        
        for node in self._nodes:
            node.unlink()

    ####################################################################
    
    def __init__(self, fname=None, string=None):
        
        """
        The instance can be initialized as an empty tree (with only a
        root node), or from a newick-formatted string. By default, the
        string is read from the file name *fname*, but it can be passed
        directly through the argument *string*. It is not allowed to
        set both *fname* and *string* at the same time. The
        newick parser expects a well-formed newick string (including
        the trailing semicolon).
        
        .. versionchanged:: 2.0.1
            Imports directly from a file. If a string is passed, it is
            interpreted as a file name by default.        
        """
        
        self._nodes = [ TreeNode() ]
        
        if fname!=None and string!=None:
            raise ValueError, 'cannot use arguments of Tree\'s constructor simulateneously'

        if fname:
            f= open(fname)
            string = f.read()
            f.close()

        if string:
            string= string.translate(None, ' \n\r\t')
            if (string[0]!='(' or string[-2:]!=');'):
                raise IOError, string+"""
----------------------
Invalid newick string!"""
            self._parse(string[:-1], self._nodes[0])

    ####################################################################

    def _parse(self, string, cur):

        """
        This parser expects a string starting with ( and ending with ).
        """

        if not len(string):
            raise IOError, 'Invalid newick string (empty string)'
        
        string= string[1:-1]
        i=0

        while(True):
            # case of a subtree
            subtree=None
            if (string[i]=='('):
                # we gather everything until the matching closing bracket
                subtree=string[i]
                acc=1
                while(acc):
                    i+=1
                    if (i==len(string)):
                        raise IOError, string+"""
----------------------
Invalid newick string!"""
                    subtree+=string[i]
                    if (string[i]=='('): acc+=1
                    if (string[i]==')'): acc-=1
                i+=1

            # gets all until the next comma or end of string
            buff=''
            while(i<len(string) and string[i]!=','):
                buff+= string[i]
                i+=1

            # converts to label and brlen
            if not len(buff):
                label= None
                brlen= None
            else:
                tbuff= buff.split(':')
                if (len(tbuff)==1):
                    label= buff
                    brlen= None
                elif (len(tbuff)==2):
                    label= tbuff[0]
                    brlen= float(tbuff[1])
                else:
                    raise IOError, string+"""
----------------------
Invalid newick string!
Cannot read label in """+buff

            # adds the node (whichever internal or terminal)
            self._nodes.append( cur.add_son(label=label,brlen=brlen) )

            # in case there is a subtree, it is recursively parsed
            if (subtree): self._parse(subtree, self._nodes[-1])

            # if we didn't reach the end of the string, we continue parsing
            if (len(string)==i):
                break
            else:
                i+=1
                continue

    ####################################################################

    def __str__(self):
        return self._nodes[0].str()+';'

    ####################################################################

    def newick(self, labels=True, brlens=True):
        
        """
        Returns the newick-formatted string representing the instance.
        If *labels* is ``False``, omits the internal branch labels. If
        *brlens* is ``False``, omits the branch lengths. Doesn't support
        closed paths.
        """
        
        return self._nodes[0].str(labels=labels, brlens=brlens)+';'

    ####################################################################

    def write(self, fname, labels=True, brlens=True):
        
        """
        Write the newick-formatted string representing the instance to
        a file named *fname*. If *labels* is ``False``, omits the
        internal branch labels. If *brlens* is ``False``, omits the
        branch lengths. Doesn't support closed paths.
        """
        
        f = open(fname, 'w')
        f.write( self._nodes[0].str(labels=labels, brlens=brlens)+';' )
        f.write('\n')
        f.close()

    ####################################################################
    
    def all_leaves(self):
        
        """
        Returns all leaves of the tree (nodes without descendant), as a
        list of :class:`~egglib.TreeNode` instances. If the tree is empty (that
        is, contains only a root node), this method returns an empty
        list.
        """

        return [node.get_label() for node in self._nodes[1:]
                                       if node.numberOfDescendants()==0]
                                       
    ####################################################################
    
    def get_terminal_nodes(self):
        
        """
        Returns the list of all :class:`~egglib.TreeNode` instances of this tree
        that don't have descendants. In case of an empty tree, an empty
        list is returned (ie the root is never returned).
        """

        return [node for node in self._nodes[1:] if node.numberOfDescendants()==0]

    ####################################################################
    
    def get_node(self, name):
        
        """
        Returns the first node of the tree bearing the given label. The
        returned object is a :class:`~egglib.TreeNode`. If no nodes of the tree
        match the passed name, ``None`` is returned. The order in which
        the nodes are examined is not defined. *name* can be of any
        type, including ``None`` (comparison is performed without
        conversion).
        """
        
        for node in self._nodes:
            if node.get_label()==name:
                return node

    ####################################################################

    def get_nodes(self, name):

        """
        Returns all nodes of the tree that bear the given label. The
        returned object is always a list of zero or more
        :class:`~egglib.TreeNode` instances.  The order in which nodes are
        sorted is not defined. *name* can be of any type, including
        ``None`` (comparison is performed without conversion).
        """

        return [node for node in self._nodes if node.get_label()==name]

    ####################################################################
    
    def get_node_re(self, regex):
        
        """
        Returns the first node of the tree matching the regular
        expression *regex*. *regex* should be a valid regular expression
        (refer to the documentation of the :mod:`re` module of the 
        standard Python library). If no nodes of the tree match the
        regular expression, ``None`` is returned. The order in which
        the nodes are examined is not defined.
        """
        
        for node in self._nodes:
            label=node.get_label()
            if label!=None and re.search(regex,label):
                return node

    ####################################################################

    def get_nodes_re(self, regex):

        """
        Returns all nodes of the tree matching the regular expression
        *regex*. *regex* should be a valid regular expression (refer to
        the documentation of the :mod:`re` module of the standard Python
        library). The returned object is always a list of zero or more
        :class:`~egglib.TreeNode` instances.  The order in which nodes are
        sorted is not defined.
        """

        return [node for node in self._nodes
                       if node.get_label()!=None and re.search(regex,node.get_label())]

    ####################################################################

    def root_node(self):
        
        """
        Returns the root node as a :class:`~egglib.TreeNode` instance.
        """
        
        return self._nodes[0]

    ####################################################################

    def last_node(self):
        
        """
        Returns the last loaded node (as a :class:`~egglib.TreeNode` instance).
        If no nodes were loaded, the root is returned.
        """
        
        return self._nodes[-1]

    ####################################################################

    def add_node(self, parent, label=None, brlen=None):
        
        """
        Adds a node to the tree. *parent* must be a :class:`~egglib.TreeNode`
        instance already present in the instance; *label* is the label
        to apply to the tree (or the taxon name if the node is intended
        to be terminal); *brlen* the length of the edge connecting
        *parent* to the new node. Their is no formal difference between
        introducing a new internal node and a terminal node (or leaf).
        The new node has initially no descendant and is therefore a
        leaf until it is itself connected to a new node. The newly
        created node can be accessed through :meth:`.last_node`.
        """
        
        self._nodes.append( parent.add_son(label, brlen) )

    ####################################################################

    def __iter__(self):
        
        for node in self._nodes:
            yield node
            
    ####################################################################

    def number_of_nodes(self):
        
        """
        Gives the number of nodes of the tree (including leaves and root).
        Returns 1 if the tree contains the root only.
        """
        
        return len(self._nodes)
        
    ####################################################################

    def number_of_leaves(self):
        
        """
        Gives the number of leaves (terminal nodes) of the tree. Returns
        0 if the tree contains the root only.
        """
        
        c = 0
        for node in self._nodes[1:]:
            if node.numberOfDescendants()==0:
                c+=1
        return c
        
    ####################################################################
    
    def total_length(self):
        """
        Returns the sum of all branch lengths across all nodes. All
        branch lengths must be defined, otherwise a :class:`ValueError`
        will be raised.
        """
        
        L = 0
        for node in self._nodes:
            for descendant in node.descendants():
                l = node.branch_to(descendant)
                if l==None:
                    raise ValueError, 'cannot compute tree\'s length: at least one edge has no defined length'
                L += l
        return L
        
    ####################################################################

    def number_of_nodes(self):
        
        """
        Gives the number of nodes of the tree (including leaves and root).
        """
        
        return len(self._nodes)

    ####################################################################

    def findMonophyleticGroup(self, taxa):
        
        """
        Checks whether a group is one of the monophyletic groups defined
        by the tree. If so, returns the first such node found as a
        :class:`~egglib.TreeNode` instance. Returns ``None`` if no such group
        is found. *taxa* is an iterable of leaf label strings. It is not
        required that all labels are unique. This method returns the
        first node encountered whose list of descending leaves matches
        exactly the passed list. This method assumes that the tree is
        rooted, ie the orientation of branches is irrelevant: for a tree
        represented by ``((A,B),(C,(D,E)),((F,G),(H,I))))``, the call
        ``findMonophyleticGroup(['A','B','C','D','E'])`` will not
        succeed (technically because the group is overlapping the root).
        If the group is searched regardless of the orientation of the
        tree, typically for unrooted trees, consider using :meth:`.findGroup`
        instead. The order of leaves is irrelevant. This method returns
        the root node if *taxa* is the list of all tree's leaves. If
        *taxa* contains a single label matching a leaf of this tree,
        then the result will be the same as with :meth:`.get_node`.
        """
        
        check_list = self.all_leaves()
        for taxon in taxa:
            if taxon not in check_list:
                raise ValueError, 'cannot find monophyletic group: %s is not leaf list' %str(taxon)
        query = sorted(taxa)

        for node in self._nodes:
            if query == sorted(node.leaves_down()):
                return node

        return None

    ####################################################################

    def findGroup(self, taxa):
        
        """
        Checks whether a group is one of the groups defined by the tree,
        regardless of the orientation of the tree. If so, returns the
        first node found as a :class:`~egglib.TreeNode` instance. Returns
        ``None`` if no such group is found. *taxa* is an iterable of
        leaf label strings. It is not required that all  labels are
        unique. This method returns the first node encountered whose
        list of descending leaves matches exactly the list *taxa* or
        the whose list of ascending leaves (that is all leaves of the
        tree that are not among the descending leaves) matches exactly
        the list *taxa*. This method disregards the tree orientation;
        for a tree represented by ``((A,B),(C,(D,E)),((F,G),(H,I))))``,
        the call ``findGroup(['A','B','C','D','E'])`` will succeed and
        return the node placed at the root of ``((F,G),(H,I))``. If a
        monophyletic group must be explicitely searched for, consider
        using :meth:`.findMonophyleticGroup` instead. The order of
        leaves is irrelevant. This method returns  the root node if
        *taxa* is the list of all tree's leaves. If *taxa* contains a
        single label matching a leaf of this tree, then the result will
        be the same as with :meth:`.get_node`.
        """
        
        check_list = self.all_leaves()
        for taxon in taxa:
            if taxon not in check_list:
                raise ValueError, 'cannot find group: %s is not leaf list' %str(taxon)
        query = sorted(taxa)

        for node in self._nodes:
            if (query == sorted(node.leaves_down()) or
                query == sorted(node.leaves_up())):
                    return node

        return None
        
    ####################################################################

    def smallest_monophyleticGroup(self, taxa, threshold=None, minimum=1):
        
        """
        Returns the most recent common ancestor of a set of leaves, as a
        :class:`~egglib.TreeNode` instance. The node returned corresponds to the
        smallest clade fulfilling the criteria. *taxa* must be a list of
        leaf labels. All labels must be found within the clade,
        including duplicates whenever appropriate. *threshold* is the
        minimum numerical label the node must exhibit to be returned.
        If threshold is ``None``, this criterion is not applied.
        Otherwise, nodes that have a label not convertible to float or
        whose label is inferior than *threshold* are not returned.
        *minimum* is the smallest number of descending leaf a clade must
        have to be returned. The root is never returned. Returns ``None``
        if no valid node can be found.
        
        .. versionchanged:: 2.0.1
            The root is never returned; duplicates are supported; the
            *minimum* argument is not checked; and nodes that don't have
            a numeric label are supported when *threshold* is not
            ``None`` (but they are excluded).
        """

        check_list = self.all_leaves()
        for taxon in taxa:
            if taxon not in check_list:
                raise ValueError, 'cannot find monophyletic group: %s is not leaf list' %str(taxon)

        query = sorted(taxa)
        smallest = None
        smallestn= 0

        for node in self._nodes[1:]:   # doesn't consider the root

            # if threshold is not None, label must be >= threshold
            if threshold!=None:
                try:
                    label = float(node.get_label())
                except TypeError:
                    continue
                if label<threshold:
                    continue

            descendants = node.leaves_down()
            n = len(descendants)
            
            # check that the node is at least smallest than the currently stored
            if n>=smallestn and smallestn!=0:
                continue

            # there must be enough descendants
            if n<minimum:
                continue

            # check that all of query (incl. duplicates) are in descendants
            flag = True
            for item in query:
                try:
                    descendants.remove(item)
                except ValueError:
                    flag = False
                    break
            
            if not flag:
                continue
                
            # if reached herem the node satisfies all conditions
            smallest = node
            smallestn= n
                    
        return smallest

    ####################################################################

    def smallest_group(self, taxa, threshold=None, minimum=1):
        
        """
        Returns the smallest clade containing a set of leaves, as a
        :class:`~egglib.TreeNode` instance, without regard to the orientation of
        the tree. The node returned corresponds to the smallest clade
        fulfilling the criteria. *taxa* must be a list of leaf labels.
        All labels must be found either within the clade, or in the rest
        of the tree (all tree leaves not in this clade). Duplicates are
        included whenever appropriate. *threshold* is the minimum
        numerical label the node must exhibit to be returned. If
        threshold is ``None``, this criterion is not applied. Otherwise,
        nodes that have a label not convertible to float or whose label
        is inferior than *threshold* are not returned. *minimum* is the
        smallest number of descending leaf a clade must have to be
        returned. The root is never returned. Returns ``None`` if no
        valid node can be found.
        
        .. versionchanged:: 2.0.1
            The root is never returned; duplicates are supported; the
            *minimum* argument is not checked; and nodes that don't have
            a numeric label are supported when *threshold* is not
            ``None`` (but they are excluded).
        """

        check_list = self.all_leaves()
        for taxon in taxa:
            if taxon not in check_list:
                raise ValueError, 'cannot find group: %s is not leaf list' %str(taxon)
        query = sorted(taxa)
            
        smallest = None
        smallestn= 0

        for node in self._nodes[1:]:   # doesn't consider the root

            # if threshold is not None, label must be >= threshold
            if threshold!=None:
                try:
                    label = float(node.get_label())
                except TypeError:
                    continue
                if label<threshold:
                    continue

            descendantsD = node.leaves_down()
            descendantsU = node.leaves_up()

            # checks descendants DOWN

            nD = len(descendantsD)
            flagD1 = (nD<smallestn or smallestn==0) and nD>=minimum
            nU = len(descendantsU)
            flagU1 = (nU<smallestn or smallestn==0) and nU>=minimum

            # check that all of query (incl. duplicates) are in descendants
            flagD2 = True
            for item in query:
                try:
                    descendantsD.remove(item)
                except ValueError:
                    flagD2 = False
                    break
            
            # the same for the complement
            flagU2 = True
            for item in query:
                try:
                    descendantsU.remove(item)
                except ValueError:
                    flagU2 = False
                    break

            if flagD1 and flagD2:
                smallest = node
                smallestn= nD
                
            if flagU1 and flagU2:
                smallest = node
                smallestn= nU
                    
        return smallest

    ####################################################################

    def collapse(self, node):
        
        """
        Collapses a branch. *node* must be one of the nodes contained in
        the tree (as a :class:`~egglib.TreeNode` instance). It must have a
        unique ascendant. If not, a ValueError is raised. Obviously, the
        tree's root cannot be collapsed. The destruction of the node
        might discard the information of its label. This information
        will be transfered to the ascending node. The ascending node's
        label will carry the information of either one whichever is not
        ``None``. If the two labels are identical, nothing will be done.
        If the two labels are different and different from ``None``,
        they will be concatenated (ascending node first), separated by a
        semicolon like in ``oldlabel;newlabel``. The length of the
        removed edge will be spread equally among all its descendants
        (see example below).
        
        Collapsing node ``[4]`` on the following tree::
        
             /------------------------------------------->[1]
             |
             |             /----------------------------->[3]
             |             |
             |----------->[2]             /-------------->[5]
             |             |              |
             |             \------------>[4]
            [0]                           |
             |                            \-------------->[6]
             |
             |              /---------------------------->[8]
             |              | 
             \------------>[7]            /------------->[10]
                            |             |
                            \----------->[9]
                                          |
                                          \------------->[11]
                                          
        will generate the following tree, with the correction of edge
        lengths as  depicted::

             /------------------------------------------->[1]
             |
             |             /----------------------------->[3]
             |             |
             |----------->[2]
             |             |
             |             |-------------------->[5]        L5 = L5+L4/2
            [0]            |
             |             \-------------------->[6]        L6 = L6+L4/2
             |
             |              /---------------------------->[8]
             |              | 
             \------------>[7]            /------------->[10]
                            |             |
                            \----------->[9]
                                          |
                                          \------------->[11]
        
        Although the total edge length of the tree is not modified, the
        relationships will be altered: the distance between the
        descendants of the collapsed node (nodes 5 and 6 in the example
        above) will be artificially increased.
        """

        # sanity checks
        
        if node not in self._nodes[1:]:
            raise ValueError, 'cannot collapse branch: passed node is not one of this tree\'s nodes or is the root'

        if node.numberOfAscendants()!=1:
            raise ValueError, 'cannot collapse branch: passed node have no or two many ascendants'

        if node.numberOfDescendants()<2:
            raise ValueError, 'cannot collapse branch: passed node doesn\'t have enough descendants'

        for son in node.descendants():
            if son.numberOfAscendants()!=1:
                raise ValueError, 'cannot collapse branch: the collapsed node\'s descendants must all have exactly one ascendant'

        if node.branch_from( node.ascendants()[0] )!=None:
            for son in node.descendants():
                if node.branch_to(son)==None:
                    raise ValueError, 'cannot collapse node: some edges have lengths, but not all!'

        # stores variables
        # the figule below indicates how variables are defined
        #
        #            (L1)
        #          /------S1
        #          |
        #    (L0)  |
        #  /------node
        #  |       |
        #  |       |
        #  F       \------S2
        #  |         (L2)
        #  |
        #  \------B
        #    (L3)
        #
        #
        # node:   the node to collapse
        # F:      its ascendant
        # L0:     the length of the ascending edge
        # S1, S2: its two descendants
        # L1, L2: the lengths of the two descending edges.
        # B:      the other descendant of F
        # L3:     the length of the other edge descending from B
        #
        
        # collection of variables
        F     = node.ascendants()[0]
        L0    = node.branch_from(F)
        sons  = node.descendants()
        lens  = map(node.branch_to, sons)
        
        # connects the sons their grandfather
        for S,L in zip(sons,lens):
            if L0!=None:
                # increases all edge lengths equally
                L += L0/len(sons)
            F.connect(S, L)
        
        # disconnects and remove the node to be collapsed
        node.remove_ascendant(F)
        for son in sons:
            node.remove_descendant(son)
        self._nodes.remove(node)
        
        # saves the label information
        Lf = F.get_label()
        Ln = node.get_label()
        if Lf!=None and Ln!=Lf:
            if Ln==None:
                node.set_label(Ln)
            else:
                node.set_label('%s;%s' %(Lf,Ln))
        F.set_label(None)

    ####################################################################

    def root(self, outgroup, branchsplit=0.5):
        
        """
        Roots the tree. *outgroup* must be a :class:`~egglib.TreeNode` instance
        contained in this tree. The root will be placed somewhere on the
        branch that leads to this node (between the current root and the
        node). If this branch doesn't have a branch length, the
        *branchsplit* argument is ignored. Otherwise, *branchsplit* must
        be a real number between 0. and 1. and gives the proportion of
        the branch that must be allocated to the basal branch leading to
        the outgroup, the complement being allocated to the branch
        leading to the rest of the tree. It is illegal to call this
        method ons tree that are already rooted (have a trifurcation at
        the root) or on trees that that have a closed path between the
        current root (or base of the tree) and the intended root. A
        :class:`ValueError` is raised whenever the tree cannot be
        rooted.
        
        If the original tree has this structure::
        
             /------------------------------------------->[1]
             |
             |             /----------------------------->[3]
             |             |
             |----------->[2]             /-------------->[5]
             |             |              |
             |             \------------>[4]
            [0]                           |
             |                            \-------------->[6]
             |
             |              /---------------------------->[8]
             |              | 
             \------------>[7]            /------------->[10]
                            |             |
                            \----------->[9]
                                          |
                                          \------------->[11]
                                          
        And rooting is requested at node ``[9]``, the root will be
        placed on the edge marked by ``[ROOT]`` below::
                         
             /------------------------------------------->[1]
             |
             |             /----------------------------->[3]
             |             |
             |----------->[2]             /-------------->[5]
             |             |              |
             |             \------------>[4]
            [0]                           |
             |                            \-------------->[6]
             |
             |              /---------------------------->[8]
             |              | 
             \------------>[7]            /------------->[10]
                            |             |
                            \---[ROOT]-->[9]
                                          |
                                          \------------->[11]

        And the outcome will be as depicted below, with the introduction
        of a new node (which would be ``[12]`` here) at the root::

                                /-------------------------[1]
                                |
                         /-----[0]     /------------------[3]
                         |      |      |
                         |      \-----[2]      /----------[5]
                         |             |       |
             /----------[7]            \------[4]
             |           |                     |
             |           |                     \----------[6]
             |           |
           [ROOT]        \--------------------------------[8]
             |
             |                     /---------------------[10]
             |                     |
             \--------[E1]--------[9]
                                   |
                                   \---------------------[11]
            
        If the length of the branch on which the root is placed (*L0*)
        is not ``None``, the length of the edge ``[E1]`` will be
        *branshsplit* * *L0* and the length of ``[E2]`` will be
        (1 - *branshsplit* ) * *L0* .
        
        Note that the label of outgroup node is copied to the node at
        the other side of the root where the root is placed. The
        rationale is that information attached to the root edge might
        have to be applied to both basal edges. In case the original
        root (or basal node) had a specified label, it will be retained
        and the outgroup label will not be copied.
        """
        
        # sanity checking
        
        if outgroup not in self._nodes[1:]:
            raise ValueError, 'cannot root tree: invalid outgroup'
            
        if not isinstance(branchsplit, (float,int)) or branchsplit<0. or branchsplit>1.:
            raise ValueError, 'cannot root tree: invalid branchsplit argument'
        
        if self._nodes[0].numberOfDescendants()!=3:
            raise ValueError, 'cannot root tree: the tree is already rooted; use collapse() first'
            
        if self._nodes[0].numberOfAscendants()!=0:
            raise ValueError, 'cannot root tree: the root has an ascendant (which is illegal)'

        # collects all nodes on the path from the new to the old root (both included)
        path = [outgroup]
        while True:
            if path[-1].numberOfAscendants()!=1:
                raise ValueError, 'cannot root tree: the involved subtree has network features'
            curr = path[-1].ascendants()[0]
            if curr in path:
                raise ValueError, 'cannot root tree: the involved subtree has network features'
            path.append(curr)
            if curr==self._nodes[0]:
                break
            # path[1] is the node at the other side of the root edge

        # makes new list of nodes with new root
        self._nodes = [TreeNode()] + self._nodes
        
        # determines the lengths of basal branches
        L = outgroup.branch_from(path[1])
        if L != None:
            L1 = L * branchsplit
            L2 = L * (1-branchsplit)
        else:
            L1 = None
            L2 = None
            
        # removes the branch where the root is placed
        outgroup.remove_ascendant( path[1] )
    
        # connects the new root to the two subtrees
        self._nodes[0].connect( outgroup, L1 )
        self._nodes[0].connect( path[1],  L2 )
    
        # reverses the ascendant/descendant relationships on the path to the root
        for i in range(1,len(path)-1):

            # takes the nodes from the old root to the new to preserve all labels
            path[-i-1].reverse(path[-i], True)

        # apply the outgroup label to the new edge (the erased label should be the root's)
        if path[1].get_label()==None:
            if outgroup.numberOfDescendants() > 0:
                path[1].set_label( outgroup.get_label() )

    ####################################################################

    def reoriente(self, new_root):
        
        """
        Moves the root location of the tree. This method is solely
        intended to alter the representation of unrooted trees
        (trees that have a trifurcation at the root). *new_root* must be
        a :class:`~egglib.TreeNode` instance contained in this tree and
        representing the position of the new root. It might be the
        current root.It is illegal to call this method on trees that
        have a closed path between the current root and the new root. A
        :class:`ValueError` is raised whenever the tree cannot be
        reoriented.
        
        If the original tree has this structure::
        
             /------------------------------------------>[1]
             |
             |             /---------------------------->[3]
             |             |
             |----------->[2]             /------------->[5]
             |             |              |
             |             \------------>[4]
            [0]                           |
             |                            \------------->[6]
             |
             |              /--------------------------->[8]
             |              | 
             \------------>[7]            /------------>[10]
                            |             |
                            \----------->[9]
                                          |
                                          \------------>[11]
                                          
        And rooting is requested at node ``[7]``, the outcome will be as
        depicted below, the edge lengths being ignored::

                         /--------------------------------[1]
                         |
             /----------[0]         /---------------------[3]
             |           |          |
             |           \---------[2]        /-----------[5]
             |                      |         |       
             |                      \--------[4]
             |                                |
            [7]                               \-----------[6]
             |
             |--------------------------------------------[8]
             |
             |                     /---------------------[10]
             |                     |
             \--------------------[9]
                                   |
                                   \---------------------[11]
            
        If the length of the branch on which the root is placed (*L0*)
        is not ``None``, the length of the edge ``[E1]`` will be
        *branshsplit* * *L0* and the length of ``[E2]`` will be
        (1 - *branshsplit* ) * *L0* .
        
        Note that the label of outgroup node is copied to the node at
        the other side of the root where the root is placed. The
        rationale is that information attached to the root edge might
        have to be applied to both basal edges. In case the original
        root (or basal node) had a specified label, it will be retained
        and the outgroup label will not be copied.
        """
        
        # sanity checking
        
        if new_root not in self._nodes[1:]:
            raise ValueError, 'cannot reoriente tree: invalid new root'
            
        if self._nodes[0].numberOfDescendants()!=3:
            raise ValueError, 'cannot reoriente tree: the tree is rooted; use collapse() first'
            
        if self._nodes[0].numberOfAscendants()!=0:
            raise ValueError, 'cannot reoriente tree: the current root has an ascendant (which is illegal)'

        if new_root.numberOfDescendants()<2:
            raise ValueError, 'cannot reoriente tree: the new root must have at least two descendants'

        # collects all nodes on the path from the new to the old root (both included)
        path = [new_root]
        while True:
            if path[-1].numberOfAscendants()!=1:
                raise ValueError, 'cannot reoriente tree: the involved subtree has network features'
            node = path[-1].ascendants()[0]
            if node in path:
                raise ValueError, 'cannot reoriente tree: the involved subtree has network features'
            path.append(node)
            if node==self._nodes[0]:
                break

        # makes new list of nodes with new root
        self._nodes = [new_root] + [node for node in self._nodes if node!=new_root]
        
        # reverses the ascendant/descendant relationships on the path to the root
        for i in range(len(path)-1):
            path[-i-2].reverse(path[-i-1], True)

    ####################################################################

    def frequency_nodes(self, trees, relative=False):
        
        """
        Labels all nodes of the current instances by integers counting
        the number of trees where the same node exists among the trees
        in the iterable *trees*. Each item must be a :class:`~egglib.Tree`
        instance defining exactly the same set of leaf lables. In case
        *relative* is ``True``, the numbers are expressed as fractions.
        The label is converted to a string in both cases.
        """
        
        # list of leaves for checking
        my_leaves = sorted(self.all_leaves())

        # sanity check
        if len(self._nodes)<1:
            raise ValueError, 'cannot compute statistics on a tree without nodes'
        
        # initializes a cache where to get the data
        cache = []
        for node in self._nodes[1:]:
            if len(node.descendants())>1:
                cache.append([node, sorted(node.leaves_down()), 0 ])
                        
        # processes the trees
        for tree in trees:
            
            # sanity check
            try:
                if my_leaves != sorted(tree.all_leaves()):
                    raise ValueError, 'an invalid item was passed to Tree.frequency_node'
            except AttributeError:
                raise ValueError, 'an invalid item was passed to Tree.frequency_node'

            tree_clades = (
                [ sorted(node.leaves_up()) for node in tree._nodes[1:]
                if node.descendants()>1]
             +  [ sorted(node.leaves_down()) for node in tree._nodes[1:]
                if node.descendants()>1])

            # checks all nodes
            for i in range(len(cache)):
                # counts the matching nodes
                if  cache[i][1] in tree_clades :
                    cache[i][2]+=1

        # annotates the nodes
        for node, leaves, count in cache:
            
            if relative:
                # make it relative
                if not len(trees):
                    raise ValueError, 'cannot express relative node support: no trees were passed'
                count = 1. * count / len(trees)
            
            node.set_label(str(count))
                        
    ####################################################################

    def copy(self):
        
        """
        Returns a deep copy of self.
        """

        clone = Tree()
        
        if self._nodes[0].numberOfAscendants()!=0:
            raise ValueError, 'this tree\'s root has ascendants: this is illegal and precludes copying'
        
        clone._nodes = [TreeNode() for i in self._nodes]
        
        for i, father in enumerate( self._nodes ):

            clone._nodes[i].set_label( father.get_label() )

            for son in father.descendants():
                L = father.branch_to(son)
                try:
                    j = self._nodes.index(son)
                except ValueError:
                    raise ValueError, 'this tree is corrupted: an edge point outside the tree'
                clone._nodes[i].connect( clone._nodes[j], L )

        return clone

    ####################################################################

    def clean_internal_labels(self):
        
        """
        Ensures that all internal labels are not different than None.
        """
        
        for node in self._nodes:
            if node.numberOfDescendants() > 0:
                node.set_label(None)

    ####################################################################

    def clean_edge_lengths(self):
        
        """
        Ensures that all edge lengths are not different than None.
        """
        
        if self._nodes[0].numberOfAscendants()!=0:
            raise ValueError, 'this tree\'s root has ascendants: this is illegal and precludes copying'

        for node in self._nodes:
            for son in node.descendants():
                node.set_branch_to( son, None )
    
    ####################################################################

    def remove_node(self, node):
        
        """
        Removes the node from the tree and all its descendants. Any
        node can be removed, including nodes without descendants,
        provided that the root is not among the nodes removed. The node
        in question must have only one ascendant. In case its ascendant
        had previously only two descendants and only one ascendant, it
        will be automatically removed.
        """

        # sanity check
        if node not in self._nodes[1:]:
            raise ValueError, 'cannot remove tree\'s node: invalid node'

        # makes the list of nodes to be removed
        def gather(node, death_list):
            if node.numberOfAscendants()>1:
                raise ValueError, 'cannot remove tree\'s node: network feature detected'
            death_list.add(node)
            for descendant in node.descendants():
                if descendant in death_list:
                    raise ValueError, 'cannot remove tree\'s node: network feature detected'
                gather(descendant, death_list)
        death_list = set()
        gather(node, death_list)
        if self._nodes[0] in death_list:
            raise ValueError, 'cannot remove tree\'s node: network feature detected'
        
        # disconnects the node from the rest of the tree
        father = node.ascendants()[0]
        node.remove_ascendant(father)
        
        # cleans the subtree
        for victim in death_list:
            victim.unlink()
            self._nodes.remove(victim)
        del death_list
        
        # if needed collapse the ascendant
        if father.numberOfDescendants()==1 and father.numberOfAscendants()==1:
            
            # gets the relatives
            grandfather = father.ascendants()[0]
            son = father.descendants()[0]
            
            # determines the branch length
            branch_up = father.branch_from(grandfather)
            branch_down = father.branch_to(son)
            if branch_up!=None and branch_down!=None:
                L = branch_up+branch_down
            elif branch_up!=None:
                L = branch_up
            elif branch_down!=None:
                L = branch_down
            else:
                L = None
            
            # disconnect the father
            father.remove_ascendant(grandfather)
            father.remove_descendant(son)
            
            # remove the father
            self._nodes.remove(father)
            
            # reconnects the link
            grandfather.connect(son, L)
                        
    ####################################################################

    def midroot(self):

        """
        Automatic rooting of the tree using the midpoint method. The
        tree must not be previously rooted, there must be not closed
        path or network-like structures, and must edges must have an
        available length value.
        """

        # sanity check
        if self._nodes[0].numberOfDescendants()<3:
            raise ValueError, 'cannot perform automatic rooting: the current root of the three must have at least three descendants'
        
        for node in self._nodes[1:]:
            if node.numberOfAscendants()!=1:
                raise ValueError, 'cannot perform automatic rooting: the tree has network features'
    
        # collects leaf nodes
        leaves  = [node for node in self._nodes[1:] if node.numberOfDescendants()==0]

        # new sanity check
        if len(leaves)<3:
            raise ValueError, 'cannot perform automatic rooting: the tree has less than three leaves'

        # collects the ancestral path (as ancestor nodes  and edge lengths

        ancestors = {}  # list of ancestors, from itself to root (excluded)
        
        for leaf in leaves:
            ancestors[leaf]=[]

            # walks the path up excluded the root
            current = leaf
            while current!=self._nodes[0]:
                ancestors[leaf].append(current)
                father = current.ascendants()[0]
                if father in ancestors[leaf]:
                    raise ValueError, 'cannot perform automatic rooting: the tree has network features (closed path)'
                current = father

        # finds the most distant pair
        longest_path = None, 0   # longest_sub_path, total_distance

        n = len(leaves)
        for i in range(0, n-1):
            for j in range(i+1, n):

                leaf1 = leaves[i]
                leaf2 = leaves[j]

                # determines the two paths
                path1= [k for k in ancestors[leaf1] if k not in ancestors[leaf2]]
                path2= [k for k in ancestors[leaf2] if k not in ancestors[leaf1]]

                # computes the total distance
                D1 = 0
                for node in path1:
                    father = node.ascendants()[0]
                    d = node.branch_from( father )
                    if d==None:
                        raise ValueError, 'cannot perform automatic rooting: all tree edges must have a valid length'
                    D1 += d
                D2 = 0
                for node in path2:
                    father = node.ascendants()[0]
                    d = node.branch_from( father )
                    if d==None:
                        raise ValueError, 'cannot perform automatic rooting: all tree edges must have a valid length'
                    D2 += d
                D = D1+D2
                
                # stores the paths if longest than current
                if D > longest_path[1]:
                    if D1>=D2:
                        longest_path = path1, D
                    else:
                        longest_path = path2, D

        # finds the midpoint
        path, D = longest_path
        limit = D/2.

        root=None
        acc=0
        for node in path:
            father = node.ascendants()[0]
            L = node.branch_from( father )
            acc += L
            
            if acc > limit:
                # midpoint found
                root = node
                offset = (L - (acc - limit)) / L
                break
        
        # sanity check
        if not root:
            raise RuntimeError, 'midpoint rooting: unexpected case (no root found)'

        # root itself
        self.root(root, offset)

    ###################################################################
    
    def lateralize(self):
        
        """
        At each node of the tree, sorts the descendants based on the
        number of leaves that descend from them. The result is a tree
        where the richest branches are pushed to the back.
        """

        for node in self._nodes:
            node.sort()



########################################################################
