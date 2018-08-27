"""
Approximate Bayesian Computation components. Part of the underlying C++
utilities are (currently) without Python wrapper, so should be used
directly throught the binding. Please refer to the C++ library
documentation for the following class:

    * :class:`~egglib.egglib_binding.ABC`: the class to perform the
      rejection-regression operations.
"""

__license__ = """
    Copyright 2009-2012 Stephane De Mita, Mathieu Siol

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

import egglib_binding, re, data, tools, simul, os, sys


# UTILITY CLASSES ######################################################

class ParamSample(object):
    
    """
    Holds a list of float values (parameters) of fixed length.
    Supports :meth:`len`, use of the subscript (``[]``) operator for
    accessing or modifying values (but not deleting) and :meth:`str`
    (``str(paramSample)`` or with :meth:`print`). String formatting
    returns a space-delimited string of the values.
    """
    
    ####################################################################
    
    def __init__(self, length):

        """
        The constructor expects a *length* argument to fix the number
        of parameters.
        """
        
        self._values = [0] * length

    ####################################################################

    def __len__(self):

        """ Number of parameters (as in ``len(paramSample)``) """

        return len(self._values)

    ####################################################################

    def values(self):
        
        """ Gets a deep copy of the values contained in the instance """
        
        return list(self._values)
        
    ####################################################################

    def __getitem__(self, index):
        
        return self._values[index]
        
    ####################################################################

    def __setitem__(self, index, value):
        
        self._values[index] = value

    ####################################################################

    def __str__(self):
        
        return ' '.join(map(str, self._values))

########################################################################

class Dataset(object):
    
    """
    Manages a set of read or simulated alignments. Supports ``len()``.
    Note that, on observed data and when there are several populations
    and/or outgroups , :meth:`.sort_aligns` should be called to handle
    the case where the different populations are mixed in alignments,    
    """

    ####################################################################

    def __init__(self):

        self._aligns = []
        self._populations = set()
        
    ####################################################################

    def __len__(self):
        
        return len(self._aligns)
        
    ####################################################################

    def add(self, align):
        
        """
        Add an :class:`~egglib.Align` instance (which will be copied by
        reference).
        """
        
        self._aligns.append(align)
        
    ####################################################################

    def pops(self):
        
        """
        Returns the set of distinct populations in the data set (note
        that all populations are not required to be represented in each
        locus). The outgroup is not considered.
        """

        pops = set()
        for align in self._aligns:
            for i in align:
                if i.group != 999:
                    pops.add(i.group)
        return pops

    ####################################################################
        
    def config(self):
        
        """
        Returns a list with - for each alignment - a tuple containing
        4 items: total sample size, list of sample size per population
        (excluding outgroups), number of outgroups, alignment length.
        The alignment length is excluding not usable sites (corresponding
        to ``lseff``). Populations are sorted by their label.
        """
        
        pops = sorted(self.pops())
        res = []
        
        for align in self._aligns:
            pol = align.polymorphism()
            ns = int(pol['nseff'])
            ls = pol['lseff']
            assert ns>0 and ls>0, 'invalid alignment (no exploitable data)'
            nsi = [0] * len(pops)
            nso = 0
            for i in align:
                if i.group == 999: nso+=1
                else: nsi[pops.index(i.group)] += 1
            if ns != sum(nsi): raise ValueError, 'invalid alignment (e.g. sequence without valid data)'
            res.append((ns, nsi, nso, ls))
        
        return res

    ####################################################################

    def sort_aligns(self):
        
        """
        Sort each alignment such as they match the exported config
        (all populations appear grouped and in the increasing order,
        with outgroups at the end).
        """

        config = self.config()
        pops = sorted(self.pops()) + [999]
        for i,align in enumerate(self._aligns):
            align2 = data.Align()
            for pop in pops:
                for seq in align:
                    if seq.group == pop:
                        align2.append(*seq)
            self._aligns[i] = align2

    ####################################################################

    def iterator(self, config=None):
        
        """
        This iterator zips the object passed as *config* the alignments
        stored in the instance. The user should ensure that the object
        passed as *config* is an iterable and has the same length as
        the current instance (it will not be done automatically). Each
        iteration round returns a ``(align, configItem)`` tuple, where
        ``configItem`` an item of *config*. If *config* is ``None``, the
        result of :meth:`.config` will be used.
        """
        
        if config==None: config=self.config()
        
        for i,j in zip(self._aligns, config): yield i,j
        
# PRIOR CLASSES ########################################################

class PriorParseError(IOError):
    
    """
    Raised by Prior parse() methods when the format is found to be
    incorrect. It can be caught to auto-detect prior types.
    """
    
    pass

########################################################################

class PriorDumb(object):
    
    """
    This prior doesn't allow covariation between parameters or discrete
    categories. The probability distribution for each parameter is
    specified has a uniform or continuous statistical distribution.
    The list below presents the available distribution types, with the
    one-letter code and the list of expected parameters, expected by the
    method :meth:`.add`:

        - **U**: uniform probability between *down* and *up*.
        - **E**: exponential distribution of mean *mean*.
        - **P**: Poisson distribution of parameter *p*.
        - **G**: gamma distribution of parameter *p*.
        - **N**: normal distribution of mean *m* and standard deviation *s*.
        - **F**: parameter fixed to the value *v*.
    """

    name = "PriorDumb"
    
    ####################################################################

    _distributions = {
    # uniform distribution between down and up
    'U': (2, lambda random,down,up: random.uniform() * (up-down)+down ),
        
    # exponential distribution of mean mean
    'E': (1, lambda random,mean:    random.erand(mean) ),
        
    # Poisson distribution of parameter p
    'P': (1, lambda random,p:       random.prand(p) ),
        
    # gamma distribution of parameter p
    'G': (1, lambda random,p:       random.grand(p) ),

    # normal distribution of mean m and standard deviation s
    'N': (2, lambda random,m,s:     random.nrand() * s + m ),
    
    # parameter fixed to a given value v
    'F': (1, lambda random,v:       v )
    }
    
    ####################################################################
    
    def __init__(self, random=None):

        """
        Constructor argument *random* must be a :class:`~egglib.egglib_binding.Random`
        instance and will be used to generate pseudorandom numbers.
        """
        
        self._random = random
        self._force_positive = False
        self.clear()

    ####################################################################

    def force_positive(self):
        
        """
        Enforces that drawn parameter values are >=0 (values <0 will be
        ignored). This flag is not cancelled if :meth:`.clear` is called.

        .. versionadded:: 2.0.2

        """
        
        self._force_positive = True

    ####################################################################
    
    def clear(self):
        
        """
        Resets the instance.
        """
        
        self._data = []
    
    ####################################################################

    def str(self):
        
        """
        Generates a string representation of the instance.
        """
            
        ret= ''
        for i,j in self._data:
            ret+= '%s(%s)\n' %(i, ';'.join(map(str,j)))
        return ret

    ####################################################################

    def parse(self, string):
        
        """
        Imports data from the string *string*. The data format is: one
        token per parameter. The token can be arranged
        as one per line or separated by any white space characters.
        Each token must follow the format ``X(...)`` where *X* is the
        one-letter code specifying the type of
        distribution (``F``, ``U``, ``N``, ``G``, ``P`` or ``E``) and
        ``...``
        represents the needed parameter values in the appropriate order,
        separated by commas or semi-colons. The brackets can be replaced
        by square brackets and the model specification is
        case-independent. The function raises a
        :class:`~egglib.fitmodel.PriorParseError` in case of invalid format.
        """

        self.clear()
        items = string.split()

        for i in items:
            
            # gets the information

            r= re.match('(.)[\(|\[](.+?)[\)|\]]', i)
            if not r:
                self.clear()
                raise PriorParseError

            model = r.group(1).upper()
            if model not in self._distributions:
                self.clear()
                raise PriorParseError

            # gets parameter values

            parameters = re.split(',|;', r.group(2))
            try:
                parameters = map(float, parameters)
            except ValueError:
                self.clear()
                raise PriorParseError
                
            if len(parameters)!=self._distributions[model][0]:
                self.clear()
                raise self.PriorParseError
               
            # adds the parameter
            
            self._data.append( (model, parameters) )
                
    ####################################################################

    def draw(self):
        
        """
        Draws a :class:`~egglib.fitmodel.ParamSample` from the instance.
        """
        
        if not len(self._data):
            raise RuntimeError, 'PriorDumb cannot draw parameters (an empty distribution)'

        if self._random == None:
            raise RuntimeError, 'PriorDumb cannot draw parameters (no random number generator)'

        ps = ParamSample(len(self._data))

        for i,(m,p) in enumerate(self._data):
            if not self._force_positive:
                ps[i] = self._distributions[m][1](self._random, *p)
            else:
                c=0
                while True:
                    X = self._distributions[m][1](self._random, *p)
                    if X>=0:
                        ps[i] = X
                        break
                    c+=1
                    if c==10000:
                        raise RuntimeError, 'unable to draw parameters >=0 in prior'

        return ps

    ####################################################################

    def add(self, type, *parameters):
        
        """
        Adds a parameter to the distribution. *type* a one-letter code
        identifying the type of statistical distribution and
        *parameters* are the distribution's parameters, given in the
        appropriate order.
        """
        
        if type not in self._distributions:
            raise ValueError, '%s is not a valid specification for a prior distribution' %str(type)

        if len(parameters)!=self._distributions[type][0]:
            raise ValueError, 'invalid number of parameters'
                
        self._data.append( (type, parameters) )

    ####################################################################

    def number_of_params(self):
        
        """
        Returns the number of parameters.
        """
        
        return len(self._data)


########################################################################

class PriorParser(object):
    
    """
    This prior is used when parameters have already been drawn by any
    other means and are placed in a file. The prior itself must be of
    the form:
        
        @fname@
        
    where *fname* is the name of an external file containing parameter
    sets.  The external file must not move while the prior is in used.
    The file must have a header with parameter names. The names
    themselves are ignored (there are recognized by their order) but
    their number constrains the rest of the file. Parameters can't
    contain spaces or tabs within their name. From the second line on,
    the file must contain only parameter values, with one parameter set
    per line with values separated by spaces or tabs. The number of
    parameters must be consistent over all file.
    
    Note that the prior argument of command line tools will be, as a
    result of this syntax, of the form (for a file "prior.txt"):
    
        prior="%@prior.txt@"
        
        .. versionadded:: 2.1.4
    """

    name = "PriorParser"
    
    ####################################################################
    
    def __init__(self, random=None):
        
        """
        The "random" argument of the constructor is here for consistency
        but is ignored (this prior is not stochastic).
        """

        self.file = None
        self.nparams = 0

    ####################################################################
    
    def parse(self, string):
        
        self.clear()

        if string[0]!='@' or string[-1]!= '@':
            raise PriorParseError
        
        fname = string[1:-1]
        
        self.fname = fname
        self.file = open(fname)
        self.nparams = len(self.file.readline().strip().split())
        
    ####################################################################

    def __del__(self):
        
        self.clear()

    ####################################################################

    def clear(self):

        if self.file: self.file.close()
        self.file = None
   
    ####################################################################

    def str(self):
        
        """
        Generates a string representation of the instance.
        """
        
        return '@{0}@'.format(self.fname)

    ####################################################################

    def draw(self):
        
        """
        Draws a :class:`~egglib.fitmodel.ParamSample` from the instance.
        """
        
        if not self.file:
            raise RuntimeError, 'PriorParser cannot draw parameters (no file specified)'

        ps = ParamSample(self.nparams)
        
        bits = map(float, self.file.readline().split())
        if len(bits) != self.nparams:
            raise RuntimeError, 'PriorParser cannot draw parameters (invalid number of parameters)'
        
        ps._values = bits

        return ps

    ####################################################################

    def number_of_params(self):
        
        """
        Returns the number of parameters (of the last parsed file if
        clear() was called; 0 if no file was parsed).
        """
        
        return self.nparams
        
        
########################################################################

class PriorDiscrete(object):
    
    """
    This prior is based on discrete categories. It consists in a set of
    weighted categories with free boundaries. Within a category, the
    probability density is uniform. It allows using uniform distribution
    with fixed bounds, discretized empirical distribution and
    theoretical laws of distribution as priors.
    
    :class:`~egglib.fitmodel.PriorDiscrete` instances have a length (the number of
    categories) and are iterable. Each iteration yields a ``(p, bounds)``
    tuple with *p* the frequency of a class and *bounds* a list giving
    the bound values (themselves as a tuple) for all parameters.
    """
    
    name = "PriorDiscrete"

    ####################################################################
    
    def __init__(self, random=None):
        
        """
        Constructor argument *random* must be a :class:`~egglib.egglib_binding.Random`
        instance and will be used to generate pseudorandom numbers.
        """
        
        self.clear()
        self._random = random
        self._force_positive = False

    ####################################################################

    def force_positive(self):
        
        """
        Enforces that drawn parameter values are >=0 (values <0 will be
        ignored). This flag is not cancelled if :meth:`.clear` is called.

        .. versionadded:: 2.0.2

        """

        self._force_positive = True

    ####################################################################
    
    def clear(self):
        
        """
        Clears the instance.
        """
        
        self._data = []
        self._total = 0

    ####################################################################

    def str(self):
        
        """
        Formats the content of the instance as a string, in a format
        appropriate for passing to :meth:`.parse`.
        """
    
        ret= ''
        for i in self._data:
            ret+= '%.12f' %i[0]
            for j,k in i[1]: ret+= ' %f;%f' %(j,k)
            ret+= '\n'

        return ret

    ####################################################################

    def parse(self, string):
        
        """
        Imports data from the string *string*. The data format is: one
        line per category (in any order), each line following the format
        ``freq down1;up1 down2;up2 ...`` where *freq* is the
        frequency of the category (needs not to be relative), *down* the
        lower bound value and *up* the upper bound value for a given
        parameter. The function raises a :class:`~egglib.fitmodel.PriorParseError`
        in case of format error.
        """

        lines = string.split('\n')
        self.clear()

        for i in lines: 

            i = i.split()
            
            if not len(i):
                continue
            
            if (len(i)<2):
                self.clear()
                raise PriorParseError
            try:
                P = float(i[0])
            except ValueError:
                self.clear()
                raise PriorParseError

            ranges = []
            for j in i[1:]:

                r= re.match('(-?[\.\d]+)[;|,](-?[\.\d]+)',j)
                if not r:
                    self.clear()
                    raise PriorParseError

                bottom= float(r.group(1))
                top= float(r.group(2))
                ranges.append( (bottom, top) )

            self._total += P
            self._data.append((P,ranges))

        # checks consistency

        if not len(self._data):
            self.clear()
            raise PriorParseError

        n = len(self._data[0][1])
        for i in self._data:
            if n != len(i[1]):
                raise PriorParseError


    ####################################################################

    def draw(self):
        
        """
        Generates a set of random values for all parameters.
        """
        
        if not len(self._data):
            raise RuntimeError, 'PriorDiscrete cannot draw parameters (an empty distribution)'
            
        if self._random == None:
            raise RuntimeError, 'PriorDiscrete cannot draw parameters (no random number generator)'

        cat=0

        if (len(self._data)>1):
            x = self._random.uniform() * self._total
            acc=self._data[0][0]
            while (x>acc):
                acc+=self._data[cat+1][0]
                cat+=1
        
        ps = ParamSample(len(self._data[0][1]))
        for i, (low,high) in enumerate(self._data[cat][1]):
            if not self._force_positive:
                ps[i] = self._random.uniform() *(high-low) + low
            else:
                c=0
                while True:
                    X = self._random.uniform() *(high-low) + low
                    if X>=0:
                        ps[i] = X
                        break
                    c+=1
                    if c==10000:
                        raise RuntimeError, 'unable to draw parameters >=0 in prior'

        return ps

    ####################################################################

    def __len__(self):
        
        return len(self._data)

    ####################################################################
    
    def number_of_params(self):
        
        """
        Returns the number of parameters, (0 if no data loaded).
        """
        
        if len(self._data)==0:
            return 0
            
        else:
            return len(self._data[0][1])

    ####################################################################

    def __iter__(self):

        for i in self._data:
            yield i

    ####################################################################

    def add(self, freq, *bounds):
        
        """
        Adds a category to the distribution. *freq* gives the frequency
        of the category. The frequencies don't need to be relative.
        *bounds* must be separate 2-item lists or tuples giving the lower
        and upper bound values for each parameter.
        """
        
        # checks that the list of parameter matches

        if len(self) and len(self._data[0][1])!=len(bounds):
                raise ValueError, 'error in PriorDiscrete: inconsistent number of parameters between categories'

        # adds the category (manual deep copy - and type checking)

        P = float(freq)
        self._total += P
        
        ranges = []
        try:
            for x,y in bounds:
                ranges.append(( float(x), float(y) ))
        except TypeError, IndexError:
            raise ValueError, 'error in PriorDiscrete: invalid bound specification'

        self._data.append((P, ranges))

# DEMOGRAPHIC MODEL CLASSES ############################################

class SNM(object):

    """
    Standard Neutral Model: constant-sized single population. Allows
    optional recombination. Parameters: THETA, RHO (optional).
    """

    ####################################################################

    name = "SNM"
    parameters = ["THETA", "[RHO]"]

    ####################################################################

    def __init__(self, recombination):
        
        """
        The constructor expects a boolean to indicate whether
        recombination must be implemented.
        """
        
        if not isinstance(recombination, bool): raise TypeError, 'invalid argument for model class constructor'
        self.recombination = recombination
        self.parameters = ['THETA']
        if recombination: self.parameters.append('RHO')

    ####################################################################

    def generate(self, cfg, ps, random):

        """
        Generates a simulated dataset based on the passed sample
        configuration and the parameter sample.
        """

        if len(ps)!=(1+self.recombination):
            raise ValueError, 'invalid number of parameters for %s' %self.name
            
        if self.recombination: rho = ps[1]
        else: rho = 0.
            
        ds = Dataset()
        for ns, nsi, nso, ls in cfg:
            mutator = simul.CoalesceFiniteAlleleMutator(ps[0]*ls)
            coalParamSet = simul.CoalesceParamSet(ns, rho=rho*ls, nsites=ls)
            align = simul.coalesce(coalParamSet, mutator, 1, random)[0]
            ds.add(align)

        return ds

########################################################################

class PEM(object):

    """
    Population Expansion Model (exponential growth), with optional
    recombination. Parameters: THETA, ALPHA, RHO (optional).
    """

    ####################################################################

    name = "PEM"
    parameters = ["THETA", "ALPHA", "[RHO]"]

    ####################################################################

    def __init__(self, recombination):
        
        """
        The constructor expects a boolean to indicate whether
        recombination must be implemented.
        """

        if not isinstance(recombination, bool): raise TypeError, 'invalid argument for model class constructor'
        self.recombination = recombination
        self.parameters = ['THETA', 'ALPHA']
        if recombination: self.parameters.append('RHO')
        
    ####################################################################

    def generate(self, cfg, ps, random):

        """
        Generates a simulated dataset based on the passed sample
        configuration and the parameter sample.
        """
      
        if len(ps)!=(2+self.recombination):
            raise ValueError, 'invalid number of parameters for %s' %self.name

        if self.recombination: rho = ps[2]
        else: rho = 0.

        ds = Dataset()

        for ns, nsi, nso, ls in cfg:
            mutator = simul.CoalesceFiniteAlleleMutator(ps[0]*ls)
            coalParamSet = simul.CoalesceParamSet(ns, rho=rho*ls, nsites=ls, alpha=ps[1])
            align = simul.coalesce(coalParamSet, mutator, 1, random)[0]
            ds.add(align)

        return ds

########################################################################

class BNM(object):

    """
    Bottleneck Model, with optional recombination. Parameters:

        - THETA
        - DATE (date of the end of the bottleneck)
        - DUR (bottleneck duration)
        - BOTZISE (size of the population during the bottleneck)
        - ANCSIZE (size of the ancestral population)
        - RHO (optional)

    Note that if botsize is >1, the model can be generalized to a
    double instant change model.
    """

    ####################################################################

    name = "BNM"
    parameters = ["THETA", "DATE", "DUR", "BOTSIZE", "ANCSIZE", "[RHO]"]

    ####################################################################

    def __init__(self, recombination):
        
        """
        The constructor expects a boolean to indicate whether
        recombination must be implemented.
        """
        
        if not isinstance(recombination, bool): raise TypeError, 'invalid argument for model class constructor'
        self.recombination = recombination
        self.parameters = ['THETA', "DATE", "DUR", "BOTSIZE", "ANCSIZE"]
        if recombination: self.parameters.append('RHO')

    ####################################################################

    def generate(self, cfg, ps, random):

        """
        Generates a simulated dataset based on the passed sample
        configuration and the parameter sample.
        """
        
        if len(ps)!=(5+self.recombination):
            raise ValueError, 'invalid number of parameters for %s' %self.name
            
        if self.recombination: rho = ps[5]
        else: rho = 0.
            
        ds = Dataset()
        for ns, nsi, nso, ls in cfg:
            mutator = simul.CoalesceFiniteAlleleMutator(ps[0]*ls)
            coalParamSet = simul.CoalesceParamSet(ns, rho=rho*ls, nsites=ls)
            coalParamSet.changeAllPopulationSizes(ps[1], ps[3])
            coalParamSet.changeAllPopulationSizes(ps[1]+ps[2], ps[4])
            align = simul.coalesce(coalParamSet, mutator, 1, random)[0]
            ds.add(align)

        return ds

########################################################################

class GDB(object):

    """
    Composite-parameter bottleneck, after the formalization of Galtier,
    Depaulis and Barton bottleneck model, with optional recombination.
    The bottleneck is implemented as a number of coalescent events
    occurring precisely at the time given by the DATE parameter. The
    STRENGTH is expressed as an amount of time of the normal coalescent
    process during which only coalescent occur (no migraton, not
    mutation) and during which the global time counter doesn't change.
    Ref: Galtier *et al.* *Genetics* **155**:981-987, 2000.
    
    Parameters: THETA, DATE, STRENGTH, RHO (optional).
    """

    ####################################################################

    name = "GDB"
    parameters = ["THETA", "DATE", "STRENGTH" "[RHO]"]

    ####################################################################

    def __init__(self, recombination):
        
        """
        The constructor expects a boolean to indicate whether
        recombination must be implemented .
        """
        
        if not isinstance(recombination, bool): raise TypeError, 'invalid argument for model class constructor'
        self.recombination = recombination
        self.parameters = ['THETA', "DATE", 'STRENGTH']
        if recombination: self.parameters.append('RHO')

    ####################################################################

    def generate(self, cfg, ps, random):

        """
        Generates a simulated dataset based on the passed sample
        configuration and the parameter sample.
        """
        
        if len(ps)!=(3+self.recombination):
            raise ValueError, 'invalid number of parameters for %s' %self.name
            
        if self.recombination: rho = ps[3]
        else: rho = 0.
            
        ds = Dataset()
        for ns, nsi, nso, ls in cfg:
            mutator = simul.CoalesceFiniteAlleleMutator(ps[0]*ls)
            coalParamSet = simul.CoalesceParamSet(ns, rho=rho*ls, nsites=ls)
            coalParamSet.bottleneck(ps[1],ps[2])
            align = simul.coalesce(coalParamSet, mutator, 1, random)[0]
            ds.add(align)

        return ds

########################################################################

class IM(object):

    """
    Island Model, with optional recombination. The number of populations
    is automatically detected from the observed dataset. Parameters:
    THETA, MIGR, RHO (optional).
    """

    ####################################################################

    name = "IM"
    parameters = ["THETA", "MIGR", "[RHO]"]

    ####################################################################

    def __init__(self, recombination):
        
        """
        The constructor expects a boolean to indicate whether
        recombination must be implemented.
        """
        
        if not isinstance(recombination, bool): raise TypeError, 'invalid argument for model class constructor'
        self.recombination = recombination
        self.parameters = ['THETA', 'MIGR']
        if recombination: self.parameters.append('RHO')

    ####################################################################

    def generate(self, cfg, ps, random):

        """
        Generates a simulated dataset based on the passed sample
        configuration and the parameter sample.
        """
        
        if len(ps)!=(2+self.recombination):
            raise ValueError, 'invalid number of parameters for %s' %self.name
            
        if self.recombination: rho = ps[2]
        else: rho = 0.
            
        ds = Dataset()
        for ns, nsi, nso, ls in cfg:
            mutator = simul.CoalesceFiniteAlleleMutator(ps[0]*ls)
            coalParamSet = simul.CoalesceParamSet(nsi, rho=rho*ls, nsites=ls, M=ps[1])
            align = simul.coalesce(coalParamSet, mutator, 1, random)[0]
            ds.add(align)

        return ds

########################################################################

class IMn(object):

    """
    Island Model with different population sizes, with optional
    recombination. The size of the first population is fixed to 1,
    therefore the size of all populations with index >1 must be
    specified as parameter. Parameters: THETA, MIGR, population sizes,
    RHO (optional).
    """

    ####################################################################

    name = "IMn"
    parameters = ["THETA", "MIGR", "SIZE1", "...", "[RHO]"]

    ####################################################################

    def __init__(self, recombination):
        
        """
        The constructor expects a boolean to indicate whether
        recombination must be implemented .
        """
        
        if not isinstance(recombination, bool): raise TypeError, 'invalid argument for model class constructor'
        self.recombination = recombination
        self.parameters = None

    ####################################################################
        
    def generate(self, cfg, ps, random):

        """
        Generates a simulated dataset based on the passed sample
        configuration and the parameter sample.
        """
        
        if not len(cfg): raise ValueError, 'empty configuration! not data?'
        if len(ps)!=(2+self.recombination+len(cfg[0][1])-1):
            raise ValueError, 'invalid number of parameters for %s' %self.name
            
        if self.recombination: rho = ps[2+len(cfg[0][1])-1]
        else: rho = 0.

        if self.parameters==None:
            self.parameters = ['THETA', 'MIGR'] + ['SIZE%d' %(i+1) for i in range(1, len(cfg[0][1]))]
            if self.recombination: self.parameters.append('RHO')
            
        ds = Dataset()
        for ns, nsi, nso, ls in cfg:
            mutator = simul.CoalesceFiniteAlleleMutator(ps[0]*ls)
            sizes = [1] + ps[2:2+len(nsi)-1]
            coalParamSet = simul.CoalesceParamSet(nsi, rho=rho*ls, nsites=ls, M=ps[1], N=sizes)
            align = simul.coalesce(coalParamSet, mutator, 1, random)[0]
            ds.add(align)

        return ds

########################################################################

class IMG(object):

    """
    Island Model with exponential Growth, with optional recombination.
    Parameters: THETA, MIGR, ALPHA, 
    """

    ####################################################################

    name = "IMG"
    parameters = ["THETA", "MIGR", "ALPHA", "[RHO]"]

    ####################################################################

    def __init__(self, recombination):
        
        """
        The constructor expects a boolean to indicate whether
        recombination must be implemented.
        """
        
        if not isinstance(recombination, bool): raise TypeError, 'invalid argument for model class constructor'
        self.recombination = recombination
        self.parameters = ['THETA', 'MIGR', 'ALPHA']
        if recombination: self.parameters.append('RHO')

    ####################################################################

    def generate(self, cfg, ps, random):

        """
        Generates a simulated dataset based on the passed sample
        configuration and the parameter sample.
        """
        
        if len(ps)!=(3+self.recombination):
            raise ValueError, 'invalid number of parameters for %s' %self.name
            
        if self.recombination: rho = ps[3]
        else: rho = 0.
            
        ds = Dataset()
        for ns, nsi, nso, ls in cfg:
            mutator = simul.CoalesceFiniteAlleleMutator(ps[0]*ls)
            coalParamSet = simul.CoalesceParamSet(nsi, rho=rho*ls, nsites=ls, M=ps[1], alpha=ps[2])
            align = simul.coalesce(coalParamSet, mutator, 1, random)[0]
            ds.add(align)

        return ds

########################################################################

class IMiG(object):

    """
    Island Model with Independent exponential Growth in each population,
    with optional recombination. The growth rate of each population must
    be provided. Parameters: THETA, MIGR, ALPHA for all populations,
    RHO (optional).
    """

    ####################################################################

    name = "IMiG"
    parameters = ["THETA", "MIGR", "ALPHA1", "ALPHA2", "...", "[RHO]"]

    ####################################################################

    def __init__(self, recombination):
        
        """
        The constructor expects a boolean to indicate whether
        recombination must be implemented.
        """
        
        if not isinstance(recombination, bool): raise TypeError, 'invalid argument for model class constructor'
        self.recombination = recombination
        self.parameters = None

    ####################################################################

    def generate(self, cfg, ps, random):

        """
        Generates a simulated dataset based on the passed sample
        configuration and the parameter sample.
        """
        
        if not len(cfg): raise ValueError, 'empty configuration! not data?'
        if len(ps)!=(2+len(cfg[0][1])+self.recombination):
            raise ValueError, 'invalid number of parameters for %s' %self.name
            
        if self.recombination:
            rho = ps[2+len(cfg[0][1])]
        else: rho = 0.
            
        if self.parameters==None:
            self.parameters = ['THETA', 'MIGR'] + ['ALPHA%d' %(i+1) for i in range(len(cfg[0][1]))]
            if self.recombination: self.parameters.append('RHO')

        ds = Dataset()
        for ns, nsi, nso, ls in cfg:
            mutator = simul.CoalesceFiniteAlleleMutator(ps[0]*ls)
            alphas = ps[2:2+len(nsi)]
            coalParamSet = simul.CoalesceParamSet(nsi,
                            rho=rho*ls, nsites=ls, M=ps[1],
                            alpha=alphas)
            align = simul.coalesce(coalParamSet, mutator, 1, random)[0]
            ds.add(align)

        return ds

########################################################################

class IMiGn(object):

    """
    Island Model with Independent exponential Growth in each population,
    different population sizes and with optional recombination (the
    size of the first population is fixed to 1). The growth rate of each
    population must be provided, and the size of all populations save
    for the first one as well. Parameters: THETA, MIGR, growth rates,
    population sizes, RHO (optional).
    """

    ####################################################################

    name = "IMiGn"
    parameters = ["THETA", "MIGR", "ALPHA1", "ALPHA1", "...", "SIZE2", "...",  "[RHO]"]

    ####################################################################

    def __init__(self, recombination):
        
        """
        The constructor expects a boolean to indicate whether
        recombination must be implemented.
        """
        
        if not isinstance(recombination, bool): raise TypeError, 'invalid argument for model class constructor'
        self.recombination = recombination
        self.parameters = None

    ####################################################################
        
    def generate(self, cfg, ps, random):

        """
        Generates a simulated dataset based on the passed sample
        configuration and the parameter sample.
        """
        
        if not len(cfg): raise ValueError, 'empty configuration! not data?'
        if len(ps)!=(2+len(cfg[0][1])+(len(cfg[0][1])-1)+self.recombination):
            raise ValueError, 'invalid number of parameters for %s' %self.name
            
        if self.recombination:
            rho = ps[2+len(cfg[0][1])+len(cfg[0][1])-1]
        else: rho = 0.
            
        if self.parameters==None:
            self.parameters = (['THETA', 'MIGR'] +
                               ['ALPHA%d' %(i+1) for i in range(len(cfg[0][1]))] +
                               ['SIZE%d' %(i+1) for i in range(1, len(cfg[0][1]))])
            if self.recombination: self.parameters.append('RHO')

        ds = Dataset()
        for ns, nsi, nso, ls in cfg:
            mutator = simul.CoalesceFiniteAlleleMutator(ps[0]*ls)
            sizes = [1] + ps[2:2+len(nsi)-1]
            alphas = ps[2:2+len(nsi)]
            coalParamSet = simul.CoalesceParamSet(nsi, N=sizes,
                            rho=rho*ls, nsites=ls, M=ps[1],
                            alpha=alphas)
            align = simul.coalesce(coalParamSet, mutator, 1, random)[0]
            ds.add(align)

        return ds

########################################################################

class SM(object):

    """
    Split Model (thinking forward), with optional recombination. The
    DATE parameter sets the split date and MIGR the migration rate
    after the split. Parameters: THETA, MIGR, DATE, RHO (optional).
    """

    ####################################################################

    name = "SM"
    parameters = ["THETA", "MIGR", "DATE", "[RHO]"]

    ####################################################################

    def __init__(self, recombination):
        
        """
        The constructor expects a boolean to indicate whether
        recombination must be implemented.
        """
        
        if not isinstance(recombination, bool): raise TypeError, 'invalid argument for model class constructor'
        self.recombination = recombination
        self.parameters = ['THETA', 'MIGR', 'DATE']
        if recombination: self.parameters.append('RHO')

    ####################################################################

    def generate(self, cfg, ps, random):

        """
        Generates a simulated dataset based on the passed sample
        configuration and the parameter sample.
        """
        
        if len(ps)!=(3+self.recombination):
            raise ValueError, 'invalid number of parameters for %s' %self.name
            
        if self.recombination: rho = ps[3]
        else: rho = 0.
            
        ds = Dataset()
        for ns, nsi, nso, ls in cfg:
            mutator = simul.CoalesceFiniteAlleleMutator(ps[0]*ls)
            coalParamSet = simul.CoalesceParamSet(nsi, rho=rho*ls, nsites=ls, M=ps[1])
            for i in range(1,len(nsi)): coalParamSet.populationFusion(ps[2], 0, i)                
            align = simul.coalesce(coalParamSet, mutator, 1, random)[0]
            ds.add(align)

        return ds

########################################################################

class AM(object):

    """
    Admixture Model, with optional recombination. The DATE argument sets
    the time when ancestral populations joined and MIGR the migration
    rate that occurred between these populations. Note that the 
    migration rate must not be 0 because coalescent time might be
    infinite. Present-day samples are not structured. Parameters: THETA,
    DATE, MIGR, RHO (optional). In `abc_sample`, specify this model as
    `AM:k` where `k` is the number of ancestral populations.
    """

    ####################################################################

    name = "AM"
    parameters = ["THETA", "DATE", "MIGR", "[RHO]"]

    ####################################################################

    def __init__(self, npop, recombination):
        
        """
        The constructor expects a boolean to indicate whether
        recombination must be implemented, and an integer giving the
        number of ancestral populations.
        """
        
        if not isinstance(npop, int) or npop<2: raise TypeError, 'invalid npop argument for AM constructor'
        self.npop = npop
        if not isinstance(recombination, bool): raise TypeError, 'invalid argument for model class constructor'
        self.recombination = recombination
        self.parameters = ['THETA', 'DATE', 'MIGR']
        if recombination: self.parameters.append('RHO')

    ####################################################################

    def generate(self, cfg, ps, random):

        """
        Generates a simulated dataset based on the passed sample
        configuration and the parameter sample.
        """
        
        if len(ps)!=(3+self.recombination):
            raise ValueError, 'invalid number of parameters for %s' %self.name
            
        if self.recombination: rho = ps[3]
        else: rho = 0.
            
        ds = Dataset()
        for ns, nsi, nso, ls in cfg:
            mutator = simul.CoalesceFiniteAlleleMutator(ps[0]*ls)
            coalParamSet = simul.CoalesceParamSet(ns, rho=rho*ls, nsites=ls)
            for i in range(self.npop-1): coalParamSet.populationSplit(ps[1], 0, 1./self.npop)
            coalParamSet.changeAllMigrationRates(ps[1], ps[2])
            align = simul.coalesce(coalParamSet, mutator, 1, random)[0]
            ds.add(align)

        return ds

########################################################################

class GGDB(object):

    """
    Generalized Galtier, Depaulis and Barton with optional recombination.
    See GDB model. ANCSIZE gives the ancestral population size.
    Parameters: THETA, DATE, STRENGTH, ANCSIZE, RHO (optional).
    """

    ####################################################################

    name = "GGDB"
    parameters = ["THETA", "DATE", "STRENGTH", "ANCSIZE", "[RHO]"]

    ####################################################################

    def __init__(self, recombination):
        
        """
        The constructor expects a boolean to indicate whether
        recombination must be implemented,.
        """
        
        if not isinstance(recombination, bool): raise TypeError, 'invalid argument for model class constructor'
        self.recombination = recombination
        self.parameters = ['THETA', 'DATE', 'STRENGTH', 'ANCSIZE']
        if recombination: self.parameters.append('RHO')

    ####################################################################

    def generate(self, cfg, ps, random):

        """
        Generates a simulated dataset based on the passed sample
        configuration and the parameter sample.
        """
        
        if len(ps)!=(4+self.recombination):
            raise ValueError, 'invalid number of parameters for %s' %self.name
            
        if self.recombination: rho = ps[4]
        else: rho = 0.

        ds = Dataset()
        for ns, nsi, nso, ls in cfg:
            mutator = simul.CoalesceFiniteAlleleMutator(ps[0]*ls)
            coalParamSet = simul.CoalesceParamSet(ns, rho=rho*ls, nsites=ls)
            coalParamSet.bottleneck(ps[1], ps[2])
            coalParamSet.changeAllPopulationSizes(ps[1]+ps[2], ps[3])
            align = simul.coalesce(coalParamSet, mutator, 1, random)[0]
            ds.add(align)

        return ds

########################################################################

class MRC(object):

    """
    Migration Rate Change, with optional recombination. MIGR0 is the
    current migration rate and MIGR1 the ancestral migration rate.
    Parameters, THETA, DATE, MIGR0, MIGR1, RHO (optional).
    """

    ####################################################################

    name = "MRC"
    parameters = ["THETA", "DATE", "MIGR0", "MIGR1", "[RHO]"]

    ####################################################################

    def __init__(self, recombination):
        
        """
        The constructor expects a boolean to indicate whether
        recombination must be implemented.
        """
        
        if not isinstance(recombination, bool): raise TypeError, 'invalid argument for model class constructor'
        self.recombination = recombination
        self.parameters = ['THETA', 'DATE', 'MIGR0', 'MIGR1']
        if recombination: self.parameters.append('RHO')

    ####################################################################

    def generate(self, cfg, ps, random):

        """
        Generates a simulated dataset based on the passed sample
        configuration and the parameter sample.
        """
        
        if len(ps)!=(4+self.recombination):
            raise ValueError, 'invalid number of parameters for %s' %self.name
            
        if self.recombination: rho = ps[4]
        else: rho = 0.
            
        ds = Dataset()
        for ns, nsi, nso, ls in cfg:
            mutator = simul.CoalesceFiniteAlleleMutator(ps[0]*ls)
            coalParamSet = simul.CoalesceParamSet(ns, rho=rho*ls, nsites=ls, M=ps[2])
            coalParamSet.changeAllMigrationRates(ps[1], ps[3])
            align = simul.coalesce(coalParamSet, mutator, 1, random)[0]
            ds.add(align)

        return ds

########################################################################

class DOM(object):

    """
    Domestication model, with optional recombination. Parameters:

        - THETA
        - SIZE (size of the cultivated population)
        - DATE (date of the bottleneck)
        - DUR (duration of the bottleneck)
        - STRENGTH (size of the bottleneck population)
        - MIGR (bidirectional migration rate)
        - RHO (optional)

    The size of the wild population is 1. The domestication date is
    DATE+DUR.
    """

    ####################################################################

    name = "DOM"
    parameters = ["THETA", "SIZE", "DATE", "DUR", "STRENGTH", "MIGR", "[RHO]"]

    ####################################################################

    def __init__(self, recombination):
        
        """
        The constructor expects a boolean to indicate whether
        recombination must be implemented.
        """
        
        if not isinstance(recombination, bool): raise TypeError, 'invalid argument for model class constructor'
        self.recombination = recombination
        self.parameters = ['THETA', 'SIZE', 'DATE', 'DUR', 'STRENGTH', 'MIGR']
        if recombination: self.parameters.append('RHO')

    ####################################################################

    def generate(self, cfg, ps, random):

        """
        Generates a simulated dataset based on the passed sample
        configuration and the parameter sample.
        """
        
        if len(ps)!=(6+self.recombination):
            raise ValueError, 'invalid number of parameters for %s' %self.name
            
        if self.recombination: rho = ps[6]
        else: rho = 0.

        ds = Dataset()
        for ns, nsi, nso, ls in cfg:
            if len(nsi)!=2:
                raise ValueError, 'invalid number of populations for DOM model'
            mutator = simul.CoalesceFiniteAlleleMutator(ps[0]*ls)
            coalParamSet = simul.CoalesceParamSet(nsi,
                                        rho=rho*ls, nsites=ls,            
                                        N=[1., ps[1]], M=ps[5])
            coalParamSet.changeSinglePopulationSize(ps[2], 1, ps[4])
            coalParamSet.populationFusion(ps[2]+ps[3], 0, 1)
            align = simul.coalesce(coalParamSet, mutator, 1, random)[0]
            ds.add(align)

        return ds
        

# SUMMARY STATISTICS ###################################################


class SDZ(object):

    """
    Computes the following statistics: S, D, H (averaged over all loci
    excluding, for D and Z (standardized H), loci without polymorphism).
    Warning: when alignments have a tMRCA member, it will be assumed
    that they are simulated and that A is always the ancestral allele.
    In that case, they should not have any outgroup sequence. Alignments
    created from fasta file don't have a tMRCA member.
    """

    ####################################################################

    name = 'SDZ'

    ####################################################################

    def compute(self, ds, cfg):

        if len(ds)==0: raise ValueError, 'empty dataset (no fasta files?)'
        self.stats = [0., 0., 0.]
        nS = 0
        nD = 0
        nZ = 0
        
        for align, (ns, nsi, nso, ls) in ds.iterator(cfg):
            
            try: align.tMRCA
            except AttributeError: pass  # assume this is real data
            else: align.append('', 'A'*align.ls(), 999) # assume this is simulated data
            
            pol = align.polymorphism(skipHaplotypeDifferentiationStats=True,
                                     skipAllHaplotypeStats=True,
                                     skipDifferentiationStats=True,
                                     skipOutgroupBasedStats=False)

            if pol['S']!=None:
                nS+=1
                self.stats[0] += pol['S']
            
            if pol['D']!=None:
                nD+=1
                self.stats[1] += pol['D']

            if pol['Z']!=None:
                nZ+=1
                self.stats[2] += pol['Z']
 
        if nS>0: self.stats[0]/=nS
        else:    self.stats[0] = None
        if nD>0: self.stats[1]/=nD
        else:    self.stats[1] = None
        if nZ>0: self.stats[2]/=nZ
        else:    self.stats[2] = None

########################################################################


class TPH(object):

    """
    Computes the following statistics: thetaW, Pi, He (averaged over
    all loci).
    """

    ####################################################################

    name = 'TPH'

    ####################################################################

    def compute(self, ds, cfg):

        if len(ds)==0: raise ValueError, 'empty dataset (no fasta files?)'
        self.stats = [0., 0., 0.]
        
        for align, (ns, nsi, nso, ls) in ds.iterator(cfg):
            pol = align.polymorphism(skipHaplotypeDifferentiationStats=True,
                                     skipDifferentiationStats=True,
                                     skipOutgroupBasedStats=True)

            # arbitrarily sets per-site stats to 0
            if pol['lseff']==0: pol['thetaW']=0; pol['Pi']=0

            self.stats[0] += pol['thetaW'] * pol['lseff'] / ls
            self.stats[1] += pol['Pi']     * pol['lseff'] / ls
            self.stats[2] += pol['He']

        self.stats = map(lambda x: x/len(ds), self.stats)

########################################################################

class TPS(object):

    """
    Computes the following statistics: total thetaW, Pi for each
    populationm and Hudson's Snn (nearest neighbor statistic). The number
    of statistics will be 2 + the number of populations. Statistics are
    averaged over all loci.
    """

    ####################################################################

    name = 'TPS'

    ####################################################################

    def compute(self, ds, cfg):

        if len(ds)==0: raise ValueError, 'empty dataset (no fasta files?)'
        assert len(cfg[0][1])>1, 'needs more than 1 population'
        self.stats = [0.] * (2+len(cfg[0][1]))

        for align, (ns, nsi, nso, ls) in ds.iterator(cfg):
            c=0
            for i,n in enumerate(nsi):
                for j in range(n):
                    align.group(c, i)
                    c+=1
            pol = align.polymorphism(skipOutgroupBasedStats=True)

            # arbitrarily sets per-site stats to 0
            if pol['lseff']==0: pol['thetaW']=0; pol['pop_Pi']=[0]*len(cfg[0][1])

            self.stats[0] += pol['thetaW'] * pol['lseff'] / ls
            for i,Pi in enumerate(pol['pop_Pi']):
                self.stats[1+i] = Pi * pol['lseff'] / ls
            self.stats[-1] += pol['Snn']

        self.stats = map(lambda x: x/len(ds), self.stats)

########################################################################

class TPF(object):

    """
    Computes the following statistics: total thetaW, Pi for each
    population, and Fst. The number of statistics will be 2 + the
    number of populations. Statistics are averaged over all loci.
    """

    ####################################################################

    name = 'TPF'

    ####################################################################

    def compute(self, ds, cfg):

        if len(ds)==0: raise ValueError, 'empty dataset (no fasta files?)'
        assert len(cfg[0][1])>1, 'needs more than 1 population'
        self.stats = [0.] * (2+len(cfg[0][1]))

        cnt=0
        for align, (ns, nsi, nso, ls) in ds.iterator(cfg):

            pol = align.polymorphism(skipOutgroupBasedStats=True)

            # arbitrarily sets per-site stats to 0
            if pol['lseff']==0:
                pol['S'] = 0
                pol['thetaW'] = 0
                pol['pop_Pi'] = [0]*len(cfg[0][1])

            self.stats[0] += pol['thetaW'] * pol['lseff'] / ls
            for i,Pi in enumerate(pol['pop_Pi']):
                self.stats[1+i] = Pi * pol['lseff'] / ls
            if pol['S']>0:
                self.stats[-1] += pol['Fst']
                cnt += 1
            
        self.stats[:-1] = map(lambda x: x/len(ds), self.stats)[:-1]
        if cnt>0: self.stats[-1] /= cnt
        
########################################################################

class TPK(object):

    """
    Computes the following statistics: total thetaW, Pi for each
    population, and Kst. The number of statistics will be 2 + the
    number of populations. Statistics are averaged over all loci.
    """

    ####################################################################

    name = 'TPK'

    ####################################################################

    def compute(self, ds, cfg):

        if len(ds)==0: raise ValueError, 'empty dataset (no fasta files?)'
        assert len(cfg[0][1])>1, 'needs more than 1 population'
        self.stats = [0.] * (2+len(cfg[0][1]))

        cnt=0
        for align, (ns, nsi, nso, ls) in ds.iterator(cfg):

            pol = align.polymorphism(skipOutgroupBasedStats=True)

            # arbitrarily sets per-site stats to 0
            if pol['lseff']==0:
                pol['S'] = 0
                pol['thetaW'] = 0
                pol['pop_Pi'] = [0]*len(cfg[0][1])

            self.stats[0] += pol['thetaW'] * pol['lseff'] / ls
            for i,Pi in enumerate(pol['pop_Pi']):
                self.stats[1+i] = Pi * pol['lseff'] / ls
            if pol['S']>0:
                self.stats[-1] += pol['Kst']
                cnt += 1
            
        self.stats[:-1] = map(lambda x: x/len(ds), self.stats)[:-1]
        if cnt>0: self.stats[-1] /= cnt
        
########################################################################

class DIV(object):

    """
    Computes the following statistics: total thetaW, total Pi, total He,
    Fst, Gst, Snn, and, for each population, thetaW, Pi and He. The
    number of statistics will be 6 + 3 * the number of populations.
    Statistics are averaged over all loci.
    """

    ####################################################################

    name = 'DIV'

    ####################################################################

    def compute(self, ds, cfg):

        if len(ds)==0: raise ValueError, 'empty dataset (no fasta files?)'
        assert len(cfg[0][1])>1, 'needs more than 1 population'
        self.stats = [0.] * (6+3*len(cfg[0][1]))

        cnt=0
        for align, (ns, nsi, nso, ls) in ds.iterator(cfg):

            pol = align.polymorphism(skipOutgroupBasedStats=True)

            # arbitrarily sets per-site stats to 0
            if pol['lseff']==0:
                pol['thetaW'] = 0
                pol['Pi'] = 0
                pol['He'] = 0

            self.stats[0] += pol['thetaW'] * pol['lseff'] / ls
            self.stats[1] += pol['Pi'] * pol['lseff'] / ls
            self.stats[2] += pol['He']
            if pol['S'] > 0:
                self.stats[3] += pol['Fst']
                self.stats[4] += pol['Gst']
                self.stats[5] += pol['Snn']
                cnt += 1
            
            c = 0
            for i, n in enumerate(nsi):
                sub = align.slice(c, c+n)
                c+=n
                pol = sub.polymorphism(skipOutgroupBasedStats=True)
                if pol['lseff']==0:
                    pol['thetaW'] = 0
                    pol['Pi'] = 0
                    pol['He'] = 0
                self.stats[6+3*i] = pol['thetaW'] * pol['lseff'] / ls
                self.stats[7+3*i] = pol['Pi'] * pol['lseff'] / ls
                self.stats[8+3*i] = pol['He']


        self.stats[:3] = map(lambda x: x/len(ds), self.stats)[:3]
        self.stats[6:] = map(lambda x: x/len(ds), self.stats)[6:]
        if cnt>0: self.stats[3:6] = map(lambda x: x/cnt, self.stats)[3:6]
        
########################################################################
        
class SFS(object):

    """
    Compute the site frequency spectrum. The statistics are the average
    thetaW over all loci, and then the relative frequency of a
    user-defined number of bins of allele minor frequencies. For example,
    if the number of bins if 4, the 5 statistics will be: average thetaW,
    and then proportion of all polymorphic sites from all loci with minor
    allele <=0.125, >0.125 and <=0.25, >0.25 and <=0.375, and >0.375 and
    <=0.5. Expected argument: number of categories in the spectrum.
    """
    
    ####################################################################

    name = 'SFS'

    ####################################################################
    
    def __init__(self, number):
        
        if number<2: raise ValueError, 'not enough categories for the SFS'
        self.number = number

    ####################################################################

    def compute(self, ds, cfg):

        if len(ds)==0: raise ValueError, 'empty dataset (no fasta files?)'
        self.stats = [0.] * (self.number+1)
        Stot = 0
        bins = [(i+1)*0.5/self.number for i in range(self.number)]
       
        for align, (ns, nsi, nso, ls) in ds.iterator(cfg):

            pol = align.polymorphism(skipAllHaplotypeStats=True,
                                     skipDifferentiationStats=True,
                                     skipOutgroupBasedStats=True)

            # arbitrarily sets per-site stats to 0
            if pol['lseff']==0: pol['thetaW']=0

            self.stats[0] += pol['thetaW'] * pol['lseff'] / ls
            Stot += pol['S']
            for site in pol['sites']:
                p,q = site.alleleFrequency(0), site.alleleFrequency(1)
                p = 1.*min((p,q))/(p+q)
                for i,v in enumerate(bins):
                    if p<=v:
                        self.stats[1+i]+=1
                        break

        self.stats[0] /= len(ds)
        if Stot>0:
            self.stats[1:] = map(lambda x: x/Stot, self.stats[1:])

########################################################################

class JFS(object):

    """
    Compute the joint frequency spectrum. This set of summary statistics
    requires two populations The first two statistics are the average
    thetaW over all loci in both populations, and then the relative
    frequency of a user-defined number of bins of the frequencies of
    the minor allele in both populations. If the number of bins if 4,
    there will be 2+4**4 = 18 statistics: average thetaW in the first
    populations, in the second populations, and then the proportion of
    mutations with the minor allele at frequency <=0.125 in both
    populations, and then at frequency <=0.125 in the first population
    but at frequency >0.125 and <=0.25 in the second population, and so
    on. Expected argument: number of categories in one dimension of the
    joint spectrum.
    
    There are some restrictions when using this summary statistics set:
    there must be exactly two populations; sequences for the first
    population must be consecutive; there must be exactly two alleles at
    each site and there cannot be any missing data.
    """
    
    ####################################################################

    name = 'JFS'

    ####################################################################
    
    def __init__(self, number):
        
        if number<2: raise ValueError, 'not enough categories for the SFS'
        self.number = number

    ####################################################################

    def compute(self, ds, cfg):
        
        if len(ds)==0: raise ValueError, 'empty dataset (no fasta files?)'
        assert len(cfg[0][1])==2, 'needs exactly 2 populations'
        self.stats = [0.] * (2+self.number*self.number)
        bins = [(i+1)*0.5/self.number for i in range(self.number)]
        Stot = 0
        
        for align, (ns, [ns1,ns2], nso, ls) in ds.iterator(cfg):

            # extracts the two populations (assuming the are consecutive)

            subs = [align.slice(0, ns1), align.slice(ns1, ns1+ns2)]
            
            # perform analysis
            
            pol = align.polymorphism(skipAllHaplotypeStats=True,
                                     skipDifferentiationStats=True,
                                     skipOutgroupBasedStats=True)

            pol1 = subs[0].polymorphism(skipAllHaplotypeStats=True,
                                     skipDifferentiationStats=True,
                                     skipOutgroupBasedStats=True)

            pol2 = subs[1].polymorphism(skipAllHaplotypeStats=True,
                                     skipDifferentiationStats=True,
                                     skipOutgroupBasedStats=True)
                        
            # computes theta per pop
            # arbitrarily sets per-site stats to 0
            if pol1['lseff']==0: pol1['thetaW']=0
            if pol2['lseff']==0: pol2['thetaW']=0

            self.stats[0] += pol1['thetaW'] * pol1['lseff'] / ls
            self.stats[1] += pol2['thetaW'] * pol2['lseff'] / ls

            # gets spectrum

            sites1 = dict(zip(pol1['siteIndices'], pol1['sites']))
            sites2 = dict(zip(pol2['siteIndices'], pol2['sites']))

            for siteIndex, site in zip(pol['siteIndices'], pol['sites']):
                
                A1 = site.allele(0)
                A2 = site.allele(1)
                
                p=[0,0]
                
                for i, sitei in enumerate([sites1,sites2]):
                    
                    if siteIndex in sitei:
                        if sitei[siteIndex].allele(0)==A1:
                            i1 = 0
                            i2 = 1
                        else:
                            i1 = 1
                            i2 = 0

                        p[i] = 1. * sitei[siteIndex].alleleFrequency(i1) / sitei[siteIndex].ns()

                    else:
                        if A1 in subs[i].column(siteIndex):
                            p[i] = 1.
                        else:
                            p[i] = 0.
    
                Stot += 1
                                
                k1 = 0
                for k in range(self.number):
                    if p[0]<=bins[k]:
                        k1=k
                        break
                k2 = 0
                for k in range(self.number):
                    if p[1]<=bins[k]:
                        k2=k
                        break
                self.stats[2+k1*self.number+k2] += 1

        self.stats[0] /= len(ds)
        self.stats[1] /= len(ds)
        if Stot>0:
            self.stats[2:] = map(lambda x: x/Stot, self.stats[2:])


# LISTS ################################################################

priors = [PriorDumb, PriorDiscrete, PriorParser]
"""
This list contains the class objects (different from class *instances*,
they are the classes themselves) corresponding to priors. They must
define :meth:`parse`, :meth:`draw` and :meth:`str` methods, and class-level string
*name* and an informative docstring but this is not (yet) enforced. This
list is designed to help interactive commands to detect automatically
available priors.
"""

summstats = [SDZ, TPH, TPS, TPF, TPK, SFS, JFS, DIV]
"""
This list contains the class objects (different from class *instances*,
they are the classes themselves) corresponding to sets of summary
statistics. They must define a :meth:`compute` method taking a dataset and
a sample configuration. This method must create a :attr:`stats` member
containing a defined number of number (statistics). The constructor
might (or might not) take integer arguments. A :meth:`name` class member and
an informative docstring are also required. All this is still not (yet)
enforced. This list is designed to help interactive commands to detect
automatically available sets of summary statistics.
"""

models = [SNM, PEM, BNM, GDB, GGDB, IM, IMn, IMG, IMiG, IMiGn, MRC, AM, SM, DOM]
"""
This list contains the class objects (different from class *instances*,
they are the classes themselves) corresponding to demographic models.
They must define a :meth:`generate` method taking a configuration list and
a param sample, their constructor must take 0 or more integer arguments
and then a boolean indicating whether recombination occurs. They must
define a class-level string :attr:`name` and a class-level list of strings
:attr:`parameters` and an informative docstring. All this is not (yet)
enforced. This list is designed to help interactive commands to detect
automatically available models.
"""

def add_model(name):
    
    """
    Adds a `name` model contained in the file `name`.py. The model will
    be accessible in the **fitmodel.models** list.
    """
    if not os.path.isfile(name+'.py'):
        raise ValueError, 'cannot add model: no such file: "%s.py"' %(name)
    try:
        sys.path.append(os.getcwd())
        mod = __import__(name)
        del sys.path[-1]
    except ImportError, SyntaxError:
        raise ValueError, 'cannot add model: cannot import %s' %name

    try:
        CustomModel = mod.__getattribute__(name)
        CustomModel.name # check that the required members are there
        CustomModel.__init__
        CustomModel.generate
    except AttributeError:
        raise ValueError, 'cannot add model: invalid specification: %s' %name
        
    models.append(CustomModel)


########################################################################

def import_posterior(fname):
    
    """
    Imports a posterior file
    
    *fname* must be the name of a file containing fitted ABC data, one
    sample per line and one parameter per column. Header line is
    optional and is automatically detected. If present, no parameter
    name can be provided as a number.
    
    Returns a tuple ``(params, data)`` where ``params`` is the list of
    parameter names (automatic if header is not present) and ``data`` is
    a list of lists, one list per parameter (note that the returned list
    is transposed with respect to the input file). Beware that headers
    with number-only parameters will be mistaken with values.
    """

    # imports the data

    f = open(fname)
        
    # reads first line 
        
    line = f.readline().split()

    try:
        line = map(float, line)
    except ValueError:
        params = line
        data = [[] for i in line]
    else:
        params = ['param%d' %(i+1) for i in range(len(line))]
        data = [[i] for i in line]

    # reads rest of file

    for line in f:
        line = line.split()
        if len(line) != len(data):
            raise IOError, 'number of lines varies in data from %s' %fname
        try:
            line = map(float, line)
        except ValueError:
            raise IOError, 'cannot import data from %s' %fname
        map(lambda (x,y): x.append(y), zip(data,  line))

    f.close()

    return params, data
