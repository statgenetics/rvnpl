"""
This module wraps the C++-implemented coalescent simulator. To perform
coalescent simulations, the user must create a :class:`~egglib.simul.CoalesceParamSet`
instance with the desired parameter values (and call the appropriate
methods to set optional demographic changes), and an instance of one of
the subclasses of :class:`~egglib.simul.CoalesceMutator`, and then pass these two
objects to the simulation function. The function will return a list of
:class:`~egglib.Align` or :class:`~egglib.SSR` instances.

.. Note:: Only polymorphic sites are returned. This can be significant
          when computed polymorphism statistics be expressed per site.
          Such statistics should then be multiplied by the length of the
          corresponding alignments (that is, the number of segregating
          sites) to obtain a gene-wise value.
"""


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

import egglib_binding, data


# thi Random is stored at the module level to ensure that consecutive calls to coalesce() use different seeds
RANDOM = egglib_binding.Random()


########################################################################

class CoalesceParamSet:
    
    """ Holder for coalescent parameters related to the reconstruction
    of genealogical trees. The value of different demographic parameters
    are passed to the constructor, and the user can subsequently add
    demographic changes. The order in which demographic changes are
    entered is not important (they are sorted automatically by date).
    
    .. versionadded:: 2.0.1 """
    
    
    def __init__(self, singleSamples, doubleSamples=None, s=0., rho=0.,
                                    nsites=1, alpha=0., M=0.1, N=1.):
        
        """ Constructor arguments:
        
        *singleSamples* and *doubleSamples* specify the number of
        sampled genes. Since the underlying model is diploid, it is
        possible to sample either one or both chromosomes of a given
        individual. *singleSamples* specifies the number of individuals
        of which one (random) chromosome was sampled. *doubleSamples*
        specifies the number of individuals of which both chromosomes
        were sampled. For each of these options, an :class:`int`, a
        sequence can be passed, as well as ``None``, which means that no
        samples of this type were collected. The number of populations
        in the model is implied by these two options. The number of
        populations is given by the length of the passed list (one
        population if integers are passed). The number of populations
        implied by these two options must be consistent, except if one
        of them is ``None`` (in what a list of ``0`` is assumed). An
        example is ``singleSamples=[10, 20], doubleSamples=[5, 0]`` 
        which sums up to 20 genes in both populations.
        
        *s* gives the selfing rate. ``s=0.`` means total panmixia and
        ``s=1.`` total autogamy.
        
        *rho* is the recombination rate, expressed as ``4Nc`` where
        ``N`` is the population size (number of diploid individuals)
        and ``c`` the per-gene instantaneous recombination rate. If
        *rho* is >0, it is required to provide a value of *nsites*
        >1.
        
        *nsites* is the number of recombining sites (ignored if *rho* is
        zero).
        
        *alpha* is the exponential growth/decline rate. A value of
        *alpha* greater than ``0`` means that the population size was
        smaller in the past. If a :class:`float` is passed, the same
        value will be applied to all populations. If a sequence is
        passed, is specifies individual, population-specific, rates and
        the length of the sequence must match the length of
        *singleSamples*  and/or *doubleSamples*.
        
        *M* is the migration rate, expressed as ``4Nm``  where ``N`` is
        the population size (number of diploid individuals) and ``m``
        the instantaneous migration rate (the probability that a given
        individual changes deme). If *M* is a :class:`float`, all non-
        diagonal values of the migration matrix will be set to *M/(k-1)*
        where ``k`` is the number of populations in the system (the
        proportion of migrants will always be *M*). Otherwise *M* must be
        a sequence of sequences of dimensions ``kxk`` where ``k`` is the
        number of populations. It is illegal to set the diagonal of the
        matrix, and therefore it is enforced that all values of the
        diagonal are ``None``. Here is an example for 2 populations with
        asymetric migration rates: ``M=[[None, 0.1], [0.2, None]]``, and
        an example for 3 populations with a "stepping stone" model of
        migration: ``M=[[None, 0.1, 0.0], [0.1, None, 0.1], [0., 0.1, 
        None]]``.
        
        *N* is the relative size of all populations. It is possible to
        pass a :class:`float`, what sets all populations to the same
        size and corresponds to rescaling the time scale. Otherwise a
        sequence of :class:`float` must be passed, and the length of
        this sequence must match the length of *singleSamples*  and/or
        *doubleSamples*. """
        
        self.paramSet = egglib_binding.ParamSet()
        
        # process sample arguments
        if isinstance(singleSamples, int):
            singleSamples = [singleSamples]
        if isinstance(doubleSamples, int):
            doubleSamples = [doubleSamples]
        if singleSamples!=None and ('__len__' not in dir(singleSamples) or '__getitem__' not in dir(singleSamples)):
            raise ValueError, 'invalid value for singleSamples'
        if doubleSamples!=None and ('__len__' not in dir(doubleSamples) or '__getitem__' not in dir(doubleSamples)):
            raise ValueError, 'invalid value for doubleSamples'

        # gets the number of populations
        if singleSamples!=None:
            k=len(singleSamples)
        elif doubleSamples!=None:
            k=len(doubleSamples)
        else:
            raise ValueError, 'at least one of singleSamples and doubleSamples must be passed'
        if singleSamples!=None and doubleSamples!=None:
            if len(singleSamples)!=len(doubleSamples):
                raise ValueError, 'inconsistent number of populations inferred from different arguments'
        if k<1:
            raise ValueError, 'invalid number of populations'

        # creates the populations
        for i in range(k-1):
            self.paramSet.addPopulation(0.1)

        # loads the sample sizes
        if singleSamples!=None:
            for i in range(k):
                if not isinstance(singleSamples[i], int):
                    raise ValueError, 'invalid number of haploid samples (population %d)' %i
                self.paramSet.singles(i, singleSamples[i])

        if doubleSamples!=None:
            for i in range(k):
                if not isinstance(doubleSamples[i], int):
                    raise ValueError, 'invalid number of diploid samples (population %d)' %i
                self.paramSet.doubles(i, doubleSamples[i])
                
        # loads the selfing rate
        if not isinstance(s, (float, int)) or s<0. or s>1.:
            raise ValueError, 'invalid selfing rate'
        self.paramSet.selfingRate(s)
                
        # loads the recombination rate
        if not isinstance(rho, (float, int)) or rho<0.:
            raise ValueError, 'invalid recombination rate'
        if (rho>0):
            if nsites<1: raise ValueError, 'invalid number of recombining segments'
            self.paramSet.numberOfSegments(nsites)
            self.paramSet.recombinationRate(rho)
        
        # loads the exponential growth/decline parameter
        if isinstance(alpha, (float, int)):
            for i in range(k):
                self.paramSet.growthRate(i, alpha)
        elif '__len__' in dir(alpha) and '__getitem__' in dir(alpha):
            if len(alpha)!=k:
                raise ValueError, 'the number of exponential growth rates doesn\'t match the number of populations'
            for i in range(k):
                self.paramSet.growthRate(i, alpha[i])
        else:
            raise ValueError, 'invalid exponential growth/decline parameter'
        
        # loads the relative population size
        if isinstance(N, (float, int)):
            for i in range(k):
                self.paramSet.populationSize(i, N)
        elif '__len__' in dir(N) and '__getitem__' in dir(N):
            if len(N)!=k:
                raise ValueError, 'the number of population sizes doesn\'t match the number of populations'
            for i in range(k):
                self.paramSet.populationSize(i, N[i])
        else:
            raise ValueError, 'invalid population size'
        
        # loads the migration rates
        if isinstance(M, (float, int)):
            if M<0.:
                raise ValueError, 'invalid migration rate (must be >=0)'
            self.paramSet.migrationRate(M)
        elif '__len__' in dir(M) and '__getitem__' in dir(M):
            if len(M)!=k:
                raise ValueError, 'migration rate matrix size doesn\'t have the expected size (%dx%d)' %(k,k)
            for i in range(k):

                if not('__len__' in dir(M[i]) and '__getitem__' in dir(M[i])):
                    raise ValueError, 'invalid migraton rate argument'
                if len(M[i])!=k:
                    raise ValueError, 'migration rate matrix size doesn\'t have the expected size (%dx%d)' %(k,k)
                for j in range(k):
                    if i==j: # diagonal
                        if M[i][j]!=None:
                            raise ValueError, 'miration rate matrix must have None at all diagonal values'
                    else:
                        if not isinstance(M[i][j], (float, int)):
                            raise ValueError, 'invalid migraton rate argument'
                        if M[i][j]<0.:
                            raise ValueError, 'invalid migration rate (must be >=0)'
                        self.paramSet.pairwiseMigrationRate(i, j, M[i][j])
        else:
            raise ValueError, 'invalid migraton rate argument'

        # list to store all Change objects
        # because they must live during the simulation (as clearly stated
        # by the documentation of SeqTools)
        self.ChangeCache = []

    ####################################################################
    
    def populationFusion(self, date, mother, daughter):
        
        """ At time *date*, all the lineages from the population
        *daughter* are moved to the population *mother* and all mutation
        rates to the population *daughter* are cancelled. """
        
        x = egglib_binding.PopulationFusion(date, mother, daughter)
        self.ChangeCache.append(x)
        self.paramSet.addChange(x)


    ####################################################################
    
    def populationSplit(self, date, population, probability):
        
        """ A the time given by *date*, the population *population* is
        split in two. An additional population (incremented from the
        current total number of populations) is created and lineages are
        randomly picked from population *population* and moved to the
        new population, with probability *probability*. If *probability*
        is ``0.``, the simulator creates an empty population (thinking
        forward in time, this corresponds to a population extinction)."""
        
        x = egglib_binding.PopulationSplit(date, population, probability)
        self.ChangeCache.append(x)
        self.paramSet.addChange(x)


    ####################################################################

    def changeAllMigrationRates(self, date, value):
        
        """ At time *date*, change all pairwise migration rates to
        *value* """
        
        x = egglib_binding.AllMigrationRateChange(date, value)
        self.ChangeCache.append(x)
        self.paramSet.addChange(x)
        

    ####################################################################
    
    def changeAllPopulationSizes(self, date, value):
        
        """ At time *date*, change all populations sizes to *value* """
        
        x = egglib_binding.AllPopulationSizeChange(date, value)
        self.ChangeCache.append(x)
        self.paramSet.addChange(x)
        

    ####################################################################
    
    def bottleneck(self, date, strength):
        
        """ At time *date*, apply a bottleneck of strength *strength* to
        all populations. The bottleneck strength corresponds to an
        amount of time where the time counted is blocked and only
        coalescences are allowed (resulting in a given - and random -
        number of instantaneous coalescence with branches of length 0). """
         
        x = egglib_binding.Bottleneck(date, strength)
        self.ChangeCache.append(x)
        self.paramSet.addChange(x)
        

    ####################################################################
    
    def changeAllGrowthRates(self, date, value):
        
        """ At time *date*, change all growth rates to *value*. """
         
        x = egglib_binding.GrowthRateChange(date, value)
        self.ChangeCache.append(x)
        self.paramSet.addChange(x)
        

    ####################################################################

    def changeSelfingRate(self, date, value):
        
        """ At time *date*, change all selfing rate to *value*. """
         
        x = egglib_binding.SelfingRateChange(date, value)
        self.ChangeCache.append(x)
        self.paramSet.addChange(x)
        

    ####################################################################

    def changePairwiseMigrationRate(self, date, source, dest, migr):
        
        """ At time *date*, change all the migration rate from
        population *source* to the population *dest* to the value
        *migr*. It is illegal to change a value of the diagonal. """
         
        x = egglib_binding.SingleMigrationRateChange(date, source, dest, migr)
        self.ChangeCache.append(x)
        self.paramSet.addChange(x)
        

    ####################################################################

    def singlePopulationBottleneck(self, date, population, strength):
        
        """ At time *date*, apply a bottleneck of strength *strength* to
        population *population*. The bottleneck strength corresponds to
        an amount of time where the time counted is blocked and only
        coalescences are allowed (resulting in a given - and random -
        number of instantaneous coalescence with branches of length 0). """
         
        x = egglib_binding.PopulationBottleneck(date, population, strength)
        self.ChangeCache.append(x)
        self.paramSet.addChange(x)
        

    ####################################################################

    def changeSinglePopulationGrowthRate(self, date, population, value):
        
        """ At time *date*, change growth rate of population
        *population* to *value*. """
         
        x = egglib_binding.PopulationGrowthRateChange(date, population, value)
        self.ChangeCache.append(x)
        self.paramSet.addChange(x)
        

    ####################################################################

    def changeSinglePopulationSize(self, date, population, value):
        
        """ At time *date*, change size of population *population* to
        *value* """
        
        x = egglib_binding.SinglePopulationSizeChange(date, population, value)
        self.ChangeCache.append(x)
        self.paramSet.addChange(x)
        
    ####################################################################
    
    def samples(self):
        
        """
        Returns a list of ``(d,s)`` tuples (one per population), with
        *d* the number of diploid samples and *s* the number of haploid
        samples from a given population
        """
        
        return [
            (self.paramSet.doubles(i), self.paramSet.singles(i))
                  for i in range(self.paramSet.numberOfPopulations()) ]


########################################################################

class CoalesceMutator:
    
    """ Base class for mutator objects.
    
    .. versionadded:: 2.0.1 """
    
    
    def __init__(self, theta):
        
        """"""
        
        self.mutator = egglib_binding.Mutator()
        self.mutator.mutationRate( theta )
        
    
    ####################################################################
        
    def fixedNumberOfMutation(self, number):
        
        """ Fix the final number of mutations to *number*. It is not
        guaranteed that it equals to the number of polymorphic sites
        unless the number of mutable sites is infinite. It is required
        that the object was created with *theta* = 0 to use this
        method. """
        
        self.mutator.fixedNumberOfMutations( number )
    
    
    ####################################################################
    
    def setSites(self, sites):    
        
        """ *sites* gives the number of mutable sites. If it equals to
        ``0`` or ``False`` or has length ``0``, an infinite number of
        sites will be assumed and all mutations will hit a different
        site. If *sites* is an integer, it specifies the number of
        mutable sites. They all will be placed at equal distance along
        a virtual chromosome bound by 0 and 1 and have equal mutation
        weight. Otherwise, *sites* must be a sequence. The length of the
        sequence gives the number of sites. Each site must be
        represented by a sequence of two items: the site position
        (between 0. and 1.) and its weight. The weight can be any
        positive number and gives the relative probability that a given
        mutation hits that particular site. """
    
        if not sites or ('__len__' in dir(sites) and len(sites)==0):
            return
            
        if isinstance(sites, int):
            if sites<0:
                raise ValueError, 'mutable sites specification cannot be a negative integer'
            self.mutator.numberOfSites(sites)
            return
            
        if not ('__len__' in dir(sites) and '__iter__' in dir(sites)):
            raise ValueError, 'mutable sites specification must be an iterable'

        self.mutator.numberOfSites(len(sites))
        for i,v in enumerate(sites):
            if '__len__' not in dir(v) or len(v)!=2 or not (isinstance(v[0], (float, int)) and isinstance(v[1], (float, int))):
                raise ValueError, 'mutable sites specification must contain pairs of floats'
            self.mutator.sitePosition(i, v[0])
            self.mutator.siteWeight(i, v[1])

    ####################################################################
    
    def numberOfMutations(self):

        """
        Returns the number of mutations yielded by the last simulation
        (returns 0 by default).
        """
    
        return self.mutator.numberOfMutations()



########################################################################

class CoalesceFiniteAlleleMutator(CoalesceMutator):
    
    """ Represents a mutation model with fixed number of alleles. At
    each mutation, alleles are drawn from a finite set of possible
    values. It is possible to set no equiprobable transitions to the
    different transition using the method :meth:`transitionWeights`. 
    This model sets by default an infinite number of mutable sites. """

    
    def __init__(self, theta=0, alleles=2, randomAncestralState=False):
        
        """ Constructors arguments: *theta* is the mutation rate,
        expressed as ``4Nu`` where ``N`` is the number of diploid
        individuals of a population and ``u`` is the per-gene mutation
        rate. *alleles* is the number of possible alleles. If
        *randomAncestralState*, the ancestral state is drawn randomly
        from possible states (otherwise, it will be zero)."""
        
        CoalesceMutator.__init__(self, theta)
        self.mutator.mutationModel('F')
        
        if isinstance(alleles, (int, long)):
            numberOfAlleles = alleles
            self.alleleValues = None
        else:
            try:
                numberOfAlleles = len(alleles)
                for i in alleles:
                    if not isinstance(i, str) or len(i)!=1:
                        raise ValueError, 'invalid `allele` argument for `CoalesceFiniteAlleleMutator`'
                self.alleleValues = str([i for i in alleles])
            except ValueError:
                raise ValueError, 'invalid `allele` argument for `CoalesceFiniteAlleleMutator`'
            
        
        if numberOfAlleles<2:
            raise ValueError, 'CoalesceFiniteAlleleMutator requires at leat 2 alleles'
        self.mutator.numberOfAlleles(numberOfAlleles)
        
        if randomAncestralState:
            self.mutator.randomAncestralAllele(True)


    ####################################################################
    
    def transitionWeights(self, matrix):
        
        """ Sets the weights to apply to each possible transition. Here,
        a transition is taken as any mutation from one character to an
        other. If the number of alleles is *k*, *matrix* must be a
        sequence of *k* sequences of *k* weights (strictly positive
        values). The higher the weight, the more likely the transition.
        The values on the diagonal are required to be ``None``. """
        
        if '__len__' not in dir(matrix) or '__getitem__' not in dir(matrix):
            raise TypeError, 'invalid matrix for transitionWeights'
        if len(matrix)!=self.mutator.numberOfAlleles():
            raise ValueError, 'invalid matrix for transitionWeights: expects dimensions %dx%d' %(self.mutator.numberOfAlleles(), self.mutator.numberOfAlleles())
        
        for i in range(self.mutator.numberOfAlleles()):
 
            if '__len__' not in dir(matrix[i]) or '__getitem__' not in dir(matrix[i]):
                raise TypeError, 'invalid matrix for transitionWeights'
            if len(matrix[i])!=self.mutator.numberOfAlleles():
                raise ValueError, 'invalid matrix for transitionWeights: expects dimensions %dx%d' %(self.mutator.numberOfAlleles(), self.mutator.numberOfAlleles())
     
            for j in range(self.mutator.numberOfAlleles()):
                
                if i==j:
                    if matrix[i][j]!=None:
                        raise ValueError, 'matrix for transitionWeights must have None at all values of the diagonal'
                    
                else:
                    if not isinstance(matrix[i][j], (float, int)):
                        raise TypeError, 'transition weights must be numbers'
                    self.mutator.transitionWeight(i, j, matrix[i][j])
    

########################################################################

class CoalesceInfiniteAlleleMutator(CoalesceMutator):
    
    """ Represents a mutation model with an infinite number of alleles.
    Each mutation creates a different new allele. This model sets by
    default one possible mutable site. """

    
    def __init__(self, theta=0):
        
        """ *theta* is the mutation rate, expressed as ``4Nu`` where
        ``N`` is the number of diploid individuals of a population and
        ``u`` is the per-gene mutation rate. """
        
        CoalesceMutator.__init__(self, theta)
        self.mutator.numberOfSites(1)
        self.mutator.mutationModel('I')


########################################################################

class CoalesceStepwiseMutator(CoalesceMutator):
    
    """ Represents the stepwise mutation model. Each mutation randomly
    increments or decrements the current allele value with a step of one
    unit. This model sets by default one possible mutable site. """


    def __init__(self, theta=0):
        
        """ *theta* is the mutation rate, expressed as ``4Nu`` where
        ``N`` is the number of diploid individuals of a population and
        ``u`` is the per-gene mutation rate. """
        
        CoalesceMutator.__init__(self, theta)
        self.mutator.numberOfSites(1)
        self.mutator.mutationModel('S')


########################################################################

class CoalesceTwoPhaseMutator(CoalesceMutator):
    
    """ Represents the two-phase mutation model. Each mutation randomly
    increments or decrements the current allele value with a variable
    step. With a given probability, the step is 1; otherwise it is drawn
    from a gamma distribution. This model sets by default one possible
    mutable site. """


    def __init__(self, theta=0, proba=0.5, param=0.5):
        
        """ *theta* is the mutation rate, expressed as ``4Nu`` where
        ``N`` is the number of diploid individuals of a population and
        ``u`` is the per-gene mutation rate. *proba* gives the
        probability of drawing the mutation step from the gamma
        distribution of parameter *param*. Both *proba* and *param* must
        be in range ``[0,1]``. """
        
        CoalesceMutator.__init__(self, theta)
        self.mutator.numberOfSites(1)
        self.mutator.mutationModel('T')
        self.mutator.TPMparam(param)
        self.mutator.TPMproba(proba)


########################################################################

def coalesce(paramSet, mutator, repets=1, random=None, maxNumberOfIterations=1000000, convert=True, forceSSR=False):
    
    """
    Generates simulated datasets, using the incorporated coalescent
    simulator.
    
    *paramSet* must be a :class:`~egglib.simul.CoalesceParamSet` instance; it holds
    parameters controlling the demographic process (and therefore the
    shape of generated genealogical trees).
    
    *mutator* must be an instance of a subclass of
    :class:`~egglib.simul.CoalesceMutator` (but not :class:`~egglib.simul.CoalesceMutator` itself);
    it holds parameters controlling the generation of genetic data
    (sequence or microsatellites).
    
    *repets* gives the number of repetitions to perform.
    
    *random* controls the pseudo-random number generator. This argument
    is polymorph: if ``None``, the generator is automatically set (using
    the system clock to determine seed values); if *random* is a
    sequence of two numbers, they will be used a seeds of the random
    number generator (therefore allowing to replicate results with a
    given set of parameters); alternatively, *seeds* can be a
    :class:`~egglib.egglib_binding.Random` instance.
    
    *maxNumberOfIterations* sets the maximum number of steps in the
    coalescent algorithm (provided as a safeguard).
    
    *convert* if ``True``, converts the generated :class:`~egglib.egglib_binding.DataMatrix`
    objects to high-level types (Align or SSR). Otherwise, returns them
    as is.

    *forceSSR* if ``True`` and if *convert* is ``True``, return
    :class:`~align.SSR` instances whatever options were passed regarding
    the mutation model (that is, even if a finite-allele model with 4 or
    less alleles was used). The value of this option is ignored if
    *convert* is ``True``.
    
    Returns a list of simulated data sets.  Data sets are instances of
    :class:`~egglib.egglib_binding.DataMatrix` if *convert* is ``False``.
    Otherwise, they are automatically converted to :class:`~egglib.SSR`
    instances or to :class:`~egglib.Align` instances. Conversion to
    :class:`~egglib.Align` is only performed if the mutator is of type
    :class:`~egglib.simul.CoalesceFiniteAlleleMutator`, the number of
    alleles is smaller than or equal to 4 and if *forceSSR* is ``False``.
    
    All returned instances have non-default members providing simulation
    statistics:
        - `tMRCA` (the time to the  most recent common ancestor),
        - `totLength` (sum of tree branch lengths),
        - `nMutations` (number of mutations that occurred),
        - `nRecomb` (number of recombination events).
    
    .. versionadded:: 2.0.1

    .. versionchanged:: 2.1.0
       *forceSSR* is added, and :class:`~egglib.Align` instances are no
       longer returned if the number of alleles is above 4.
    
    """
    
    # checks arguments
    if not isinstance(paramSet, CoalesceParamSet):
        raise TypeError, 'argument `paramSet` of coalesce must be a CoalesceParamSet instance'
    if not isinstance(mutator, (CoalesceFiniteAlleleMutator, 
            CoalesceInfiniteAlleleMutator, CoalesceStepwiseMutator,
            CoalesceTwoPhaseMutator)):
        raise TypeError, 'argument `mutator` of coalesce is not of one of the expected types'
    if not isinstance(repets, int) or repets<1:
        raise ValueError, 'argument `repets` of coalesce must be an integer > 0'

    # initializes the random number generator
    if random==None:
        random = RANDOM
    elif isinstance(random, egglib_binding.Random):
        pass
    else:
        try:
            seed1, seed2 = map(int, random)
        except (ValueError, TypeError):
            raise TypeError, 'argument `random` of `coalesce` must be `None`, a Random instance or a sequence of two numbers'
        random = egglib_binding.Random(seed1, seed2)
    
    # performs the iteration over repetitions
    results = []
    
    for c in range(repets):
        
        # creates the controller
        controller = egglib_binding.Controller(paramSet.paramSet, random)

        # reconstructs the tree
        iter=0
        while controller.step()>1:
            iter+=1
            if iter>=maxNumberOfIterations:
                raise RuntimeError, 'reached maximum number of coalescence iterations (there might be an unconnected population)'

        # applies mutations
        arg = controller.getArg()
        simul = mutator.mutator.mute(arg, random)

        # sets group labels
        c=0
        p=0
        for d,s in paramSet.samples():
            for i in range(d+s):
                simul.populationLabel(c, p)
                c+=1
            p+=1
        
        # transforms in appropriate data type
        if not convert:
            simul.tMRCA = arg.ageUltimateMRCA()
            #simul.vMRCA = [arg.ageMRCA(i) for i in range(arg.numberOfSegments)]
            simul.totLength = arg.totalLength
            simul.nRecomb = arg.numberOfRecombinationEvents
            simul.nMutations = mutator.numberOfMutations()
            results.append(simul)
            
        elif (isinstance(mutator, (CoalesceInfiniteAlleleMutator,
                    CoalesceStepwiseMutator, CoalesceTwoPhaseMutator))
              or forceSSR or mutator.mutator.numberOfAlleles()>4):
            
            ssr = data.SSR()
            config = [(paramSet.paramSet.doubles(i), paramSet.paramSet.singles(i))
                        for i in range(paramSet.paramSet.numberOfPopulations())]
            
            simul.shift(0)
            
            ssr.load(simul, config)
            ssr.tMRCA = arg.ageUltimateMRCA()
            #ssr.vMRCA = [arg.ageMRCA(i) for i in range(arg.numberOfSegments)]
            ssr.totLength = arg.totalLength
            ssr.nRecomb = arg.numberOfRecombinationEvents
            ssr.nMutations = mutator.numberOfMutations()
            results.append(ssr)
        
        else:
            align = data.Align()
            align._object = egglib_binding.Convert.align(simul, 0, random)
            align.tMRCA = arg.ageUltimateMRCA()
            #align.vMRCA = [arg.ageMRCA(i) for i in range(arg.numberOfSegments)]
            align.totLength = arg.totalLength
            align.nRecomb = arg.numberOfRecombinationEvents
            align.nMutations = mutator.numberOfMutations()
            results.append(align)
                
    # done
    return results

########################################################################
