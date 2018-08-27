"""
Executable tools. This module defines the classes :class:`Option` and
:class:`BaseCommand` that allow to develop automatically executable
commands. A utilitary script ``egglib`` will load this module at runtime
and let the user run any of this commands as if they were independent
programs.

"""



__license__ = """
    Copyright 2010-2012 Stephane De Mita, Mathieu Siol

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


import abc, sys, os, math, re, StringIO, glob, inspect, tempfile, multiprocessing, random, time, string
import egglib_binding, tools, fitmodel, data, wrappers, simul


########################################################################

class Option:
    
    """
    Defines the type, default value and requirements of all program
    options (except flags that are dealt with otherwise).
    """

    ####################################################################

    def __init__(self, name, doc, convert=str, default=None, requirements=[]):
        
        """
        *name* is the name of the option, *doc* must be a string
        providing documentation. *convert* should be a function
        taking a string and returning a value of the appropriate type
        (can also be classes such as :class:`int`, :class:`float` or
        *lambda* expressions, providing that they take a string and
        process it the appropriate way), *default* is the default value
        (``None`` means that the option must be specified) and
        *requirements* is a list of requirements, each of them under the
        form of a function taking a possible option value and returning
        ``True`` is the option is valid.
        """
        
        self.name = name
        self.doc = doc
        self.convert = convert
        self.default = default
        self.requirements= requirements

########################################################################

class BaseCommand:
    
    """
    Abstract base class for executable commands. Several members and
    function have to be overriden to provide all information needed by
    the ``egglib`` script.
    """

    ####################################################################
    
    __metaclass__ = abc.ABCMeta
    
    ####################################################################
    
    def __new__(cls):
        
        o = object.__new__(cls)
        o._debug = False
        o._quiet = False
    
        return o
    
    ####################################################################

    options = []
    """
    This member must be overriden. List of :class:`~Option` instances.
    """
    
    flags = []
    """
    This member must be overriden. List of flag, each given as a tuple
    of two strings: ``(label, description)``.
    
    """
    
    brief = ""
    """
    This member must be overriden. One-line summary.
    """
    
    description = ""
    """
    This member must be overriden. Full description.
    """

    ####################################################################
    
    def get_debug(self): return self._debug
    def set_debug(self, value): self._debug = value
    debug = property(get_debug, set_debug, doc='Flag indicating whether \
        full error messages should be returned.')

    ####################################################################
    
    def get_quiet(self): return self._quiet
    def set_quiet(self, value): self._quiet = value
    quiet = property(get_quiet, set_quiet, doc='Flag indicating whether \
        information should be displayed in the standard output stream \
        (if *quiet* is ``True``, information is not displayed).')

    ####################################################################
    
    @classmethod
    def doc(cls):

        try:
            name = re.match('<class \'egglib\.utils\.(.+)\'>', str(cls)).group(1)
        except AttributeError:
            try:
                name = re.match('<class \'__main__\.(.+)\'>', str(cls)).group(1)
            except AttributeError:
                raise

        LINE = 72
        brief = cls.brief
        descr = tools.wrap(cls.description, 72, 0)

        string = """%s: %s

%s

General usage:

    egglib %s OPTION1=VALUE OPTION2=VALUE ... FLAG1 FLAG2 ...

Options:

""" %(name, brief, descr, cls.__name__)

        L = max([len(option.name) for option in cls.options]) + 2
        for option in cls.options:
            substring = '    '
            substring += (option.name + ' ').ljust(L, '.')
            substring += ' '
            substring += option.doc

            if option.default != None:
                substring += ' (default: `%s`)' %option.default
            else:
                substring += ' (required)'
            
            string += tools.wrap(substring, LINE, L+5)
            string += '\n'
            
        string += '\nFlags (inactive by default):\n\n'
        flags = cls.flags + [('quiet', 'Runs without console output'),
                             ('debug', 'Show complete error messages')]
        L = max([len(name) for name,doc in flags]) + 2
        for name,doc in flags:
            name += ' '
            string += tools.wrap('    %s %s' %(name.ljust(L,'.'), doc), LINE, L+5)
            string += '\n'
        
        return string


    ####################################################################

    def run(self, *fargs, **kwargs):
    
        """
        Execute commands. *fargs* are flag arguments and *kwargs* are
        keyword arguments.
        """
        
        fargs = list(fargs)
        
        self.quiet = False
        self.debug = False
        
        if 'quiet' in fargs:
            self.quiet = True
            fargs.remove('quiet')
        
        if 'debug' in fargs:
            self.debug = True
            fargs.remove('debug')
        
        kwargs = self._process_arguments(*fargs, **kwargs)
        
        self._run(*fargs, **kwargs)

    ####################################################################
    @abc.abstractmethod
    def _run(self, kwargs, fargs):

        pass
    
    ####################################################################

    def process_cmdline_arguments(self, arguments):
        
        """
        Processes *arguments*; returns a ``(fargs, kwargs)`` tuple.
        Don't change anything of the instance (i.e. don't set anything).
        Don't check anything.
        """
        
        # imports arguments
        kwargs = {}
        fargs = []

        for argument in arguments:
            
            c = argument.count('=')

            # flag (options without keyword)
            
            if c==0:
                
                if argument in fargs:
                    raise ValueError, 'duplicated flag: `%s`' %argument
                
                fargs.append(argument)

            # normal options KEYWORD=VALUE
            
            elif c==1:

                key, val = argument.split('=')
                
                if key in kwargs:
                    raise ValueError, 'duplicated option: `%s`' %key
                
                kwargs[key] = val
                
            else:
                raise ValueError, 'invalid argument: `%s`' %argument
                
        return (fargs, kwargs)

    ####################################################################

    def _process_arguments(self, *fargs, **kwargs):
        
        """
        Processes *arguments*. Checks that arguments correspond to 
        expected arguments and returns a dictionary with converted
        and value-checked *kwargs* and adds to it all other arguments
        (with default values).
        """
        
        for farg in fargs:
            if farg not in [name for name,doc in self.flags]:
                raise ValueError, 'invalid flag: `%s`' %farg

        for key, val in kwargs.iteritems():

            try:
                option = [option for option in self.options
                                                if option.name==key][0]
            except IndexError:
                raise ValueError, 'invalid option: `%s`' %key

            try:
                val = option.convert(val)
                
            except (ValueError, AttributeError):
                raise ValueError, 'invalid value for option `%s`: `%s`' %(key, val)

            for requirement in option.requirements:
                if not requirement(val):
                    raise ValueError, 'invalid value for option `%s`: `%s`' %(key, val)
            
            kwargs[key] = val
                
                
        # applies default values to keyword arguments (options)
        
        for option in self.options:
            
            if option.name not in kwargs:
                
                if option.default==None:
                    raise IOError, 'argument `%s` must be specified' %option.name
                
                kwargs[option.name] = option.default

        return kwargs

########################################################################

commands = []

########################################################################

class abc_sample(BaseCommand):
    
    """
    Documentation
    """

    brief = 'Generates samples to fit Approximate Bayesian\
 Computation models'

    description = 'This command draws a given number of random\
 sample from the prior distribution and generates associated set of\
 summary statistics. Note that the output file is overwritten without\
 prompting.'

    options = [
        Option('dir',
               'Directory containing fasta files',
                str, '.', []),
        Option('ext',
               'Extension of files to import. If an empty string is\
 passed (as in `ext=`), only files without extension (without any dot\
 in their name) are processed',
                str, 'fas', []),
        Option('params',
               'Name of report file.',
                str, 'abc_sample.txt', []),
        Option('data',
               'name of main output file',
               str, 'abc_sample.out', []),
        Option('model',
               'Demographic model (use option `model?` for more\
 information) (must be specified) If an argument is needed, it must be\
 given as in the following example: `AM:2` (for the model AM)',
               str, '', []),
        Option('prior',
               'Prior distribution file (use option `prior?` for more\
 information) (must be specified)',
                str, '', []),
        Option('stats',
               'Set of summary statistics (use option `stats?` for more\
 information) (must be specified). If an argument is needed, it must be\
 given as in the following example: `SFS:4` (for the statistic set SFS)',
                str, '', []),
        Option('post',
               'Number of points to sample',
               int, 10000, [lambda x: x>0]),
        Option('seeds',
               'Seeds of the random number generator. They must be\
 given as two integer separated by a comma, as in `seeds=12345,67890`.\
 By default, the random number generator is seeded from the system\
 time',
               lambda x: map(int, x.split(',')), 0,
               [lambda x: len(x)==2]),
        Option('restart',
               'Complete an interrupted run. The arguments are read from\
 the file and all other command line arguments are ignored. The argument\
 must be the name of a `params` file (or an empty string to disable this\
 function). Note that it is currently impossible to restore the random\
 number generator status (meaning that the seeds will be lost and that\
 the new run will be based on seeds based from system time)',
               str, '', []),
        Option('add_model',
               'The name of a Python module containing a model\
 definition. Pass a module name (without dots or dashes), such as "MyModel"\
 and create a file "MyModel.py" (with a py extension in addition of the\
 module name. The class defining the model must have the same name\
 ("MyModel")',
               str, '',
               []),
        Option('max_threads',
               'Maximum number of threads to start for parallel\
 computations. The maximum number of threads is the number of CPUs\
 available. By default (max_threads=0), all CPUs are used',
               int, 0,
               [lambda x: x>=0])
    ]

    ####################################################################
    
    flags = [
        ('prior?', 'Show instructions for specifying priors'),
        ('model?', 'Show the list of available demographic models'),
        ('stats?', 'Show the list of available sets of summary stats'),
        ('force_positive',  'Forces all values drawn from priors to be >=0')
    ]

    ####################################################################
    
    def thread(self, howmany, pipe, index):
        try:
            for i in range(howmany):
                ps = self.prior[index].draw()
                ds = self.model.generate(self.config, ps, self.prior[index]._random) # ugly way to get random[i]
                self.stats.compute(ds, self.config)
                pipe.send((str(ps), ' '.join(map(str, self.stats.stats))))
        except BaseException, e:
            pipe.send([e])
        finally:    
            pipe.close()

    ####################################################################

    def _run(self, *fargs, **kwargs):
        
        ### auto-reinitializes from restart file
        
        if kwargs['restart']!= '':
            
            if not os.path.isfile(kwargs['restart']):
                raise ValueError, 'cannot restart: %s parameter file not found' %kwargs['restart']
            
            # import arguments
            
            arguments_in = [ line.strip() for line in open(kwargs['restart']) ]
            arguments = []
            nloci = None
            fnames = []
            self.config = []
            obs = None
            
            # import sample configuration data (will not read again fasta files)
            
            for i in arguments_in:

                # gets number of loci
                
                m = re.match('number_of_loci=(\d+)', i)
                if m:
                    nloci = int(m.group(1))
                    continue

                # gets a locus configuration

                m = re.match('locus:(.+)=ls:(\d+),ns:(\d+),pop_ns:(.+),outgroup:(\d+)', i)
                if m:
                    fname = m.group(1)
                    ls = int(m.group(2))
                    ns = int(m.group(3))
                    try:
                        ns_pop = map(int, m.group(4).split(','))
                    except ValueError:
                        raise IOError, 'invalid line in %s: %s' %(kwargs['restart'], i.strip())
                    nso = int(m.group(5))
                    fnames.append(fname)
                    self.config.append((ns, ns_pop, nso, ls))
                    continue

                # gets observed values

                m = re.match('observed=(.+)', i)
                if m:
                    try:
                        obs = map(float, m.group(1).split(','))
                    except ValueError:
                        raise IOError, 'invalid line in %s: %s' %(kwargs['restart'], i.strip())
                    obs = ','.join(map(str, obs))
                    continue
                    
                # gets prior
                m = re.match('prior=[\'"](.+)[\'"]', i)
                if m:
                    string = m.group(1)
                    string.replace('\r\n', '\n')
                    string.replace('\r', '\n')
                    arguments.append('prior=%%%s' %string)
                    continue
                    
                # ignores the seeds
                if re.match('seeds=(.+)', i): continue
                
                # ignores the model parameters
                if re.match('parameters=(.+)', i): continue

                # ignores the number of threads
                if re.match('number_of_threads=(.+)', i): continue

                # uses other entries are command line arguments
                arguments.append(i)

            if nloci==None: raise IOError, 'no number of loci given in %s' %kwargs['restart']
            if nloci!=len(self.config): raise IOError, 'incorrect number of loci specified in %s' %kwargs['restart']
            if len(set([i[2] for i in self.config]))!=1: raise IOError, 'inconsistent number of populations in loci specified in %s' %kwargs['restart']
            if obs==None: raise IOError, 'observed values not found in %s' %kwargs['restart']

            # process non-locus-related arguments
            
            fargs, kwargs = self.process_cmdline_arguments(arguments)
            
            self.quiet = False
            self.debug = False
            
            if 'quiet' in fargs:
                self.quiet = True
                fargs.remove('quiet')
            
            if 'debug' in fargs:
                self.debug = True
                fargs.remove('debug')
                
            kwargs = self._process_arguments(*fargs, **kwargs)
            kwargs['restart'] = ''

            if not os.path.isfile(kwargs['data']):
                raise ValueError, 'cannot restart: %s results file not found' %kwargs['data']
        
            restart = True
        
        else: restart = False
        
        ### processes the helper functionalities ###
        
        help_requested = False
        
        for arg in fargs:
            
            if arg=='prior?':
                help_requested = True
                self.help_prior()

            elif arg=='model?':
                help_requested = True
                self.help_model()

            elif arg=='stats?':
                help_requested = True
                self.help_stats()

            elif arg=='force_positive':
                pass
                
            else:
                raise ValueError, 'invalid argument for `abc_sample`: `%s`' %arg

        # if at least one helper was requested, don't actually do anything
        if help_requested:
            return None

        ### computes the number of threads ###

        try:
            np = multiprocessing.cpu_count()
        except NotImplementedError:
            nthreads = 1
        else:
            if kwargs['max_threads']!=0 and kwargs['max_threads']<np:
                nthreads = kwargs['max_threads']
            else:
                nthreads = np

        ### creates the random objects (one per thread) ###
        
        if kwargs['seeds'] == 0:
            random = [egglib_binding.Random()]
        else:
            random = [egglib_binding.Random( *kwargs['seeds'] )]
        for i in range(1, nthreads):
            random.append(egglib_binding.Random(
                    random[-1].irand(999999999), random[-1].irand(999999999)))

        ### processes prior ###

        if kwargs['prior'] == '':
            raise ValueError, '`abc_sample`: prior must be specified'
        
        self.prior = None
        
        # if file, imports the file
        if kwargs['prior'][0]=='%':
            priorstring = kwargs['prior'][1:]
            # note: sys takes the argument as `raw`
            # the following is an ugly workaround (fixme)
            priorstring = priorstring.replace(r'\n', '\n')
            priorstring = priorstring.replace(r'\t', '\t')
        
        else:
            f = open(kwargs['prior'])
            priorstring = f.read()
            f.close()

        # tries all possible prior type
        
        for Prior in fitmodel.priors:
            
            tentativeprior = map(Prior, random)
            try:
                [i.parse(priorstring) for i in tentativeprior]
                self.prior = tentativeprior
            except IOError:
                pass

        # checks that one of the prior types fit the file or string
        
        if self.prior == None:
            raise ValueError, 'invalid prior string'
            
        # bounds prior if needed
        
        if 'force_positive' in fargs:
            [i.force_positive() for i in self.prior]

        ### process model ###

        if kwargs['add_model'] != '':
            fitmodel.add_model(kwargs['add_model'])

        if kwargs['model'] == '':
            raise ValueError, '`abc_sample`: model must be specified'

        # checks if recombination is to be included and options

        if ':' in kwargs['model']:
            label, options = kwargs['model'].split(':')
            options = options.split(',')
            try:
                options = map(int, options)
            except:
                raise ValueError, '`abc_sample`: invalid options in model: `%s`' %kwargs['model']
        else:
            label = kwargs['model']
            options = []
        
        if label[-1]=='R':
            label = label[:-1]
            recomb = True
        else:
            recomb = False

        # identifies the model

        self.model = None
        for Model in fitmodel.models:
            if Model.name == label:
                self.model = Model
                break
        if self.model == None:
            raise ValueError, '`abc_sample`: invalid model: `%s`' %kwargs['model']

        options.append(recomb)
        try:
            self.model = self.model(*options)
        except TypeError:
            raise ValueError, '`abc_sample`: invalid options for model: `%s`' %kwargs['model']

        #### process stats rule ###
        
        if kwargs['stats'] == '':
            raise ValueError, '`abc_sample`: stats must be specified'
        
        if ':' in kwargs['stats']:
            label, options = kwargs['stats'].split(':')
            options = options.split(',')
            try:
                options = map(int, options)
            except:
                raise ValueError, '`abc_sample`: invalid options in stats: `%s`' %kwargs['stats']
        else:
            label = kwargs['stats']
            options = []

        self.stats = None
        for Stats in fitmodel.summstats:
            if Stats.name == label:
                self.stats = Stats
                break
        if self.stats == None:
            raise ValueError, '`abc_sample`: invalid stats: `%s`' %kwargs['stats']

        try:
            self.stats = self.stats(*options)
        except TypeError:
            raise ValueError, '`abc_sample`: invalid options for stats: `%s`' %kwargs['stats']


        ### imports data (only if not restart) ###

        if not restart:

            # checking of dir argument:
            
            if not os.path.isdir(kwargs['dir']):
                raise ValueError, '`abc_sample`: invalid directory path: `%s`' %kwargs['dir']

            # gets list of file names

            pattern = kwargs['dir']+'/*'
            if len(kwargs['ext']): pattern += '.'+kwargs['ext']
            fnames = sorted(glob.glob(pattern))
            
            # import alignments
            
            dataset = fitmodel.Dataset()
            for i in fnames:
                align = data.Align(i, groups=True)
                dataset.add(align)

            # sorts alignments
            
            dataset.sort_aligns()

            # computes  observed statistics
            
            self.config = dataset.config()
            self.stats.compute(dataset, self.config)
            obs = ','.join(map(str, self.stats.stats))
            del dataset

        # needs to draw a parameter to ensure the list of parameters if generated 
        # the exported seeds will be after this draw
        
        ps = self.prior[0].draw()
        ds = self.model.generate(self.config, ps, random[0])

        ### saves parameters ###

        params = open(kwargs['params'], 'w')

        kwargs.update({'quiet': self.quiet, 'debug': self.debug, 'number_of_threads': str(nthreads)})

        for arg, value in kwargs.iteritems():
                    
            if arg=='seeds':
                seeds = int(random[0].seed1()), int(random[0].seed2())
                value = '%d,%d' %seeds
                        

            if arg=='prior':
                value = repr(self.prior[0].str())

            if isinstance(value, bool):

                if value==True:
                    params.write(arg + '\n')
                        
            else:
                params.write( '%s=%s\n' %(arg, value) )
                
            if arg=='model':
                params.write('parameters=%s\n' %(','.join(self.model.parameters)))


        params.write( 'number_of_loci=%d\n' %len(self.config))
        for fname, cfg in zip(fnames, self.config):
            params.write( 'locus:%s=ls:%d,ns:%d,pop_ns:%s,outgroup:%d\n' %(fname, cfg[3], cfg[0], ','.join(map(str,cfg[1])),cfg[2]))
            
        params.write( 'observed=%s\n' %obs )
        params.close()

        ### imports already sampled data

        done = 0
        
        if restart:
            
            # draw a paramSample (doesnt matter since the seeds are lost anyway)            
            check = len(self.prior[0].draw())
            
            # counts points and check integrity
            f = open(kwargs['data'])
            cache = ''
            for line in f:
                try:
                    paramsi, statsi = line.split('#')
                except ValueError:
                    raise IOError, 'invalid data file (invalid line): %s' %kwargs['data']
                if len(statsi.split()) != len(obs.split(',')):
                    raise IOError, 'invalid data file (incorrect number of statistics): %s' %kwargs['data']
                if len(paramsi.split()) != check:
                    raise IOError, 'invalid data file (incorrect number of statistics): %s' %kwargs['data']
                done += 1
                cache = line
            f.close()
            if len(line) and line[-1] != '\n':
                    raise IOError, 'invalid data file (last line miss a newline): %s' %kwargs['data']

            fstream = open(kwargs['data'], 'a')

        else:
            fstream = open(kwargs['data'], 'w')

        if not self.quiet:
            updater = tools.Updater(kwargs['post']-done)
            updater.refresh('%d points sampled of %d (%d%%) - $REMAINING left' %(done, kwargs['post'],100*done/kwargs['post']))

        # defines the thread lengths
            
        min_thread_size = 100
            
        sizes = []
        cake = kwargs['post']-done
        while True:
            share = int(0.75*cake/nthreads)
            if share<=min_thread_size: break
            sizes += [share]*nthreads
            cake -= share*nthreads
        if not len(sizes):
            sizes = [min_thread_size] * nthreads

        try:
            
            # initializes and starts the threads
            
            threads = []
            pipes = []
            for i in range(nthreads):
                conn1, conn2 = multiprocessing.Pipe(False)
                pipes.append(conn1)
                thread = multiprocessing.Process(
                         target=self.thread, args=(sizes.pop(0), conn2, i))
                thread.start()
                threads.append(thread)

            # main loop

            while True:
                time.sleep(0.01)
                for i in range(nthreads):
                    res = pipes[i].poll(0.1)
                    if not res and not threads[i].is_alive():
                        if len(sizes): size = sizes.pop(0)
                        else: size = min_thread_size
                        conn1, conn2 = multiprocessing.Pipe(False)
                        pipes[i] = conn1
                        threads[i] = multiprocessing.Process(
                                     target=self.thread, args=(size, conn2, i))
                        threads[i].start()
                    else:
                        if done==kwargs['post']: break
                        X = pipes[i].recv()
                        if len(X)==1: raise X[0]
                        fstream.write('%s # %s\n' %X)
                        fstream.flush()
                        done+=1
                        if not self.quiet:
                            updater.refresh('%d points sampled of %d (%d%%) - $REMAINING left' %(
                                    done, kwargs['post'], 100*done/kwargs['post']), 1)
                if done==kwargs['post']: break            

        except KeyboardInterrupt:
            sys.exit()

        except Exception, e:
            raise

        finally:

            if not self.quiet:
                updater.wipe()
                updater.refresh('%d points sampled - $ELAPSED elapsed' %done, grain=0)
                updater.close()

            map(multiprocessing.Process.terminate, threads)


    ####################################################################
    
    def help_prior(self):
        
        # gets the list of prior types
        
        priors = [i.name for i in fitmodel.priors]
                        
        print """Prior specification for `abc_sample`

There are two ways of specifying priors: by passing the name of a file
containing a prior specification string, and by passing this string
itself. The prior specification format depends on the prior type and is
given in the documentation of the `fitmodel` module of the EggLib python
package, and examples are given later in this document. Note that the
prior type is automatically detected from the string.

<INSERT HERE>

An example of prior specification for `PriorDiscrete` is:

    0.8 0.00;0.05 0.0;0.5
    0.1 0.05;0.10 0.0;0.5
    0.1 0.00;0.05 0.5;5.0

It specifies an almost flat uniform prior from 0. to 0.1 on the first
axis and from 0. to 5.0 on the second axis, with an increased
probability for values with THETA lesser than 0.05 and ALPHA lesser than
0.5.

An example of prior specification for `PriorDumb` is:

    U(0.,0.5) E(0.1)

This prior specifies a flat uniform prior distribution for the first
parameter and an exponential distribution with mean 0.1 for the second
parameter . Note that it is also possible to write the specification for
individual parameters on separated lines.

To pass a file name, use the `prior` option normally, as in:

    egglib abc_sample prior=filename
    
To pass a raw string and avoid that it is mistaken for a file name, use
a % character as below:

    egglib abc_sample prior="%0.9 0.00;0.10"

For prior specifications that require more than one line, use the line
separator `\\n` as below:

    egglib abc_sample prior="%0.9 0.00;0.05\\n0.1 0.05;0.10"
""".replace('<INSERT HERE>', tools.wrap('Currently available prior types: ' + ', '.join(priors), 72, 0))

    ####################################################################
    
    def help_model(self):
        
        # gets the list of models
        
        print """
Demographic models (with list of parameters) for `abc_sample`:
"""

        for model in fitmodel.models:
            
            print '=' * (len(model.name)+2)
            print ' %s' %model.name
            print '=' * (len(model.name)+2)

            print
            params = ' '.join(model.parameters)
            print '-' * (len(params)+2)
            print ' %s' %params
            print '-' * (len(params)+2)

            print model.__doc__

    ####################################################################
    
    def help_stats(self):
        
        # gets the list of stats functions
        
        print """
Sets of summary statistics for `abc_sample`:
"""

        for stats in fitmodel.summstats:
            
            print '=' * (len(stats.name)+2)
            print ' %s' %stats.name
            print '=' * (len(stats.name)+2)
            print stats.__doc__

commands.append(abc_sample)

########################################################################

#class abc_fix(BaseCommand):
#    
#    """
#    Documentation
#    """
#
#    ####################################################################
#
#    brief = 'Fix the last line of a corrupted ABC sample file.'
#
#    description = """This command generates a new file and removes the \
# last line in case it is invalid (truncated)."""
#
#    options = [
#        Option('input', 'Input ABC sample file', str, None, []),
#        Option('output', 'Output ABC sample file', str, None, [])
#    ]
#        
#    ####################################################################
#
#    flags = []
#
#    ####################################################################
#
#    def _run(self, *fargs, **kwargs):
#        
#        fin = open(kwargs['input'])
#        fout = open(kwargs['output'], 'w')
#        cache = ''
#        while True:
#            line = fin.readline()
#            if line=='': break   # more fixing needed? doesnt it just remove the last line???
#            fout.write(cache)
#            cache = line
#        fout.close()
#        fin.close()
        
########################################################################

class abc_fit(BaseCommand):
    
    """
    Documentation
    """

    brief = 'Uses samples to fit models using Approximate Bayesian\
 Computation'

    description = 'Performs rejection-regression method of Beaumont et\
 al. Genetics 2002. Note: ensure that enough samples will pass the\
 tolerance threshold.'

    options = [
        Option('input',
               'Name of data file to analyze. The file must be the\
 parameter file generated by `abc_sample` (by default:\
 `abc_sample.txt`)',
                str, None, []),
        Option('tolerance',
               'Proportion of samples to include in the local region\
 (example: a value of 0.05 specifies that the 5% closest samples should\
 be used).',
                float, None, [lambda x: x>=0]),
        Option('transform',
                'Data transformation to apply. Accepted values are `none`, `log` and `tan`',
                str, 'none', [lambda x: x in set(['none', 'tan', 'log'])]),
        Option('output',
               'Name of the output file.',
               str, 'abc_fit.out', [])
    ]

    ####################################################################
    
    flags = []

    ####################################################################
    
    def _run(self, *fargs, **kwargs):
        
        # imports the data
        
        params = {}

        f = open(kwargs['input'])

        for line in f:
            if line.count('=') == 1:
                k,v = line.strip().split('=')
                params[k] = v

        if 'prior' not in params:
            raise IOError, 'invalid ABC config file: %s' %kwargs['input']
        if 'observed' not in params:
            raise IOError, 'invalid ABC config file: %s' %kwargs['input']
        if 'data' not in params:
            raise IOError, 'invalid ABC config file: %s' %kwargs['input']
        if 'parameters' not in params:
            raise IOError, 'invalid ABC config file: %s' %kwargs['input']
        observed = map(float, params['observed'].split(',')) 

        # gets number of parameters from the prior string
        
        prior = None
        for Prior in fitmodel.priors:
            xprior = Prior()
            try:
                sprior = re.search('[\'"](.+)[\'"]', params['prior']).group(1)
                sprior = sprior.replace(r'\r\n', '\n')
                sprior = sprior.replace(r'\n', '\n')
                sprior = sprior.replace(r'\r', '\n')
                sprior = sprior.replace(r'\t', '\t')
            except AttributeError:
                raise IOError, 'invalid prior in this ABC config file: %s' %kwargs['input']
            try: xprior.parse(sprior)
            except fitmodel.PriorParseError:
                continue
            else:
                prior = xprior
                break
        if prior==None: raise IOError, 'invalid prior in this ABC config file: %s' %kwargs['input']
        np = prior.number_of_params()
        
        if not self.quiet:
            print 'Number of parameters: %d' %np

        param_names = params['parameters'].split(',')
        if len(param_names)!=np:
            raise ValueError, 'abc_fit: number of parameter names given in {0} is incorrect'.format(kwargs['input'])

        # initalizes the ABC instance

        ABC = egglib_binding.ABC()
        ABC.number_of_statistics(len(observed))
        for i,v in enumerate(observed):
            ABC.obs(i,v)

        if not self.quiet:
            print 'Number of statistics: %d' %len(observed)

        # 1st step

        ABC.add_fname(params['data'], np)
        ABC.get_threshold(kwargs['tolerance'])
        
        if not self.quiet:
            for i,v in enumerate(observed):
                print 'Statistic %d: %f (%f)' %(i+1, v, ABC.sd(i))
            print 'Number of points: %d' %ABC.number_of_samples()
            print 'Euclidean distance threshold: %f' %ABC.threshold()
        
        # 2nd step (rejection)
        
        handle, tmpfile = tempfile.mkstemp()
        
        try:
            os.close(handle)
        
            acc = ABC.rejection(tmpfile, )

            if not self.quiet:
                print 'Number of accepted points: %d' %acc

            if kwargs['transform']=='none':
                transform = ABC.NONE
            if kwargs['transform']=='tan':
                transform = ABC.TAN
            if kwargs['transform']=='log':
                transform = ABC.LOG

        # 3rd step (regression)

            ABC.regression(tmpfile, kwargs['output'], transform, ' '.join(params['parameters'].split(',')))

            if not self.quiet:
                print 'Posterior file %s generated' %kwargs['output']

        finally:
            if os.path.isfile(tmpfile): os.remove(tmpfile)

commands.append(abc_fit)
            
########################################################################

class abc_compare(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Compares several models.'

    description = """The same set of summary statistics must have been\
 used during simulations. This command expects a list of config files\
 that must all present the same statistics but may have been generated\
 under different models, or models with differing constraints.\
 This command will display the proportion of accepted points from each\
 file in the console (and ignore the `quiet` arguments). Ref.: Fagundes\
 et al. PNAS 2007."""

    options = [
        Option('input', 'One or several ABC config files, separated by\
 commas when more than one.', lambda x: x.split(','), None, []),
        Option('tolerance',
               'Proportion of samples to include in the local region\
 (example: a value of 0.05 specifies that the 5% closest samples should\
 be used).',
                float, None, [lambda x: x>=0])
    ]
        
    ####################################################################

    flags = []

    ####################################################################

    def _run(self, *fargs, **kwargs):
        
        observed = set()
        models = []
        for fname in kwargs['input']:
            params = {}
            f = open(fname)
            for line in f:
                if line.count('=') == 1:
                    k,v = line.strip().split('=')
                    params[k] = v
            if 'observed' not in params:
                raise IOError, 'invalid ABC config file: %s' %fname
            observed.add( tuple(map(float, params['observed'].split(','))) )
            if 'data' not in params:
                raise IOError, 'invalid ABC config file: %s' %fname
            if 'parameters' not in params:
                raise IOError, 'invalid ABC config file: %s' %fname
            models.append((params['data'], len(params['parameters'].split(','))))

        if len(observed) != 1:
            raise IOError, 'inconsistent observed statistics over posterior files - aborting'
        observed = observed.pop()

        print 'Number of files to process: %d' %len(models)

        handle, tmpfile = tempfile.mkstemp()
        
        try:
            ABC = egglib_binding.ABC()
            ABC.number_of_statistics(len(observed))
            for i,v in enumerate(observed): ABC.obs(i,v)
            for fname,np in models: ABC.add_fname(fname, np)
            ABC.get_threshold(kwargs['tolerance'])
            
            print 'Total number of points: %d' %ABC.number_of_samples()

            accept = ABC.rejection(tmpfile, True)
            
            print 'Number of accepted points: %d' %accept
            if accept==0:
                print '    (You should increase the tolerance.)'

            else:
                print 'Model probabilities:'
                res = [0] * len(models)

                f = open(tmpfile)
                for line in f:
                    match = re.match('\[(\d+)\]', line)
                    if not match:
                        raise IOError, 'problem with ABC rejection output'
                    index = int(match.group(1))-1
                    if index <0 or index>=len(res):
                        raise IOError, 'problem with ABC rejection output'
                    res[index] += 1
                    
                T = sum(res)
                if T!=accept:
                    raise IOError, 'problem with ABC rejection output'
                for i, v in enumerate(kwargs['input']):
                    print '    %s\t%f' %(v, 1.*res[i]/T)

        finally:
            if os.path.isfile(tmpfile): os.remove(tmpfile)

commands.append(abc_compare)

########################################################################

class abc_bin(BaseCommand):
    
    """
    Documentation
    """
    
    ####################################################################
    
    class _ranges:
        
        ################################################################

        def __init__(self, string):
            
            if not isinstance(string, basestring):
                raise ValueError
            
            self._items = []
            
            if len(string):

                bits = string.split(',')
                
                for bit in bits:
                    value1, value2 = bit.split(':')
                    self._items.append((float(value1), float(value2)))
    
        ################################################################

        def __str__(self):
            
            return ','.join(['%f:%f' %(a,b) for a,b in self._items])
            
        ################################################################
        
        def __iter__(self):
            
            for x,y in self._items:
                yield x,y

        ################################################################
        
        def __len__(self):
            
            return len(self._items)

    ####################################################################

    brief = 'Binarizes a posterior distribution'

    description = 'Uses the output file of the `abc_fit` command to\
 binarize the posterior and generate a "PriorDiscrete"-compatible file.\
 The `quiet` argument is ignored.'

    options = [
        Option('input',
               'Name of data file to analyze. The file must be the\
 output file generated by `abc_fit` (by default: `abc_fit.out`)',
                str, None, []),
        Option('bins',
               'Number of categories for all parameters. If specified,\
 the argument `parambins` overwrites this argument',
                int, 12, [lambda x: x>0]),
        Option('parambins',
               'Specifies specific number of categories for one or more\
 parameters. The argument must be a list of integers (separated by\
 commas) giving the number of categories for all parameters',
                lambda x: map(int,x.split(',')), [], []),    
        Option('ranges',
               'Specifies the prior ranges for one or more parameters.\
 The argument must be a list of ranges separated by commas (such as\
 `min:max,min:max,min:max,min:max`) giving minimum and maximum value for\
 all parameters (if values lie outside of ranges, and error will be\
 generated)',
                _ranges, _ranges(''), []),               
        Option('output',
               'Name of the output file',
               str, 'abc_bin.out', [])
    ]

    ####################################################################
    
    flags = [
    ]

    ####################################################################
    
    def _run(self, *fargs, **kwargs):
        
        params, data = fitmodel.import_posterior(kwargs['input'])
        
        # creates the binner

        bin = tools.Bin(data)

        # defines numbers of categories

        ncat = [kwargs['bins']] * len(data)

        if len(kwargs['parambins'])!=0:
            if len(kwargs['parambins'])!=len(ncat):
                raise IOError, 'invalid number of items in parambins argument'
            ncat = kwargs['parambins']

        # defines the specified ranges

        ranges = [(None,None)] * len(data)
        
        if len(kwargs['ranges'])!=0:
            if len(kwargs['ranges'])!=len(ranges):
                raise IOError, 'invalid number of items in ranges argument'

            for i,(value1,value2) in enumerate(kwargs['ranges']):
                if value2<=value1:
                    raise ValueError, 'invalid range for parameter %d' %(i+1)

                if value1 > min(data[i]) or value2 < max(data[i]):
                    raise ValueError, 'invalid range for parameter %d [empirical range: %f->%f]' %(i+1, min(data[i]), max(data[i]))

                ranges[i] = (value1, value2)
            
        # performs binarization
        
        posterior = [bin]
        
        for i in range(len(data)):

            if not self.quiet:
                print '-- binarizing dimension %d --' %(i+1)
    
            posterior = [
                tools.Bin.slice(bin, i, ncat=ncat[i], bot=ranges[i][0],
                                top=ranges[i][1]) for bin in posterior ]

            posterior = reduce(lambda a,b: a+b, posterior, [])
        

        # prepares the prior formalization
        
        posterior.sort(lambda x,y: cmp(len(y), len(x)))
        prior = fitmodel.PriorDiscrete()
        
        N = sum(map(len, posterior))
        
        acc = 0
        c = 0
        
        for bin in posterior:
        
            n = 1. * len(bin) / N
            acc += n
            if n!=0: c+=1

            prior.add(n, *bin.ranges())
            
        if not self.quiet:
            print 'total frequency: %f (expects 1.)' %acc
            print 'total number of bins: %d' %len(prior)
            print 'number of non-empty bins: %d' %c

        # saves the prior
        
        f = open(kwargs['output'], 'w')
        f.write(str(prior.str()))
        f.close()

commands.append(abc_bin)

########################################################################

class abc_statsmarg(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Marginal properties of a posterior distribution'

    description = 'Uses the output file of the `abc_fit` command and\
 computes properties of the marginal distribution of each parameter.\
 Results are displayed in the console. The argument `quiet` is ignored.'

    options = [
        Option('input',
               'Name of data file to analyze. The file must be the\
 output file generated by `abc_fit` (by default: `abc_fit.out`)',
                str, None, [])
    ]

    ####################################################################
    
    flags = [
    ]

    ####################################################################
    
    def _run(self, *fargs, **kwargs):
        
        # imports the data
        
        params, data = fitmodel.import_posterior(kwargs['input'])

        print 'Analyzing data in %s' %kwargs['input']

        print '------------------------'
        for i in range(len(data)):
            print params[i]
            print '------------------------'
            print '  observations: %d' %len(data[i])
            m, V, sd, se = tools.stats(data[i])
            print '       average: %f' %m
            print ' std deviation: %f' %sd

            values = sorted(data[i])

            def quantile(p):
                return values[int(math.ceil(p*len(values)))-1]
            
            print '--- quantiles ----------'
            print '           min: %f' %values[0]
            print '            1%%: %f' %quantile(0.01)
            print '            5%%: %f' %quantile(0.05)
            print '           10%%: %f' %quantile(0.10)
            print '           50%%: %f' %quantile(0.50)
            print '           90%%: %f' %quantile(0.90)
            print '           95%%: %f' %quantile(0.95)
            print '           99%%: %f' %quantile(0.99)
            print '           max: %f' %values[-1]
            print '------------------------'
        
        if len(data)>1:
        
            print 'Parameter correlation'
            print '------------------------'

            for i in range(len(data)):
                for j in range(i+1, len(data)):

                    print '%s - %s ' %(params[i],params[j])
                    r, r2, a = tools.correl(data[i], data[j])
                    print '     r : %f' %r
                    print '     r2: %f' %r2
                    print '     a : %f' %a
                    print '------------------------'

commands.append(abc_statsmarg)

########################################################################

class abc_statsdisc(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Properties of a discretized posterior distribution'

    description = 'The posterior must have been discretized using the\
 command `abc_bin`. The joint properties of distribution are computed\
 and displayed in the console. The argument `quiet` is ignored.'

    options = [
        Option('input',
               'Name of data file to analyze. The file must be of the\
 `PriorDiscrete` form as the output file of `abc_bin` (by default:\
 `abc_bin.out`), but any `PriorDiscrete` data is supported',
                str, None, []),
        Option('q',
               'Which credible interval to output (by default, the\
 bounds of the 95% density set are presented)',
                float, 0.95, [lambda x:x>=0, lambda x:x<=1])
    ]

    ####################################################################
    
    flags = [
    ]

    ####################################################################
    
    def _run(self, *fargs, **kwargs):
        
        # imports the data
        
        data = fitmodel.PriorDiscrete()
        
        f = open(kwargs['input'])
        string = f.read()
        f.close()

        try:
            data.parse(string)
            
        except IOError:
            raise IOError, 'cannot import discrete data from `%s`' %kwargs['input']

        del string

        # computes the number of non-empty classes and of points
        
        P = 0
        K = 0
        
        for p, bounds in data:
            P+=p
            if p>0:
                K+=1

        # output main statistics
        
        print 'Analyzing data in %s' %kwargs['input']
        
        print 'Sum of class frequencies: %f' %P
        print 'Number of categories: %d' %len(data)
        print 'Number of non-empty categories: %d' %K
        
        # finds the best category (MAP)
        # doesn't assume the prior is sorted
        # takes care that no ex-aequo

        best = 0., [None]
        for p, bounds in data:
            if p==best[0]:
                best[1].append(bounds)
            if p>best[0]:
                best = p, [bounds]
            
        if best[0]==0:
            raise ValueError, 'best posterior category has 0 frequency'
        
        if len(best[1])>1:
            raise ValueError, 'several categories of the posterior have the maximum frequency %f' %best[0]

        print 'Midpoint of the best category (MAP estimates):'

        for i,(b,t) in enumerate(best[1][0]):
            m = b + 0.5 * (t-b)
            print '    parameter %d: %f' %(i+1, m)

        # finds the bounds of the q density set
        
        bins = sorted([i for i in data], lambda a,b: cmp(b[0], a[0]))

        acc = 0
        lims = []
        for i in bins[0][1]:
            lims.append(list(i))
        
        for p, bounds in bins:

            for i,(bot,top) in enumerate(bounds):
    
                if bot < lims[i][0]:
                    lims[i][0] = bot

                if top > lims[i][1]:
                    lims[i][1] = top
                    
            acc += p
            if acc >= kwargs['q'] * P:
                break
        
        print 'Limits of the %f credible interval:' %kwargs['q']
        for i,(bot,top) in enumerate(lims):
            print '    parameter %d: %f, %f' %(i+1, bot, top)
            
commands.append(abc_statsdisc)

########################################################################

class abc_plot1D(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Plots marginal distributions from a discretized posterior'

    description = 'The posterior must have been discretized using the\
 command `abc_bin`. The command will plot the marginal distribution of\
 of either one (specified) or all parameters as png (portable network\
 graphics) files. The graphics will be histogram, when the class limits\
 will be fully defined by the discretization step accomplished\
 previously. The Python module matplotlib is needed to use this command.'

    options = [
        Option('input',
               'Name of data file to analyze. The file must be of the\
 `PriorDiscrete` form as the output file of `abc_bin` (by default:\
 `abc_bin.out`), but any `PriorDiscrete` data is supported',
                str, None, []),
        Option('index',
               'Which parameter to plot. By default (and with an empty\
 argument), all parameters are plotted',
                lambda x: int(x), '', [lambda x: x>0]),
        Option('params',
               'Name(s) to use in the graphic axes. It must match the\
 number of parameter(s) to be plotted. By default, the index of\
 parameters will be used',
                lambda x: x.split(','), []),
        Option('root',
               'Root name of output files. The template output file\
 name is `<root>_<param>.png` where <root> is the value of this argument\
 and <param> is the name of the parameter being plotted',
                str, 'abc_plot', [])
    ]

    ####################################################################
    
    flags = [
    ]

    ####################################################################
    
    def _run(self, *fargs, **kwargs):
        
        # checks
        
        try:
            from matplotlib import pyplot
        
        except ImportError:
            raise ImportError, 'this command requires that your Python installation contains the matplotlib module'
        
        # imports the data
        
        posterior = fitmodel.PriorDiscrete()
        
        f = open(kwargs['input'])
        string = f.read()
        f.close()

        try:
            posterior.parse(string)
            
        except IOError:
            raise IOError, 'plot1D: cannot import discrete data from `%s`' %kwargs['input']

        del string
        
        # determines which parameter(s) to plot

        if kwargs['index'] != '':

            if kwargs['index'] > posterior.number_of_params():
                raise ValueError, 'abc_plot1D: invalid parameter index: %d (max: %d)' %(kwargs['index'], posterior.number_of_params())

            indices = [kwargs['index']-1]
            
        else:
            
            indices = range(posterior.number_of_params())

        # iterates

        if len(kwargs['params'])!=0:
            if len(kwargs['params'])!=len(indices):
                raise ValueError, 'abc_plot1D: invalid number of parameter labels'

        for index in indices:
            
            # collects the data

            classes = {}
            
            for p, bounds in posterior:

                if p==0.: continue

                key = bounds[index]
                
                if key not in classes:
                    classes[key] = 0.

                classes[key] += p
            
            # transforms the data to a list
            # sorts and checks that no overlap
            
            classes = [(b,p) for b,p in classes.iteritems()]
            classes.sort(lambda x,y: cmp(x[0][0], y[0][0]))
            
            for i in range(1, len(classes)):
                if classes[i-1][0][1] > classes[i][0][0]:
                    raise ValueError, 'posterior in `%s` has overlapping classes' %kwargs['input']
            
            # prepares the plot data

            left = [b[0] for b,p in classes]
            height = [p for b,p in classes]
            width = [b[1]-b[0] for b,p in classes]

            # gets label
            
            if len(kwargs['params']):
                label = kwargs['params'][index]
            else:
                label = 'Parameter %d' %(index+1)
            
            # plots and saves

            pyplot.bar(left, height, width=width, color='0.2')
            pyplot.xlabel(label)
            pyplot.ylabel('frequency')
            fname = '%s_%s.png' %(kwargs['root'], index+1)
            pyplot.savefig(fname)
            pyplot.clf()

            if not self.quiet:
                print ' picture `%s` saved' %fname
            
commands.append(abc_plot1D)

########################################################################

class abc_psimuls(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Performs posterior simulations'

    description = 'This command generates a defined number of a user-\
 defined list of statistics for one locus. A different set of parameter\
 values is randomly drawn for each repetition. Simulations are\
 conditioned on the number(s) of sequences and alignment length(s)\
 passed as arguments. The command generates a comma-separated table\
 without header that is displayed in the console. `None` denote\
 unavailable statistics (such as those that are undefined because of the\
 lack of polymorphism). The argument `quiet` is ignored.'

    options = [
        Option('model',
               'Model to use for simulation. This argument corresponds\
 to the model specification in the `abc_sample` command',
               str, None, []),
        Option('prior',
               'Distribution of parameters. This argument corresponds to\
 the prior specification in the `abc_sample` command. Note that\
 binarized posterior files generated by the `abc_bin` command are\
 compatible.',
                str, None, []),
        Option('ns',
               'Sample configuration: gives the number of sequence\
 sampled in one or more subpopulations. Each value must be an integer\
 and, when more than one, values must be separated by commas. Each locus\
 must contain at least two samples (in any subpopulation)',
               lambda x: map(int, x.split(',')), None,
               [lambda x: sum(x)>=2]),
        Option('ls',
               'Sample configuration: gives the number of sites to\
 simulate. The argument must be an integer',
               lambda x: int(x), None, [lambda x: x>=1]),
        Option('nrepets',
               'Number of repetitions to perform',
               int, None, [lambda x: x>0]),
        Option('stats',
               'Labels of the statistics to compute. The statistic names\
 correspond to the arguments of the EggLib function `polymorphism` (note\
 that some statistics are only available when more than one population\
 is defined and/or when EggLib\'s core was linked to the Bio++ libraries\
 at compile-time). The statistics are printed to the console in the\
 order given by this option, one line per simulation',
                lambda x: x.split(','), None, []),
        Option('seeds',
               'Seeds of the random number generator. They must be\
 given as two integer separated by a comma, as in `seeds=12345,67890`.\
 By default, the random number generator is seeded from the system\
 time',
               lambda x: map(int, x.split(',')), 0,
               [lambda x: len(x)==2]),
        Option('add_model',
               'The name of a file containing a model definition',
               str, '',
               [])
    ]

    ####################################################################
    
    flags = [
    ]

    ####################################################################
    
    def _run(self, *fargs, **kwargs):

        # defines the random object

        if kwargs['seeds'] == 0:
            random = egglib_binding.Random()
            
        else:
            random = egglib_binding.Random( *kwargs['seeds'] )

        ### processes prior ###

        if kwargs['prior'] == '':
            raise ValueError, '`abc_sample`: prior must be specified'
        
        prior = None
        
        # if file, imports the file
        if kwargs['prior'][0]=='%':
            priorstring = kwargs['prior'][1:]
            # note: sys takes the argument as `raw`
            # the following is an ugly workaround (fixme)
            priorstring = priorstring.replace(r'\n', '\n')
            priorstring = priorstring.replace(r'\t', '\t')
        
        else:
            f = open(kwargs['prior'])
            priorstring = f.read()
            f.close()

        # tries all possible prior type
        
        for Prior in fitmodel.priors:
            
            tentativeprior = Prior(random)
            try:
                tentativeprior.parse(priorstring)
                prior = tentativeprior
            except IOError:
                pass

        # checks that one of the prior types fit the file or string
        
        if prior == None:
            raise ValueError, 'invalid prior string'

        # checks that one of the prior types fit the file or string
        if prior == None:
            raise ValueError, 'invalid prior string'

        ### process model ###
        
        if kwargs['add_model'] != '':
            fitmodel.add_model(kwargs['add_model'])

        if kwargs['model'] == '':
            raise ValueError, '`abc_psimuls`: model must be specified'

        # checks if recombination is to be included and options

        if ':' in kwargs['model']:
            label, options = kwargs['model'].split(':')
            options = options.split(',')
            try:
                options = map(int, options)
            except:
                raise ValueError, '`abc_psimuls`: invalid options in model: `%s`' %kwargs['model']
        else:
            label = kwargs['model']
            options = []
        
        if label[-1]=='R':
            label = label[:-1]
            recomb = True
        else:
            recomb = False

        # identifies the model

        model = None
        for Model in fitmodel.models:
            if Model.name == label:
                model = Model
                break
        if model == None:
            raise ValueError, '`abc_psimuls`: invalid model: `%s`' %kwargs['model']

        options.append(recomb)
        try:
            model = model(*options)
        except TypeError:
            raise ValueError, '`abc_psimuls`: invalid options for model: `%s`' %kwargs['model']
                           
        # makes a pseudodataset of the configuration list

        config = [(sum(kwargs['ns']), kwargs['ns'], 0, kwargs['ls'])]
                        
        # performs the required number of simulations

        for i in range(kwargs['nrepets']):

            # generates an alignment

            dataset = model.generate( config, prior.draw(), random )

            data,bidon = list(dataset.iterator(config))[0]

            # adds the group labels
            c=0
            for p, ns in enumerate(kwargs['ns']):
                for i in range(ns):
                    data.group(c, p)
                    c+=1

            if len(data)==0:
                stats = [None] * len(kwargs['stats'])
                
            else:
                pol = data.polymorphism()
                        
                stats = []
                for stat in kwargs['stats']:
                    
                    if stat not in pol:
                        raise ValueError, 'abc_psimuls: invalid statistic label: %s' %stat

                    X = pol[stat]
                    if stat in ('thetaW', 'Pi') and pol['S']:
                        X = X*pol['S']/kwargs['ls']
                    stats.append(X)

            stats = map(str, stats)
            print ','.join(stats)
            
commands.append(abc_psimuls)

########################################################################

class abc_plot2D(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Plots discretized posterior on a two-dimensional plan'

    description = 'The posterior must have been discretized using the\
 command `abc_bin`. The command will plot the (marginal) distribution of\
 of two specified parameters as a png (portable network\ graphics) file.\
 The graphics will be a two-dimensional density plot, where the class\
 limits will be fully defined by the discretization step accomplished\
 previously. The distribution should be called `marginal` if the model\
 has more than two parameters and will be the full posterior\
 distribution (with all information visible) if the model has two\
 parameters.'

    options = [
        Option('input',
               'Name of data file to analyze. The file must be of the\
 `PriorDiscrete` form as the output file of `abc_bin` (by default:\
 `abc_bin.out`), but any `PriorDiscrete` data is supported',
                str, None, []),
        Option('index1',
               'Index of the parameter to plot on the first axis',
                int, None, [lambda x:x>0]),
        Option('index2',
               'Index of the parameter to plot on the second axis',
                int, None, [lambda x:x>0]),
        Option('param1',
               'Name of the parameter to use as first axis label',
                str, '', []),
        Option('param2',
               'Name of the parameter to use as second axis label',
                str, '', []),
        Option('output',
               'Name of the output file. The default corresponds to\
 `abc_plot_PARAM1-PARAM2.png`',
                str, '', [])
    ]

    ####################################################################
    
    flags = [
        ('CI', 'displays the 95% credible interval as colored region'),
    ]

    ####################################################################
    
    def _run(self, *fargs, **kwargs):
        
        # checks
        
        try:
            from matplotlib import pyplot
        
        except ImportError:
            raise ImportError, 'this command requires that your Python installation contains the matplotlib module'
        
        # imports the data
        
        posterior = fitmodel.PriorDiscrete()
        
        f = open(kwargs['input'])
        string = f.read()
        f.close()

        try:
            posterior.parse(string)
            
        except IOError:
            raise IOError, 'cannot import discrete data from `%s`' %kwargs['input']

        del string
        
        # determines which parameter(s) to plot

        if kwargs['index1'] > posterior.number_of_params():
            raise ValueError, 'abc_plot2D: invalid parameter index: %d (max: %d)' %(kwargs['index1'], posterior.number_of_params())
        index1 = kwargs['index1']-1

        if kwargs['index2'] > posterior.number_of_params():
            raise ValueError, 'abc_plot2D: invalid parameter index: %d (max: %d)' %(kwargs['index2'], posterior.number_of_params())
        index2 = kwargs['index2']-1

        # collects the data

        classes = {}
            
        for p, bounds in posterior:
            
            key =  ((bounds[index1], bounds[index2]))
                
            if key not in classes:
                classes[key] = 0.

            classes[key] += p

        # extracts limits and checks that joint classes

        xbins = sorted(set([i[0] for i in classes]), lambda x,y: cmp(x[0],y[0]))
        for i in range(1, len(xbins)):
            if xbins[i-1][1]!=xbins[i][0]:
                raise ValueError, 'missing categories in `%s`' %kwargs['input']
            
        ybins = sorted(set([i[1] for i in classes]), lambda x,y: cmp(x[0],y[0]))
        for i in range(1, len(ybins)):
            if ybins[i-1][1]!=ybins[i][0]:
                raise ValueError, 'missing categories in `%s`' %kwargs['input']

        ## creates matrix
        # matrix = [[ classes[(x,y)] for y in ybins ] for x in xbins]
        
        # determines color bins

        if 'CI' in fargs:
            colors = zip([ 0.95, 0.99,  1.01],
                         [ 'r',  'y',   'w'])
        else:
            colors = zip([ 0.05, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99,  1.01],
                         ['0.',  '.2', '.4', '.6', '.8', '.9', '.95', '1.0'])
        
        lst = [i for i in classes.items()]
        lst.sort(lambda x,y: cmp(y[1], x[1]))
        total = sum([v for k,v in lst])
        lims = [(i*total,j) for
                        i,j in colors]
        acc = 0
        for (k,v) in lst:
            acc+= (1.*v/total)
            if acc >= lims[0][0]:
                del lims[0]
            classes[k] = lims[0][1]

        # plots

        for x1,x2 in xbins:
            for y1,y2 in ybins:
                
                k=0
                color = classes[((x1,x2),(y1,y2))]
                pyplot.broken_barh([(x1, x2-x1)], (y1, y2-y1), color=color)

        # saves

        if len(kwargs['param1']):
            label1 = kwargs['param1']
        else:
            label1 = 'Parameter %d' %(index1+1)

        if len(kwargs['param2']):
            label2 = kwargs['param2']
        else:
            label2 = 'Parameter %d' %(index2+1)


        if kwargs['output'] == '':
            output = 'abc_plot_%s_%s.png' %(index1+1, index2+1)
        else:
            output = kwargs['output']

        pyplot.xlabel(label1)
        pyplot.ylabel(label2)
        pyplot.savefig(output)

commands.append(abc_plot2D)

########################################################################

class analyzer(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Extended port of samplestat'

    description = 'This command reads `ms` output and computes several\
 statistics. Results are presented in the console in a format similar to\
 the `samplestat` output, one simulation per line (although the number,\
 and idenity and order of statistics are different). To analyse data\
 from standard input with default options, you have to type:\
 `egglib analyzer input=`. This command always displays results in the\
 standard output stream; the `quiet` option is ignored'

    options = [
        Option('input',
               'Name of the `ms` output file to read. By default (empty\
 string), data are read from standard input)',
                str, '', ''),
        Option('config',
               'Sample configuration. In case of a structured sample,\
 this option gives the number of samples from each population, each\
 separated by a comma, as  in `config=20,20,18`. For a unique and\
 non-subdivised population, a single integer should be passed',
                lambda x: map(int, x.split(',')), None,
                [lambda x: min(x)>=0]),
        Option('mis',
               'Misorientation rate (if >0, reverse randomly the\
 assignation ancestral/derived with the probability)',
                float, 0., [lambda x: x>=0, lambda x: x<=1]),
        Option('stats',
               'Specifies the list of stats (and the order) to compute.\
 The list must be comma-separated and contain only names of valid\
 statistics that can be computed from the `ms` data passed. Still,\
 invalid statistic will be silently skipped. Refer to the documentation\
 of EggLib\'s `Align.polymorphism` method for details about the\
 statistics. By default, a pre-defined list of statistics is used',
                lambda x: x.split(','), [], [])
    ]

    ####################################################################
    
    flags = [
    ]

    ####################################################################
    
    def _run(self, *fargs, **kwargs):
        
        # defines the list of statistics
        
        if kwargs['stats'] != []:

            stats = kwargs['stats']

        else:

            stats= ['D', 'tHnew_gene', 'tW_gene', 'H', 'Z', 'He', 'K',
                    'tH_gene', 'S', 'Pi_gene', 'Fst', 'Gst', 'Hst',
                    'Kst', 'Snn']
        
        # imports the data
        
        if kwargs['input'] == '':

            string = sys.stdin.read()
        
        else:
            
            f = open(kwargs['input'])
            string = f.read()
            f.close()

        # parses the individual simulations

        items = string.split('//')

        for item in items[1:]:
            
            # converts to Align
            
            item = '//\n' + item.lstrip()
            dataMatrix = egglib_binding.Ms.get(item, sum(kwargs['config']))
            align = egglib_binding.Convert.align(dataMatrix)
            Align = data.Align()
            Align._object = align

            # applies group labels
            
            if sum(kwargs['config']) != len(Align):
                raise ValueError, 'sum of `config` argument doesn\'t match with the lenght of simulations'
                
            acc = 0
            for i, ns in enumerate(kwargs['config']):
                Align.group(acc, i)
                acc+=1

            # adds fake outgroup
            
            Align.append('outgroup', 'A'*Align.ls(), 999)

            # computes statistics

            pol = Align.polymorphism()
            results = []
            for stat in stats:
                if stat in pol:
                    results.append((stat, str(pol[stat])))
            
            # writes to standard output
            
            print ' '.join(['%s: %s' %i for i in results])

commands.append(analyzer)

########################################################################

class blastgb(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Blasts all coding sequences from a GenBank file.'

    description = 'Imports a GenBank record and performs a BLAST search\
 of all `CDS` features against a given (local) database. Generates\
 another GenBank record with the name of the best hit(s) (separated by\
 the // string when more than one) appended to the `note` field.'

    options = [
        Option('input', 'Name of a GenBank file', str, None, []),
        Option('output', 'Name of the output file', str, None, []),
        Option('db',
               'Path of the target database. By default, the database\
 should be a fasta-formatted file of nucleotide sequences but flags\
 `prot` and `formatted` can control this',
              str, None, []),
        Option('evalue',
               'Expectaction value: expected number of random hits by\
 chance alone, depending on the database size. The default value is e^-6\
 (therefore much less - and more stringent - than `blastn`\'s default\
 value which is 10)',
              float, math.exp(-6), [lambda x: x>=0]),
        Option('nresults',
               'Maximum number of hits to output',
               int, 1, [lambda x: x>=0])
    ]

    ####################################################################
    
    flags = [
        ('prot', 'Performs protein-against-protein BLAST searches. With\
 this flag activated, the database passed through `db` must contain\
 protein sequences'),
        ('formatted', 'Pass this flag is the file named by the `db`\
 option is a pre-formatted BLAST database (using the `formatdb` command)\
 rather than a fasta file. In this case, the base name of the database\
 should be passed')
    ]

    ####################################################################
    
    def _run(self, *fargs, **kwargs):
        
        # import the genbank file

        genBank = data.GenBank(kwargs['input'])

        # prepares the blast database and blast objects
        
        if 'prot' in fargs:
            mode = 'prot'
        else:
            mode = 'nucl'

        if 'formatted' not in fargs:
            db = wrappers.BLASTdb(data.Container(kwargs['db']), mode)
        else:
            db = kwargs['db']

        blast = wrappers.BLAST()

        # processes all CDS features
        
        features = [
            feature for feature in genBank if feature.type() == 'CDS' ]

        if not self.quiet:
            updater = tools.Updater(len(features))

        for i,feature in enumerate(features):
                                
            sequence = feature.get_sequence()

            # translates if necessary

            if 'prot' in fargs:
                sequence = tools.translate(sequence)
            
            # makes blast
            
            if 'prot' in fargs:
                results = blast.blastp(sequence, db, evalue=kwargs['evalue'])

            else:
                results = blast.blastn(sequence, db, evalue=kwargs['evalue'])

            # crops results
            
            results = results.values()[0][:kwargs['nresults']]

            if not len(results):
                results = '"No homology found"'
            
            else:
                results = '"%s"' %' // '.join([j['subject'] for j in results])

            # adds results to the "note" field
            
            feature.add_qualifier('note', results)
        
            if not self.quiet:
                updater.refresh(
                        'feature %d of %d remaining: $REMAINING' %(
                        i+1, len(features)))

        updater.refresh('%d features processed' %len(features), grain=0)
        updater.close()

        # saves the output file
    
        genBank.write( kwargs['output'] )

commands.append(blastgb)

########################################################################

class clean_seq(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Removes ambiguity characters from nucleotide sequences.'

    description = 'The `quiet` option is ignored.'

    options = [
        Option('input', 'Name of a the input fasta file', str, None, []),
        Option('output', 'Name of the output file', str, None, []),
        Option('chars',
                'A string listing all valid characters. Note that the\
 comparisons are case-insensitive.',
              str, 'ACGTN-', [])
    ]

    ####################################################################
    
    flags = [
    ]

    ####################################################################
    
    def _run(self, *fargs, **kwargs):

        self._chars = kwargs['chars'].upper()
            
        container = data.Container(kwargs['input'])
        
        for i, (n,s,g) in enumerate(container):
            container.sequence(i, self._filter(s))

        container.write(kwargs['output'])

    ####################################################################

    def _filter(self, seq):

        result = ''

        for i in seq:

            if i.upper() in self._chars:
                result += i

            else:
                result += 'N'

        return result

commands.append(clean_seq)

########################################################################

class clean_tree(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Cleans a newick tree.'

    description = 'This command removes internal branch labels and/or\
 branch lengths from a newick tree. The `quiet` option is ignored.'

    options = [
        Option('input', 'Name of an input newick file', str, None, []),
        Option('output', 'Name of the output file', str, None, [])
    ]

    ####################################################################
    
    flags = [
        ('keep_labels', 'Don\'t remove internal branch labels'),
        ('keep_brlens', 'Don\'t remove branch lengths')
    ]

    ####################################################################
    
    def _run(self, *fargs, **kwargs):

        data.Tree(kwargs['input']).write(
                                fname=kwargs['output'],
                                labels=('keep_labels' in fargs),
                                brlens=('keep_brlens' in fargs))

commands.append(clean_tree)

########################################################################

class codalign(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Protein-based alignment of coding sequences.'

    description = 'This command accepts codings (nucleotide) sequences\
 and perform alignment at the protein level (or accept the corresponding\
 protein alignment) and generates a coding sequence alignment guaranteed\
 to fit the reading frame (gaps are multiple of three and don\'t split\
 codons apart). Note that errors can be  generated by the presence of\
 stop codons in sequences. By default, this command crops the final stop\
 codon of coding sequences. Use the `keepstop` flag to prevent this.'

    options = [
        Option('input', 'Name of an input fasta file', str, None, []),
        Option('output', 'Name of the output file', str, None, []),
        Option('prot',
               'Name of a fasta file containing aligned proteins. The\
 proteins sequences should match exactly the conceptual traduction of\
 coding sequences. If an empty string is passed (the default), the\
 option is ignored and the alignment is performed automatically on\
 conceptual translations', str, '', [])
    ]

    ####################################################################
    
    flags = [
        ('muscle', 'Uses the program `muscle` (default is `clustalw`)'),
        ('keepstop', 'Don\'t crop final stop codon of coding sequences')
    ]

    ####################################################################
    
    def _run(self, *fargs, **kwargs):

        # imports cds sequences
        cds = data.Container(kwargs['input'])
        
        # uses clustal-friendly names
        mapping = cds.encode()
        
        # removes trailing stop codons
        if 'keepstop' not in fargs:
            for i in range(len(cds)):
                if cds.sequence(i)[-3:].upper() in ('TGA', 'TAA', 'TAG'):
                    cds.sequence(i, cds.sequence(i)[:-3])
                
        # translates and aligns

        if kwargs['prot']=='':
            prot = tools.translate(cds)
            
            if 'muscle' not in fargs:
                prot = wrappers.clustal(prot, quiet=self.quiet)

            else:
                prot = wrappers.muscle(prot, quiet=self.quiet)
            
            prot.rename(mapping)
        
        # or imports already-aligned sequences
        else:
            prot = data.Align(kwargs['prot'])
            
        cds.rename(mapping)

        # calls backalign function
        align = tools.backalign(cds, prot)
        
        # saves
        align.write(kwargs['output'])
        
commands.append(codalign)

########################################################################

class concat(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Concatenation of sequence alignments.'

    description = 'Combines sequence information from fasta-formatted\
 sequence alignments. Sequences are concatenated when their names match\
 (either exact or partial matches), regardless of the order of\
 sequences. When sequences are missing in one of the alignment, they are\
 replaced by a stretch of missing data of appropriate length. Spacers of\
 missing data can be placed between concatenated alignments, depending\
 on option values. The `quiet` option is ignored. By default, the full\
 name of sequences is used for comparison. It is possible to restrict\
 the comparison to the beginning of the sequence (at a fixed length) or\
 using a specified separator character, but not both.'

    options = [
        Option('input',
                'A list of fasta file names. The names must be separated\
 by commands, as in `file1,file2,file3`. It is possible to use UNIX\
 wild cards (*, ~). File names might be duplicated',
                lambda x: x.split(','), None, []),
        Option('output', 'Name of the output file', str, None, []),
        Option('spacer',
               'Gives the length of stretches of missing data to be\
 introduced between concatenated alignments. If the argument is an\
 integer, the same spacer in introduced between all pairs of consecutive\
 alignments. To introduce variable-length spacers, a list of\
 comma-separated integers must be passed, and the number of values must\
 be equal to the number of alignments minus 1. By default, no spacers\
 are inserted.',
               lambda x: map(int, str(x).split(',')), 0,
               [lambda x: min(x)>=0]),
        Option('character',
               'Character to use for spacer stretches and for missing\
 segments. This argument should be changed to `X` when dealing with\
 protein sequences', str, '?', [lambda x: len(x)==1]),
        Option('sep',
               'Character to use as separator (only characters before\
 the first occurrence of the character are considered; the whole string\
 is considered if the character is not present in a sequence name)', str,
               '', [lambda x: len(x)==1]),
        Option('len',
               'Maximum number of characters to considered (the rest of\
 the string is discarded)', int, -1, [lambda x: x>0])
    ]

    ####################################################################
    
    flags = [
        ('partial', 'The comparison of sequence names is performed only\
 over the length of the shortest name, such as `anaco` and `anaconda`\
 are held as identical, and concatenated under the name `anaconda`'),
        ('case',    'Ignore case for comparison (all names are converted\
 to lower case)')
    ]

    ####################################################################
    
    def _run(self, *fargs, **kwargs):

        fnames = []
        for name in kwargs['input']:
            fnames += glob.glob(name)
        
        if kwargs['len']!=-1 and kwargs['sep']!='':
            raise ValueError, 'cannot specified both `len` and `sep` arguments'

        aligns = [data.Align(i) for i in fnames]

        for align in aligns:
            for i in align:
                if 'case' in fargs:
                    i.name = i.name.lower()
                if kwargs['len']!=-1:
                    i.name = i.name[:kwargs['len']]
                if kwargs['sep']!='':
                    i.name = i.name.split(kwargs['sep'])[0]

        spacer = kwargs['spacer']
        if spacer==0: spacer = [0]
        if len(spacer)==1: spacer *= (len(aligns)-1)  

        if not self.quiet:
            print 'Processing %d files...' %len(aligns)

        result = tools.concat(aligns=aligns, spacer=spacer,
                ch=kwargs['character'], strict='partial' not in fargs)
        
        result.write(kwargs['output'])

        if not self.quiet:
            print 'Number of input alignments: %d' %len(aligns)
            print 'Minimum number of sequences in input: %d' %min(map(len, aligns))
            print 'Maximum number of sequences in input: %d' %max(map(len, aligns))
            print 'Final number of sequences: %d' %len(result)
            print 'List of names in the final output file:\n    %s' %tools.wrap(', '.join(result.names()), 70, 4)
            print 'Minimum alignment length in input: %d' %(min(map(data.Align.ls, aligns)))
            print 'Minimum alignment length in input: %d' %(max(map(data.Align.ls, aligns)))
            print 'Final alignment length: %d' %(result.ls())
            print 'Result went to `%s`' %(kwargs['output'])


commands.append(concat)

########################################################################

class concatgb(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Concatenation of GenBank records.'

    description = 'The `quiet` option is ignored.'

    options = [
        Option('file1', 'First GenBank record', str, None, []),
        Option('file2', 'First GenBank record', str, None, []),
        Option('output', 'Name of the output file', str, None, []),
        Option('spacer',
               'Number of characters to insert between records',
               int, 0, [lambda x: x>=0]),
        Option('character',
               'Character to use for the spacer ',
               str, 'N', [lambda x: len(x)==1]),
    ]

    ####################################################################
    
    flags = [
    ]

    ####################################################################
    
    def _run(self, *fargs, **kwargs):

        # imports the input sequences

        gb1 = data.GenBank(kwargs['file1'])
        gb2 = data.GenBank(kwargs['file2'])
        
        # connects and shifts all gb2's features
        
        for feature in gb2:
            copy = feature.copy(gb1)
            copy.shift( len(gb1) + kwargs['spacer'] )
            gb1.add_feature(copy)

        # adds spacer and gb2's sequence to gb1's sequence
        
        gb1.set_sequence( gb1.get_sequence()
                + ''.join([kwargs['character']] * kwargs['spacer'])
                + gb2.get_sequence())

        # outputs
        
        gb1.write(kwargs['output'])

commands.append(concatgb)

########################################################################

class consensus(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Builds consensus of sequences with matching names.'

    description = 'From a nucleotide sequence alignment, the consensus\
 of all pairs of sequences that share the same prefix is computed, and\
 only unique names are exported. By default, names `spam_a001`,\
 `spam_b145`, as well as `spam` are considered as unique and merged.\
 The resulting sequence will be named `spam`. More information is\
 available in the documentation of the C++ class `Consensus`.'

    options = [
        Option('input', 'Nucleotide sequence alignment file', str, None, []),
        Option('output', 'File name for results', str, None, []),
        Option('separator', 'Character used to separate the common\
 prefix from variable part of sequence names', str, '_',
                [lambda x: len(x)==1]),
        Option('missing',
               'Character intepreted as missing data (always ignored)',
               str, '?', [lambda x: len(x)==1]),
        Option('inconsistency',
               'Character used to identify inconsistencies in the\
 `conservative` mode',
               str, 'Z', [lambda x: len(x)==1])
    ]

    ####################################################################
    
    flags = [
        ('conservative', 'Conservative mode of consensus: all\
 differences between two sequences are considered as inconsistencies and\
 are marked by the `Z` character (by default)')
    ]

    ####################################################################
    
    def _run(self, *fargs, **kwargs):

        align = data.Align(kwargs['input'])
        consensus = egglib_binding.Consensus()
        consensus.setDisagreement(kwargs['inconsistency'])
        consensus.setMissing(kwargs['missing'])
        result = data.Align()
        result._object = consensus.consensus(align._object, 
                    kwargs['separator'], ('conservative' not in fargs))
        
        if not self.quiet:
            print 'number of sequences: %d -> %d' %(len(align),len(result))
            print 'number of ambiguous positions: %d' %sum(consensus.ambiguousPositions())
            print 'number of inconsistent positions: %d' %sum(map(len,consensus.inconsistentPositions()))

        result.write(kwargs['output'])

commands.append(consensus)

########################################################################

class cprimers(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Finds consensus primers.'

    description = 'Generates consensus primers (degenerated if needed)\
 from a nucleotide sequence alignment. Ideally, expects a coding\
 sequence alignment as `input` from which the primers will be designed\
 and an annotated sequence as `gbin` containing the full sequence with\
 introns which will be used to select only primers contained in exons\
 and filter the primers overlapping splicing sites out. Generates three\
 output files where `output` is an optional base name passed as option:\
 `output.list.txt`, `output.pairs.txt` and `output.primers.gb`. The\
 first file contains a list of the generated primers, the second\
 contains a list of the generated pairs and the last one present the\
 reference sequence with annotations showing the position of all\
 primers.'

    options = [
        Option('input', 'Nucleotide sequence alignment file', str, None, []),
        Option('output', 'Base file name for results', str, 'cprimers', []),
        Option('gbin', 'Reference genbank file (if empty, the first\
 sequence of the alignment will be used', str, '', []),
        Option('ndeg',
               'Maximum number of degenerate positions allowed, per pair',
               int, 3, [lambda x: x>=0]),
        Option('liml',
               'Left limit of the selected region (based on the\
 reference sequence, not the alignment)', int, 1, [lambda x: x>=0]),
        Option('limr',
               'Right limit of the selected region  (based on the\
 reference sequence, not the alignment) (`-1` means the end of the\
 sequence)', int, -1, [lambda x: x>=0 or x==-1]),
        Option('clean_ends',
               'Number of clean positions (without degenerated bases)\
 at the end of primers', int, 3, [lambda x: x>=0]),
        Option('nseq',
               'Number of sequences to include (the default, 0,\
 corresponds to all)', int, 0, [lambda x: x>=0])
    ]

    ####################################################################
    
    flags = [
        ('no_check', 'Don\'t check for primer dimerization and other\
 primer pair problems')
    ]

    ####################################################################
    
    def _run(self, *fargs, **kwargs):

        # imports fasta

        align= data.Align(kwargs['input'])
        if not self.quiet:
            print '# file: '+kwargs['input']
            print '# %d sequences' %len(align)
            
        # selects subset of sequences

        if kwargs['nseq']==0:
            nseq= len(align)
        else:
            nseq = kwargs['nseq']
        align= align.slice(0, nseq)

        # generates consensus

        consensus= align.consensus()
        if not self.quiet:
            L = sum([consensus.count(i) for i in 'ACGT'])
            print '# consensus length: %d' %len(consensus)
            print '# non-ambiguous bases: %d' %L
            
        consensus = consensus.replace('-', 'N')

        # imports the reference sequence

        if kwargs['gbin']!='':
            ref= data.GenBank(kwargs['gbin'])
        else:
            ref= data.GenBank()
            sequence = align.sequence(0).replace('-','').replace('?','')
            ref.set_sequence( sequence )
            if not self.quiet:
                print '# used %s as reference sequence' %align.name(0)
                print '# warning: primers may overlap splicing sites'
            
        # generates primers

        params = {}
        params['PRIMER_MAX_NS_ACCEPTED'] = kwargs['ndeg']
        params['PRIMER_MIN_TM']= 54.
        params['PRIMER_OPT_TM']= 55.
        params['PRIMER_MAX_TM']= 55.
        primers= wrappers.Primer3(consensus, **params)
        a = primers.find_primers()
        if not self.quiet:
            print '# primer3 found %d forward and %d reverse primers' %tuple(a)
        primers.find_pairs()
        if 'no_check' not in fargs:
            primers.check_pairs()
        if not self.quiet:
            print '# pairing yields %d pairs' %len(primers.pairs())

        # cleaning ambiguity at ends

        primers.clean_primer_ends(kwargs['clean_ends'])
        a = primers.clean_pair_ends(kwargs['clean_ends'])
        if not self.quiet:
            print '# drops to %d after end cleaning' %a

        # removes primers with too many amb characters

        pairs1 = []
        for i in primers.pairs():
            concat = i['F']['seq']+i['R']['seq']
            valid = sum([concat.count(j) for j in 'ACGT'])
            amb = len(concat)-valid
            if amb <= kwargs['ndeg']:
                pairs1.append(i)
        if not self.quiet:
            print '# %d pairs with maximum %d degenerated positions' %(len(pairs1), kwargs['ndeg'])

        # locates primers

        sequence = ref.get_sequence()
        liml = kwargs['liml'] -1
        limr = kwargs['limr']
        if limr == -1:
            limr = len(sequence)-1
        else:
            limr -= 1
        fprimers={}
        rprimers={}
        pairs2 = []
        for i in pairs1:

            if i['F']['seq'] not in fprimers:
                a = tools.locate(sequence, i['F']['seq'], start=liml, stop=limr)
                fprimers[i['F']['seq']] = a
            else:
                a = fprimers[i['F']['seq']]

            if i['R']['seq'] not in rprimers:
                seq = tools.rc(i['R']['seq'])
                b = tools.locate(sequence, seq, start=liml, stop=limr)
                rprimers[i['R']['seq']] = b
            else:
                b = rprimers[i['R']['seq']]

            if a != None and b != None:
                i['startpos'] = a
                i['endpos'] = b + len(i['R']['seq']) - 1
                pairs2.append(i)

        if not self.quiet:
            print '# %d pairs could be mapped on the gene (between\
 positions %d and %d)' %(len(pairs2), liml+1, limr+1)

        # output individual primers

        s = StringIO.StringIO()
 
        pairs2.sort(lambda x,y: cmp(x['startpos'], y['startpos']))
        s.write('label\tstart\tend\tsize\tsequence\tQ\tTm\tGC%\tANY\tEND\n')

        fprimers = {}
        rprimers = {}
        for pair in pairs2:
            fprimers[ pair['F']['seq'] ] = pair
            rprimers[ pair['R']['seq'] ] = pair

        n = len(str(len(fprimers)))

        for i,k in enumerate(sorted(fprimers, lambda x,y: cmp(
                    fprimers[x]['startpos'], fprimers[y]['startpos']))):

            v = fprimers[k]

            fprimers[k]['Flabel'] = '%sF' %(str(i+1).rjust(n+1,'0'))
            s.write(fprimers[k]['Flabel'] + '\t')
            s.write('%d\t' %(v['startpos']))
            s.write(str(v['startpos']+len(v['F']['seq'])-1) + '\t')
            s.write(str(len(v['F']['seq'])) + '\t')
            s.write(v['F']['seq'] + '\t')
            s.write(str(v['F']['Q']) + '\t')
            s.write(str(v['F']['Tm']) + '\t')
            s.write(str(v['F']['GC%']) + '\t')
            s.write(str(v['F']['ANY']) + '\t')
            s.write(str(v['F']['END']) + '\n')

        n = len(str(len(fprimers)))

        for i,k in enumerate(sorted(rprimers, lambda x,y: cmp(
                        rprimers[x]['endpos'], rprimers[y]['endpos']))):

            v = rprimers[k]

            rprimers[k]['Rlabel'] = '%sR' %(str(i+1).rjust(n+1,'0'))
            s.write(rprimers[k]['Rlabel'] + '\t')
            s.write(str(v['endpos']-len(v['R']['seq'])+1) + '\t')
            s.write('%d\t' %(v['endpos']))
            s.write(str(len(v['R']['seq'])) + '\t')
            s.write(v['R']['seq'] + '\t')
            s.write(str(v['R']['Q']) + '\t')
            s.write(str(v['R']['Tm']) + '\t')
            s.write(str(v['R']['GC%']) + '\t')
            s.write(str(v['R']['ANY']) + '\t')
            s.write(str(v['R']['END']) + '\n')

        f = open(kwargs['output']+'.list.txt', 'w')
        f.write(s.getvalue())
        f.close()

        # output pairs
 
        s = StringIO.StringIO()

        s.write('forward\treverse\tstart\tend\tsize\t(mean)Q\tTm\t(mean)GC%\t(mean)ANY\t(mean)END\n')

        for i in sorted(pairs2, lambda x,y: cmp((x['F']['Q']+x['R']['Q'])/2.,
                                                (y['F']['Q']+y['R']['Q'])/2.)):

            s.write(fprimers[i['F']['seq']]['Flabel']+'\t')
            s.write(rprimers[i['R']['seq']]['Rlabel']+'\t')
            s.write(str(i['startpos'])+'\t')
            s.write(str(i['endpos'])+'\t')
            s.write(str(i['endpos']-i['startpos']+1)+'\t')
            s.write(str((i['F']['Q']+i['R']['Q'])/2.)+'\t')
            s.write(str(i['F']['Tm'])+' '+str(i['R']['Tm'])+'\t')
            s.write(str((i['F']['GC%']+i['R']['GC%'])/2.)+'\t')
            s.write(str((i['F']['ANY']+i['R']['ANY'])/2.)+'\t')
            s.write(str((i['F']['END']+i['R']['END'])/2.)+'\n')

        f = open(kwargs['output']+'.pairs.txt', 'w')
        f.write(s.getvalue())
        f.close()

        # adds the primers as features in the gb

        for i in sorted(pairs2, lambda x,y: cmp((x['F']['Q']+x['R']['Q'])/2.,
                                                (y['F']['Q']+y['R']['Q'])/2.)):

            locationF = data.GenBankFeatureLocation()
            locationF.addBaseRange(i['startpos'],
                                   i['startpos']+len(i['F']['seq'])-1)

            locationR = data.GenBankFeatureLocation()
            locationR.addBaseRange(i['endpos']-len(i['R']['seq'])+1,
                                   i['endpos'])
            locationR.setComplement()

            featureF = data.GenBankFeature(ref)
            featureF.set(
                type='primer_bind',
                location=locationF,
                label=fprimers[i['F']['seq']]['Flabel'],
                note=i['F']['seq'],
                color='2'
            )
            
            featureR = data.GenBankFeature(ref)
            featureR.set(
                type='primer_bind',
                location=locationR,
                label=rprimers[i['R']['seq']]['Rlabel'],
                note=i['R']['seq'],
                color='2'
            )

            ref.add_feature(featureF)
            ref.add_feature(featureR)
             
        # writes the gb
        ref.write(kwargs['output']+'.primers.gb')

commands.append(cprimers)

########################################################################

class extract(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Extract specified ranges of an alignment.'

    description = 'The command reads a fasta alignment, and generates\
 another fasta alignment consisting of one or several ranges of\
 positions.'

    options = [
        Option('input', 'Nucleotide sequence alignment file', str, None, []),
        Option('output', 'File name for results', str, None, []),
        Option('ranges', 'List of positions or ranges to extract,\
 separated by commas. Each item of the list can come as a unique integer\
 (for a unique position) or as an expression `start-stop` to extract the\
 positions `start` to `stop` (both included). It is possible to mix both\
 forms, as in `ranges=1-200,225,250,280,300-800` where `225`, for\
 example, is strictly equivalent to `225-225`',
        lambda x: map(lambda y: map(int, y.split('-')), x.split(',')),
            None, [])   ]

    ####################################################################
    
    flags = [ ]

    ####################################################################
    
    def _run(self, *fargs, **kwargs):

        # imports alignment

        align = data.Align(kwargs['input'])
        
        # the ranges must be converted manually
        
        positions = []

        for pos in kwargs['ranges']:

            if len(pos)==1:
                positions.append(pos[0]-1)
                bounds = pos[0]-1, pos[0]
            elif len(pos)==2:
                positions  += range(pos[0]-1, pos[1])
                bounds = pos[0]-1, pos[1]
            else:
                raise ValueError, 'invalid `ranges` argument (expect a single position or two bounds)'

            if bounds[0] < 0:
                raise ValueError, 'negative bound value in `ranges` argument'

            if bounds[1] > align.ls():
                raise ValueError, 'a bound value of `ranges` is over the end of the alignment'

        # extracts sequences

        extracted = align.extract(positions)
        
        # exports the resulting alignment
        
        extracted.write(kwargs['output'])

commands.append(extract)

########################################################################

class extract_clade(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Extracts the sequences corresponding to a tree clade.'

    description = 'The command takes a phylogenetic tree and a fasta\
 file containing the corresponding sequences (aligned or not). The\
 smallest clade containing all specified names will be extracted as\
 another fasta file. By default, clades encompassing the root (which\
 would be paraphyletic groups under the assumption that the tree is\
 rooted) are exported as well; use the flag `monophyletic` to prevent\
 this behaviour. Note that the root (or base of the tree) is never\
 returned.'

    options = [
        Option('sequences', 'Fasta file containing sequences', str, None, []),
        Option('tree', 'Newick file containing tree', str, None, []),
        Option('output', 'Name of resulting fasta file', str, None, []),
        Option('names', 'Name of a least one leaf of the tree (separated\
 by commas when more than one). The command supports lists containing\
 repeated leaves',
            lambda x: x.split(','), None, [lambda x: len(x)>=1]),
        Option('threshold', 'Minimum value the node must have as label\
 to be returned (only positive values are supproted). Nodes that have a\
 label not convertible to float and those whose label is inferior than\
 threshold are not returned. By default (-1), this criterion is not\
 applied at all (all nodes are returned when they contain the requested\
 names). This is different than 0 (then, only nodes that have a number\
 as label can be returned.',
            float, -1, [lambda x: x>=0 or x==-1]),
        Option('minimum', 'Smallest number of descending leaves a clade\
 must have to be returned. Clades with less nodes are ignored',
            int, 2, [lambda x: x>=2])  ]

    ####################################################################
    
    flags = [
        ('monophyletic', 'Consider only monophyletic clades (assuming\
        the tree is rooted'),
        ('exact', 'The clade must contain exactly (rather than `at\
 least`) the number of leaves given by the `minimum` option')
    ]

    ####################################################################
    
    def _run(self, *fargs, **kwargs):

        # imports data

        fasta = data.Container(kwargs['sequences'])
        tre = data.Tree(kwargs['tree'])
        kwargs['names']

        # finds the good node
        
        leaves = None
        
        for node in tre:
            
            # checks that label is sufficient
            
            if kwargs['threshold']!=-1:
                try:
                    label = int(node.get_label())
                except (TypeError, ValueError):
                    continue
                if label < kwargs['threshold']:
                    continue
              
            # collect the list(s) of leaves
            
            leaves1 = node.leaves_down()
            bleaves1 = list(leaves1)
            
            if 'monophyletic' in fargs:
                leaves2 = None
                bleaves2 = None
            else:
                leaves2 = node.leaves_up()
                bleaves2 = list(leaves2)
                
            # checks all the leaves are in either of the leaves list
            
            for name in kwargs['names']:

                if bleaves1 != None:
                    if name not in bleaves1:
                        bleaves1 = None
                    else:
                        bleaves1.remove(name)

                if bleaves2 != None:
                    if name not in bleaves2:
                        bleaves2 = None
                    else:
                        bleaves2.remove(name)
                        
                if bleaves1==None and bleaves2==None:
                    break
            
            # checks one or both lists are good

            if bleaves1 != None:
                if leaves == None or len(leaves1)<len(leaves):

                    if 'exact' in fargs:
                        good = len(leaves1)==kwargs['minimum']
                    else:
                        good = len(leaves1)>=kwargs['minimum']

                    if good:
                        leaves = leaves1

            if bleaves2 != None:
                if leaves == None or len(leaves2)<len(leaves):

                    if 'exact' in fargs:
                        good = len(leaves2)==kwargs['minimum']
                    else:
                        good = len(leaves2)>=kwargs['minimum']

                    if good:
                        leaves = leaves2

        # checks
        
        if leaves==None:
            raise ValueError, 'cannot find the specified node'

        # imports the sequences

        result = data.Container()

        for i in sorted(leaves):

            x = fasta.find(i)

            if x==None:
                raise ValueError, '%s not in %s' %(i, kwargs['sequences'])
                
            result.append( *fasta[x] )

        # exports data

        result.write(kwargs['output'])

        if not self.quiet:
            print '%d sequences exported' %len(result)

commands.append(extract_clade)

########################################################################

class family(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Finds homologs of a gene family using BLAST.'

    description = 'This command uses all sequences from a fasta file of\
 source sequences to blast against a database and reports (in a fasta\
 file) all sequences of the target database that produce a significant\
 hit with any of the source sequences. To use this command, you need to\
 have the NCBI BLAST+ package installed. You need a fasta file of\
 protein of nucleotide sequences. You need a target database (from which\
 sequences should be extracted) as a fasta files.'

    options = [
        Option('input', 'Fasta file containing source sequences', str, None, []),
        Option('target', 'Fasta file containing target database', str, None, []),
        Option('output', 'Name of resulting fasta file', str, None, []),
        Option('mode', 'Program to use: `blastn` for nucleotide source\
 against nucleotide database, `blastp` for protein source against\
 protein database, `blastx` for (translated) nucleotide source against\
 protein database, `tblastn` for protein source against (translated)\
 nucleotide database, `tblastx` for (translated) nucleotide source\
 against (translated) nucleotide database',
            str, 'blastn', [lambda x: x in set(['blastn', 'blastp', 
            'blastx', 'tblastn', 'tblastx'])]),
        Option('evalue', 'Maximum threshold to report hits. The\
 parameter used is E-value, that is for a given BLAST hit the\
 theoretical probability of obtaining such hit by chance alogn, given\
 the length of the database. It can be necessary to decrease this\
 parameter to obtain results',
            float, math.exp(-6), [lambda x: x>=0])  ]
    
    ##############################################################
    
    flags = [    ]

    ####################################################################
    
    def _run(self, *fargs, **kwargs):


        # configuring blast dabase

        if not self.quiet:
            print '# importing database from `%s`' %kwargs['target']
        
        target = data.Container(kwargs['target'])
        if not self.quiet:
            print '  >`%s` contains %d sequences' %(kwargs['target'], len(target))
        
        if kwargs['mode'] in set(['blastp', 'blastx']):
            type = 'prot'
        else:
            type = 'nucl'
        db = wrappers.BLASTdb(target, type)

        if not self.quiet:
            print '  >target database ready'

        # imports source data

        if not self.quiet:
            print '# importing source sequences from `%s`' %kwargs['input']

        input = data.Container(kwargs['input'])
        if not self.quiet:
            print '  >`%s` contains %d sequences' %(kwargs['input'], len(input))

        # replace all invalid characters by Ns
        
        if kwargs['mode'] in ['blastn', 'blastx', 'tblastx']:
            valid = set('ACGT')
            invalid = 'N'
        elif kwargs['mode'] in ['blastp', 'tblastn']:
            valid = set('ARNDCEQGHILKMFPSTWYV')
            invalid = 'X'
            for i in range(len(input)):
                for j in range(input.ls(i)):
                    if input.get(i,j).upper() not in valid:
                        input.set(i,j,invalid)
        else: raise ValueError, 'invalid mode: {0}'.format(kwargs['mode'])

        # blast phase

        if not self.quiet:
            print '# blasting...'
        
        blast = wrappers.BLAST()
        results = blast._search(kwargs['mode'], input, db, evalue=kwargs['evalue'])

        if not self.quiet:
            print '  >blast finished'
            print (tools.wrap('  >number of hits per sequence: '+
                    ', '.join(map(str,[len(results[n]) for n,s,g in
                    input])), 72, 3))

        hits = set()
        for i in results:
            hits.update( [j['subject'] for j in results[i]] )

        if not self.quiet:
            print '  >total number of hits: %d' %len(hits)

        # creates the final file
        
        final = data.Container()

        for i in sorted(hits):
            final.append( i, target.sequenceByName(i) )

        final.write(kwargs['output'])
        
        if not self.quiet:
            print '# %s successfully saved' %(kwargs['output'])

commands.append(family)

########################################################################

class fasta2mase(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Converts a fasta alignment to the mase format.'

    description = 'The `quiet` argument is ignored.'

    options = [
        Option('input', 'Name of a fasta-formatted alignment', str, None, []),
        Option('output', 'Name of resulting mase file', str, None, []) ]
    
    ##############################################################
    
    flags = [ ]

    ####################################################################
    
    def _run(self, *fargs, **kwargs):

        align = data.Align(kwargs['input'])
        mase = tools.Mase(align)
        f = open(kwargs['output'], 'w')
        f.write(str(mase))
        f.close()

commands.append(fasta2mase)

########################################################################

class fasta2nexus(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Converts a fasta alignment to the NEXUS format.'

    description = 'The `quiet` argument is ignored.'

    options = [
        Option('input', 'Name of a fasta-formatted alignment', str, None, []),
        Option('output', 'Name of resulting NEXUS file', str, None, []),
        Option('type', '`nucl` for nucleotides or `prot` for proteins',
                    str, None, [lambda x: x in set(['nucl', 'prot'])]) ]

    ##############################################################

    flags = [ ]

    ####################################################################
    
    def _run(self, *fargs, **kwargs):

        align = data.Align(kwargs['input'])
        f = open(kwargs['output'], 'w')
        f.write(align.nexus(kwargs['type']=='prot'))
        f.close()

commands.append(fasta2nexus)

########################################################################

class fasta2phyml(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Converts a fasta alignment to the `phyml` format.'

    description = 'The so-called `phyml` format is a modification of the\
 PHYLIP file format suitable for importing data to the programs PAML\
 and PHYML. The `quiet` argument is ignored.'

    options = [
        Option('input', 'Name of a fasta-formatted alignment', str, None, []),
        Option('output', 'Name of resulting `phyml` file', str, None, []) ]
        
    ####################################################################

    flags = [ ]

    ####################################################################
    
    def _run(self, *fargs, **kwargs):

        align = data.Align(kwargs['input'])
        f = open(kwargs['output'], 'w')
        f.write(align.phyml())
        f.close()

commands.append(fasta2phyml)

########################################################################

class fg2gb(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Generates a GenBank record from fgenesh output.'

    description = 'The command requires the sequence of the annotated\
 regions as a fasta file and the fgenesh output as a separate text file.\
 Obviously, all features must fit in the sequence length. A GenBank file\
 incorporate the information of predicted genes as `gene`, and `CDS`\
 annotation features.'

    options = [
        Option('seq', 'File with fasta-formatted sequence', str, None, []),
        Option('ann', 'File with fgenesh output', str, None, []),
        Option('output', 'Name of the resulting GenBank file', str, None, []) ]
        
    ####################################################################

    flags = [ ]

    ####################################################################
    
    def _run(self, *fargs, **kwargs):

        # imports sequence

        sequence = data.Container(kwargs['seq'])

        if len(sequence)!=1:
            raise ValueError, '%s should contain exactly one fasta-formatted sequence' %kwargs['seq']

        name, sequence, group = sequence[0]

        if not self.quiet:
            print 'imported sequence `%s`' %name
        
        # creates gb
        
        gb = data.GenBank()
        gb.title = name
        gb.set_sequence(sequence)
        
        if not self.quiet:
            print 'sequence of length %d imported' %len(sequence)
        
        # imports annotation
        
        annotations = tools.get_fgenesh(kwargs['ann'])
        
        # imports annotations
        
        cds = 0
        gene = 0
        
        for ann in annotations:

            feature = data.GenBankFeature(gb)
            location = data.GenBankFeatureLocation()
            
            for range in ann['pos']:
                location.addBaseRange(*range)
            del ann['pos']
            
            if ann['strand']=='minus':
                location.setComplement()
            del ann['strand']
            
            if ann['type'] == 'CDS':
                cds += 1
                type = 'CDS'
            elif ann['type'] == 'gene':
                gene += 1
                type = 'gene'
            else:
                raise ValueError, 'invalid feature type: `%s`' %ann['type']

            del ann['type']

            feature.set(type, location, **ann)
            gb.add_feature(feature)

        if not self.quiet:
            print '%d gene and %d CDS features imported' %(gene, cds)
        
        # export gb
        
        gb.write(kwargs['output'])
        
        if not self.quiet:
            print '`%s` created' %kwargs['output']
        
commands.append(fg2gb)

########################################################################

class gb2fas(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Converts GenBank records to fasta.'

    description = 'The command takes one or more GenBank records and\
 generates a single fasta file. Each GenBank record can be multiple\
 (contain multiple entries). Each sequence is named after the title of\
 the GenBank entry (disregarding the file name).'

    options = [
        Option('input', 'One or more GenBank file names, separated by\
 commas when more than one', lambda x: x.split(','), None, []),
        Option('output', 'Name of the fasta file to generate', str, None, []) ]
        
    ####################################################################

    flags = [ ]

    ####################################################################
    
    def _run(self, *fargs, **kwargs):

        fasta = data.Container()

        for i in kwargs['input']:

            if not self.quiet:
                print 'processing: `%s`' %i

            f = open(i)
            s = f.read()
            f.close()
            s = s.split('//')
            if len(s)<2:
                raise ValueError, 'invalid GenBank file: `%s`' %i
            for j in s[:-1]:
                gb = data.GenBank(string=j+'//')
                fasta.append(gb.title, gb.get_sequence())
        
        if not self.quiet:
            print 'number of sequences imported: %d' %len(fasta)
            
        fasta.write(kwargs['output'])
        
commands.append(gb2fas)

########################################################################

class infos(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Displays basic information from fasta files.'

    description = 'The commands displays the number of sequences and\
 alignment length (length of the longest sequences for unaligned sets\
 of sequences) for all fasta files passed. The `quiet` option is\
 ignored.'

    options = [
        Option('input', 'One or more fasta file names, separated by\
 commas when more than one', lambda x: x.split(','), None, []) ]
        
    ####################################################################

    flags = [ ]

    ####################################################################
    
    def _run(self, *fargs, **kwargs):

        for fname in kwargs['input']:

            container = data.Container(fname)

            print fname
            print ' ... %d sequence%s' %(len(container),
                                  {True:'s',False:''}[len(container)>1])
            if len(container)==0:
                print ' ... empty file'
            elif len(container)==1:
                print ' ... single sequence'
                print ' ... length: %d' %container.ls(0)
            elif container.isEqual():
                print ' ... alignment'
                print ' ... length: %d' %container.ls(0)
            else :
                print ' ... not an alignment'
                print ' ... max length: %d' %max([container.ls(i)
                                        for i in range(len(container))])
        
commands.append(infos)

########################################################################

class interLD(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Computes linkage disequilibrium statistics between two loci.'

    description = 'Computes association statistics between two\
 alignments. It is required that both sequence alignments contain the\
 exact same list of sequence names (duplicates are not supported) and\
 there should be at least four sequences in each alignment. In the\
 definition of statistics, an allele is a haplotype as determnined by\
 the method Align.polymorphism(). The frequency of allele i at one locus\
 is Pi, the frequency of the combination i,j (i at locus 1 and j at\
 locus 2) is Pij. For a given pair of alleles i,j (i at locus 1 and j at\
 locus j), Dij is Pij - PiPj. D\'ij is Dij/Dijmax if Dij>=0 and\
 Dij/Dijmin if Dij<0, where Dijmax is min(Pi(1-Pj), (1-Pi)Pj) and Dijmin\
 is min(PiPj, (1-Pi)(1-Pj)). To obtain the complete LD estimates both\
 measures are averaged over all allele pairs as Dijtot = sum(PiPj|Dij|)\
 of all i,j pairs).'

    options = [
        Option('align1', 'First alignment', str, None, []),
        Option('align2', 'Second alignment', str, None, []),
        Option('permus', 'Number of permutations to perform. If the\
 value is larger than 0, the distribution of linkage statistics is\
 computed by randomly shuffle the sequences of one of the alignments.',
 int, 0, [lambda x: x>=0]),
        Option('output', 'Name of output file', str, 'interLD.txt', [])]
        
    ####################################################################

    flags = [ ]

    ####################################################################

    def _run(self, *fargs, **kwargs):

        # imports alignments

        align1 = data.Align(kwargs['align1'], groups=True)
        align2 = data.Align(kwargs['align2'], groups=True)

        # checks that names are correct

        if align1.contains_duplicates():
            raise ValueError, 'error: `%s` contains duplicates' %kwargs['align1']

        if align2.contains_duplicates():
            raise ValueError, 'error: `%s` contains duplicates' %kwargs['align2']

        if set(align1.names()) != set(align2.names()):
            raise ValueError, 'error: sequence names from `%s` and `%s` don\'t match' %(kwargs['align1'], kwargs['align2'])
            
        if len(align1)<4:
            raise ValueError, 'error: alignments must contain at least 4 sequences'

        # computes LD
        
        ns1, ns2, S1, S2, K1, K2, D, Dp = tools.LD(align1, align2, False)

        output = open(kwargs['output'], 'w')
        output.write('\tLocus 1\tLocus2\n')
        output.write('File:\t%s\t%s\n' %(kwargs['align1'], kwargs['align2']))
        output.write('ns:\t%d\t%d\n' %(ns1, ns2))
        output.write('S:\t%d\t%d\n' %(K1, S2))
        output.write('K:\t%d\t%d\n' %(S1, K2))
        output.write('\n')
        output.write('D:\t%f\n' %D)
        output.write('D\':\t%f\n' %Dp)

        # performs the test
        
        if kwargs['permus']>0:

            P1 = 0.
            P2 = 0.
            
            if not self.quiet:
                print 'performing permutations'
                updater=tools.Updater(kwargs['permus'])

            for i in range(kwargs['permus']):
                
                a, b, c, d, e, f, Di, Dpi = tools.LD(align1, align2, True)
                if Di>=D: P1+=1
                if Dpi>=Dp: P2+=1
                
                if not self.quiet:
                    updater.refresh()

            if not self.quiet:
                updater.close()
                
            P1 /= kwargs['permus']
            P2 /= kwargs['permus']

            output.write('\nP-values from %d iterations:\n')
            output.write('D:\t%f\n' %P1)
            output.write('D\':\t%f\n' %P2)

        output.close()

commands.append(interLD)

########################################################################

class matcher(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Finds homologous regions between two sequences.'

    description = 'This command performs a `bl2seq` search using the\
 first sequence as query and the second sequence as target. It then\
 produces a genbank record containing the first sequence with\
 annotation features indicating the positions of the hits with the\
 second sequence. The *long* sequence should not contain gaps.'

    options = [
        Option('long', 'Fasta file containing the first sequence (this\
 sequence must be a nucleotide sequence)', str, None, []),
        Option('short', 'Fasta file containing the sequence sequence\
 depending on the `mode` option value, this sequence must be either a\
 nucleotide sequence or a protein sequence)', str, None, []),
        Option('output', 'Name of resulting GenBank file', str, 'matcher.gb', []),
        Option('mode', 'Program to use: `blastn` for nucleotide source\
 against nucleotide database, `blastx` for (translated) nucleotide\
 source against protein database, `tblastx` for (translated) nucleotide\
 source against (translated) nucleotide database',
            str, 'blastn', [lambda x: x in set(['blastn', 'blastx',
            'tblastn', 'tblastx'])]),
        Option('evalue',
               'Expectaction value: expected number of random hits by\
 chance alone, depending on the database size. The default value is e^-6\
 (therefore much less - and more stringent - than `blastn`\'s default\
 value which is 10)',
              float, math.exp(-6), [lambda x: x>=0]) ]
        
    ####################################################################

    def _run(self, *fargs, **kwargs):

        # imports sequences

        fas1 = data.Container(kwargs['long'], groups=True)
        fas2 = data.Container(kwargs['short'], groups=True)
        
        # checks that files look normal
        
        if len(fas1) != 1:
            raise ValueError, '%s must contain exactly one sequence' %kwargs['long']

        if len(fas2) != 1:
            raise ValueError, '%s must contain exactly one sequence' %kwargs['sort']
        
        # performs the blast search
        
        bl2seq = wrappers.BL2SEQ()
        res = bl2seq._search(kwargs['mode'], fas1.sequence(0),
                            fas2.sequence(0), evalue=kwargs['evalue'])

        # collects hit positions
        
        pos = []
        
        for i in res:
            
            loc = sorted([i['qstart'], i['qend']])

            minus =  ((i['qend']>i['qstart'] and i['send']<i['sstart']) or
                      (i['qend']<i['qstart'] and i['send']>i['sstart']))
            
            pos.append((loc,minus))
            
        pos.sort(lambda x,y: cmp(x[0][0], y[0][0]))
        
        # removes redundant hits
        
        i = 0
        
        while i<len(pos):

            flag = False

            for j in range(len(pos)):
                
                if (i!=j and pos[i][0][0]>=pos[j][0][0] 
                         and pos[i][0][1]<=pos[j][0][1]):

                    flag = True
                    break
            
            if flag:
                del pos[i]
            
            else:
                i+=1

        # small feedbacks

        if not self.quiet:
            print 'number of hits: %d' %len(pos)

        # initializes GenBank record
        
        gb = data.GenBank()
        gb.set_sequence(fas1.sequence(0))
        gb.title = fas1.name(0)

        # adds hits as features
        
        for (start,stop),minus in pos:

            loc = data.GenBankFeatureLocation()
            loc.addBaseRange(start,stop)
            if minus:
                loc.setComplement()
            
            feat = data.GenBankFeature(gb)
            feat.set('misc_feature', loc, note='bl2seq hit', color='2')
            
            gb.add_feature(feat)

        # saves
        
        gb.write(kwargs['output'])

commands.append(matcher)

########################################################################

class names(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Lists sequence names from a fasta file.'

    description = 'The order of names is preserved. The `quiet` flag is\
 ignored. By default, one name is displayed per line.'

    options = [
        Option('input', 'Name of fasta-formatted sequence file', str, None, []) ]
        
    ####################################################################

    flags = [ ('wrap', 'Displays several sequence names per line.\
 Activate this flag when sequence names are short and don\'t contain spaces!') ]

    ####################################################################

    def _run(self, *fargs, **kwargs):

        fas = data.Container(kwargs['input'])
        
        if 'wrap' in fargs:
            print tools.wrap(', '.join(fas.names()), 72, 0)

        else:
            for name in fas.names():
                print name

commands.append(names)

########################################################################

class rename(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Rename sequences according to a replacement list.'

    description = 'The replacements must be given in a text file. It is\
 not necessary to specify all names of the fasta file. The command does\
 not require either that all replacements of the list are performed. If\
 present, group labels are preserved and are not considered (they should\
 not be included in the replacement list). If leading or trailing spaces\
 are present in either old or new names, they will be removed.'

    options = [
        Option('input', 'Name of fasta-formatted sequence file', str, None, []),
        Option('list', 'Name of a text file giving the list of\
 replacements to perform. Each replacement must take one line and give\
 the old name and the new name, in that order, separated by a tabulation',
                str, None, []),
        Option('output', 'Name of output file', str, None, []) ]
        
    ####################################################################

    flags = [ ]

    ####################################################################

    def _run(self, *fargs, **kwargs):

        # imports data

        input = data.Align(kwargs['input'], groups=True)
        
        rule = {}
        
        f = open(kwargs['list'])

        for line in f:
            
            bits = line.split('\t')
            if len(bits) != 2:
                raise IOError, 'error: invalid line in %s: %s' %(kwargs['list'], line.strip())
        
            old = bits[0].strip()
            new = bits[1].strip()
        
            if bits[0] in rule:
                raise IOError, 'error: name %s is duplicated in %s' %(old, kwargs['list'])
            
            rule[old] = new
        
        f.close()
                
        # performing replacements
        
        output = data.Align()
        
        cnt = 0
        for n,s,g in input:
            
            if n in rule:
                output.append(rule[n], s, g)
                cnt+=1
                
                if not self.quiet:
                    print '%s -> %s' %(n, rule[n])
        
            else:
                output.append(n, s, g)
        
        if not self.quiet:
            print '%d replacements performed' %cnt
        
        # saves
        
        output.write(kwargs['output'])

commands.append(rename)

########################################################################

class truncate(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Truncates sequence names.'

    description = 'The user can specify a separator or a number of\
 characters (`length`) or both. By default (if neither argument\
 `separator` or `length` is specified), nothing is done. If both actions\
 are requested, they are always performed in the order: first\
 `separator`, then `length`. If present, group labels are preserved and\
 are not considered.'

    options = [
        Option('input', 'Name of fasta-formatted sequence file', str, None, []),
        Option('output', 'Name of output file', str, None, []),
        Option('separator', 'The separator can be a single character or\
 a string. Whenever it occurs in a sequence name, everything right of\
 its first occurrence (as well as the separator itself) will be\
 removed. The default (an empty string) means that this criterion is not\
 applied', str, '', []),
        Option('length', 'Maximum length of names. The default (an\
 empty string) means that this criterion is not applied',
                int, 0, [lambda x: x>=0]) ]
        
    ####################################################################

    flags = [ ]

    ####################################################################

    def _run(self, *fargs, **kwargs):

        # imports data

        input = data.Align(kwargs['input'], groups=True)
                        
        # performing replacements
        
        output = data.Align()
        
        for n,s,g in input:
            
            new = n
            
            if kwargs['separator'] != '':

                pos = n.find(kwargs['separator'])

                if pos != None:
                    new = n[:pos]

            if kwargs['length'] != 0:
                new = new[:kwargs['length']]

            if not self.quiet:
                sys.stdout.write('%s -> %s\n' %(n,new))

            output.append(new, s, g)
        
        # saves
        
        output.write(kwargs['output'])

commands.append(truncate)

########################################################################

class reroot(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Changes the orientation of a newick tree.'

    description = 'This command doesn\'t actually root (or reroot) the\
 tree; the original tree must not be rooted (it must have a trifurcation\
 at the root, and resulting tree will be likewise (only the\
 representation will be altered to present the outgroup as one of the\
 basal groups. A list of leaves representing a monophyletic group of the\
 current tree (without encompassing the root) must be passed. The\
 `quiet` argument is ignored. By default, the command uses the midpoint\
 method.'

    options = [
        Option('input', 'Name of newick-formatted tree file', str, None, []),
        Option('output', 'Name of output file', str, None, []),
        Option('outgroup', 'List of leaves constituting the outgroup,\
 separated by commas when more than one. It is possible to place the\
 list in a file (one per line) and pass the name of the file (say,\
 `fname`) using the `@` prefix, as in `outgroup=@fname` (there must be\
 exactly one item and no comma separator in that case). By default\
 (empty string) the command uses the midpoint method',
            lambda x: x.split(','), '', []) ]
        
    ####################################################################

    flags = [ ]

    ####################################################################

    def _run(self, *fargs, **kwargs):

        if len(kwargs['outgroup'])==1 and kwargs['outgroup'][0][:1] == '@':

            f = open(kwargs['outgroup'][0][1:])
            kwargs['outgroup'] = [i.strip() for i in f.readlines() if len(i.strip())]
            f.close()
            
        tree = data.Tree(kwargs['input'])
        
        if len(kwargs['outgroup']):
        
            root = tree.findMonophyleticGroup(kwargs['outgroup'])

            if root == None or len(root.ascendants()) != 1:
                raise ValueError, 'cannot reroot: invalid outgroup'

            root = root.ascendants()[0]

            tree.reoriente(root)

        else:
            
            tree.midroot()

        tree.write(kwargs['output'])

commands.append(reroot)

########################################################################

class select(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Select a given list of sequences from a fasta file.'

    description = 'The names should not include any group label (`@0`,\
 `@1`, `@999` etc. tags) is they are present in the file (group labels\
 are ignored). When a name is duplicated in the file (whether the\
 different duplicates bear different group label or not), they are all\
 exported to the output file. It is required that all names passed are\
 found at least once. Sequences are exported in the order as they appear\
 in the passed list.'

    options = [
        Option('input', 'Name of fasta-formatted file', str, None, []),
        Option('output', 'Name of the output file', str, None, []),
        Option('list', 'List of names of sequences that should be\
 selected, separated by commas when more than one. It is possible to\
 place the list of names in a file (one per line) and pass the name of\
 the file (say, `fname`) using the `@` prefix, as in `list=@fname`\
 (there must be exactly one item and no comma separator in that case)',
            lambda x: x.split(','), '', []) ]
        
    ####################################################################

    flags = [ ]

    ####################################################################

    def _run(self, *fargs, **kwargs):

        if len(kwargs['list'])==1 and kwargs['list'][0][:1] == '@':

            f = open(kwargs['list'][0][1:])
            kwargs['list'] = [i.strip() for i in f.readlines() if len(i.strip())]
            f.close()
            
        input = data.Container(kwargs['input'], groups=True)
        output = data.Container()
        
        for name in kwargs['list']:
            
            found = False
            for n,s,g in input:
                if n==name:
                    found = True
                    output.append(n,s,g)

            if not found:
                raise ValueError, 'sequence %s not found in %s' %(name, kwargs['input'])

        if not self.quiet:
            print '%d names requested, %s sequences found' %(len(kwargs['list']), len(output))

        output.write(kwargs['output'])

commands.append(select)

########################################################################

class sprimers(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Design copy-specific PCR primers from an alignment.'

    description = 'This command designs PCR primers that a specific to\
 genes from a sequence alignment (and are as unlikely as possible to\
 amplify other genes from the alignment). Primers are generated using\
 PRIMER3. Next, they are filtered according to several criteria. The\
 preferred primers must be close to the end of sequences (by default),\
 with low homology to other sequences of the alignment. A BLAST search\
 is performed and primers whose 3\' end matches any other sequence and\
 excluded. Finally, a pair check is performed using PRIMER3. The\
 corresponding programs must be available.'

    options = [
    Option    ('input', 'Name of the input fasta alignment', str, None, []),
    Option    ('output', 'Name of the output file', str, 'sprimers.csv', []),
    Option    ('sizemin', 'Mininal product size', int, 70, [lambda x: x>50]),
    Option    ('sizemax', 'Maximal product size', int, 150, [lambda x: x>50]),
    Option    ('minTm', 'Minimal annealing temperature', float, 58, [lambda x: x>20]),
    Option    ('optTm', 'Optimal annealing temperature', float, 60, [lambda x: x>20]),
    Option    ('maxTm', 'Maximal annealing temperature', float, 62, [lambda x: x>20]),
    Option    ('minGc', 'Minimal GC content percentage', float, 30, [lambda x: x>0]),
    Option    ('optGc', 'Optimal GC content percentage', float, 50, [lambda x: x>0]),
    Option    ('maxGc', 'Maximal GC content percentage', float, 80, [lambda x: x>0]),
    Option    ('numAmb', 'Maximal number of degenerate bases in primers', int, 0,[lambda x: x>=0]),
    Option    ('filter1', 'Pre-selection filter (before BLAST) as a maximal number of pairs to process', int, 5000, [lambda x: x>0]),
    Option    ('filter2', 'Pre-selection filter (after BLAST) as a maximal number of pairs to process', int, 100, [lambda x: x>0]),
    Option    ('threshold1', 'Maximum homology score to other genes (a real number between 0. and 1.)', float, 0.75, [lambda x: x<=1 and x>=0]),
    Option    ('threshold2', 'Maximum homology score to other regions of the same gene (a real number between 0. and 1.)', float, 0.50, [lambda x: x<=1 and x>=0]),
    Option    ('show', 'How many pairs to export in the output file', int, 10,[lambda x: x>0]),
    ]
        
    ####################################################################

    flags = [
        ('selection', 'Restrict the primer search to one or more\
 sequences of the alignment. The user should tag the names of selected\
 sequences with labels such as @1 (any number larger to or equal to 1\
 is allowed)'),
        ('prefer_end', 'Prefer pairs closer the end of genes')

    ]

    ####################################################################

    def _run(self, *fargs, **kwargs):

        align = data.Align(kwargs['input'], groups=True)
        output = open(kwargs['output'], 'w')

        self.db = None

        try:

            # write output file header

            output.write('gene,start,end,size,Fname,Fseq,FTm,FANY,FEND,FQ,FGC%,FBlastScore,FAutoBlastScore,Rname,Rseq,RTm,RANY,REND,RQ,RGC%,RBlastScore,RAutoBlastScore\n')

            c=0

            # prepare blast database

            self.container = data.Container()
            for n,s,g in align:
                s = s.translate(string.maketrans('?', 'N'), '-')
                self.container.append(n, s, g)
            
            self.db = wrappers.BLASTdb(self.container, 'nucl')
            
            for n,s,g in self.container:
                
                
                # applies to a subset of sample only (if required)
                
                if 'selection' in fargs and g==0:
                    continue

                if not self.quiet:
                    print n

                primers= wrappers.Primer3(
                    s,
                    PRIMER_MIN_TM=kwargs['minTm'],
                    PRIMER_OPT_TM=kwargs['optTm'],
                    PRIMER_MAX_TM=kwargs['maxTm'],
                    PRIMER_MIN_GC=kwargs['minGc'],
                    PRIMER_OPT_GC_PERCENT=kwargs['optGc'],
                    PRIMER_MAX_GC=kwargs['maxGc'],
                    PRIMER_MAX_NS_ACCEPTED=kwargs['numAmb']
                )
    
                # searches primers

                nf, nr = primers.find_primers()
                np = primers.find_pairs(kwargs['sizemin'], kwargs['sizemax'])

                if not self.quiet:
                    print '    %d forward, %d reverse primers' %(nf, nr)
                    print '    %d pairs' %np,
                    sys.stdout.flush()

                # score primers and sort

                pairs = primers.pairs()

                for i,v in enumerate(pairs):
                                        
                    a = tools.locate(s, v['F']['seq'])
                    b = tools.locate(s, tools.rc(v['R']['seq']))

                    if (a==None):
                        raise RuntimeError, 'sprimers: this primer was not found: {0}'.format(v['F']['seq'])
                    if (b==None):
                        print tools.rc(v['R']['seq'])
                        print s
                        raise RuntimeError, 'sprimers: this primer was not found: {0}'.format(v['R']['seq'])

                    b += len(v['R']['seq'])-1
                    pairs[i]['coord'] = (a,b)
                    pairs[i]['pscore'] = ('prefer_end' in fargs)*self.pscore(len(s), a+(b-a)/2.)

                # the pairs closest of the end are ranked first

                pairs.sort(cmp= lambda x,y: cmp(y['pscore'], x['pscore']))

                # pre-selects primers

                pairs.sort(lambda x,y: cmp(
                    x['F']['Q']+x['R']['Q']+(1-x['pscore']),
                    y['F']['Q']+y['R']['Q']+(1-y['pscore']))
                )
                
                pairs = pairs[:kwargs['filter1']]
                
                if not self.quiet:
                    print '-> %d preselected (max %d)' %(len(pairs), kwargs['filter1'])

                # selects based on blast
                
                blast_results= {}
                accepted= []
                
                if not self.quiet:
                    updater = tools.Updater(len(pairs))

                for i in pairs:

                    if not i['F']['seq'] in blast_results:
                        blast_results[i['F']['seq']] = self.blprimer(i['F']['seq'], n)
                        
                    if not i['R']['seq'] in blast_results:
                        blast_results[i['R']['seq']] = self.blprimer(i['R']['seq'], n)

                    fscores = blast_results[i['F']['seq']]
                    rscores = blast_results[i['R']['seq']]
                    
                    if ((fscores[0]+rscores[0])/2. <= kwargs['threshold1'] and
                        (fscores[1]+rscores[1])/2. <= kwargs['threshold2']):
                        i['fscores']= fscores
                        i['rscores']= rscores
                        accepted.append(i)
                        
                        if len(accepted) == kwargs['filter2']:
                            break

                    if not self.quiet:
                        updater.refresh('    blasting %d/%d $ELAPSED (unknown; max: $REMAINING)' %(len(accepted),kwargs['filter2']))
                        
                if not self.quiet:
                    updater.refresh('    blasting done')
                    updater.close()

                pairs = accepted
                pairs.sort(lambda x,y: cmp((x['fscores'][1]+x['rscores'][1]), (y['fscores'][1]+y['rscores'][1])))
                pairs.sort(lambda x,y: cmp((x['fscores'][0]+x['rscores'][0]), (y['fscores'][0]+y['rscores'][0])))
                
                if not self.quiet:
                    print '    %d primer pairs validated (max %d)' %(len(accepted), kwargs['filter2'])

                # check
                
                del primers.pairs()[:]
                primers.pairs().extend(pairs)
                
                primers.check_pairs()

                if not self.quiet:
                    print '    Final pairs: %d reduced to %d' %(len(primers.pairs()), min(kwargs['show'], len(primers.pairs())))

                pairs = primers.pairs()[:kwargs['show']]

                # outputting
                
                d=0
                for i in pairs:
                    
                    d += 1

                    output.write('%s,%d,%d,%d' %(n,i['start'],i['end'],i['size']))
                    output.write(',%s,%s,%f,%f,%f,%f,%f,%f,%f' %(
                        '%s_%sF' %(n,str(d).rjust(2,'0')),
                        i['F']['seq'],
                        i['F']['Tm'],
                        i['F']['ANY'],
                        i['F']['END'],
                        i['F']['Q'],
                        i['F']['GC%'],
                        i['fscores'][0],
                        i['fscores'][1]))
                    output.write(',%s,%s,%f,%f,%f,%f,%f,%f,%f\n' %(
                        '%s_%sR' %(n,str(d).rjust(2,'0')),
                        i['R']['seq'],
                        i['R']['Tm'],
                        i['R']['ANY'],
                        i['R']['END'],
                        i['R']['Q'],
                        i['R']['GC%'],
                        i['rscores'][0],
                        i['rscores'][1]))
                    output.flush()
                c+=1

        finally:

            output.close()

            if self.db: del self.db

        if not self.quiet:
            print 'finished'

    ####################################################################

    def pscore(self, length, position):
 
        """
        Computes a position score
        """
 
        return 1/(1+math.exp(5-10.*position/length))

    ####################################################################

    def blprimer(self, primer, current):

        """
        Blasts a primer
        """

        bl = wrappers.BLAST()

        hits = bl.blastn(primer, self.db, evalue=1000,
                         word_size=7, penalty=-1, gapopen=1000)

        allo_hits = []
        auto_hits = []

        for i in hits['']:
            
            if i['subject'] == current:
                auto_hits.append(i)
            else:
                allo_hits.append(i)

        for i in auto_hits+allo_hits:

            (a,b)= sorted([i['qstart'], i['qend']])
            (c,d)= sorted([i['hstart'], i['hend']])

            c -= a
            d += (len(primer)-1)-b

            seq = self.container.sequenceByName(i['subject'])

            hit = seq[max(0,c):d+1]
            if ((i['qstart']>i['qend'] or i['hstart']>i['hend']) and
                not (i['qstart']>i['qend'] and i['hstart']>i['hend'])):
                hit= tools.rc(hit)

            for j in range(c,0):
                hit='N'+hit
            for j in range(len(seq)-1,d):
                hit=hit+'N'

            score=0.
            mscore=0.

            for j in range(len(primer)):
                if (primer[j]==hit[j]):
                    score += self.pscore(len(primer), j)
                mscore += self.pscore(len(primer), j)
            score /= mscore
            i['score'] = score

        auto_hits.sort(lambda x,y: cmp(y['score'], x['score']))
        allo_hits.sort(lambda x,y: cmp(y['score'], x['score']))

        try:
            X = allo_hits[0]['score']
        except IndexError:
            X = 0.
        try:
            Y = auto_hits[1]['score']
        except IndexError:
            Y = 0.

        return X,Y

commands.append(sprimers)

########################################################################

class staden2fasta(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Converts a STADEN GAP4 dump file to fasta.'

    description = 'The file must have been generated using the command\
 `dump contig to file` of the GAP4 contig editor. This command will\
 generate a fasta alignment file, padding sequences with `?` whenever\
 necessary.'

    options = [
        Option('input', 'Input dump file', str, []),
        Option('output', 'Output fasta alignment file', str, []),
        Option('consensus', 'Defines what should be done with the\
 sequence named `CONSENSUS`. Three values are possible: `remove`:\
 removes the `CONSENSUS` sequence (if it is present); `keep`: keeps the\
 `CONSENSUS` sequence (if it is present) and `only`: removes all other\
 sequences and keeps only the `CONSENSUS` sequence (it must be present)',
                str, 'remove', [lambda x: x in set(['remove', 'keep', 'only'])])
    ]
        
    ####################################################################

    flags = []

    ####################################################################

    def _run(self, *fargs, **kwargs):

        f = open(kwargs['input'])
        string = f.read()
        f.close()

        try:
            staden = egglib_binding.Staden.parse(string, False)
        except IOError:
            raise IOError, 'cannot read GAP4 dump file from %s' %kwargs['input']

        align = data.Align()

        for i in range(staden.ns()):
            align.append(staden.name(i), staden.sequence(i))
        
        if kwargs['consensus'] == 'remove':
            i = align.find('CONSENSUS')
            if i != None:
                del align[i]
                
        if kwargs['consensus'] == 'only':
            i = align.find('CONSENSUS')
            if i == None:
                raise ValueError, '`CONSENSUS` is not in %s' %kwargs['input']
            CONSENSUS = align[i]
            align.clear()
            align.append(*CONSENSUS)

        if not self.quiet:
            print 'number of sequences: %d' %align.ns()
            print 'length of alignment: %d' %align.ls()

        align.write(kwargs['output'])
        
commands.append(staden2fasta)

########################################################################

class translate(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Translates coding sequences to protein sequences.'

    description = 'The `quiet` argument is ignored.'

    options = [
        Option('input', 'Input fasta file', str, '', []),
        Option('output', 'Output fasta file', str, '', []),
        Option('code', 'Genetic code specification. Should be an integer\
 among the available values. Use the flag `codes` to display the list of\
 avalaible genetic codes. The default corresponds to the standard code',
                int, 1, [])
    ]
        
    ####################################################################

    flags = [('codes', 'displays available genetic codes')]

    ####################################################################

    def _run(self, *fargs, **kwargs):
        
        if 'codes' in fargs:
            
            print 'Available codes:'

            for i, short, long in tools.GeneticCodes.codes():
                print str(i).rjust(3)+':', tools.wrap(long, 72, 5)
            
        else:
            
            if kwargs['input']=='':
                raise ValueError, 'error: input argument must be given'

            if kwargs['output']=='':
                raise ValueError, 'error: output argument must be given'

            fasta = data.Container(kwargs['input'])
            fasta = tools.translate(fasta, kwargs['code'])
            fasta.write(kwargs['output'])
        
commands.append(translate)

########################################################################

class ungap(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Removes gaps from a sequence alignment.'

    description = 'This command either removes all gaps from an\
 alignment (break the alignment) or removes alignment positions (column)\
 where the frequency of gaps is larger than a given threshold.'

    options = [
        Option('input', 'Input fasta file', str, None, []),
        Option('output', 'Output fasta file', str, None, []),
        Option('threshold', 'Proportion giving the threshold frequency\
 for removing gaps. All sites for which the frequency of gaps is equal\
 to or larger than the specified values will be removed. A value of 0\
 will remove all sites and a value of 1 will remove only columns\
 consisting of gaps only. If the flag `all` is used, the value of this\
 option is ignored',
                float, 1, [lambda x: x<=1. and x>=0.])
    ]
        
    ####################################################################

    flags = [('all', 'Removes all gaps of the alignment - regardless of\
 the value of the `threshold` argument. The output sequences will not\
 be aligned anymore (save for special cases)'),
    ('triplets', 'Removes only complete triplets (the alignment length\
 must be a multiple of 3)')]

    ####################################################################

    def _run(self, *fargs, **kwargs):
        
        align = data.Align(kwargs['input'])

        if not self.quiet:
            print 'alignment length: %d' %align.ls()

        if 'all' in fargs:
            
            result = tools.ungap_all(align)

            if not self.quiet:
                print 'longest sequence: %d' %max([result.ls(i) for i in range(len(result))])
            
        else:

            if 'triplets' in fargs:
                result = tools.ungap_triplets(align, kwargs['threshold'])

            else:
                result = tools.ungap(align, kwargs['threshold'])

            if not self.quiet:
                print 'new alignment length: %d' %result.ls()


        result.write(kwargs['output'])

commands.append(ungap)

########################################################################

class winphyml(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Computes tree likelihood along a sliding window.'

    description = 'This command runs a sliding window along a sequence\
 alignment. For each window, it computes the likelihood of the\
 maximum-likelihood tree along as well as the likelihood of a given set\
 of trees. It can detect regions of a sequence that support a given tree\
 rather than an other. The command expects nucleotide sequences.'

    options = [
        Option('input', 'Input fasta file', str, None, []),
        Option('trees', 'Input newick file containing one or more trees',
                str, None, []),
        Option('output', 'Main output file name', str, 'winphyml.csv', []),
        Option('wsize', 'Window size', int, 500, [lambda x: x>=10]),
        Option('wstep', 'Window step length', int, 100, [lambda x: x>=1]),
    ]
        
    ####################################################################

    flags = [('savetrees', 'Saves the maximum-likelihood tree for each\
 windows. Each window tree will be saved as a file named\
 `<base>_<start>_<end>.tre` when `base` is the name of the main output\
 file minus the extension if there is one, and `start` and `end` are the\
 limits of the window. With default values, the trees will be saved as\
 `winphyml_1_200`, `winphyml_21_220`, etc.')]

    ####################################################################

    def _run(self, *fargs, **kwargs):
        
        # import data

        align = data.Align(kwargs['input'], groups=True)
                
        f = open(kwargs['trees'])
        string = f.read().strip()
        f.close()
        trees = [ data.Tree(string=bit+';') for bit in string.split(';')[:-1] ]

        if not self.quiet:
            print 'sequence data: %d sequences - %d sites' %(align.ns(), align.ls())
            print 'number of input trees: %d' %len(trees)
            
        # check that names match
        
        names = sorted(align.names())
        for tree in trees:
            if sorted(tree.all_leaves()) != names:
                print
                print tree.all_leaves()
                print names
                raise ValueError, 'error: the names of the alignment don\'t match those of all trees'

        # initializes output file

        output = open(kwargs['output'], 'w')
        output.write(','.join(['Wstart', 'Wend','Best tree']+
                           ['Tree %i' %(i+1) for i in range(len(trees))])+'\n')

        # defines bounds
        
        bounds = []
        i=0
        while i<align.ls():
            bounds.append((i, i+kwargs['wsize']))
            i += kwargs['wstep']
        
        print '%d windows will be processed' %len(bounds)
        
        if not self.quiet:
            updater = tools.Updater(len(bounds))
        
        # loops
        
        pos = kwargs['output'].rfind('.')
        root = kwargs['output'][:pos]

        c = 0
        
        for start, end in bounds:
        
            c+=1
            
            if not self.quiet:
                updater.refresh('window %d: %d - %d [->%d]  $ELAPSED\
 [$REMAINING]' %(c,start,end,align.ls()), increment=0, grain=-1)

            window = align.extract(start, end)
            
            ml_tree, lk0 = wrappers.phyml(window, model='GTR', rates=1,
                                          boot=0, quiet=True)

            if 'savetrees' in fargs:
                ml_tree.write('%s_%d_%d.tre' %(root, start+1, end+1))

            lk1 = []

            for tree in trees:
                tree, lki = wrappers.phyml(window, topo=tree, model='GTR',
                                           rates=1, boot=0, quiet=True)
                lk1.append(lki)

            output.write( ','.join(map(str,[start,end,lk0] + lk1)) + '\n' )
            output.flush()

            if not self.quiet:
                updater.increment(1)
            
        # finishes

        if not self.quiet:
            updater.refresh('done', grain=0)
            updater.close()

        output.close()
        
commands.append(winphyml)

########################################################################

class phyml(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Performs maximum-likelihood phylogenetic reconstruction.'

    description = 'This command reconstructs a phylogenetic tree and,\
 optionnally, performs bootstrap repetitions. Crashes occurring during\
 the bootstrap procedure due to estimation problems are ignored,\
 allowing to complete the run. The substitution model name implies data\
 type (HKY85, JC69, K80, F81, F84, TN93 and GTR imply nucleotides,\
 others imply amino acids). `LG` is default for amino acids in the\
 stand-alone `phyml` software.'

    options = [
        Option('input', 'Input alignment fasta file', str, None, []),
        Option('output', 'Output newick tree file', str, None, []),
        Option('boot_trees', 'Output newick tree file for bootstrap\
 trees (by default, no file)', str, '', []),
        Option('boot', 'Number of bootstrap repetions', int, 0,
               [lambda x: x>=0]),
        Option('model', 'substitution model name', str, 'HKY85', []),
        Option('rates', 'number of gamma substitution rate categories',
                int, 1, [lambda x: x>0]),
        Option('search', 'tree topology search operation option (`NNI`,\
 `SPR` (slower) or `BEST`', str, 'NNI', [lambda x: x in ['NNI', 'SPR',
  'BEST']]),
        Option('start', 'name of a file containing the tree topology to\
 use as starting tree (by default, use a distance tree)', str, '', []),
        Option('recover', 'name of a file containing already generated\
 bootstrap trees - must contain the same taxa and can be the same as\
 `boot_trees` (ignored if `boot` is 0)', str, '', [])
    ]
        
    ####################################################################

    flags = []

    ####################################################################

    def _run(self, *fargs, **kwargs):
        
        # imports data
        align = data.Align(kwargs['input'], groups=True)

        if kwargs['start']!='':
            start = data.Tree(kwargs['start'])
            start.clean_edge_lengths()
            start.clean_internal_labels()
            if sorted(start.all_leaves())!=sorted(align.names()):
                raise ValueError, 'the starting tree must have the same taxon names as the input alignment'
        else:
            start = None

        # builds main tree

        if self.quiet == False and kwargs['boot']==0:
            quiet = False
        else:
            quiet = True


        main_tree, lk = wrappers.phyml(align,
                                       model=kwargs['model'],
                                       search=kwargs['search'],
                                       rates=kwargs['rates'],
                                       boot=0,
                                       start=start,
                                       quiet=quiet)


        # performs bootstrap

        if kwargs['boot'] > 0:
            
            trees = []
            
            if kwargs['recover'] != '':
                f = open(kwargs['recover'])
                for line in f:
                    tree = data.Tree(string=line)
                    if sorted(tree.all_leaves()) != sorted(main_tree.all_leaves()):
                        raise ValueError, 'the leaves in the recovered bootstrap tree do not correspond to the alignment!'
                    trees.append(tree)
                f.close()
                
            if kwargs['boot_trees'] != '':
                f = open(kwargs['boot_trees'], 'w')
                for tree in trees:
                    f.write(str(tree) + '\n')
                f.flush()
            
            random = egglib_binding.Random()

            while(len(trees)<kwargs['boot']):

                if not self.quiet:
                    print 'bootstrap %d of %d...' %(len(trees)+1, kwargs['boot'])

                l = align.ls()
                
                sites = [random.irand(l-1) for i in range(l)]
                sites.sort()
                
                boot = data.Align()
                for seq in align:
                    boot.append(seq.name, ''.join([seq.sequence[i] for i in sites]))

                try:
                    tree, lk = wrappers.phyml(boot,
                                              model=kwargs['model'],
                                              search=kwargs['search'],
                                              rates=kwargs['rates'],
                                              boot=0,
                                              quiet=quiet)

                except RuntimeError:
                    if not self.quiet:
                        print '  .. failed! ignoring...'
                    continue

                if kwargs['boot_trees'] != '':
                    f.write(str(tree) + '\n')
                    f.flush()
                trees.append(tree)

            if kwargs['boot_trees'] != '':
                f.close()

            main_tree.frequency_nodes(trees)

        main_tree.write(kwargs['output'])
        if not self.quiet:
            print 'done'

commands.append(phyml)

########################################################################

class diffscan(BaseCommand):
    
    """
    Documentation
    """

    ####################################################################

    brief = 'Scans loci for signatures of adaptive differentiation.'

    description = 'This command applies the method of Beaumont and\
 Nichols 1996 (Proceedings R Soc. Biol. Sci. 263:1619-1626). It uses a\
 large number of loci to estimate the genome-wide between-population\
 index of fixation based on Weir and Cockerham 1984 (Evolution\
 18:1358-1370; equation 10). The fixation index is called theta in Weir\
 and Cockerham but we retain the Fst terminology here.\
\n\
The coalescent migration rate 4Nm (here called M) is estimated as:\n\
\n\
M = (1-Fst)*(d-1)/(Fst*d),\n\
\n\
assuming an island model where `d` is the number of populations in the\
 system. Coalescent simulations are performed assuming the number of\
 populations, the actual set of samples and assuming a single mutation\
 per locus.\n\
\n\
Input file:\n\
\n\
The input file is a string of one or more loci. Each locus is\
 represented by populations  (demes in Weir and Cockerham). There must\
 be at least two populations. The number of populations must be\
 consistent over loci. Note that white lines are ignored throughout the\
 file and can be used as separators but are not required and need not to\
 be homogeneously used. Spaces and tabs can also be used to align the\
 file and are ignored when they occur at the beginning and end of lines.\n\
\n\
Comments:\n\
\n\
Comments are lines starting with with a hash (`#`) symbol. White\
 spaces before the hash are ignored. Comments cannot be placed at the\
 end of lines containing data.\n\
\n\
Loci:\n\
\n\
Loci take a single line each. The type of the locus is given by\
 reserved symbols. `$` denotes reference loci (they will be used for\
 computing genome-wide parameters and tested) and `%` denotes candidate\
 loci (which will be skipped for estimating genome-wide parameters).\
 Type symbol must appear before data. An optional label can be placed\
 before the symbol. Labels are used to name the locus (by default, an\
 automatic label based on its rank is be applied). The same label might\
 be used several times and labelled and unlabelled loci might be mixed.\
 Labels cannot start by a hash (`#`) symbol, otherwise the line is taken\
 as a comment. Labels cannot contain the dollar (`$`) and percent (`%`)\
 symbols. The general structure is therefore: `[label]$ data` for\
 reference loci and `[label]% data` for test loci. See definition of\
 data and example below.\n\
\n\
Locus data:\n\
\n\
Locus data is given by pairs of allele counts, one for each\
 population. The number of populations must be the same across loci.\
 The alleles are provided in an arbitrary order.\
 Counts for both alleles must be provided as two integer values\
 separated by a comma (`,`). Population (ie, pairs counts) are\
 separated from each other and from the type symbol by any number of\
 spaces (tabs are also supported). Frequencies must be 0 or more.\
 Unsampled populations (`0,0`) are allowed. (If a population is\
 missing for all loci, better use the `k` option to specify the real\
 number of populations). The first and second alleles are equivalent (no\
 orientation) but they must be the same across all populations of a\
 given locus.\n\
\n\
An example of the input file is provided below. This data set comprises\
 two reference loci and one locus to be used for testing only, a total\
 of five populations with varying sample size.\n\
\n\
    # A comment\n\
    Reference locus 1 $  10,4  4,1  5,0  12,1   6,3\n\
    Reference locus 2 $   4,5  3,2  2,4   4,6   3,5\n\
    Test locus 1      %  15,1  8,0  3,1   2,12  0,10\n\
\n\
The command will first perform the number of requested simulations,\
 unless the argument `load_simuls` is set. In this case, simulations will\
 be imported from a text files containing He, Fst values (one pair per\
 line). In this case, simulations parameters (`simuls` and `step`) will\
 be ignored. The option `save_simuls` (inactive by default) allows to\
 save simulations and to import them in a following run, eg for trying\
 out different values of the binarization factor.\n\
\n\
The final file contains, for all loci, its He and Fst value, and the\
 p-value derived from the distribution.'

    options = [
        Option('input', 'Genotypes iput file (see description above)', str, None, []),
        Option('output', 'Test results output file', str, None, []),
        Option('plot', 'Graphical output file (requires matplotlib)', str, '', []),
        Option('k', 'Number of populations. By default (0), the number\
 of populations in the data set is used. If specified, the value must be\
 at least equal to the number of populations in the data set', int, 0, [lambda x: x>=0]),
        Option('simuls', 'Number of simulation rounds', int, 10000, []),
        Option('step', 'Number of simulation between screen update', int, 100, []),
        Option('bins', 'Number of bins in the [0,0.5] heterozygosity range', int, 12, [lambda x: x>1]),
        Option('save_simuls', 'Name of file where to save simulations (default: no save)', str, '', []),
        Option('load_simuls', 'Name of file where to import simulations (this will skip simulations; default: perform simulations)', str, '', [])
    ]
        
    ####################################################################

    flags = []

    ####################################################################

    def _run(self, *fargs, **kwargs):

        if not self.quiet: print '<diffscan>'
        
        # imports data
        self.parse(kwargs['input'])
        
        # checks number of populations
        if not self.quiet:
            print '> Number of populations: {0}'.format(self.k)
        if kwargs['k']!=0:
            if kwargs['k'] < self.k:
                raise ValueError, 'Argument `k` is too small!'
            self.k = kwargs['k']
            print '> Actual number of populations: {0}'.format(self.k)

        # computes fstatistics
        self.stats()

        # determines all available sampling schemes
        
        self.sam = []
        for locus, ref in zip(self.loci, self.types):
            if ref==False:
                continue
            self.sam.append(map(sum, locus))
            
        if len(self.sam) != self.n1:
            raise RuntimeError, 'Internal error'

        # performs simulations if needed
        
        if kwargs['load_simuls'] == '':
            self.simulate(kwargs['simuls'], kwargs['step'], kwargs['save_simuls'])
        else:
            self.import_simuls(kwargs['load_simuls'])

        # discretizes
        self.discr(kwargs['bins'])

        # plot if needed
        if kwargs['plot'] != '': self.plot(kwargs['plot'])
           
        # perform tests
        self.test(kwargs['output'])
               
    ####################################################################
        
    def parse(self, fname):
        
        self.data = []
        self.k = set()
        self.nb_comments = 0
        self.line_c = 0
        self.n1 = 0
        self.n2 = 0
        
        if not self.quiet:
            print '> Reading data from: {0}'.format(fname)
        
        f = open(fname)
        
        for line in f:
            
            self.line_c += 1
        
            # remove leading and trailing white spaces
            line = line.strip()
            
            # skip comments
            if line[0] == '#':
                self.nb_comments += 1
                continue

            # determine type, extract data
            
            dollar = line.count('$')
            percen = line.count('%')
            
            if dollar>1: raise ValueError, 'Invalid format: line {0} contains more than one $'.format(self.line_c)
            if percen>1: raise ValueError, 'Invalid format: line {0} contains more than one $'.format(self.line_c)
            if dollar+percen>1: raise ValueError, 'Invalid format: line {0} contains both $ and %'.format(self.line_c)
            if dollar+percen==0: raise ValueError, 'Invalid format: line {0} contains none of $ and % and is not comment'.format(self.line_c)

            if dollar==1:
                symb = '$'
                reference = True
                self.n1 += 1
            elif percen==1:
                symb = '%'
                reference = False
                self.n2 += 1
            else: symb = None

            if symb == line[0]:
                label = None
                data = line[1:]

            else:
                line = line.split(symb)
                if len(line) != 2:
                    raise ValueError, 'Invalid format: line {0} does not have a valid structure'.format(self.line_c)
                label = line[0].strip()
                data= line[1].strip()

            self.data.append([reference, label, []])
            
            # processes populations
            
            data = data.split()
            k = 0
            for pop in data:
                k += 1
                pop = pop.split(',')
                if len(pop) != 2:
                    raise ValueError, 'Invalid format: line {0}, population {1} does not contain two items'.format(self.line_c, k)
                try:
                    freq = map(int, pop)
                except ValueError:
                    raise ValueError, 'Invalid format: line {0}, population {1} - does not contain integers'.format(self.line_c, k)

                self.data[-1][2].append(freq)
            
            self.k.add(k)

        f.close()

        # general checking

        if len(self.k) == 0:
            raise ValueError, 'Invalid format: no locus defined'
        if len(self.k) != 1:
            raise ValueError, 'Invalid format: inconsistent number of populations'
        self.k = self.k.pop()
        if self.k < 2:
            raise ValueError, 'Invalid format: unsufficient number of populations'
        
        if len(self.data)!=self.n1+self.n2:
            raise ValueError, 'Invalid format: internal error'

        if len(self.data)==0:
            raise ValueError, 'Invalid format: empty dataset'

        if self.n1 == 0:
            raise ValueError, 'Not enough reference loci: cannot continue'

        # feedback
        
        if not self.quiet:
            print '> Number of read lines (including comments): {0}'.format(self.line_c)
            print '> Number of reference loci: {0}'.format(self.n1)
            print '> Number of test loci: {0}'.format(self.n2)
            print '> Total number of loci: {0}'.format(self.n1 + self.n2)

        # reformat data
        
        self.types, self.labels, self.loci = map(list, zip(*self.data))
        del self.data

        # generate default names
        
        L = len(str(len(self.loci)))
        for i in range(len(self.loci)):
            if self.labels[i]==None:
                self.labels[i] = 'locus{0}'.format(str(i+1).rjust(L, '0'))

    ####################################################################
    
    def stats(self):
        
        self.locus_Fst = []
        self.locus_He = []
        self.locus_p = []
        
        T1 = 0.
        T2 = 0.
        
        for reference, locus in zip(self.types, self.loci):

            hfs = egglib_binding.HFStatistics()
            n = sum(map(sum,locus))
            p = 0
            hfs.reserve( n )
            for pop, (A, B) in enumerate(locus):
                p += B
                for j in range(A):
                    hfs.loadIndividual(0, pop)
                for j in range(B):
                    hfs.loadIndividual(1, pop)

            p = 1.*p/n
            self.locus_p.append(p)
            self.locus_He.append(2*p*(1-p))

            t1 = hfs.T1()
            t2 = hfs.T2()
                
            if t2>0: self.locus_Fst.append(t1/t2)
            else: self.locus_Fst.append(None)

            if reference==True:
                T1 += t1
                T2 += t2

        # computes F-stat
        
        if T2==0:
            raise ValueError, 'Cannot estimate Fst: no variance'
            
        Fst = T1/T2
        
        # deduce M
        
        if Fst == 0:
            raise ValueError, 'Cannot estimate M: Fst is null'
            
        self.M = (1.-Fst)*(self.k-1)/(Fst*self.k)

        # feedback
        if not self.quiet:
            print '> Estimated Fst: {0}'.format(Fst)
            print '> Estimated M: {0}'.format(self.M)

        if Fst < 0:
            raise ValueError, 'Fst is negative - cannot continue'
        
    ####################################################################

    def simulate(self, nrepet, step, fname):

        self.distri = []
        
        if not self.quiet:
            print '> Performing {0} simulations now'.format(nrepet)
        if fname != '':
            print '> Saving simulations to {0}'.format(fname)
            f = open(fname, 'w')

        while len(self.distri)<nrepet:

            # pick a sampling scheme and complete with unsample populations
            
            sam = random.choice(self.sam)
            sam.extend([0] * (self.k - len(sam)))
            if len(sam) != self.k:
                raise RuntimeError, 'Internal error'
            
            # makes <step> simulations with these parameters
            
            ps = simul.CoalesceParamSet(sam, M=self.M)
            mu = simul.CoalesceFiniteAlleleMutator(0)
            mu.fixedNumberOfMutation(1)
            dms = simul.coalesce(ps, mu, step, convert=False)

            for dm in dms:

                # generate the simulated table and analyze
                
                if sum(sam) != dm.numberOfSequences():
                    raise RuntimeError, 'Internal error'
                    
                hfs = egglib_binding.HFStatistics()
                hfs.reserve(sum(sam))

                p = 0
                for i in range(dm.numberOfSequences()):
                    allele = dm.get(i, 0)
                    label = dm.populationLabel(i)
                    hfs.loadIndividual(allele, label)
                    p += allele

                p = 1.*p/sum(sam)
                He = 2 * p * (1-p)
                Fst = hfs.theta()

                self.distri.append((He, Fst))

                if fname != '':
                    f.write('{0} {1}\n'.format(He, Fst))

                if len(self.distri)==nrepet:
                    break

            if not self.quiet:
                print '  --- {0} simulations of {1} ---'.format(len(self.distri), nrepet)
        if fname != '':
            f.close()
        if not self.quiet:
            print '> Simulation finished'

    ####################################################################

    def import_simuls(self, fname):

        if not self.quiet:
            print '> Loading simulations from {0}'.format(fname)

        self.distri = []
        f = open(fname)
        c = 0
        for line in f:
            c += 1
            bits = line.split()
            if len(bits) != 2:
                raise IOError, 'Cannot read {0}: invalid number of tokens line {1}'.format(fname, c)
            try:
                a, b = map(float, bits)
            except ValueError:
                raise IOError, 'Cannot read {0}: invalid tokens are line {1} are not integers'.format(fname, c)
            self.distri.append((a,b))
        f.close()

        if not self.quiet:
            print '> {0} simulated imported'.format(c)

    ####################################################################

    def discr(self, bins):
        
        p = 0.01, 0.05, 0.5, 0.95, 0.99
        
        if not self.quiet:
            print '> Binarizing'

        # binarizes the distribution
        
        self.distri.sort(lambda x,y: cmp(x[0],y[0]))

        self.bins = []
        c = 0
        step = 0.5/bins
        while len(self.distri)>0:
            self.bins.append(((c, c + step/2., c+step), [None, []]))
            while len(self.distri)>0 and self.distri[0][0] <= c + step:
                self.bins[-1][1][1].append(self.distri.pop(0)[1])

            # finds median Fst
            L = len(self.bins[-1][1][1])
            if L >= 100:

                if L%2==1:
                    self.bins[-1][1][0] = self.bins[-1][1][1][L//2]

                else:
                    a = self.bins[-1][1][1][L//2]
                    b = self.bins[-1][1][1][L//2+1]
                    self.bins[-1][1][0] = a+(b-a)/2
                
            c += step

        # computes quantiles of each vector
        
        self.L1 = []
        self.L2 = []
        self.L3 = []
        self.L4 = []
        self.L5 = []
        for (beg, mid, end), (median,dist) in self.bins:
            dist.sort()
            L = len(dist)
            q = [i*L for i in p]
            q[0] = int(math.floor(q[0]))
            q[1] = int(math.floor(q[1]))
            q[2] = int(math.floor(q[2]))
            q[3] = int(math.ceil(q[3]))
            q[4] = int(math.ceil(q[4]))
            
            if L >= 1000: self.L1.append((mid, dist[q[0]]))
            if L >= 200:  self.L2.append((mid, dist[q[1]]))
            if L >= 20:   self.L3.append((mid, dist[q[2]]))
            if L >= 200:  self.L4.append((mid, dist[q[3]]))
            if L >= 1000: self.L5.append((mid, dist[q[4]]))
            
    ####################################################################
  
    def plot(self, fname):

        if not self.quiet:
            print '> Plotting to {0}'.format(fname)

        try:
            from matplotlib import pyplot
        except ImportError:
            raise ValueError, 'Cannot plot: the matplotlib module is not available'

        # plot quantile lines

        if len(self.L1):
            x, y = zip(*self.L1)
            pyplot.plot(x, y, 'k:')

        if len(self.L2):
            x, y = zip(*self.L2)
            pyplot.plot(x, y, 'k--')

        if len(self.L3):
            x, y = zip(*self.L3)
            pyplot.plot(x, y, 'k-')

        if len(self.L4):
            x, y = zip(*self.L4)
            pyplot.plot(x, y, 'k--')

        if len(self.L5):
            x, y = zip(*self.L5)
            pyplot.plot(x, y, 'k:')

        # plot points
        
        x0 = [i for (i,j) in zip(self.locus_He, self.types) if j]
        x1 = [i for (i,j) in zip(self.locus_He, self.types) if not j]
        y0 = [i for (i,j) in zip(self.locus_Fst, self.types) if j]
        y1 = [i for (i,j) in zip(self.locus_Fst, self.types) if not j]

        pyplot.plot(x0, y0, 's', mfc='None', mec='b')
        pyplot.plot(x1, y1, 'o', mfc='None', mec='r')
        
        pyplot.xlabel('He')
        pyplot.ylabel('Fst')
        pyplot.savefig(fname)
        pyplot.clf()
        
    ####################################################################

    def test(self, fname):
        
        if not self.quiet:
            print '> Performing tests'
        

        f = open(fname, 'w')
        f.write('Locus\tp\tHe\tFst\tp-value\n')
        c1 = 0
        c2 = 0
        t = 0
        
        for typ, label, p, He, Fst in zip(self.types, self.labels,
                        self.locus_p, self.locus_He, self.locus_Fst):
            
            if Fst==None:
                f.write('{0}\t{1}\t{2}\tNA\tNA\n'.format(label, p, He))
                continue
            
            # spot right category
            for (beg, mid, end), (median, data) in self.bins:

                if len(data)<100:
                    P = 'NA'
                    break

                if He >= beg and He <= end:

                    lims = [None,None]

                    if Fst < median:
                        lims[0] = Fst
                        lims[1] = median + (median-Fst)
                    else:
                        lims[0] = median - (Fst-median)
                        lims[1] = Fst
                    
                    P = 0
                    c = 0
                    while data[c]<=lims[0]:
                        c+=1
                        P+=1
                    c = -1
                    while data[c]>=lims[1]:
                        c-=1
                        P+=1
                    
                    P = 1. * P / len(data)
                    t += 1
                    if typ:
                        if P<0.05: c1+=1
                        if P<0.01: c2+=1

                    break
            
            else:
                raise ValueError, 'Internal error - cannot identify bin for {0}'.format(He)
        
            f.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(label, p, He, Fst, P))
      
        
        if not self.quiet:
            print '> Number of loci: {0}'.format(len(self.labels))
            print '> Number of loci tested: {0}'.format(t)
            if t>0:
                print '> Reference loci significant at 5%: {0} ({1:.2%})'.format(c1,1.*c1/t)
                print '> Reference loci significant at 1%: {0} ({1:.2%})'.format(c2,1.*c2/t)
        
        f.close()
       
    ####################################################################

commands.append(diffscan)


########################################################################

dcommands = dict(zip([i.__name__ for i in commands], commands))


########################################################################

def execute(*fargs, **fkwargs):

    """
    Execute utils commands. This functions takes arguments to specify
    the command name and its options. There must be at least one
    positional argument: the first positional argument gives the
    function name and other (optional) positional arguments give the
    command flags to be activated. The keyword arguments give the
    command options. Flag, option names but also option values should be
    string (option values will be converted automatically). In case
    options are of a simple type (int, float), they can be passed as
    such. But options that are described as a string presenting a list
    of values separated by commas CANNOT be passed as a list using the
    function. If there is no argument whatsoever, arguments will be read
    from :data:`sys.argv`. In this case, the first argument will be
    ignored; the second argument must be the  command name; and keyword
    arguments must be passed under the form ``key=value`` where *key* is
    the option name.
    
    For example, running the command::

        egglib ungap input=file1.fas output=file2.fas threshold=0.5 triplets
        
    is the equivalent of calling the function::

        >>> egglib.utils.execute('ungap', 'triplets', input='file1.fas', output='file2.fas', threshold=0.5)

    where *threshold* can also take the string ``"0.5"``.
    
    .. versionchanged 2.1.0:: Accepts arguments in order to be run from
        python scripts.
    """

    global commands
    
    if len(fargs):
        command = fargs[0]
        args = list(fargs[1:]) + [str(i)+'='+str(fkwargs[i]) for i in fkwargs]
    else:
        if len(sys.argv)==1: command = None
        else:  command = sys.argv[1]
        args = sys.argv[2:]

    
    # if no argument at all, print script doc
    if command==None:

        from . import version
        commandstring = 'Available commands: %s' %', '.join(sorted(dcommands))
        commandstring = tools.wrap(commandstring, 72, 8)

        print __license__
        print """Version number: %s

Usage:
        egglib <command> [<options> [debug] [quiet]]
        
    %s""" %(version, commandstring)
        return

    else:
        
        if command not in dcommands:

            sys.stderr.write('error - command not available\n')

        else:

            cmd = dcommands[command]
        
            # runs a particular command

            # if command requested without argument, prints specific doc
            if len(args)==0:
                print cmd.doc()

            # runs the command
            else:
                
                obj = cmd()
                flags, options = obj.process_cmdline_arguments(args)
                
                try:
                    obj.run(*flags, **options)
                    
                except Exception, e:
                    a= re.match('<type \'exceptions.(.+)\'>',str(e.__class__))
                    if a: s= '[%s] ' %a.group(1)
                    else: s=''
                    sys.stderr.write( 'An error occurred: %s%s\n' %(s,str(e)))
                    if obj.debug: raise
                    sys.exit(-1)


########################################################################

if __name__ == '__main__':
    execute()

