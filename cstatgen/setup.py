# $File: setup.py $
# $LastChangedDate:  $
# $Rev:  $
# Copyright (c) 2014, Gao Wang <ewanggao@gmail.com>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)

from distutils.core import setup, Extension
# from distutils.dep_util import newer
try:
   from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:
   from distutils.command.build_py import build_py

# monkey-patch for parallel compilation
import multiprocessing, multiprocessing.pool

def compile_parallel(
        self,
        sources,
        output_dir=None,
        macros=None,
        include_dirs=None,
        debug=0,
        extra_preargs=None,
        extra_postargs=None,
        depends=None):

    # Copied from distutils.ccompiler.CCompiler
    macros, objects, extra_postargs, pp_opts, build = self._setup_compile(
        output_dir, macros, include_dirs, sources, depends, extra_postargs)
    cc_args = self._get_cc_args(pp_opts, debug, extra_preargs)
    #
    def _single_compile(obj):

        try:
            src, ext = build[obj]
        except KeyError:
            return
        self._compile(obj, src, ext, cc_args, extra_postargs, pp_opts)
    # convert to list, imap is evaluated on-demand
    list(multiprocessing.pool.ThreadPool(multiprocessing.cpu_count()).imap(_single_compile, objects))
    return objects

import distutils.ccompiler
distutils.ccompiler.CCompiler.compile=compile_parallel
   
import sys, os, subprocess, platform
from glob import glob
from src import NAME, VERSION

if sys.platform == "win32":
    sys.exit('Windows OS is not supported.')


# use ccache to speed up build
try:
    if subprocess.call(['ccache'], stderr = open(os.devnull, "w")):
        os.environ['CC'] = 'ccache clang -Qunused-arguments' if platform.system() == 'Darwin' else 'ccache gcc'
except OSError:
    pass
    
SWIG_OPTS = ['-c++', '-python', '-O', '-shadow', '-keyword',
             '-w-511', '-w-509', '-outdir', '.']

if sys.version_info.major == 2:
    PYVERSION = 'py2'
else:
    SWIG_OPTS.append('-py3')
    PYVERSION = 'py3'
#
def getfn(fn, prefix = "src/umich"):
    if type(fn) is list:
        return sum([glob(os.path.join(prefix, x)) if "*" in x else [os.path.join(prefix, x)] for x in fn], [])
    else:
        return glob(os.path.join(prefix, fn)) if "*" in fn else os.path.join(prefix, fn)
#
HEADER = getfn("*.hpp", "src")
CPP = getfn("*.cpp", "src")
WRAPPER_CPP = getfn('cstatgen_{0}.cxx'.format(PYVERSION), "src")
WRAPPER_PY = getfn('cstatgen_{0}.py'.format(PYVERSION), "src")
WRAPPER_I = getfn('cstatgen.i', "src")

# generate wrapper files
try:
    ret = subprocess.call(['swig', '-python', '-external-runtime', 'swigpyrun.h'], shell=False)
    if ret != 0:
        sys.exit('Failed to generate swig runtime header file. Is "swig" installed?')
    #
    if (not os.path.isfile(WRAPPER_PY) or not os.path.isfile(WRAPPER_CPP) or \
        os.path.getmtime(WRAPPER_CPP) < max([os.path.getmtime(x) for x in [WRAPPER_I] + HEADER + CPP])):
        ret = subprocess.call(['swig'] + SWIG_OPTS + ['-o', WRAPPER_CPP, WRAPPER_I], shell=False)
        if ret != 0:
            sys.exit('Failed to generate C++ extension.')
        os.rename('cstatgen.py', WRAPPER_PY)
    os.remove('swigpyrun.h')
except OSError as e:
    sys.exit('Failed to generate wrapper file: {0}'.format(e))

# Under linux/gcc, lib stdc++ is needed for C++ based extension.
libs = ['stdc++'] if sys.platform == 'linux2' else []
compile_args_umich = ["-O3", "-std=c++11", "-D_FILE_OFFSET_BITS=64", "-D__ZLIB_AVAILABLE__"]
link_args = ["-lm", "-lz"]
#
UMICH_FILES = getfn(["clusters/*.cpp", "libsrc/*.cpp", "merlin/*.cpp",
                     "pdf/*.cpp", "klib/*.c", "general/*.cpp", "vcf/*.cpp"])
CSTATGEN_MODULE = [
    Extension('{}._cstatgen'.format(NAME),
              sources = [WRAPPER_CPP] + CPP + UMICH_FILES,
              extra_compile_args = compile_args_umich,
    	      extra_link_args = link_args,
              libraries = libs,
              library_dirs = [],
              include_dirs = getfn(["general", "klib", "vcf", "clusters", "libsrc", "merlin", "pdf"]) + ["src"]
              )
]
#
compile_args_egglib = ["-O3", "-std=c++11", "-UHAVE_LIBBPP_SEQ", "-UHAVE_LIBBPP_CORE", "-UHAVE_LIBGSLCBLAS"]
# exclude two modules due to lack of gsl and bio++; egglib should have used macro to control for it, though
EGGLIB_FILES = [x for x in getfn("egglib-cpp/*.cpp", prefix = "src/egglib") if not (x.endswith("ABC.cpp") or x.endswith("BppDiversity.cpp"))]
EGGLIB_MODULE = [
    Extension('_egglib_binding',
              sources = ["src/egglib/egglib_binding.cpp"] + EGGLIB_FILES,
              extra_compile_args = compile_args_egglib,
    	      extra_link_args = link_args,
              libraries = libs,
              library_dirs = [],
              include_dirs = getfn(["egglib/egglib-cpp", "egglib"], prefix = "src")
              )
]

setup(name = NAME,
    version = VERSION,
    description = "Gao Wang's statgen library",
    author = "Gao Wang",
    packages = [NAME, NAME + ".egglib"],
    package_dir = {NAME:'src', NAME + ".egglib":'src/egglib'},
    package_data = {NAME + ".egglib":['apps.conf']},
    cmdclass = {'build_py': build_py},
    ext_modules = CSTATGEN_MODULE + EGGLIB_MODULE
)
