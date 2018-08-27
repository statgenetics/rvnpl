from setuptools import setup, Extension
from src import NAME, VERSION
import sys,os

try:
    __import__('cstatgen')
except ImportError:
    sys.stderr.write('Installing cstatgen ...\n')
    cwdPath = os.getcwd()
    os.chdir('./cstatgen/')
    cmd = "python2.7 setup.py install {}".format(" ".join(sys.argv[2:]))
    os.system("{} > /dev/null".format(cmd))
    os.chdir(cwdPath)

setup(
	name = NAME,
	version = VERSION,
      	description = "A novel approach to use sequence data for nonparametric linkage analysis",
      	author = "Linhai Zhao",
      	packages = [NAME,'RVNPLcpp','SEQLinkage'],
      	scripts = ['src/rvnpl'],
      	package_dir = {NAME:'src','RVNPLcpp':'cppextend','SEQLinkage':'seqlink'},
      	install_requires = ['scipy', 'numpy','sympy','prettyplotlib', 'tornado', 
			    'brewer2mpl','faulthandler', 'matplotlib'],
      	ext_modules = [
		Extension('RVNPLcpp.cconv2',
			['cppextend/cconv2.cpp'],
			libraries=['boost_python'],
			extra_compile_args=['-std=c++11']
		),
		Extension('RVNPLcpp.cmissingparents',
			['cppextend/missingparents.cpp'],
			libraries=['boost_python'],
			extra_compile_args=['-std=c++11']
		),
		Extension('RVNPLcpp.cmissing_infer',
			['cppextend/missing_infer.cpp'],
			libraries=['boost_python'],
			extra_compile_args=['-std=c++11']
		),
		Extension('RVNPLcpp.cpostInv',
			['cppextend/postInv.cpp'],
			libraries=['boost_python'],
			extra_compile_args=['-std=c++11']
		),
		Extension('RVNPLcpp.ibd_m_cpp',
			['cppextend/ibd_m.cpp'],
			libraries=['boost_python'],
			extra_compile_args=['-std=c++11']
		),
		Extension('RVNPLcpp.ibd_rv_cpp',
			['cppextend/ibd_m_rv.cpp'],
			libraries=['boost_python'],
			extra_compile_args=['-std=c++11']
		),
		Extension('RVNPLcpp.sall_cpp',
			['cppextend/sall_rv_complete.cpp'],
			libraries=['boost_python'],
			extra_compile_args=['-std=c++11']
		),
	]
)
