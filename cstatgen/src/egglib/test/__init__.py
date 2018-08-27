"""
This subpackage contains test functions organized in modules matching
the module of the EggLib Python package. Use the top-level
:meth:`test_all` function to launch all tests, as in the following::

    >>> from egglib import test
    >>> test.test_all()

Note that all optional
dependencies are required to run the tests, and that test procedures
make use of unsecure local temporary files that might result of user
files if they accidently have the same name as one of the test files
("test.fas", for example). If you wish to run the test suite, we
strongly advice you to run them in an empty directory.
"""

import t_egglib_binding, t_data, t_tools, t_wrappers
import t_simul, t_fitmodel, t_utils


def test_all():
    """
    Call all :meth:`test_all` functions of test modules.
    """

    t_egglib_binding.test_all()
    t_data.test_all()
    t_tools.test_all()
    t_simul.test_all()
    t_wrappers.test_all()
    t_fitmodel.test_all()
    t_utils.test_all()
