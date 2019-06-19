#!/usr/bin/env python2
#
# Copyright (c) Linhai Zhao
# Distributed under the terms of the MIT License.

import unittest
import pickle

class TestParser(unittest.TestCase):
    def setUp(self):
        # things to do before running any unit test
        data = dict(a=1, b=2)
        with open('data/trivial.pkl', 'wb') as f:
            pickle.dump(data, f)

    def tearDown(self):
        # things to do after all tests are complete
        pass

    def testExample(self):
        result = dict(a=1,b=2)
        with open('data/trivial.pkl', 'rb') as f:
            saved_result = pickle.load(f)
        self.assertEqual(result, saved_result)

if __name__ == '__main__':
    #suite = unittest.defaultTestLoader.loadTestsFromTestCase(TestParser)
    # unittest.TextTestRunner(, suite).run()
    unittest.main()
