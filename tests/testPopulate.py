#!/usr/bin/python

import sys
import os

import math
import unittest

# get parent directory name to import module to be tested
parentdir = os.path.abspath(os.pardir)
sys.path.append(parentdir)

# import the module to be tested
from populate import generate

class testPopulate(unittest.TestCase):
    def setUp(self):
        self.populate = generate(1038)

    def test_double_sided_exp_zero(self):
        result = self.populate. _double_sided_exp(0.0, origin=0.0)
        self.assertEqual(result, 0.0)

    def test_double_sided_exp_origin(self):
        result = self.populate. _double_sided_exp(0.0, origin=10.0)
        self.assertEqual(result, 10.0)

    def test_double_sided_exp_scale(self):
        # this one is more tricky since there are randoms involved
        # simple test - check in right range
        result =  self.populate. _double_sided_exp(10.0, origin=0.0)
        self.assertTrue(math.fabs(result) < 10**9  and math.fabs(result) > 0)
