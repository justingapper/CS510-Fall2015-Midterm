from attractor import Attractor
from random import uniform, randint
from math import sqrt
from nose import with_setup
from random import randint

###
# Test Suite for specified attractor interface
#
# Run with the command: "nosetests test_attractor.py"
###

#import unittest

class attr_test:
    """Test setup and default/overide initial parameters"""
    def setup(self):
        s = randint(1,10)
        p = randint(1,10)
        b = randint(1,10)
        attr = Attractor(s, p, b)
        attr.s_t = s
        attr.p_t = p
        attr.b_t = b
        attr = Attractor(attr.s_t, attr.p_t, attr.b_t)
        assert attr.params[0] == s
        assert attr.params[1] == p
        assert attr.params[2] == b
        assert attr.start == 0
        assert attr.end == 80
        assert attr.points == 10000

    def test_euler(self):
        attr = Attractor()
        attr.evolve([10,10,10],1)
        assert attr.solution['x'].count() > 0
        assert attr.solution['y'].count() > 0
        assert attr.solution['z'].count() > 0

    def test_rk2(self):
        attr = Attractor()
        attr.evolve([10,10,10],2)
        assert attr.solution['x'].count() > 0
        assert attr.solution['y'].count() > 0
        assert attr.solution['z'].count() > 0

    def test_rk4(self):
        attr = Attractor()
        attr.evolve([10,10,10],4)
        assert attr.solution['x'].count() > 0
        assert attr.solution['y'].count() > 0
        assert attr.solution['z'].count() > 0