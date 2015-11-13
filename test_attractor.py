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
    """Test setup and default/overide initial parameters.  This test evaluates if the initial parameters can be modified to accept something other than the default values.  """
    def setup(self):
        s = randint(1,10)
        p = randint(1,10)
        b = randint(1,10)
        attr = Attractor(s, p, b)
        attr.s_t = s
        attr.p_t = p
        attr.b_t = b
        attr = Attractor(attr.s_t, attr.p_t, attr.b_t)
        assert attr.params[0] == s, "\nError in assignment of value for s parameter"
        assert attr.params[1] == p, "\nError in assignment of value for p parameter"
        assert attr.params[2] == b, "\nError in assignment of value for b parameter"
        assert attr.start == 0, "\nError in initial value for start parameter"
        assert attr.end == 80, "\nError in initial value for end parameter"
        assert attr.points == 10000, "\nError in initial value for points parameter"

    def test_euler(self):
        attr = Attractor()
        attr.evolve([10,10,10],1)
        assert attr.solution['x'].count() > 0, "\nError in solution for x using Euler method, Order=1"
        assert attr.solution['y'].count() > 0, "\nError in solution for y using Euler method, Order=1"
        assert attr.solution['z'].count() > 0, "\nError in solution for z using Euler method, Order=1"

    def test_rk2(self):
        attr = Attractor()
        attr.evolve([10,10,10],2)
        assert attr.solution['x'].count() > 0, "\nError in solution for x using rk2 method, Order=2"
        assert attr.solution['y'].count() > 0, "\nError in solution for y using rk2 method, Order=2"
        assert attr.solution['z'].count() > 0, "\nError in solution for z using rk2 method, Order=2"

    def test_rk4(self):
        attr = Attractor()
        attr.evolve([10,10,10],4)
        assert attr.solution['x'].count() > 0, "\nError in solution for x using rk4 method, Order=4"
        assert attr.solution['y'].count() > 0, "\nError in solution for y using rk4 method, Order=4"
        assert attr.solution['z'].count() > 0, "\nError in solution for z using rk4 method, Order=4"