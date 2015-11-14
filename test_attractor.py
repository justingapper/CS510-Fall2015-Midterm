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
        """Test that initial conditions are setup correctly for default values (start, end, points) and default over-ride (s, p, q).  I chose this test so that the initial conditions can be verified and if there is an error identify which parameter the problem is associated with."""
        s = randint(1,10)
        p = randint(1,10)
        b = randint(1,10)
        attr = Attractor(s, p, b)
        attr.s_t = s
        attr.p_t = p
        attr.b_t = b
        attr = Attractor(attr.s_t, attr.p_t, attr.b_t)
        print "Assert value for s parameter over-ride"
        assert attr.params[0] == s
        print "    PASSED!!!!"
        print "Assert value for p parameter over-ride"
        assert attr.params[1] == p
        print "    PASSED!!!!"
        print "Assert value for b parameter over-ride"
        assert attr.params[2] == b
        print "    PASSED!!!!"
        print "Assert default value for start parameter"
        assert attr.start == 0
        print "    PASSED!!!!"
        print "Assert default value for end parameter"
        assert attr.end == 80
        print "    PASSED!!!!"
        print "Assert default value for points parameter"
        assert attr.points == 10000
        print "    PASSED!!!!"

    def test_euler(self):
        """Test Euler method yields results.  I chose this test to verify that the Euler Method yields results.  If there is an error is it associated with the calculation of the x-value, y-value, or z-value."""
        attr = Attractor()
        attr.evolve([10,10,10],1)
        print "Assert output for Euler method, x parameter."
        assert attr.solution['x'].count() > 0
        print "    PASSED!!!!"
        print "Assert output for Euler method, y parameter."
        assert attr.solution['y'].count() > 0
        print "    PASSED!!!!"
        print "Assert output for Euler method z parameter."
        assert attr.solution['z'].count() > 0
        print "    PASSED!!!!"

    def test_rk2(self):
        """Test second-order RK method yields results.  I chose this test to verify that the 2nd order RK method yields results and, if not, is the error associated with the calculation of the x-value, y-value, or z-value."""
        attr = Attractor()
        attr.evolve([10,10,10],2)
        print "Test for error in solution for x using rk2 method, Order=2"
        assert attr.solution['x'].count() > 0
        print "    PASSED!!!!"
        print "Test for error in solution for y using rk2 method, Order=2"
        assert attr.solution['y'].count() > 0
        print "    PASSED!!!!"
        print "Test for error in solution for z using rk2 method, Order=2"
        assert attr.solution['z'].count() > 0
        print "    PASSED!!!!"

    def test_rk4(self):
        """Test fourth-order RK method yields results.  I chose this test to verify that the 4th order RK method yields results and, if not, is the error associated with the calculation of the x-value, y-value, or z-value."""
        attr = Attractor()
        attr.evolve([10,10,10],4)
        print "Test for error in solution for x using rk4 method, Order=4"
        assert attr.solution['x'].count() > 0
        print "    PASSED!!!!"
        print "Test for error in solution for y using rk4 method, Order=4"
        assert attr.solution['y'].count() > 0
        print "    PASSED!!!!"
        print "Test for error in solution for z using rk4 method, Order=4"
        assert attr.solution['z'].count() > 0
        print "    PASSED!!!!"