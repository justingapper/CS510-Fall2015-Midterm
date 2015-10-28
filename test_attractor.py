from attractor import Attractor
from random import uniform, randint
from math import sqrt
from nose import with_setup

###
# Test Suite for specified attractor interface
#
# Run with the command: "nosetests test_attractor.py"
###

#import unittest

class attr_test:
#   """Test setup and default/overide initial parameters"""
    def test_setup(self):
        attr = Attractor()
        assert attr.params[0] == 10
        assert attr.params[1] == 28
        assert attr.params[2] == 8.0/3.0
        assert attr.start == 0
        assert attr.end == 80
        assert attr.points == 10000
    
    """Test euler method outputs results"""
    
    def test_euler(self):
        attr = Attractor()
        attr.evolve([10,10,10],1)
        assert attr.solution['x'].count() > 0
        #assert attr.solution['y'].count() > 0
        #assert attr.solution['z'].count() > 0

    """Test RK2 method outputs results"""
    
    def test_euler(self):
        attr = Attractor()
        attr.evolve([10,10,10],2)
        assert attr.solution['x'].count() > 0
        #assert attr.solution['y'].count() > 0
        #assert attr.solution['z'].count() > 0
        
        """Test RK4 method outputs results"""
        
    def test_euler(self):
        attr = Attractor()
        attr.evolve([10,10,10],4)
        assert attr.solution['x'].count() > 0
        #assert attr.solution['y'].count() > 0
        #assert attr.solution['z'].count() > 0