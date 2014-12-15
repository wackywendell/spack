from . import *
from unittest import TestCase

import numpy as np
from numpy import sqrt, array
from numpy.testing import assert_allclose

class PackingTest3D(TestCase):
    def setUp(self):
        self.L = 2.0066668050219723
        self.diameters = array([ 0.96,  0.97,  0.98,  0.99,  1.  ,  1.01,  1.02,  1.03,  1.04])
        self.locs = array([[ 1.40776762,  1.26647724,  0.73389219],
                           [ 0.58704249,  2.11399   ,  1.52956579],
                           [ 1.75917911,  0.54290089,  1.27577478],
                           [ 2.13750384,  0.87508242,  0.21938647],
                           [ 1.07283961,  0.87692084,  1.9060841 ],
                           [ 0.09550267,  1.94404465,  0.56463369],
                           [ 1.07636871,  2.1942971 ,  0.63752152],
                           [ 0.49922725,  1.20002224,  1.13360082],
                           [-0.27724757,  1.62152603,  1.67262247]])
        
        self.pack = Packing(self.locs, self.diameters, L=self.L)
        
    def test_contacts(self):
        self.assertEqual(self.pack.contacts(), (25, 25, 0))
    
    def test_DM(self):
        dm_freqs_exp = array([  7.01815161e-01,   6.71651700e-01,   6.41648465e-01,
                                6.37510086e-01,   6.13838144e-01,   5.84421942e-01,
                                5.80378277e-01,   5.65948152e-01,   5.53149405e-01,
                                5.13210148e-01,   4.62667794e-01,   4.47085008e-01,
                                4.28671950e-01,   3.95836233e-01,   3.72709231e-01,
                                3.31390557e-01,   2.98647061e-01,   2.66650336e-01,
                                2.48403245e-01,   2.30577214e-01,   1.97602540e-01,
                                1.53048522e-01,   8.36861967e-02,   1.12936087e-01,
                                5.69943494e-09,   6.12301687e-09,   1.11078076e-09])
        dm_freqs_exp.sort()
        dm_freqs = self.pack.DM_freqs()
        dm_freqs.sort()
        
        atol, rtol = 1e-7, 1e-7
        for f1, f2 in zip(dm_freqs, dm_freqs_exp):
            diff = abs(f1 - f2)
            maxdiff = atol + rtol * abs(f2)
            if diff > maxdiff:
                print(f1, f2, atol, rtol*abs(f2))
        
        assert_allclose(dm_freqs, dm_freqs_exp, rtol=rtol, atol=atol, verbose=True)
