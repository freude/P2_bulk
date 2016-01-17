#!/usr/bin/python

"""The module of the ConfSet class."""

import numpy as np
from confset import ConfSet


class ConfSet_total(object):
    """The class determine a set of configurations for FCI method."""

    def __init__(self, a, b):
        self.str_a = a
        self.str_b = b
        self.num_conf = a.num_conf*b.num_conf
        self.conf = np.array(self.num_conf,2)

        for j1 in xrange(a.num_conf):
            for j2 in xrange(b.num_conf):
                self.conf[a.num_conf*j2+j1,0]=a.num_conf*j2+j1
                self.conf[a.num_conf*j2+j1,1]=j1
                self.conf[a.num_conf*j2+j1,2]=j2


#---------------------------------------------------------------------

if __name__ == "__main__":

    cs = ConfSet(Nel=1, Norb=6, M=0)
    cs.print_conf()
