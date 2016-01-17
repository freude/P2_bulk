#!/usr/bin/python

import itertools
import math
import random
import numpy as np

class McGint(object):

    def __init__(self,X,F):
        self.X=X
        self.F=F
        self.L=len(X[0])

    @property
    def sampler(self):
        while True:

            x = random.randint(0, self.L)
            y = random.randint(0, self.L)
            z = random.randint(0, self.L)

            yield (self.X[x],self.X[y],self.X[z])


    def integrand(self,x):

        return (x[0]**2 + x[1]**2+ x[1]**2)

    def integrate(measure=1.0, n=100):
        """Sum elements and elements squared"""
        total = 0.0
        total_sq = 0.0
        for x in itertools.islice(self.sampler, n):
            f = self.integrand(x)
            total += f
            total_sq += (f**2)
        # Return answer
        sample_mean = total/n
        sample_var = (total_sq - ((total/n)**2)/n)/(n-1.0)
        return (measure*sample_mean, measure*math.sqrt(sample_var/n))


if __name__ == "__main__":

    result, error = mcint.integrate(self.integrand, self.sampler(), measure=math.pi/4)
