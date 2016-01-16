#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
from const import MyConst

class CoordSys(object):

    def __init__(self, **kw):
        #function self = CoordSys(num_cells,arrays_sizes,units)

        self.num_cells = kw.get('num_cells', 1)
        self.arrays_sizes = kw.get('arrays_sizes', 1)
        self.units = kw.get('units', 'au')

        if self.units=='au':
            self.lattice_const=0.5431*1e-9/MyConst.ab;
        else:
            self.lattice_const=0.5431*1e-9;

        self.coord_sizes = self.num_cells*self.arrays_sizes
        dist = self.num_cells * self.lattice_const
        self.coord_stps = dist/(self.coord_sizes)

        self.origin_cells = 1

        self.origin_inds = self.arrays_sizes * (self.origin_cells - 1) + 1
        self.coord_limits = [0, 0]
        self.coord_limits[0] = -self.origin_inds * self.coord_stps
        self.coord_limits[1] = -self.origin_inds * self.coord_stps + self.num_cells * self.lattice_const


    def set_origin_cells(self,origin_cells):

        self.origin_cells = origin_cells
        self.origin_inds=self.arrays_sizes*(self.origin_cells-1)+1
        self.coord_limits[0] = -self.origin_inds*self.coord_stps+self.coord_stps
        self.coord_limits[1] = -self.origin_inds*self.coord_stps+self.num_cells*self.lattice_const


    def x(self):
        x = np.linspace(self.coord_limits[0],self.coord_limits[1],self.coord_sizes, endpoint = True)
        return x


if __name__ == "__main__":

    cs = CoordSys(num_cells=2, arrays_sizes=5, units='au')
    cs.set_origin_cells(2)
    x=cs.x()

    print cs.coord_sizes

    plt.plot(x, marker='*')
    plt.show()

