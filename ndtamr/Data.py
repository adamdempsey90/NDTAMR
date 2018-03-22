import numpy as np
import copy

class SimpleTest2D():
    cols=['x','y','d']
    def __init__(self,coords=(0,0),file=None,data=None):
        self.xc,self.yc = coords
        if data is not None:
            for c in self.cols:
                setattr(self,c,getattr(data,c))
        elif file is None:
            for c in self.cols:
                setattr(self,c,self.func())
        else:
            for c in self.cols:
                setattr(self,c,file[c][...])
    def func(self):
        
        cx = .65 + .5* 2.**(-8)
        cy = .65 + + .5* 2.**(-8)
        s = .1
        res = np.exp(-((self.xc-cx)**2+(self.yc-cy)**2)/(2*s**2))
        cx = .3 + .5* 2.**(-8)
        cy = .3 + .5* 2.**(-8)
        s = .1
        res += np.exp(-((self.xc-cx)**2+(self.yc-cy)**2)/(2*s**2))
        
    
        return res
    def get_refinement_data(self):
        return self.d
    def copy(self):
        import copy
        return copy.copy(self)
    def save(self,file):
        grp = file.create_group('Data')
        for c in self.cols:
            grp.create_dataset(c,data=getattr(self,c))
    def __rmul__(self,val):
        newdat = copy.deepcopy(self)
        newdat.d *= val
        return newdat
    def __mul__(self,val):
        newdat = copy.deepcopy(self)
        newdat.d *= val
        return newdat
    
    def __add__(self,val):
        newdat = copy.deepcopy(self)
        newdat.d += val.d
        return newdat
    def __radd__(self,val):
        newdat = copy.deepcopy(self)
        newdat.d += val.d
        return newdat
class SpiralTest2D():
    cols=['x','y','d']
    def __init__(self,coords=(0,0),file=None,data=None):
        self.xc,self.yc = coords
        if data is not None:
            for c in self.cols:
                setattr(self,c,getattr(data,c))
        elif file is None:
            for c in self.cols:
                setattr(self,c,self.func())
        else:
            for c in self.cols:
                setattr(self,c,file[c][...])
    def func(self):
        r = np.sqrt( self.xc**2 + self.yc**2)
        p = np.arctan2(self.yc,self.xc)
        
        ps = np.log(r/1)/.2
        xs = r*np.cos(ps)
        ys = r*np.sin(ps)
        res = np.exp(-((self.xc-xs)**2 + (self.yc-ys)**2)/(2*.3**2))
        if np.isnan(res) or np.isinf(res):
            res = 1
        return res
    def get_refinement_data(self):
        return self.d
    def copy(self):
        import copy
        return copy.copy(self)
    def save(self,file):
        grp = file.create_group('Data')
        for c in self.cols:
            grp.create_dataset(c,data=getattr(self,c))
    def __rmul__(self,val):
        newdat = copy.deepcopy(self)
        newdat.d *= val
        return newdat
    def __mul__(self,val):
        newdat = copy.deepcopy(self)
        newdat.d *= val
        return newdat
    
    def __add__(self,val):
        newdat = copy.deepcopy(self)
        newdat.d += val.d
        return newdat
    def __radd__(self,val):
        newdat = copy.deepcopy(self)
        newdat.d += val.d
        return newdat
class SpiralTest3D():
    cols=['x','y','z','d']
    def __init__(self,coords=(0,0,0),file=None,data=None):
        self.xc,self.yc,self.zc = coords
        self.d = 0.
        if data is not None:
            for c in self.cols:
                setattr(self,c,getattr(data,c))
        elif file is None:
            for c in self.cols:
                setattr(self,c,self.func())
        else:
            for c in self.cols:
                setattr(self,c,file[c][...])
    def func(self):
        r = np.sqrt( self.xc**2 + self.yc**2)
        p = np.arctan2(self.yc,self.xc)
        
        ps = np.log(r/1)/.2
        xs = r*np.cos(ps)
        ys = r*np.sin(ps)
        res = np.exp(-((self.xc-xs)**2 + (self.yc-ys)**2)/(2*.3**2))
        if np.isnan(res) or np.isinf(res):
            res = 1
        return res * np.exp(-self.zc**2/(2*.1**2))
    def get_refinement_data(self):
        return self.d
    def copy(self):
        import copy
        return copy.copy(self)
    def save(self,file):
        grp = file.create_group('Data')
        for c in self.cols:
            grp.create_dataset(c,data=getattr(self,c))
    def __rmul__(self,val):
        newdat = copy.deepcopy(self)
        newdat.d *= val
        return newdat
    def __mul__(self,val):
        newdat = copy.deepcopy(self)
        newdat.d *= val
        return newdat
    
    def __add__(self,val):
        newdat = copy.deepcopy(self)
        newdat.d += val.d
        return newdat
    def __radd__(self,val):
        newdat = copy.deepcopy(self)
        newdat.d += val.d
        return newdat