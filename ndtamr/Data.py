import numpy as np
import copy

class GenericData():
    def __init__(self,coords=(0,0),file=None,data=None):
        self.coords = [c for c in coords]
        if data is not None:
            for c in self.data_cols:
                setattr(self,c,getattr(data,c))
        elif file is None:
            for c in self.data_cols:
                setattr(self,c,self.func())
        else:
            for c in self.data_cols:
                setattr(self,c,file[c][...])
    def copy(self):
        import copy
        return copy.copy(self)
    def save(self,file):
        grp = file.create_group('Data')
        for c in self.data_cols:
            grp.create_dataset(c,data=getattr(self,c))
    def __rmul__(self,val):
        newdat = copy.deepcopy(self)
        for d in newdat.data_cols:
            try:
                setattr(newdat,d,getattr(newdat,d)*getattr(val,d))
            except AttributeError:
                setattr(newdat,d,getattr(newdat,d)*val)
        return newdat
    def __mul__(self,val):
        return self.__rmul__(val)
    def __radd__(self,val):
        newdat = copy.deepcopy(self)
        for d in newdat.data_cols:
            setattr(newdat,d,getattr(newdat,d) + getattr(val,d))
        return newdat
    def __add__(self,val):
        return self.__radd__(val)
    def __rsub__(self,val):
        newdat = copy.deepcopy(self)
        for d in newdat.data_cols:
            setattr(newdat,d,getattr(val,d)-getattr(newdat,d))

        return newdat
    def __sub__(self,val):
        newdat = copy.deepcopy(self)
        for d in newdat.data_cols:
            setattr(newdat,d,getattr(newdat,d) - getattr(val,d))
        return newdat
    def __rtruediv__(self,val):
        newdat = copy.deepcopy(self)
        for d in newdat.data_cols:
            try:
                setattr(newdat,d,getattr(val,d)/getattr(newdat,d))
            except AttributeError:
                 setattr(newdat,d,val/getattr(newdat,d))
               
        return newdat
    def __truediv__(self,val):
        newdat = copy.deepcopy(self)
        for d in newdat.data_cols:
            try:
                setattr(newdat,d,getattr(newdat,d)/getattr(val,d))
            except AttributeError:
                 setattr(newdat,d,getattr(newdat,d)/val)
        return newdat

class Empty(GenericData):
    data_cols = ['d']
    def __init__(self,coords=(0,0),file=None,data=None):
        GenericData.__init__(self,coords=coords,file=file,data=data)
        self.d = 0
    def get_refinement_data(self):
        return self.d


class SimpleTest2D(GenericData):
    data_cols = ['d']    
    def __init__(self,coords=(0,0),file=None,data=None):
        GenericData.__init__(self,coords=coords,file=file,data=data)
    def func(self):
        xc,yc = self.coords
        cx = .65 + .5* 2.**(-8)
        cy = .65 + + .5* 2.**(-8)
        s = .1
        res = np.exp(-((xc-cx)**2+(yc-cy)**2)/(2*s**2))
        cx = .3 + .5* 2.**(-8)
        cy = .3 + .5* 2.**(-8)
        s = .1
        res += np.exp(-((xc-cx)**2+(yc-cy)**2)/(2*s**2))
        return res
    def get_refinement_data(self):
        return self.d
      
class SpiralTest2D(GenericData):
    data_cols = ['d']
    def __init__(self,coords=(0,0),file=None,data=None):
        GenericData.__init__(self,coords=coords,file=file,data=data)

    def func(self):
        xc,yc = self.coords
        r = np.sqrt( xc**2 + yc**2)
        p = np.arctan2(yc,xc)
        
        ps = np.log(r/1)/.2
        xs = r*np.cos(ps)
        ys = r*np.sin(ps)
        res = np.exp(-((xc-xs)**2 + (yc-ys)**2)/(2*.3**2))
        if np.isnan(res) or np.isinf(res):
            res = 1
        return res
    def get_refinement_data(self):
        return self.d
    
    
class SpiralTest3D(GenericData):
    data_cols = ['d']
    def __init__(self,coords=(0,0,0),file=None,data=None):
        GenericData.__init__(self,coords=coords,file=file,data=data)
        
    def func(self):
        xc,yc,zc = self.coords
        r = np.sqrt( xc**2 + yc**2)
        p = np.arctan2(yc,xc)
        
        ps = np.log(r/1)/.2
        xs = r*np.cos(ps)
        ys = r*np.sin(ps)
        res = np.exp(-((xc-xs)**2 + (yc-ys)**2)/(2*.3**2))
        if np.isnan(res) or np.isinf(res):
            res = 1
        return res * np.exp(-zc**2/(2*.4**2))
    def get_refinement_data(self):
        return self.d
