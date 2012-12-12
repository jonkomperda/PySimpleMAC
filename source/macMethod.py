from mac import *
from numpy import *

class domain():
    """Inializes a domain within which the fluid resides"""
    def __init__(self, r=0.1):
        self.fluid = fluid
        dim.r = self.r = r
        
        self.nx = int(dim.nx)
        self.ny = int(dim.ny)
        self.dx = float(dim.dx)
        self.dy = float(dim.dy)
        
        self.x = [x*self.dx for x in range(self.nx)]
        self.y = [y*self.dy for y in range(self.ny)]
        self.z = 0
    


class fluid():
    """Initializes a fluid for use with the solver"""
    def __init__(self, domain, re=100):
        dim.re = self.re = re
        self.domain = domain
    
    
    def u(self):
        array = swapaxes(solver.u,0,1)
        return array
    
    
    def v(self):
        array = swapaxes(solver.v,0,1)
        return array
    
    
    def p(self):
        array = swapaxes(solver.p,0,1)
        return array
    
    def w(self):
        array = swapaxes(solver.vort(),0,1)
        return array
    
    def setRe(self,re):
        dim.re = re
    


class MACSolver():
    """Marker and Cell method solver for incompressible fluids"""
    def __init__(self, fluid, domain, steps=5000):
        self.domain = domain
        self.fluid = fluid
        self.steps = steps
        # zero out the timestep
        dim.t = 0
        
        # initialize the domain with zeros
        solver.initzero()
        
        # apply initial conditions
        solver.ghostcondition()
        solver.lidcondition()
        
        # calculate the timestep for initialization
        solver.calctstep()
    
    def reset(self):
        solver.initzero()
        solver.ghostcondition()
        solver.lidcondition()
        solver.calctstep()
    
    def step(self):
        solver.ghostcondition()
        solver.lidcondition()
        solver.calctstep()
        solver.calcfngn()
        solver.calcqn()
        solver.poisson()
        #solver.parpoisson()
        solver.calcvel()
    
    def run(self,pinterval=1000):
        self.pinterval = pinterval
        
        for k in range(self.steps):
            self.step()
            if k%pinterval==0:
                self.visitPrint(k)
            
    def visitPrint(self,k):
        import pyvtk
        uVel = pyvtk.Scalars(self.fluid.u().flatten(),name='u')
        vVel = pyvtk.Scalars(self.fluid.v().flatten(),name='v')
        totp = pyvtk.Scalars(self.fluid.p().flatten(),name='p')
        celldat = pyvtk.PointData(uVel,vVel,totp)
        grid = pyvtk.RectilinearGrid(self.domain.x,self.domain.y,self.domain.z)
        vtk = pyvtk.VtkData(grid,celldat)
        vtk.tofile('data%i'%k)

if __name__ == '__main__':
    import pyvtk
    square = domain(r=0.05)
    
    air = fluid(square, re=8000)
    
    method = MACSolver(air, square, steps=50001)
    
    method.run()
