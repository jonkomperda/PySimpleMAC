from mac import *
from numpy import *

class domain():
    """
    Defines an object that behaves as a computatational domain.
    """
    def __init__(self, r=0.1):
        """
        Sets convergence criteria, domain sizes, step sizes, and directional arrays upon initialization.
        """
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
    """
    Defines a fluid object. Contains fluid parameters.
    """
    def __init__(self, domain, re=100):
        """
        Sets Reynolds number upon initialization. Also inherits the domain as a subset of the fluid.
        """
        dim.re = self.re = re
        self.domain = domain
    
    
    def u(self):
        """
        Accesses the Fortran persistent data and returns U-Velocity in the form of a Numpy array.
        """
        array = swapaxes(solver.u,0,1)
        return array
    
    
    def v(self):
        """
        Accesses the Fortran persistent data and returns V-Velocity in the form of a Numpy array.
        """
        array = swapaxes(solver.v,0,1)
        return array
    
    
    def p(self):
        """
        Accesses the Fortran persistent data and returns the Pressure in the form of a Numpy array.
        """
        array = swapaxes(solver.p,0,1)
        return array
    
    def w(self):
        """
        Calls a Fortran function to calculate Vorticity and returns it's value in the form of a Numpy array.
        """
        array = swapaxes(solver.vort(),0,1)
        return array
    
    def setRe(self,re):
        """
        Sets the Reynolds number in real time for the Fortran program.
        """
        dim.re = re
    


class MACSolver():
    """
    A generalized Marker and Cell method solver for the incompressible Navier-Stokes equations.
    """
    def __init__(self, fluid, domain, steps=5000):
        """
        Given a fluid and domain, this initializes our solver. All preprocessing is performed.
        """
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
        """
        Resets the solver to the initial condition state.
        """
        solver.initzero()
        solver.ghostcondition()
        solver.lidcondition()
        solver.calctstep()
    
    def step(self):
        """
        Takes a single computational step. Data may be accessed from the fluid class.
        """
        solver.ghostcondition()
        solver.lidcondition()
        solver.calctstep()
        solver.calcfngn()
        solver.calcqn()
        solver.poisson()
        #solver.parpoisson()
        solver.calcvel()
    
    def run(self,pinterval=1000):
        """
        Runs the method as a 'black box' for a predetermined number of steps, 
        set upon initialization of the solver. Prints using PyVTK at predetermined intervals.
        """
        self.pinterval = pinterval
        
        for k in range(self.steps):
            self.step()
            if k%pinterval==0:
                self.visitPrint(k)
            
    def visitPrint(self,k):
        """
        Uses PyVTK to write out data in a VisIt Visualization Tool capable format.
        """
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
