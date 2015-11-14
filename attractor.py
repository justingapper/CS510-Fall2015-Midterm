import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

class Attractor(object):
    """The Attractor class calculates the Euler, 2nd Order Runge-Kutta, 3rd Order Runge-Kutta, and 4th Order Runge-Kutta integration method for a Lorenz Attractor System."""
    def __init__(self, s=10.0, p=28.0, b=8.0/3.0, start=0.0, end=80.0, points=10000):
        """Define initial conditions:
        s = Prandtl Number, the ratio of momentum diffusivity and thermal diffusivity. (Default value = 10.0)
        p = Rayleigh Number to determine whether the heat transfer is primarily in the form of conduction or convection. (Default value = 28.0)
        b = geometric factor (Default value = 8.0/3.0)
        start, end, and points are initial, final and time increment variables, respectively. (Default values = 0.0, 80.0, and 10,000)"""
        self.s = s
        self.p = p
        self.b = b
        self.start = start
        self.end = end
        self.points = points
        self.params = np.array([s, p, b])
        self.dt = (self.end-self.start)/self.points   
        self.x = np.zeros([])
        self.y = np.zeros([])
        self.z = np.zeros([])
        self.t = np.zeros([])
        
    def fxn1(self, yt, xt, ss, ddt):
        """Define function 1 for dx(t)/dt = s[y(t) - x(t)]"""
        return ss*(yt - xt)*ddt
    
    def fxn2(self, yt, xt, zt, pp, ddt):
        """Define function 2 for dy(t)/dt = x(t)[p - z(t)] - y(t)"""
        return (xt*(pp - zt) - yt)*ddt
    
    def fxn3(self, yt, xt, zt, bb, ddt):
        """Define function 3 for dz(t)/dt = x(t)y(t) - bz(t)"""
        return (xt*yt - bb*zt)*ddt
    
    def euler(self, xyzarray):
        """Define k1 euler method.  Create numpy array for output and iterate through the differential equations fxn1, fxn2, and fxn3"""
        self.x = np.zeros(self.points + 1)
        self.y = np.zeros(self.points + 1)
        self.z = np.zeros(self.points + 1)
        self.t = np.zeros(self.points + 1)
        
        self.x[0] = xyzarray[0]
        self.y[0] = xyzarray[1]
        self.z[0] = xyzarray[2]
                
        for i in xrange(self.points):
            self.x[i+1] = self.x[i] + self.fxn1(self.y[i],self.x[i],self.s,self.dt)
            self.y[i+1] = self.y[i] + self.fxn2(self.y[i],self.x[i],self.z[i],self.p,self.dt)
            self.z[i+1] = self.z[i] + self.fxn3(self.y[i],self.x[i],self.z[i],self.b,self.dt)
            
            self.t[i+1] = self.t[i] + self.dt
                
        return (self.x, self.y, self.z)
    
    def rk2(self, xyzarray):
        """Define 2nd order Runge-Kutta method.  Create numpy array for results and iterate through the differential equations fxn1, fxn2, and fxn3.  k2 = y(t+dt/2, x(t) + k1dt/2)"""
        self.x = np.zeros(self.points + 1)
        self.y = np.zeros(self.points + 1)
        self.z = np.zeros(self.points + 1)
        self.t = np.zeros(self.points + 1)
        
        self.x[0] = xyzarray[0]
        self.y[0] = xyzarray[1]
        self.z[0] = xyzarray[2]
        
        for i in xrange(self.points):
            k1_x = self.fxn1(self.y[i],self.x[i],self.s,self.dt)
            k2_x = self.fxn1(self.y[i],k1_x*(.5*self.dt)+self.x[i],self.s,self.dt)
            self.x[i+1] = self.x[i] + (.5*(k1_x + k2_x))
            
            k1_y = self.fxn2(self.y[i],self.x[i],self.z[i],self.p,self.dt)
            k2_y = self.fxn2(self.y[i]+k1_y*(.5*self.dt),self.x[i],self.z[i],self.p,self.dt)
            self.y[i+1] = self.y[i] + (.5*(k1_y + k2_y))
            
            k1_z = self.fxn3(self.y[i],self.x[i],self.z[i],self.b,self.dt)
            k2_z = self.fxn3(self.y[i],self.x[i],self.z[i]+k1_z*(.5*self.dt),self.b,self.dt)
            self.z[i+1] = self.z[i] + (.5*(k1_z + k2_z))
            
            self.t[i+1] = self.t[i] + self.dt
               
        return (self.x,self.y,self.z)
    
    def rk4(self, xyzarray):
        """Define 4th order Runge-Kutta method.  Create numpy array for results and iterate through the differential equations fxn1, fxn2, and fxn3.
        k3 = y(t+dt/2, x(t) + k1dt/2)
        k4 = y(t+dt/2,x(t) + k2dt/2)"""
        self.x = np.zeros(self.points + 1)
        self.y = np.zeros(self.points + 1)
        self.z = np.zeros(self.points + 1)
        self.t = np.zeros(self.points + 1)
        
        self.x[0] = xyzarray[0]
        self.y[0] = xyzarray[1]
        self.z[0] = xyzarray[2]
       
        for i in xrange(self.points):
            k1_x = self.fxn1(self.y[i],self.x[i],self.s,self.dt)
            k2_x = self.fxn1(self.y[i],k1_x*(.5*self.dt)+self.x[i],self.s,self.dt)
            k3_x = self.fxn1(self.y[i],k2_x*(.5*self.dt)+self.x[i],self.s,self.dt)
            k4_x = self.fxn1(self.y[i],k3_x*self.dt+self.x[i],self.s,self.dt)
            self.x[i+1] = self.x[i] + (1.0/6.0)*(k1_x + 2*k2_x + 2*k3_x + k4_x)

            k1_y = self.fxn2(self.y[i],self.x[i],self.z[i],self.p,self.dt)
            k2_y = self.fxn2(self.y[i]+k1_y*(.5*self.dt),self.x[i],self.z[i],self.p,self.dt)
            k3_y = self.fxn2(self.y[i]+k2_y*(.5*self.dt),self.x[i],self.z[i],self.p,self.dt)
            k4_y = self.fxn2(self.y[i]+k3_y*self.dt,self.x[i],self.z[i],self.p,self.dt)
            self.y[i+1] = self.y[i] + (1.0/6.0)*(k1_y + 2*k2_y + 2*k3_y + k4_y)

            k1_z = self.fxn3(self.y[i],self.x[i],self.z[i],self.b,self.dt)
            k2_z = self.fxn3(self.y[i],self.x[i],self.z[i]+k1_z*(.5*self.dt),self.b,self.dt)
            k3_z = self.fxn3(self.y[i],self.x[i],self.z[i]+k2_z*(.5*self.dt),self.b,self.dt)
            k4_z = self.fxn3(self.y[i],self.x[i],self.z[i]+k3_z*self.dt,self.b,self.dt)
            self.z[i+1] = self.z[i] + (1.0/6.0)*(k1_z + 2*k2_z + 2*k3_z + k4_z)

            self.t[i+1] = self.t[i] + self.dt

        return (self.x, self.y, self.z)

    def evolve(self,r0=[1.0,0.0,0.0], order=4):
        """Define evolve.  Create numpy array r0 with default values.  Allow selection of euler, RK2, or RK4."""
        x0 = r0[0]
        y0 = r0[1]
        z0 = r0[2]

        if order == 1:
            self.euler(r0)
        elif order == 2:
            self.rk2(r0)
        else:
            self.rk4(r0)

        self.solution = pd.DataFrame({'time': self.t, 'x': self.x, 'y': self.y, 'z': self.z})

        return self.solution

    def save(self):
        """Define save to csv"""
        self.solution.to_csv('data.csv')

    def plotx(self):
        """Define plot x"""
        self.solution.plot(x='time',y='x', title = 'x(t) vs Time Solution')
        plt.show()

    def ploty(self):
        """Define plot y"""
        self.solution.plot(x='time',y='y', title = 'y(t) vs Time Solution')
        plt.show()

    def plotz(self):
        """Define plot z"""
        self.solution.plot(x='time',y='z', title = 'z(t) vs Time Solution')
        plt.show()

    def plotxy(self):
        """Define plot xy"""
        plt.plot(self.solution['x'],self.solution['y'])
        plt.show()

    def plotyz(self):
        """Define plot yz"""
        plt.plot(self.solution['y'],self.solution['z'])
        plt.show()

    def plotzx(self):
        """Define plot zx"""
        plt.plot(self.solution['z'],self.solution['x'])
        plt.show()

    def plot3d(self):
        """Define 3d plot"""
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot_trisurf(self.x, self.y, self.z, cmap=cm.jet, linewidth=.2)
        ax.set_xlabel('x(t)')
        ax.set_ylabel('y(t)')
        ax.set_zlabel('z(t)')
        plt.show()
        
    def plot3dwf(self):
        """Define 3d wireframe plot"""
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot_wireframe(self.x, self.y, self.z, rstride=1, cstride=1)
        ax.set_xlabel('x(t)')
        ax.set_ylabel('y(t)')
        ax.set_zlabel('z(t)')
        plt.show()