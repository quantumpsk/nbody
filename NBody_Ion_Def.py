# -*- coding: utf-8 -*-
"""
Created on Sat Apr  1 21:45:59 2017
Saturday, April Fools Day, 01/04
I'm gonna try something here with the n body problem..
I'll start with simulating a simple solar system..
or maybe an even more simple three body system..
we shall later try to scale it for any physics problem..
lets begin with the basic class structures..

Sunday, 02/04
So ive got the basic structure down, going ahead with calculating some
of the variables..
realized i have to store all vectors as ndarrays using numpy since i cant
assign lists to a scope of a ndarray.. so thats next.

Tuesday, 04/04
Well the System and its Objects are set up and i have the complete initial state.
Next is the evolution of the system..with a time step..lets brute this shit..
well that was simple, a time.step, a calculate.state in a system.evolve and
we have what we're after.
Lets get started with plotting it..(after tonights date)

Wednesday, 05/04
So after date night, I've got to pin down the exact way im storing all the
positions of the different objects over the course of the sim to be able to plot
the full path.. and apparently the best way to dynamically grow a numpy array
is to not do it, rather grow a list and then create the array from the list !

Thursday, 06/04
Ok.. its been over 24 hours now that im stuck with what seems should be an easy
problem to fix..
The issue is i need to store all the positions of the different objects. I thought
of growing a list as each new position is calculated. But for some reason each
log of each object is storing not only its own position but also that of other
objects and i dont know why that is..weird...i think my compiler is possessed.
...
So for now im temporarily giving up on this method of arrays bull to store the
positions, instead im storing it in a file, maybe itll be easier.

Friday, 07/04
early in the morning, here i am sitting for one more attempt at a position
log array, to convert to np.array, to plot. I dont know yet why append is adding
more than it should, but it is, so lets work with and ill can learn something..
...
oh man, marathon session of coding. i solved the log and plot problem. i still
dont know why append is behaving the way it is, but some ghetto work around and
i have the positions and now the plot.
Only then did i realize i hadnt accounted for the zero-distance error where
after being attracted close to each other objects then experience immense g forces
which end up scattering them..so thats next..just had a eureka moment, elastic
collisions ! thatll prevent the zero distance issue..gotta add a radius tho.
...
I've forked this version over cos i decided i wanted to have the choice of gravity
or electrostatics as the medium of interaction. I'm going to introduce charge into
the mix, should make it fun.

Sunday, 09/04
I spent all the time up to now working on the original program to get the plotting
process to work and finally did, so this one doesnt have that, gonna have to add
that in before i get started with changing things up from gravity to electromagnetism.



@author: quantumpsk
"""

import math
import numpy as np
np.set_printoptions(precision=2)
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os


class Object:
    name = " "
    mass = 0.0
    radius = 1.0
    charge = 1.0
    position = np.zeros(shape=(1, 3), dtype=float)
    velocity = np.zeros(shape=(1, 3), dtype=float)
    momentum = np.zeros(shape=(1, 3), dtype=float)
    netforce = np.zeros(shape=(1, 3), dtype=float)
    netacc = np.zeros(shape=(1, 3), dtype=float)
    Xposlog = []
    Yposlog = []
    Zposlog = []

    # initialize samples with pre or user defined values
    def __init__(self, name="Object", mass=10.0,radius=1.0, charge=1.0,
                 position=[0.0, 0.0, 0.0], velocity=[0.0, 0.0, 0.0]):
        self.name = name
        self.mass = mass
        self.radius = radius
        self.charge = charge * 1.602*10**-19 # Coulombs
        self.position = np.ndarray(shape=(1, 3), buffer=np.array(position))
        self.velocity = np.ndarray(shape=(1, 3), buffer=np.array(velocity))
        self.momentum = mass*self.velocity
        self.netforce = np.zeros(shape=(1, 3), dtype=float)
        self.netacc = np.zeros(shape=(1, 3), dtype=float)

    # accept initial values
    def initial(self, num):
        print("Object #", num)
        self.name = input("Enter the name: ")
        try:
            self.mass = float(input("Enter the mass: "))
            self.radius = float(input("Enter the radius "))
            self.charge = 1.602*10**-19 * float(input("Enter the charge: "))
            self.position = np.ndarray(shape=(1, 3), buffer=np.array([float(val)
                            for val in input("Enter the position: (x y z):").split()]))
            self.velocity = np.ndarray(shape=(1, 3), buffer=np.array([float(val)
                            for val in input("Enter the velocity: (Vx Vy Vz):").split()]))
        except:
            print("Don't be silly !")

    # display details
    def putObj(self):
        print("\n%s, with a mass of %d kgs, radius of %d ms, at\n" % (self.name, self.mass, self.radius, self.charge),
              self.position, " ms, \nmoves at ", self.velocity, " m/s\n",
              "with a net force of ", self.netforce, " N\n and net acc of ", self.netacc)


class System:
    G = 6.67*10**-6         # should be actually e-11, but what the hell
    Eps0 = 8.854*10**-12    # epsilon naught
    e = 1.602*10**-19       # charge of one electron/proton
    ke = (4*math.pi*Eps0)**-1   # coulombs constant
    size = 0
    name = ""
    Sysdist = np.zeros(shape=(1, 1), dtype=float)
    Sys_gforces = np.zeros(shape=(1, 1), dtype=float)
    Sys_eforces = np.zeros(shape=(1, 1), dtype=float)
    System = []
    dt = 0.01       # time step
    Tt = 100        # total time
    tt = 0          # current time
    Poslog = [[]]   # position log
    Xlog = np.zeros(shape=(1, 1), dtype=float)
    Ylog = np.zeros(shape=(1, 1), dtype=float)
    Zlog = np.zeros(shape=(1, 1), dtype=float)
    path = ""
    # start with the initialization
    def __init__(self):
        self.name = "System1"
        try:
            start = 1 # **for ease of trial runs, disable later**
            #int(input("Do you want a (1)Sample System or (2)User-defined System?\n"))
            if start == 1:
                self.samplesys()
            elif start == 2:
                self.usersys()
            else:
                print("Dont be an ass!")
        except:
            print("Choose one of the above only, oh and Dont be an ass!")

    # user defined system
    def usersys(self):
        self.name = input("Enter the name of the System: ")
        self.size = int(input("Enter number of Objects: "))
        self.Sysdist = np.zeros(shape=(self.size, self.size), dtype=float)
        self.Sys_gforces = np.zeros(shape=(self.size, self.size, 1, 3), dtype=float)
        self.Sys_eforces = np.zeros(shape=(self.size, self.size, 1, 3), dtype=float)
        for i in range(self.size):
            self.System.append(Object())
            self.System[i].initial(i)
        self.Tt = int(input("Enter the duration of the sim: "))
        self.dt = float(input("Enter the time step value: "))
        self.totalsteps = self.Tt/self.dt
        self.tt = 0.0

    # set up sample system
    def samplesys(self):
        self.name = "SamSys"
        self.size = 5
        self.Sysdist = np.zeros(shape=(self.size, self.size), dtype=float)
        self.Sys_gforces = np.zeros(shape=(self.size, self.size, 1, 3), dtype=float)
        self.System.append(Object(name="Prospero", mass=1.0, radius=1.0, charge=1.0, position=[10.0, 10.0, 0.0], velocity=[0.0, 0.0, 0.0]))
        self.System.append(Object(name="Caliban", mass=1.0, radius=1.0, charge=1.0, position=[10.0, 0.0, 0.0], velocity=[0.0, 0.0, 0.0]))
        self.System.append(Object(name="Ariel", mass=1.0, radius=1.0, charge=-1.0, position=[0.0, 0.0, 0.0], velocity=[0.0, 0.0, 0.0]))
        self.System.append(Object(name="Miranda", mass=1.0, radius=1.0, charge=1.0, position=[0.0, 10.0, 0.0], velocity=[0.0, 0.0, 0.0]))
        self.System.append(Object(name="Sycorax", mass=1.0, radius=1.0, charge=1.0, position=[5.0, 5.0, 0.0], velocity=[0.0, 0.0, 0.0]))
        self.Tt = 10000
        self.dt = 1
        self.tt = 0.0
        self.totalsteps = self.Tt/self.dt
        for i in range(self.size):
            self.System[i].Xposlog.append(self.System[i].position[0, 0])
            self.System[i].Yposlog.append(self.System[i].position[0, 1])
            self.System[i].Zposlog.append(self.System[i].position[0, 2])
        self.path = os.getcwd()
        self.path = os.path.join(self.path, self.name)
        self.Box = np.array([[-500.0, 500.0], [-500.0, 500.0], [-500.0, 500.0]])

    # set up initial conditions
    def startstate(self):
        for i in range(self.size):
            for j in range(self.size):
                if j == i:
                    continue            # prevent ZeroError
                self.Sysdist[i][j] = np.linalg.norm(self.System[j].position - self.System[i].position)    # calc D(i,j)
                self.Sys_gforces[i, j, 0, :] = (self.G * self.System[i].mass * self.System[j].mass /      # calc Grav F(i,j)
                                                self.Sysdist[i][j]**2)
                self.Sys_eforces[i, j, 0, :] = (self.ke * self.System[i].charge *
                                        self.System[j].charge) / self.Sysdist[i][j]**2
        for i in range(self.size):
            temp1 = np.zeros(shape=(1, 3), dtype=float)
            temp2 = np.zeros(shape=(1, 3), dtype=float)
            for j in range(i, self.size):
                temp1 += self.Sys_gforces[i, j, 0, :]
                temp2 += self.Sys_eforces[i, j, 0, :]
            self.System[i].netforce = temp1 + temp2              # calc Fnet(i)
            self.System[i].netacc = self.System[i].netforce/self.System[i].mass     # calc Anet(i)
        print("sys state started")

    # single time step evolution
    def timestep(self):
        # calc new velocity  and position of objects
        for i in range(self.size):
            self.System[i].velocity += (self.System[i].netacc * self.dt)
            self.System[i].position += (self.System[i].velocity * self.dt) + (0.5 * self.System[i].netacc * self.dt**2)
        self.poslogupdate()

    # calculate current state values
    def calcstate(self):
        #calc new distance and force of system
        for i in range(self.size):
            for j in range(self.size):
                if j == i:
                    continue            # prevent ZeroError
                self.Sysdist[i][j] = np.linalg.norm(self.System[j].position-self.System[i].position)    # calc D(i,j)
                self.Sys_gforces[i, j, 0, :] = (self.G * self.System[i].mass * self.System[j].mass /      # calc F(i,j)
                               self.Sysdist[i][j]**3) * (self.System[j].position - self.System[i].position)
                self.Sys_eforces[i, j, 0, :] = (self.ke * (self.System[i].charge *self.e                     # calc Elec F(i,j)
                                                * self.System[j].charge * self.e)) / (self.Sysdist[i][j]**2)
        # calc new net force and acc for each object
        for i in range(self.size):
            temp1 = np.zeros(shape=(1, 3), dtype=float)
            temp2 = np.zeros(shape=(1, 3), dtype=float)
            for j in range(i, self.size):
                temp1 += self.Sys_gforces[i, j, 0, :]
                temp2 += self.Sys_gforces[i, j, 0, :]
            self.System[i].netforce = temp1 + temp2              # calc Fnet(i)
            self.System[i].netacc = self.System[i].netforce / self.System[i].mass     # calc Anet(i)

    # update the position log
    def poslogupdate(self):
        # update object position log
        for i in range(self.size):
            self.System[i].Xposlog.append(self.System[i].position[0, 0])
            self.System[i].Yposlog.append(self.System[i].position[0, 1])
            self.System[i].Zposlog.append(self.System[i].position[0, 2])

#        # write the whole system into a JSON file for later potential use
#        sysfile = open(filepath, mode='a')
#        # write the position log into file for later plotting

    # detect zero distance errors
    def testZDE(self):
        for i in range(self.size):
            for j in range(self.size):
                if i == j:
                    pass
                mindist = self.System[i].radius + self.System[j].radius
                if mindist >= self.Sysdist[i][j]:
                    okcoll = self.collision(i, j)

    # detect possible boundary event
    def testbound(self):
        for i in range(self.size):      # for each obj
            for j in range(3):          # for each axis
                if self.System[i].position[0, j] < self.Box[j, 0] or self.System[i].position[0, j] > self.Box[j, 1]:
                    self.System[i].velocity[0, j] *= -1

    # collision alert! avoid zero distance error (ZDE)!
    def collision(self, num1, num2):
        m1 = self.System[num1].mass
        m2 = self.System[num2].mass
        v1 = self.System[num1].velocity
        v2 = self.System[num2].velocity
        a = m1 * (m1 + m2)
        b = -2 * m1 * (m1 * v1 + m2 * v1)
        c = (m1 * v1**2 * (m1 - m2)) + (2 * m1 * m2 * v1 * v2)
        D = b**2 - (4 * a * c)
        if D.any() < 0:
            print("Zero Error!")
            return False
        root1 = (-b + math.sqrt(D.all())) / 2 * a
        root2 = (-b - math.sqrt(D.all())) / 2 * a
        if root1.all() == v1.all():
            v1_f = root2
        elif root2.all() == v1.all():
            v1_f = root1
        v2_f = (m1 * v1 + m2 * v2 - m1 * v1_f) / m2
        self.System[num1].velocity = v1_f
        self.System[num2].velocity = v2_f
        return True

    # evolve the system
    def sysevolve(self):
        # start at 0 to Tt
        while self.tt < self.Tt:
            self.timestep()
            self.calcstate()
 #           self.testZDE()
            self.testbound()
            self.tt += self.dt

    # plot the system simulation
    def sysplot(self):
        # start with the plot
        mpl.rcParams['legend.fontsize'] = 10
        fig = plt.figure()
        systplot = fig.gca(projection='3d')
        # obtain the orbital coords of each object
        self.Xlog = np.ndarray(shape=(len(self.System[0].Xposlog), 1), buffer=np.array(self.System[0].Xposlog))
        self.Ylog = np.ndarray(shape=(len(self.System[0].Yposlog), 1), buffer=np.array(self.System[0].Yposlog))
        self.Zlog = np.ndarray(shape=(len(self.System[0].Zposlog), 1), buffer=np.array(self.System[0].Zposlog))
        # reshape the logs acc to object index
        try:
            shape = (int(len(self.Xlog)/self.size), self.size)
            self.Xlog = self.Xlog.reshape(shape)
            self.Ylog = self.Ylog.reshape(shape)
            self.Zlog = self.Zlog.reshape(shape)
        except:
            print("\n...\n", shape, "..", len(self.Xlog), self.size)
        for i in range(self.size):
            systplot.plot(self.Xlog[:, i], self.Ylog[:, i], self.Zlog[:, i], label=self.System[i].name)
        systplot.legend()
        plt.show()

    # show off your system
    def display(self):
        for Obj in self.System:
            Obj.putObj()


# Main Code
Sys = System()
Sys.startstate()
Sys.sysevolve()
Sys.display()
Sys.sysplot()
