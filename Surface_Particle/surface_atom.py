import math
import tkinter as tk
import numpy as np


class Particle:
    """
    Parameters
    -----------
    element:['Au', 'Ag', 'Cu', 'Ni', 'Pt', 'Co', 'Rh', 'Pd, 'Ir']
            A str. Defines the element to be calculated.

    shape: ['sphere', 'rod', 'cube']
            A str. Defines the shape of the nanoparticle.

    **radius: A parameter for sphere and rod, entered in the unit of cm, eg 20nm = 2*E6 cm

    **length: A parameter for cube, entered in the unit of cm.

    Attributes
    ----------
    atomic_weight, atomic_radius, atomic_density, surf_atoms: a tuple(# of surface atoms, atom % in bulk)

    Examples
    -----------
    Pd = Particle('Pd','Sphere',radius = 1E-7)
    Pd.surf_atoms()
    """

    # Class Variable
    metal = {'Au': (196.97, 174, 19.3), 'Ag': (107.87, 165, 10.49), 'Cu': (63.55, 145, 8.98), 'Ni': (58.69, 149, 8.91), 'Pd': (
        106.42, 169, 12.02), 'Pt': (195.08, 177, 21.45), 'Co': (69.72, 136, 8.9), 'Rh': (102.91, 173, 12.45), 'Ir': (192.22, 180, 22.56)}

    # Instant Initiialization
    def __init__(self, element, shape, **kwargs):
        self.element = element
        self.shape = shape.lower()
        # Initiate arbitrary number of keyword arguments
        self.__dict__.update(kwargs)

    # Instant Method: Getting atomic weight
    def atomic_weight(self):
        self.a_weight = self.metal[self.element][0]
        return self.a_weight

    # Instant Method: Getting atomic weight
    def atomic_radius(self):
        self.a_radius = self.metal[self.element][1]
        # print('{} pm'.format(self.atomic_radius))
        return self.a_radius

        # Instant Method: Getting atomic density
    def atomic_density(self):
        self.a_density = self.metal[self.element][2]
        # print('{} g/cm^3'.format(self.atomic_density))
        return self.a_density

    # Instant Method: Getting total number of atoms in surface layer
    def surf_atoms(self):
        if self.shape == 'sphere':
            self.particle_volume = 4/3 * np.pi * (self.radius)**3

            self.new_radius = self.radius - self.atomic_radius()*2*1E-10
            self.new_particle_volume = 4/3 * np.pi * (self.new_radius)**3

        elif self.shape == 'rod':
            # Two half-spheres and a cylinder
            self.particle_volume = 4/3 * np.pi * \
                (self.radius)**3 + np.pi * (self.length) * (self.radius)**2

            self.new_radius = self.radius - self.atomic_radius()*2*1E-10
            self.new_particle_volume = 4/3 * np.pi * \
                (self.new_radius)**3 + np.pi * \
                (self.length) * (self.new_radius)**2

        elif self.shape == 'cube':
            self.particle_volume = self.length**3

            self.new_length = self.length - 2*self.atomic_radius()*2*1E-10
            self.new_particle_volume = self.new_length**3

        self.particle_mass = self.atomic_density() * self.particle_volume
        self.particle_num = self.particle_mass / self.atomic_weight() * 6.02E23

        self.new_particle_mass = self.atomic_density() * self.new_particle_volume
        self.new_particle_num = self.new_particle_mass / self.atomic_weight() * 6.02E23

        return (self.particle_num - self.new_particle_num, (self.particle_num - self.new_particle_num)/self.particle_num*100)


# Creating GUI
calc = tk.Tk()
# to rename the title of the window
calc.title("Surface Atom Calculator")
# pack is used to show the object in the window
tk.Label(calc, text="\n\nREADME: this calculator is programmed to calculate the number and percentage of surface atoms present in a metal nanoparticle. The chemical composition (['Au', 'Ag', 'Cu', 'Ni', 'Pt', 'Co', 'Rh', 'Pd, 'Ir']), and the shape (['sphere', 'rod', 'cube']) of the particle need to be specified. Sphere and rod need a radius input. Rod and cube need a length input. Set the value to 0 if not used!!!\n\n\n", wraplength=750, justify=tk.CENTER).grid(
    row=0, column=0, columnspan=5)

tk.Label(calc, text="Element: ").grid(row=1, column=0)
tk.Label(calc, text="Shape: ").grid(row=1, column=1)
tk.Label(calc, text="Radius: ").grid(row=1, column=2)
tk.Label(calc, text="cm").grid(row=2, column=3)
tk.Label(calc, text="Length: ").grid(row=1, column=4)
tk.Label(calc, text="cm").grid(row=2, column=5)


e1 = tk.Entry(calc)
e1.insert(0, 'Au')
e1.grid(row=2, column=0)

e2 = tk.Entry(calc)
e2.insert(1, 'sphere')
e2.grid(row=2, column=1)

e3 = tk.Entry(calc)
e3.insert(2, '1E-6')
e3.grid(row=2, column=2)


e4 = tk.Entry(calc)
e4.insert(2, '0')
e4.grid(row=2, column=4)


def button1():
    element = e1.get()
    shape = e2.get()
    radius = float(e3.get())
    length = float(e4.get())

    nano = Particle(element, shape, radius=radius, length=length)
    tk.Label(calc, fg='blue', text='\nThere are {} surface atoms, which account for {:.2f} percent of total atoms!\n'.format(
        math.floor(nano.surf_atoms()[0]), nano.surf_atoms()[1]), wraplength=750, justify=tk.CENTER).grid(row=6, column=0, columnspan=5)


button1 = tk.Button(calc, text='Calculate', padx=30, pady=10,
                    command=button1).grid(row=5, column=1, columnspan=2)

calc.mainloop()
