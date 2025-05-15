import numpy as np
import os

# ---------------------
# building mesh arrays 
# ---------------------

class Mesh():
    def __init__(self, directory=""):
        if len(directory) > 1:
            if directory[-1] != '/':
                directory += '/'                
        # -----
        # radius
        # -----
        try:
            domain_rad = np.loadtxt(directory+"used_rad.dat")  # radial interfaces of grid cells
        except IOError:
            print('IOError')
        self.redge = domain_rad                              # r-edge
        self.rmed = 2.0*(domain_rad[1:]*domain_rad[1:]*domain_rad[1:] - domain_rad[:-1]*domain_rad[:-1]*domain_rad[:-1]) / (domain_rad[1:]*domain_rad[1:] - domain_rad[:-1]*domain_rad[:-1]) / 3.0  # r-center 

        # -----
        # azimuth
        # -----
        if self.fargo3d == 'No':
            if self.fargo_orig == 'No':
                self.pedge = np.linspace(0.,2.*np.pi,self.nsec+1)
                # CUIDADIN!!
                '''
                domain_azi = np.loadtxt(directory+"used_azi.dat")  # azimuthal interfaces of grid cells
                if self.nsec > 1:
                    self.pedge = np.append(domain_azi[:,1],domain_azi[-1:,2][0])
                else:
                    self.pedge = np.linspace(0.,2.*np.pi,self.nsec+1)
                '''
            else:
                # case used_azi.dat does not exist, for instance with original FARGO code:
                self.pedge = np.linspace(0.,2.*np.pi,self.nsec+1)  # phi-edge   
            # cuidadin (May 2025)
            self.pmed  = np.linspace(0.,2.*np.pi,self.nsec+1)   
            self.pedge = self.pmed - 0.5*(self.pmed[1]-self.pmed[0])
        else:
            try:
                domain_azi = np.loadtxt(directory+"domain_x.dat")  # radial interfaces of grid cells
            except IOError:
                print('IOError')
            self.pedge = domain_azi
            self.pmed = 0.5*(self.pedge[:-1] + self.pedge[1:]) # phi-center

        # -----
        # colatitude / latitude
        # -----
        if self.cylindrical_grid == 'No' and self.nz > 1:
            try:
                domain_col = np.loadtxt(directory+"domain_z.dat")  # radial interfaces of grid cells
            except IOError:
                print('IOError')
            self.tedge = 0.5*np.pi - domain_col[3:-3]           # lat-edge
            self.tedge = self.tedge[::-1]
            self.tmed  = 0.5*(self.tedge[:-1] + self.tedge[1:]) # lat-center
            #
            # cylindrical grid for 3D FARGO3D runs using spherical grid
            #
            # define number of cells in vertical direction
            self.nver = self.nz
            # define an array for vertical altitude above midplane
            # (half-disc assumed)
            zbuf = self.redge.max()*np.sin(self.tedge)
            self.zedge = np.linspace(0.0,zbuf.max(),self.nver+1)
            self.zmed  = 0.5*(self.zedge[:-1] + self.zedge[1:]) # z-center

        # -----
        # altitude in 3D cylindrical runs
        # -----        
        if self.cylindrical_grid == 'Yes' and self.nz > 1:
            try:
                domain_col = np.loadtxt(directory+"domain_z.dat")  # radial interfaces of grid cells
            except IOError:
                print('IOError')
            self.zedge = domain_col[3:-3]           # altitude-edge
            self.zmed  = 0.5*(self.zedge[:-1] + self.zedge[1:]) # altitude-center

        # -----
        # cartesian grid for FARGO3D runs using cartesian coordinates
        # -----        
        if self.fargo3d == 'Yes' and self.cartesian_grid == 'Yes':
            try:
                domain_x = np.loadtxt(directory+"domain_x.dat")
            except IOError:
                print('IOError for domain_x.dat')
            self.xedge = domain_x[3:-3] 
            self.nx = int(len(self.xedge)-1)
            self.xmed  = 0.5*(self.xedge[:-1] + self.xedge[1:]) # x-center
            try:
                domain_y = np.loadtxt(directory+"domain_y.dat")
            except IOError:
                print('IOError for domain_y.dat')
            self.yedge = domain_y[3:-3] 
            self.ny = int(len(self.yedge)-1)
            self.ymed  = 0.5*(self.yedge[:-1] + self.yedge[1:]) # x-center
