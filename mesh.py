import numpy as np

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
            domain_azi = np.loadtxt(directory+"used_azi.dat")  # azimuthal interfaces of grid cells
            self.pedge = np.append(domain_azi[:,1],domain_azi[-1:,2][0])
        else:
            self.pedge = np.linspace(0.,2.*np.pi,self.nsec+1)  # phi-edge         
        self.pmed = 0.5*(self.pedge[:-1] + self.pedge[1:]) # phi-center

        # -----
        # colatitude / latitude
        # -----
        if self.ncol > 1:
            try:
                domain_col = np.loadtxt(directory+"domain_z.dat")  # radial interfaces of grid cells
            except IOError:
                print('IOError')
            self.tedge = 0.5*np.pi - domain_col[3:-3]           # lat-edge
            self.tedge = self.tedge[::-1]
            self.tmed  = 0.5*(self.tedge[:-1] + self.tedge[1:]) # lat-center

        # -----
        # cylindrical grid for 3D FARGO3D runs using spherical grid
        # -----        
        if self.fargo3d == 'Yes' and self.ncol > 1:
            # define number of cells in vertical direction
            self.nver = self.ncol 
            # define an array for vertical altitude above midplane
            # (half-disc assumed)
            zbuf = self.redge.max()*np.sin(self.tedge)
            self.zedge = np.linspace(0.0,zbuf.max(),self.nver+1)
            self.zmed  = 0.5*(self.zedge[:-1] + self.zedge[1:]) # z-center

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
