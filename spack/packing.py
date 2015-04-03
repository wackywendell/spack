# -*- coding: UTF-8 -*-

import sys, os.path
import numpy as np
from math import sqrt

def rand_sphere(d0):
    """
    Get random points within a sphere of radius 1. Returns array of shape (d0, 3).
    """
    p1 = np.random.randn(d0, 3)
    m = np.sqrt(np.sum(p1**2, axis=1))
    
    rad = pow(np.random.rand(d0), 1.0/3.0)
    return (p1.T * (rad/m)).T

def rand_disk(d0):
    """
    Get random points within a disk of radius 1. Returns array of shape (d0, 2).
    """
    phi = np.random.rand(d0)*(2*np.pi)
    rad = pow(np.random.rand(d0), 1.0/2.0)
    return (np.array([cos(phi), sin(phi)]) * rad).T

def Vec2_diff(r1, r2, shear=0.0, Lx = 1.0, Ly = 1.0):
    """
    Difference between two points in a 2D box with shear
    """
    dx,dy = (r1 - r2).T
    im = np.round(dy/Ly)
    dy = dy - im*Ly
    dx = dx-np.round(dx/Lx-im*shear)*Lx-im*shear*Lx
    return np.array((dx,dy)).T

class Packing:
    """
    A class representing a packing of spheres in a periodic box.
    """
    def __init__(self, rs, diameters, shear = 0.0, L=1.0):
        self.rs = np.array(rs) / float(L)
        self.diameters = np.array(diameters) / float(L)
        self.shear = shear
        self.L = L
        
        self.N = len(self.diameters)
        n, self.ndim = np.shape(self.rs)
        if n != self.N:
            raise ValueError("Need shape N for diameters, Nx2 or Nx3 for rs; got {} and {}x{}".format(
                self.N, n, self.ndim))
        if self.ndim == 3:
            if self.shear != 0:
                raise NotImplementedError("Can't do shearing for 3D")
        elif self.ndim != 2:
            raise ValueError("Number of dimensions must be 2 or 3; got {}".format(self.ndim))
    
    def neighbors(self, tol=1e-8):
        """
        For a set of particles at xs,ys with diameters diameters, finds the
        distance vector matrix (d x N x N) and the adjacency matrix.
        
        Assumes box size 1, returns (adjacency matrix, diffs)
        """
        diffs = np.remainder(np.array([np.subtract.outer(xs, xs) for xs in self.rs.T])+.5, 1) -.5
        
        if self.shear != 0:
            xdiff, ydiff = diffs[:2]
            im = np.round(ydiff)
            xdiff -= im*self.shear
            ydiff = ydiff - im
            xdiff -= np.round(xdiff)
            diffs[:2] = xdiff, ydiff
        
        sigmadists = np.add.outer(self.diameters, self.diameters)/2.
        dists = np.sqrt(np.sum(diffs**2, axis=0))
        
        return dists - sigmadists < tol, diffs*self.L
    
    def backbone(self, tol=1e-8):
        """Returns (backbone indices, neighbor matrix)"""
        areneighbors, _ = self.neighbors(tol)
        notfloaters = np.sum(areneighbors, axis=0) >= self.ndim + 2 # self.ndim + 1 for stability, +1 for itself
        
        oldNpack = -1
        Npack = np.sum(notfloaters)
        while Npack != oldNpack:
            areneighbors[~notfloaters] = 0
            areneighbors[:, ~notfloaters] = 0
            notfloaters = np.sum(areneighbors, axis=0) >= self.ndim + 2
            oldNpack, Npack = Npack, np.sum(notfloaters)
        
        return notfloaters, areneighbors
    
    def contacts(self, tol=1e-8):
        """Returns (number of backbone contacts, stable number, number of floaters)"""
        idx, nbor = self.backbone(tol=tol)
        return np.sum(np.triu(nbor, 1)), (np.sum(idx)-1) * self.ndim + 1, np.sum(~idx)
    
    def size_indices(self, tol=1e-8):
        """Returns [idx of sigma1, idx of sigma2, ...]"""
        sigs = np.array(np.round(self.diameters/tol), dtype=int)
        sigset = set(sigs)
        return [sigs == s for s in sorted(sigset)]
    
    def dist_tree(self, other, tol=1e-8):
        """
        Find the distance between two packings.
        
        Requires pyparm.
        """
        if self.ndim == 2:
            from pyparm import d2 as sim
        else:
            from pyparm import d3 as sim
            
        assert self.ndim == other.ndim
        
        if self.shear != 0 or other.shear != 0:
            assert np.abs(self.shear - other.shear) <= tol
            shear = (self.shear + other.shear) / 2.0
            box = sim.LeesEdwardsBox(sim.Vec(1,1), shear)
        else:
            box = sim.OriginBox(1.0)
        
        sz1 = self.size_indices()
        assert(len(sz1) == 2)
        cutoff1 = int(np.sum(sz1[0]))
        
        sz2 = other.size_indices()
        assert(len(sz2) == 2)
        cutoff2 = int(np.sum(sz2[0]))
        
        vs1 = [sim.Vec(*xy) for idx in sz1 for xy in self.rs[idx]]
        vs2 = [sim.Vec(*xy) for idx in sz2 for xy in other.rs[idx]]
        
        tree = sim.jammingtreeBD(box, sim.vecvector(vs1), sim.vecvector(vs2), cutoff1, cutoff2)
        return tree
    
    def dist(self, other, tol=1e-8, maxt = 1000000):
        tree = self.dist_tree(other, tol=tol)
        tree.expand(maxt)
        return sqrt(tree.curbest().distsq)*self.L
        
    @staticmethod
    def _cage_pts(xyz, neighbor_xyzs, sigma, neighbor_diameters, L, M, R):
        """Finds points within a distance R of point xyz that do not conflict with neigbors"""
        pts = rand_sphere(M)*R + xyz
        for nxyz, nsig in zip(neighbor_xyzs, neighbor_diameters):
            dpts = np.remainder(pts - nxyz + L/2.0, L) - L/2.0
            dists_sq = np.sum(dpts**2, axis=1)
            goodix = dists_sq >= ((nsig + sigma)/2.0)**2
            pts = pts[goodix, :]
        return pts

    def cages(self, M=10000, R=None, Rfactor = 1.2, padding=0.1, Mfactor=0.1):
        """
        Find all cages in the current "packing".
        
        The algorithm uses Monte Carlo: it finds M random points within a sphere of radius R from
        each particle, and sees if that particle could sit there without conflicting with other particles.
        Then (number of accepted points) / (number of test points) * (volume of sphere) is the
        volume of the cage.
        
        The algorithm is adaptive: if not enough test points are accepted (n < M * Mfactor), it tries
        more test points. If any test points are within `padding` of the edge, `R` is (temporarily)
        expanded.
        
        Parameters
        ----------
        M : Number of points in the sphere to test
        R : Size of sphere to test (will be expanded if necessary)
        Rfactor : How much to increase R by when the cage doesn't fit
        padding : How much larger the sphere should be than the cage (if it isn't, the sphere is
                    expanded)
        Mfactor : Mfactor * M is the minimum number of points to find per cage. If they aren't
                    found, more points are tested.
        
        Returns
        -------
        points : a list of (A x 3) lists, A indeterminate (but larger than M * Mfactor), with each
                    list corresponding to the points within one cage.
        Vs : The approximate volumes of each cage.
        """
        if R is None: R = min(self.diameters) * 0.2
        neighbordict = {}
        
        psets = []
        Vs = []
        for n, (xyz, s) in enumerate(zip(self.rs, self.diameters)):
            curR = R
            curpow = -1
            nxyzs = nsigs = None
            
            def get_pts():
                pts = self._cage_pts(xyz, nxyzs, s, nsigs, 1.0, M, curR)
                maxdist = np.max(np.sqrt(np.sum((pts-xyz)**2, axis=1))) if len(pts) > 0 else 0
                return pts, maxdist
            
            pts, maxdist = [], curR
            while maxdist * (1. + padding) > curR:
                # print(n, curpow, maxdist * (1. + padding), curR)
                curpow += 1
                curR = R * pow(Rfactor, curpow)
                curM = M
                if curpow not in neighbordict:
                    pack = Packing(self.rs, self.diameters + curR)
                    cur_neighbors, _ = pack.neighbors(tol=0)
                    cur_neighbors[np.diag_indices_from(cur_neighbors)] = False
                    neighbordict[curpow] = cur_neighbors
                cur_neighbors = neighbordict[curpow]
                nix = cur_neighbors[n]
                nxyzs = self.rs[nix, :]
                nsigs = self.diameters[nix]
                pts, maxdist = get_pts()
                if maxdist * (1. + padding) > curR: continue
            
                while len(pts) < Mfactor * M:
                    # print(curM, len(pts), Mfactor * M, len(pts) < Mfactor * M)
                    pts2, maxdist2 = get_pts()
                    maxdist = max((maxdist, maxdist2))
                    pts = np.concatenate((pts, pts2))
                    curM += M
                    if maxdist * (1. + padding) > curR: break
                if maxdist > 0.5:
                    raise ValueError('Cage size filling entire space, cannot continue.')
            
            fracgood = len(pts) / curM
            r = fracgood**(1./3.) * curR
            V = fracgood * 4 * np.pi / 3 * (curR**3)
            psets.append(pts)
            Vs.append(V)
        return psets**self.L, np.array(Vs)*(self.L**self.ndim)
    
    def scene(pack, cmap=None, rot=0, camera_height = 0.7,  camera_dist = 1.5, angle=None,
            lightstrength = 1.1, orthographic=False, pad=None, floatercolor=(.6,.6,.6),
            bgcolor = [1,1,1]):
        """
        Render a 3D scene.
        
        Requires `vapory` package, which requires the `povray` binary.
        
        Parameters
        ----------
        cmap : a colormap
        
        Returns
        -------
        scene : vapory.Scene, which can be rendered using its `.render()` method.
        """
        import vapory
        import numpy as np
        
        try:
            import matplotlib as mpl
            import matplotlib.cm as mcm
            vmin, vmax = min(pack.diameters), max(pack.diameters)
            sm = mcm.ScalarMappable(norm=mpl.colors.Normalize(vmin, vmax), cmap=cmap)
            cols = [sm.to_rgba(s) for s in pack.diameters]
        except ImportError:
            if not isinstance(cmap, list):
                raise ValueError("matplotlib could not be imported, and cmap not recognizeable as a list")
            cols = list(cmap)
        except TypeError:
            if not isinstance(cmap, list):
                raise ValueError("matplotlib could not convert cmap to a colormap, and cmap not recognizeable as a list")
            cols = list(cmap)

        if floatercolor is not None:
            ix, _ = pack.backbone()
            ns, = np.nonzero(~ix)
            for n in ns: cols[n] = floatercolor
        rs = np.remainder(pack.rs+.5, 1)-.5
        spheres = [
            vapory.Sphere(xyz, s/2., vapory.Texture( vapory.Pigment( 'color', col[:3] )))
            for xyz, s, col in zip(rs, pack.diameters, cols)
        ]

        extent = (-.5, .5)
        corners = [np.array((x,y,z)) for x in extent for y in extent for z in extent]
        pairs = [(c1, c2) for c1 in corners for c2 in corners if np.allclose(np.sum((c1-c2)**2), 1) and sum(c1-c2) > 0]

        radius = 0.01
        col = vapory.Texture( vapory.Pigment( 'color', [.5,.5,.5]))
        cyls = [vapory.Cylinder(c1, c2, 0.01, col) for c1,c2 in pairs]
        caps = [vapory.Sphere(c, radius, col) for c in corners]
        
        
        light_locs = [
            [ 8.,  5., -3.],
            [-6.,  6., -5.],
            [-6., -7., -4.],
            [ 10., -5., 7.]
        ]
        rotlocs = [[x*np.cos(rot) - z*np.sin(rot),y, z*np.cos(rot) + x*np.sin(rot)] for x,y,z in light_locs]
        lights = [
            #vapory.LightSource( [2,3,5], 'color', [1,1,1] ),
            vapory.LightSource(loc, 'color', [lightstrength]*3 ) for loc in rotlocs
        ]
        cloc = [np.cos(rot)*camera_dist,camera_dist*camera_height,np.sin(rot)*camera_dist]
        mag = sqrt(sum([d**2 for d in cloc]))
        direction = [-v *2/ mag for v in cloc]
        
        if angle is None:
            if pad is None: pad = max(pack.diameters)
            w = sqrt(2) + pad
            angle = float(np.arctan2(w, 2*camera_dist))*2*180/np.pi
        camera = vapory.Camera('location', cloc, 'look_at', [0,0,0], 'angle', angle)
        # vapory.Camera('orthographic', 'location', cloc, 'direction', direction, 'up', [0,2,0], 'right', [2,0,0])
         
        return vapory.Scene( camera, objects= lights + spheres + cyls + caps + [vapory.Background( "color", bgcolor )])
    
    def plot_disks(self, ax=None, color=None, alpha=0.4, reshape=True):
        """
        Plot the packing as a set of disks.
        
        Color can be None (uses the standard sets), 'diameter' (colors by diameter), or a list of colors.
        
        'reshape' means set axis scaled, etc.
        """
        import matplotlib as mpl
        from itertools import cycle
        import numpy as np
        
        if color == 'diameter':
            dset = sorted(set(self.diameters))
            cold = dict(zip(dset, cycle(mpl.rcParams['axes.color_cycle'])))
            color = [cold[d] for d in self.diameters]
        if not np.iterable(color):
            color = cycle((color,))
            
        if ax is None: ax = mpl.pyplot.gca()
        
        rs = np.remainder(self.rs+.5, 1)-.5
        L = self.L
        dloc = (0, 1)
        for (x0,y0),d,c in zip(self.rs, self.diameters, color):
            for x,y in [np.array((x0+dx,y0+dy)) for dx in dloc for dy in dloc]:
                if (x + d > 0 and x - d < 1 and y + d > 0 and y - d < 1):
                    circ = mpl.patches.Circle((x*L,y*L), d*L/2, axes=ax, ec='none', fc=c, alpha=alpha)
                    ax.add_patch(circ)
        
        if reshape:
            ax.axis([0,L,0,L])
            ax.set_aspect('equal')
            ax.set_xticks([],[])
            ax.set_yticks([],[])
    
    def plot_contacts(self, ax=None, tol=0, reshape=True, **kw):
        """Designed for use with plot_disks, this will plot a line between neighboring particles."""
        import matplotlib as mpl
        import numpy as np
        
        if ax is None: ax = mpl.pyplot.gca()
        
        kw.setdefault('color', 'k')
        kw.setdefault('alpha', 0.4)
        
        adj, diffs = self.neighbors(tol=tol)
        
        for i, (adjrow, drow) in enumerate(zip(adj, diffs.T)):
            for j, (a, d) in enumerate(zip(adjrow, drow)):
                if not a: continue
                
                ri = np.remainder(self.rs[i], 1)
                x,y = ri
                rid = ri + d
                
                if i < j and max(abs(rid*2 - 1)) < 1:
                    # if we're going to get to (j,i) later
                    # and both points are in the box, skip it this time
                    continue
                x2,y2 = rid
                
                x,y = ri*self.L
                x2,y2 = rid*self.L
                ax.plot([x,x2], [y,y2], **kw)
        
        if reshape:
            ax.axis([0,self.L,0,self.L])
            ax.set_aspect('equal')
            ax.set_xticks([],[])
            ax.set_yticks([],[])
                
    def DM(self, masses=None):
        """Dynamical matrix for array rs, size ds. Assumes epsilon is the
        same for all.
        
        Parameters
        ----------
        masses : an array of length N of the masses of the particles.
        """
        N = len(self.diameters)
        rs = self.rs
        d = self.ndim
        M=np.zeros((d*N,d*N))
        
        for i in range(N):
            sigi = self.diameters[i]
            for j in range(i):
                rijvec=rs[i,:]-rs[j,:]
                rijvec=rijvec-np.around(rijvec)
                rijsq=np.sum(rijvec**2)
                dij=(sigi+self.diameters[j])/2
                dijsq=dij**2
                if rijsq < dijsq:
                    rij=np.sqrt(rijsq)
                    rijouter = np.outer(rijvec,rijvec)
                    # U(r) = ½(1 - r/d)²
                    # d²U/dxdy = (dr/dx)(dr/dy)/d² - (1 - r/d)(d²r/dxdy)/d
                    # dr/dx = x/r
                    # d²r/dxdy = -(x y) / r³
                    # d²U/dxdy = -(x y)/(r² d²) + (1 - r/d)((x y)/r²)/(d r)
                    # d²U/dx² = (dr/dx)²/d² - (1 - r/d)(d²r/dx²)/d
                    # d²r/dx² = -x² / r³ + 1/r
                    # d²U/dxᵢdxⱼ = -(xᵢ xⱼ)/(r² d²) + (1 - r/d)((xᵢ xⱼ)/r² - δᵢⱼ)/(d r)
                    
                    Mij1=-rijouter/rijsq/dijsq
                    Mij2=(1-rij/dij)*(rijouter/rijsq-np.eye(d))/rij/dij
                    Mij=Mij1+Mij2
                    
                    M[d*i:d*i+d,d*j:d*j+d]=Mij
                    M[d*j:d*j+d,d*i:d*i+d]=Mij
                    M[d*i:d*i+d,d*i:d*i+d]-=Mij
                    M[d*j:d*j+d,d*j:d*j+d]-=Mij
        
        np.divide(M, self.L**2, out=M)
        if masses is None:
            return M
        
        # TODO: is the mass part of this really part of this?
        marr = np.array(masses)
        assert np.shape(masses) == np.shape(self.diameters)
        marr = np.array([masses]*d)
        marr = marr.T.flatten()
        # marr is now [m1,m1,m2,m2,...] (in 2D)
        mm = np.eye(d*N)
        np.multiply(mm, marr**-.5, out=mm)
        # mm is now M^-½, where M is the mass matrix
        
        mm.dot(M, out=M)
        M.dot(mm, out=M)
        return M
    
    def DM_freqs(self, masses=None):
        """Find the frequencies corresponding to the eigenvalues of the dynamical matrix.
        
        This is just a short wrapper around DM()."""
        ew, ev = np.linalg.eig(self.DM(masses=masses))
        # this used to be over 2pi; I don't know where the 2 went, but it seems to be gone now...
        return np.sqrt(np.abs(ew)) / (np.pi)
    
    def forces(self):
        """Find Fij on each particle, assuming a harmonic potential, U = 1/2 (1 - r/σ)^2
        
        Returns a dxNxN matrix."""
        adj, diffs = self.neighbors(tol=0)
        
        dists = np.sqrt(np.sum(diffs**2, axis=0))
        sigij = np.add.outer(self.diameters, self.diameters)/2.

        adj2 = np.triu(adj, k=1)
        overlaps = (sigij - dists)[adj2]
        
        dr = (1 - dists / sigij)
        with np.errstate(divide='ignore'):
            Fmags = dr / sigij / dists
        Fmags[~adj] = 0
        Fmags[np.diag_indices_from(Fmags)] = 0
        return diffs * Fmags / self.L
    
    def tess(self):
        """Get a `tess.Container` instance of this.
        
        Requires `tess`.
        """
        import tess
        return tess.Container(self.rs*self.L % self.L, limits=self.L, radii=self.sigmas/2., periodic=True)
