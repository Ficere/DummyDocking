import os
from math import radians, degrees, pi
import numpy as np
TWOPI = np.pi * 2.0

import rotation, nlopt
from rotation import normalize
from shoemake import cube3_to_quaternion, quaternion_to_cube3

class BFGSsearch:

    # infinitesimal rotation, i.e.: delta_radians
    # find most stable value, this might be too big/small near shoemake singularities TODO
    infinitesimal_radians = 1E-5 
        
    def __init__(self, adfr):
        self._neval = 0
        self._nbTors = 0
        self._adfr = adfr
        self._atoms = adfr.FT.atomSet
        self._elements = self.getElements(adfr.ligand._ag)
        self._rotCenter = adfr.FT.rotationCenter

        from ADFRcc.adfr import Parameters
        parameters = Parameters.getParameters()
        parameters.setHbondMode(0)# use 10-12 potential

        pop = adfr.createPopulation(1)
        ind = pop[0]
        var = ind.genomePy.getVariablesPy()
        self._nbTors = len(var)-7
        self._best = [99999999999999, ind.genomePy.getIdentityGenesPy()]
        
        self.setupGBFS()
        
    def getElements(self, atoms):
        elems = []
        for atype in atoms.getData("AD_element"):
            if atype in ['HD', 'OA', 'NA', 'SA']:
                elems.append(atype[0])
            elif atype == 'A':
                elems.append('C')
            else:
                elems.append(atype)
        return elems

    def write_xyz(fname, coords, mode='a'):
        #import pdb; pdb.set_trace()
        f = open(fname, mode)
        f.write('%d\n\n' % len(coords))
        for elem, (x, y, z) in zip(self._elements, coords):
            f.write('%5s %12.6f %12.6f %12.6f\n' % (elem, x, y, z)) 
        f.close()

    def forces_to_delta_genes(self, genes, coords, forces):
        # genes
        # 
        genes_gradient = np.zeros(len(genes))# + len(coords_obj.zmol.tors))  
        # # # # # # # # # #
        #                 #  
        #   translation   #  
        #                 #  
        # # # # # # # # # #
        genes_gradient[0:3] = np.sum(forces, 0) 
        # # # # # # # # # #
        #                 #  
        #    quaternion   #  
        #                 #  
        # # # # # # # # # #
        torque = np.zeros(3)
        max_dist_sq = 0.
        for xyz, f in zip(coords, forces):
            r = xyz - self._rotCenter
            torque += np.cross(r, f) 

        # calculate angle
        torque_length = np.sqrt(np.sum(torque**2))
        #print '------------------'
        #print 'torque: %10.5f %10.5f %10.5f' % tuple(torque)
        #print 'torque magnitude: %20.5f' % torque_length

        # quaternion that performs infinitesimal rotation around torque axis
        quaternion_torque = rotation.axisangle_to_q(torque, self.infinitesimal_radians)

        # infinitesimal rotation corresponds to a move in shoemake space
        # from "current_shoemake" to "target_shoemake"
        current_shoemake = genes[3:6]
        current_quaternion = cube3_to_quaternion(current_shoemake)
        target_quaternion = rotation.q_mult(quaternion_torque, current_quaternion)
        target_shoemake = quaternion_to_cube3(target_quaternion)

        # The infinitesimal rotation will produce an infinitesimal displacement
        # in shoemake space. This is to guarantee that the direction of
        # the displacement in shoemake space is not distorted.
        # The correct amount of displacement in shoemake space is obtained
        # by multiplying the infinitesimal displacement by shoemake_scaling:
        shoemake_scaling = torque_length / self.infinitesimal_radians

        # calculate displacement in shoemake space
        # u2 and u3 are treated as angles (i.e.:  1deg - 359deg = +2, not -358)
        genes_gradient[3:4] = shoemake_scaling * (target_shoemake[0] - current_shoemake[0]) # u1
        genes_gradient[4:5] = shoemake_scaling * ((target_shoemake[1] - current_shoemake[1] + np.pi) % TWOPI - np.pi) # u2
        genes_gradient[5:6] = shoemake_scaling * ((target_shoemake[2] - current_shoemake[2] + np.pi) % TWOPI - np.pi) # u3

        # corrections to derivatives
        if 0.0 < current_shoemake[0] < 1.0:
            genes_gradient[3:4] *= (1.0 / current_shoemake[0] + 1.0 / (1.0 - current_shoemake[0]))
        genes_gradient[4:5] *= 4.0 * (1.0 - current_shoemake[0])
        genes_gradient[5:6] *= 4.0 * current_shoemake[0]

        #print 'quaternion derivative: %17.15f %17.15f %17.15f %17.15f' % tuple(quaternion_torque)
        #print '------------------'

        # iterate over torsions
        for i, tors in enumerate(self._adfr.FT.allMotions[2:]):
            t1, t2, t3, t4 = tors.getDihedralAtomIndexesPy()

            # "a" and "b" are the two ends of the rotatable bond
            a = coords[t2]   # coords [x,y,z] of atom "a"
            b = coords[t3]   # coords [x,y,z] of atom "b"

            # axis of rotation, length 1.0
            unitvec = normalize(b - a)             # array [x, y, z]

            # collect atom indexes that rotate with this bond
            affected_by_tor = self._adfr.FT.atomsAffectedByTorsion[tors]

            # collect atom coordinates that rotate with this bond
            # N x 3 array
            points = [coords[j] for j in affected_by_tor]

            # collect forces acting on corresponding atoms
            # N x 3 array
            forces_ = [forces[j] for j in affected_by_tor]

            # number of atoms, integer
            N = len(points)

            # initialize torque
            torque = np.zeros(3)            # 3 scalars (x, y, z)

            # iterate over each atom
            for j in range(N):

                atom_coords = points[j]  # 3 scalars (x, y, z)
                atom_force = forces_[j]   # 3 scalars (dx, dy, dz)

                # calculate torque on point "a"
                # could be any other point (e.g. "b") along the axis of rotation
                torque += np.cross(atom_coords - a, atom_force)

            # project torque on axis of rotation
            torque_on_axis = np.dot(unitvec, torque)

            # assign value to gradient of the gene 
            genes_gradient[6+i] = torque_on_axis # 6+i because we are in shoemake space

        return genes_gradient

    def _eval(self, genes, grad=None):
        """ 
        Input:  genes
        Output: score
        """
        ind = self._ind
        self._neval += 1
        d,a,b,c = cube3_to_quaternion(genes[3:6])
        adfrgenes = ind.genomePy.getGenesForValuesPy(list(genes[:3])+[a,b,c,d]+[degrees(x) if x < pi else  -degrees(2*pi-x) for x in genes[6:]])
        ind.setGenes(adfrgenes)
        score = self._ind.score()
        #print score, self._best[0],
        if score > self._best[0]:
            self._best = [score, adfrgenes]
            #print 'saving'
        #else:
        #    print
        #print 'SC ##################', neval, -score
        #print list(genes[:3])+[a,b,c,d]+list(genes[6:])
        #print adfrgenes
        gridScorer = ind.genomePy.scorer.getLrrGridScorer()
        #coords = ind.genomePy.scorer._ligandAtomSet.getCoordsPy()
        #forces = ind.genomePy.scorer._ligandAtomSet.getGradientArray()
        coords = self._atoms.getCoordsPy()
        forces = self._atoms.getGradientArray()
        #write_xyz('traj.xyz', coords, adElems, mode='a')

        #forces = gradAT+gradE*charges+gradS*np.abs(charges)
        # grad[:] = forces_to_delta_genes(genes, coords, forces*charges)
        genes_gradient = self.forces_to_delta_genes(genes, coords, forces)
        #for i in range(len(grad)):
        #    grad[i] = genes_gradient[i]
        grad[:] = genes_gradient[:]
        return -score

    def setupGBFS(self):
        adfr = self._adfr
        cx, cy, cz = adfr.maps[0].getCenterPy()
        spacing = adfr.maps[0].getDistBetweenGridPoints()
        sx, sy, sz = adfr.maps[0].getNumGridPointsPy()
        mini = [ cx-((sx-1)/2)*spacing,
                 cy-((sy-1)/2)*spacing,
                 cz-((sz-1)/2)*spacing]
        maxi = [ cx+((sx-1)/2)*spacing,
                 cy+((sy-1)/2)*spacing,
                 cz+((sz-1)/2)*spacing]
        df = 6 + self._nbTors
        opt = nlopt.opt(nlopt.LD_LBFGS, df)
        lower_lim = mini + [0.,0.,0.] + [0.]*self._nbTors
        upper_lim = maxi + [1., 2. * np.pi, 2. * np.pi] + [TWOPI]*self._nbTors
        self._lowerBounds = lower_lim
        self._upperBounds = upper_lim
        opt.set_lower_bounds(lower_lim)
        opt.set_upper_bounds(upper_lim)
        opt.set_min_objective(self._eval)
        self._opt = opt
        
    def search(self, ind, steps=200):
        self._opt.set_maxeval(steps)
        self._opt.set_xtol_rel(0.01)
        self._ind = ind
        self._neval = 0
        self._best = [ind._score, ind.getGenesPy()]
        ind.score() # needed for ind.genomePy.getVariablesPy() to return correct variables
        var = ind.genomePy.getVariablesPy()
        a,b,c,d = var[3:7]
        u1 = quaternion_to_cube3([d,a,b,c])
        indopt = list(var[:3])+list(u1)+[radians(x) if x>0 else radians(360+x) for x in var[7:]]
        #import pdb; pdb.set_trace()
        try:
            xopt = self._opt.optimize(indopt)
            end = 'Steps'
        except RuntimeError:
            end = 'RuntimeError'
        return self._best

if __name__=='__main__':
    from ADFR import ADFRRigidReceptor
    from ADFRcc import adfrcc
    from MolKit2 import Read
    ligandfile = "7cpa_lig.pdbqt"
    mapsFolder = "7cpa_rec"

    ligand = Read(ligandfile)
    adfr = ADFRRigidReceptor(ligand, mapsFolder)
    pop = adfr.createPopulation(1)
    ind = pop[0]
    atset = adfr.FT.atomSet
    goodPoints = np.load(os.path.join(mapsFolder, 'translationPoints.npy'))
    adfr.FT.allMotions[0].setPreferredPoints(goodPoints)
    ind.randomize()
    bfgs = BFGSsearch(adfr)
    bfgs.search(ind)
    print bfgs._neval, bfgs._best
