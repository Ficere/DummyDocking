import openbabel as ob
import numpy, math
from numpy.linalg import norm
from mglutil.math.rotax import rotax
from mglutil.math.rigidFit import RigidfitBodyAligner

def cross(a,b):
    return numpy.array([a[1]*b[2]-a[2]*b[1],a[2]*b[0]-
                        a[0]*b[2],a[0]*b[1]-a[1]*b[0]],'g')

def unit(a):
    return a/norm(a)

def closestCanonical( t1, t2, t3, angles=None):
    if angles is None:
        from puckerCoords import HillReilly as angles        
    n = 0
    mini = 3*360.
    best = None
    for name,c1,c2,c3 in angles:
        diff = abs(c1-t1) + abs(c2-t2) + abs(c3-t3)
        #print diff, mini, best, c1, c2, c3, t1, t2, t3
        if diff < mini:
            mini = diff
            best = n
        n += 1
    return angles[best]

class RingConformers:
    """
    This class can be used to compute the data neeed to build an
    FTDiscreteConformation node for a rigid body containing a flexible 5 or 6
    ring (Only 6 rings tested so far).

    coords are the 3D coordinates of the atom that rigidBodyIndices indexes into
    rigidBodyIndices are the indices of the atoms in the rigid body (i.e. ring +
    rigid fragments attached to ring atoms.
    ringIndices is a list of 5 or 6 atom indices describing the ring. They are
    expected to be ordered to follow the ring cycle.
    attachedIndices is a list of 5 or 6 lists of indices pointing the atoms
    attached to a given ring atoms.
    subTreeInd is a list of integers ranging from 0 to len(ring)-1 indicating
    the ring atoms for which there will be children in the FT, e.g. [0,3]
    indicated that there is a rotatable bond in the subtrees of atoms attached
    to ring atom 0 and ring atom 3. 
    """
    def __init__(self, coords, rigidBodyIndices, ringIndices, attachedIndices, subTreeInd):
        assert len(ringIndices)==6 or len(ringIndices)==5
        assert len(attachedIndices)==6 or len(attachedIndices)==5
        assert len(ringIndices)==len(attachedIndices)
        self._ringLength = len(ringIndices)
        self._coords =  numpy.array(coords)
        self._rigidBodyIndices = rigidBodyIndices
        self._ringIndices = ringIndices
        self._attachedIndices = attachedIndices
        self._origAngles = self.HillReillyAngles(coords[ringIndices])
        self._matrices = [numpy.identity(4)]*6
        self._subTree = numpy.array([False]*6)
        self._subTree[subTreeInd] = True
        self._subTreeInd = subTreeInd
        self._align = RigidfitBodyAligner()
        
    def HillReillyAngles(self, coords):
        #  return Hill-Reilly angles of puckering from atomic coordinates
        # coordinates are assumes ordered according to the ring
        assert len(coords)==6 or len(coords)==5
        
        atoms = numpy.zeros((len(coords)+1,3),dtype='float64')
        atoms[:-1] = coords
        atoms[-1] = coords[0]
        lcoords = len(coords)
        hlcoords = lcoords/2
        r = numpy.zeros((lcoords,3),dtype='float64')
        a = numpy.zeros((hlcoords,3),dtype='float64')
        p = numpy.zeros((lcoords,3),dtype='float64')
        q = numpy.zeros((hlcoords,3),dtype='float64')

        for i in range(0,lcoords):
            r[i] = atoms[i+1] - atoms[i]
        for i in range(0,hlcoords):
            a[i] = atoms[2*(i+1)] - atoms[2*i]
        for i in range(1,lcoords):
            p[i] = cross(r[i-1],r[i])
        for i in range(0,hlcoords):
            q[i] = cross(a[i],p[2*i+1])
        n = cross(a[1],a[0])
        theta=[]
        for i in q:
            theta.append(90 - math.degrees(math.acos(numpy.dot(i,n) /
                                                     (norm(i)*norm(n)))))
        if lcoords==6:
            return theta[0], theta[1], theta[2]
        else:
            return theta[0], theta[1]
            
    def refFrame(self, p1, p2, p3):
        # return 4 points forming an orthogonal coordinate system with
        # p2 at the origin. 2 vectors in the (p1,p2,p3) plane and one orthogonal
        x1, y1, z1 = p1 
        x2, y2, z2 = p2
        x3, y3, z3 = p3
        v1 = unit([x1-x2, y1-y2, z1-z2])
        v2 = unit([x3-x2, y3-y2, z3-z2])
        v3 = cross(v1,v2)
        return (p2, p2+v1, p2+v2, p2+v3)

    def coordsForThetas(self, thetas):
        #self._matrices = [numpy.identity(4)]*6
        ## self._matrices = [None]*6
        thetaDiff = numpy.array(thetas)-self._origAngles
        #print 'orig', r._origAngles
        #print 'target', theta1, theta2, theta3
        #print 'diff', thetaDiff
        coords = self._coords.copy()

        # compute ref frames at 3 atoms held in plane
        refs = []
        for i in range(len(self._ringIndices)):
            #print 'REF', self._ringIndices[2*i-1], self._ringIndices[2*i], self._ringIndices[2*i+1]
            #print i-1, i, (i+1)%6
            refs.append(self.refFrame(self._coords[self._ringIndices[i-1]],
                                      self._coords[self._ringIndices[i]],
                                      self._coords[self._ringIndices[(i+1)%self._ringLength]]))
        for i in range(len(self._ringIndices)/2):
            # define atom indices
            # flap axis
            axis0 = self._ringIndices[2*i]
            axis1 = self._ringIndices[(2*(i+1))%self._ringLength]

            ## rotate flap
            # first flap tip
            indices = [self._ringIndices[2*i+1]]
            # add the atoms attached to the flap tip
            #indices.extend(self._attachedIndices[2*i+1])
            # build flap rotation matrix
            mat = rotax(self._coords[axis0], self._coords[axis1], math.radians(thetaDiff[i]))
            # transform flap tip and attached atoms 
            c = numpy.ones((len(indices),4), 'f')
            c[:,:3] = coords[indices]
            coords[indices] = numpy.dot(c, mat)[:, :3]
            #print 'AAA', 2*i+1, mat
        # transform atoms attached to atoms _ringIndices and build matrices affecting
        # attached atoms

        for i in range(len(self._ringIndices)):
            newref = self.refFrame(coords[self._ringIndices[i-1]],
                                   coords[self._ringIndices[i]],
                                   coords[self._ringIndices[(i+1)%self._ringLength]])
            self._align.setRefCoords(newref)
            self._align.rigidFit(refs[i])
            indices = self._attachedIndices[i]
            coords[indices] = self._align.transformCoords(coords[indices])[:,:3]
        ##     if self._subTree[i]: # there is an attachment we save the matrix
        ##         self._matrices[i] = self._align.getMatrix()
        ##         #print 'bbb', i, self._matrices[i]
        return coords

    def conformerCoords(self, angleList):
        coreIndices = self._rigidBodyIndices
        allCoords = numpy.zeros( (len(angleList), len(coreIndices), 3), 'f')

        # compute reference frame for each ring atom
        refs = []
        for i in range(len(self._ringIndices)):
            refs.append(self.refFrame(self._coords[self._ringIndices[i-1]],
                                      self._coords[self._ringIndices[i]],
                                      self._coords[self._ringIndices[(i+1)%self._ringLength]]))

        # find index of last, first and second atom in ring within the list of
        # indices of that rigid body to compute reference frame used to align rings
        
        # compute coordinate system at ring entry point
        # This inverse of this transformation needs to be added to the entire ring
        entryRef = self.refFrame(self._coords[self._ringIndices[-1]],
                                 self._coords[self._ringIndices[0]],
                                 self._coords[self._ringIndices[1]])
        # we need a matrix at each ring atom for each conformer
        matrices = []
        for n in range(len(self._subTreeInd)):
            matrices.append([])
        nb = 0
        #print 'RING', self._ringIndices
        #for name, angles in angleList:#[3:4]:
        for angles in angleList:#[3:4]:
            coords1 = self.coordsForThetas(angles)
            # coordinate at ring entry point in this conformer
            newref = self.refFrame(coords1[self._ringIndices[-1]],
                                  coords1[self._ringIndices[0]],
                                  coords1[self._ringIndices[1]])
            self._align.setRefCoords(entryRef)
            self._align.rigidFit(newref)
            coords2  = self._align.transformCoords(coords1)[:,:3]
            allCoords[nb] = coords2[coreIndices]

            # recompute transformation for attachment points to account
            # for transformation used to align conformers to input structure
            #nn = 0
            for n, i in enumerate(self._subTreeInd):
                newref = self.refFrame(coords2[self._ringIndices[i-1]],
                                       coords2[self._ringIndices[i]],
                                       coords2[self._ringIndices[(i+1)%self._ringLength]])
                self._align.setRefCoords(newref)
                self._align.rigidFit(refs[i])
                matrices[n].append(
                    self._align.getMatrix().transpose().tolist())
            nb+=1

        return allCoords, matrices


class FlexibleRings:
    """ given a molecule, perceive flexible rings of maxSize or within
        a list of acceptable sizes (acceptSize)

        >>> acceptedSize = [5,6,7] # list of ring sizes accepted
        >>> rigid_body_indices = [20,21,22,23,24, 25,44,45,46,47,48,49,50,51,52]
        >>> excludeDoubleBonds = True # rings containing at least one double bond will be rejected
        >>> multi = False # stop after the first ring is found
        >>> mol = ob.OBMol()
        >>> rp = RingPerception(mol)
        >>> rp.findFlexibleRing(rigid_body_indices, maxSize, acceptedSizes, excludeDoubleBonds, multi)
        Processing ring [8, 9, 10, 11, 30] FAIL: Ring size not acceptable 5 [6]
        Processing ring [2, 3, 4, 5, 6, 7] FAIL: Ring is aromatic
        Processing ring [14, 15, 16, 17, 18, 19] FAIL: Ring is aromatic
        Processing ring [20, 21, 22, 23, 24, 25] *** FOUND RING ***
        Ring found [[19, 20, 21, 22, 23, 24]]
        Ring members [[[42], [43, 44], [45, 46], [], [47, 48], [49, 50]]]
        >>> rings = rp.getRings()
        [[19, 20, 21, 22, 23, 24]]
        >>> attached = rp.getAttached()
        [[[42], [43, 44], [45, 46], [], [47, 48], [49, 50]]]
    """

    def __init__(self, obmol, perceivedRings=None): 
        """ object constructor
        the object gets built for a molecule and its findFlexibleRing can be
        called repeatedly for various rigid bodies of that molecule to identify
        flexible rings within this rigid body
        """
        assert isinstance(obmol, ob.OBMol)
        self.obmol = obmol
        if perceivedRings is None:
            self.perceivedRings = obmol.GetSSSR()
        else:
            self.perceivedRings = perceivedRings

    def _initialize(self):
        self._rings = []
        self._attachedInd = []
        self._rigidBodyIdx = []
        
    def findFlexibleRing(self, rigidBodyIdx, acceptedSizes=[6], 
                         excludeDoubleBonds=True, multi=False):
        """
        for a given subset of atoms specified through a list of atom indices in self.obmol
        this method will identify flexible rings. Rings and atoms moving together
        (as a rigid body) with each ring atom can be retrieved using as getRings() and
        getAttached()
        """
        #print 'OBRB', rigidBodyIdx
        self._initialize()
        self._rigidBodyIdx = rigidBodyIdx
        # find rings of acceptable sizes
        self._scanRings(acceptedSizes=acceptedSizes, multi=multi)
        # exclude rings containing not entirely made of single bonds
        if excludeDoubleBonds:
            self._checkBonds()
        # ?? not sure what this does
        self._completeRings()
        rings = self._getRings()
        attached = self._getAttached()
        return rings, attached
    
    def _scanRings(self, acceptedSizes, multi):
        """ find rings of size acceptedSizes in the current rigid body
            if multi is false, only the first one will be returned
        """
        found = False
        for r in self.perceivedRings:
            ringAtoms = list(r._path)
            #print "Ring", ringAtoms,
            # discard aromatics
            if r.IsAromatic():
                #print "IGNORED: Ring is aromatic"
                continue
            ringAtoms = list(r._path)
            # ring of desired size
            if not len(ringAtoms) in acceptedSizes:
                #print "IGNORED: Ring size not acceptable", len(ringAtoms), acceptedSizes
                continue
            # this ring is entirely in the rigid body
            if len( set(self._rigidBodyIdx) & set(ringAtoms) ) == len(ringAtoms):
                #print 'Flexbile Ring', r._path
                self._rings.append(r._path)
                
                if not multi:
                    #print "Multi is OFF, returning first ring"
                    break
            #else:
                #print "Ring/rigid body overlap failed"
                #print len( set(self._rigidBodyIdx) & set(ringAtoms) ) == len(ringAtoms)
                #print self._rigidBodyIdx
                #print ringAtoms
                    
    def _checkBonds(self):
        """ check if the ring contains any double bonds"""
        if self._rings == []:
            return
        remove = []
        for r in self._rings:
            size = len(r)
            for i in range(size+1):
                idx1 = i % size
                idx2 = (i+1) % size
                bond = self.obmol.GetBond(r[idx1],r[idx2])
                if not bond.GetBO() == 1:
                    #print "FAIL non-single bond found %d %d! removing"%(
                    #    r[idx1], r[idx2]),r
                    remove.append(r)
                    break
                atom = self.obmol.GetAtom(r[idx1])
                if not atom.GetHyb() == 3:
                    print "FAIL non-sp3 atom in ring found"
                    remove.append(r)
                    break
        for x in remove:
            self._rings.remove(x)
        
    def _completeRings(self):
        """ this part will be modified to walk through
            rotatable bonds in the ring and simplify
            completion when aromatic rings are attached
            to a flexible ring
        """
        if self._rings == []:
            return
        remove = []
        for r in self._rings:
            attached = []
            accepted = True
            for aIdx in r:
                complete = []
                a = self.obmol.GetAtom(aIdx)
                for neigh in ob.OBAtomAtomIter(a):
                    n = neigh.GetIdx()
                    # the neighbor is another atom in this ring
                    if n in r:
                        continue
                    # the neighbor is in another rigid body
                    if not n in self._rigidBodyIdx:
                        continue
                    # the neighbor is in a ring (which must be another ring, now)
                    if neigh.IsInRing():
                        # not accepted, marked for deletion
                        remove.append(r)
                        accepted = False
                        break
                    complete.append(n)
                if not accepted:
                    break
                else:
                    attached.append(complete)
                    
            # attached atoms are saved only if ring is accepted
            if accepted:
                self._attachedInd.append(attached)
            else:
                self._attachedInd.append(None)

        for r in remove:
            #print 'FAGA remove', r
            idx = self._rings.index(r)
            self._rings = self._rings[:idx]+self._rings[idx+1:]
            self._attachedInd = self._attachedInd[:idx] + self._attachedInd[idx+1:]

        # if self._attachedInd contains a "None", there's a problem

    ## def _completeRings(self):
    ##     """ this part will be modified to walk through 
    ##         rotatable bonds in the ring and simplify 
    ##         completion when aromatic rings are attached
    ##         to a flexible ring
    ##     """
    ##     if self._rings == []:
    ##         return
    ##     for r in self._rings:
    ##         attached = []
    ##         for aIdx in r:
    ##             complete = []
    ##             a = self.obmol.GetAtom(aIdx)
    ##             for neigh in ob.OBAtomAtomIter(a):
    ##                 n = neigh.GetIdx()
    ##                 # the neighbor is another atom in this ring
    ##                 if n in r:
    ##                     continue
    ##                 # the neighbor is in another rigid body
    ##                 if not n in self._rigidBodyIdx:
    ##                     continue
    ##                 # the neighbor is in a ring (which must be another ring, now)
    ##                 if neigh.IsInRing():
    ##                     continue
    ##                 complete.append(n)
    ##             attached.append(complete)
    ##         self._attachedInd.append(attached)

    def _getRings(self):
        """ return 0-based indices"""
        rings = []
        for r in self._rings:
            rings.append([ x-1 for x in r]) 
        return rings

    def _getAttached(self):
        """ return 0-based indices"""
        attached = []
        for aList in self._attachedInd:
            aRing = []
            for a in aList:
                aRing.append( [x-1 for x in a])
            attached.append(aRing)
        return attached
