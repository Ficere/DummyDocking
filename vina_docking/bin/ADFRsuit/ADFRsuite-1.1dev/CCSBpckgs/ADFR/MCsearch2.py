import sys
from MolKit2 import Read
from ADFR import ADFR, getLigandFromFile
from ADFR.utils.maps import MapsFile, flexResStr2flexRes
from ADFR.GA import Population
from mglutil.math.rmsd import RMSDCalculator

from geomutils import efitlib
ellipse = efitlib.ellipsoid()
ellipseInfo = efitlib.efit_info()
cov_scale = 1.75
ell_scale = 1.0

def minimizeBFGS(bfgs, individual, *args, **kw):
    before = individual._score
    score = bfgs.minimize(individual, **kw)
    individual.score() # force the Python Individual to update _score and phenotype etc...
    return score, bfgs.nevals

## def _mkADFR(targetFile, ligandFile, seed):

##     if ligandFile is None:
##         ligandFile = '%s_lig.pdbqt'%self.ligand.name
##         self.agfr.saveLigandAsPDBQT(self.ligand.name)
##         self.ligand.filename = ligandFile

##     ligand, error = getLigandFromFile(ligandFile)

##     mf = MapsFile(targetFile)
##     mf.unzipMaps()
##     mapsFolder = mf.getMapsFolder()
##     receptorFilename = os.path.join(mf.getMapsFolder(),
##                                     mf.getReceptorFilename())
##     receptor = Read(receptorFilename)
##     flexResStr = mf.getFlexResStr()
##     covalentRec = mf.getCovalentBond()

##     if covalentRec is not None:
##         covalentRec += [d.get('covalentBondTorsionAtoms')[-1]]
##     else:
##         covalentRec = None    
##     tPtsFile = os.path.join(mapsFolder, "translationPoints.npy")
##     if os.path.exists(tPtsFile):
##         tpointsFilename = os.path.join(mapsFolder, "translationPoints.npy")
##     else:
##         tpointsFilename = None
##     mapFilesRoot = 'rigidReceptor'

##     adfr = ADFR(ligand, mapsFolder, logFile=None,
##                 receptor=receptor,
##                 flexibleResidues=flexResStr,
##                 mapFilesRoot=mapFilesRoot,
##                 tpointsFilename=tpointsFilename,
##                 seedValue=seed,
##                 ## localSearchMode = args['localSearchMode'],
##                 ## covalentIndices=covalentIndices,
##                 ## FTRecsrc=args.get('FTRecsrc', None),
##                 ## fixedRoot=fixedRoot,
##                 ## neighborSearchCutoff=neighborSearchCutoff
##                 )
##     return adfr

from math import exp, sqrt, pi, ceil, cos, sin
from random import random, randint
import numpy as np
from numpy.random import choice
 
def metropolis_oracle(old_e, new_e, temperature):
    if (new_e < old_e): return True
    acceptance_probability = exp((old_e - new_e) / temperature)
    return random() < acceptance_probability

def random_direction():
    while True: 
        r1 = 2*random() - 1
        r2 = 2*random() - 1
        r3 = 2*random() - 1
        if (r1*r1 + r2*r2 + r3*r3) < 1.0:
            return (r1,r2,r3)

def gyration_radius(fw):
    coords = fw.phenotype[1] 
    n = len(coords)
    g = np.sum(coords, 0)/n
    d = (coords-g)
    dist = np.sqrt(np.sum(d*d, 1))
    gr  = sqrt(np.sum(dist)/n)
    return gr

def normalize_angle(x):
    # subtract or add enough 2*pi's to make x be in [-pi, pi]
    if x > 3*pi: # very large
        n = (x - pi) / (2*pi) # how many 2*pi's do you want to subtract?
        x -= 2*pi*ceil(n) # ceil can be very slow, but this should not be called often

    elif x < -3*pi: # very small
        n = (-x - pi) / (2*pi) # how many 2*pi's do you want to add?
        x += 2*pi*ceil(n) # ceil can be very slow, but this should not be called often

    if x > pi: # in (   pi, 3*pi]
        return x - 2*pi;

    elif x < -pi: # in [-3*pi,  -pi)
        return x + 2*pi;

def axis_angle_to_quaternion(axis, angle):
    # axis is assumed to be a unit vector
    normalize_angle(angle) # this is probably only necessary if angles can be very big
    c = cos(angle/2);
    s = sin(angle/2);
    return (c, s*axis[0], s*axis[1], s*axis[2])

def angle_to_quaternion(rotation):
    angle = np.linalg.norm(rotation)
    if angle > sys.float_info.epsilon:
        angle1 = 1.0/angle
        axis = [angle1*rotation[0], angle1*rotation[1], angle1*rotation[2]]
        return axis_angle_to_quaternion(axis, angle)
    return [0.5, 0.5, 0.5, 1.0] # identity rotation

def quaternion_normalize_approx(q, tolerance = 1e-6):
    s = q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]
    if abs(s - 1) < tolerance:
        return q
    else:
        a = 1./sqrt(s)
        return [q[0]*a, q[1]*a, q[2]*a, q[3]*a]

def quaternion_multiply(quaternion1, quaternion0):
    x0, y0, z0, w0 = quaternion0
    x1, y1, z1, w1 = quaternion1
    return np.array([x1 * w0 + y1 * z0 - z1 * y0 + w1 * x0,
                     -x1 * z0 + y1 * w0 + z1 * x0 + w1 * y0,
                     x1 * y0 - y1 * x0 + z1 * w0 + w1 * z0,
                     -x1 * x0 - y1 * y0 - z1 * z0 + w1 * w0], dtype=np.float64)

def quaternion_increment(q, rotation):
	q = quaternion_multiply(angle_to_quaternion(rotation), q)
	return quaternion_normalize_approx(q)

def flipPose(pose):
    coords = pose.phenotype[1].astype('f')
    status = efitlib.fitEllipse(coords, ell_scale, cov_scale, ellipseInfo, ellipse)
    if status==0: # OK
        center = ellipse.getPosition() 
	size = ellipse.getAxis().astype('f')
	orient = np.identity(4, 'f')
	orient[:3, :3] = ellipse.getOrientation()
	axis = randint(0,2)
	qx, qy, qz = orient[axis, :3]
        qw = 0.5 # 180 degrees ?
        genes = pose.getGenesPy()
        flipQuat = quaternion_multiply(genes[3:7], [qx, qy, qz, qw])
        genes[3:7] = quaternion_normalize_approx(flipQuat)
        pose.setGenes(genes)
    return status

def pose_mutate_orientation(pose, amplitude, flip):
    # for rotation we migh not have to go through variables, but can
    # change genes right away
    if flip:
        if random()<.1: # 10% change to do a flip
    	    status = flipPose(pose)
    	    if status==0: # ok
                #print "FLIP", pose._score
    	        return
    gr = gyration_radius(pose)
    if gr>sys.float_info.epsilon:
        amp = amplitude / gr
        v = random_direction()
        rotation = [amp*v[0], amp*v[1], amp*v[2]]
        # configure the genome with the genes of this individual
        genes = pose.getGenesPy()
        pose.genomePy.setGenes(genes)
        # get the variables
        var = pose.genomePy.getVariablesPy()
        newquat = quaternion_increment(var[3:7], rotation)
        var[3:7] = newquat
        newgenes = pose.genomePy.getGenesForValuesPy(var)
        pose.setGenes(newgenes)
        
def pose_mutate_position(pose, amplitude):
    genes = pose.getGenesPy()
    pose.genomePy.setGenes(genes)
    var = pose.genomePy.getVariablesPy()
    dt = random_direction()
    for i in xrange(3):
        var[i] += amplitude * dt[i]
    newgenes = pose.genomePy.getGenesForValuesPy(var)
    pose.setGenes(newgenes)

def pose_mutate_torsion(pose):
    # randomize torsion
    genes = pose.getGenesPy()
    if len(genes)==7: return # no torsions
    pose.genomePy.setGenes(genes)
    var = pose.genomePy.getVariablesPy()
    if len(genes)>10 and random()<.1:
        #print '2 TORSIONS', pose._score
        # 2 torsion 180
	index = randint(7,len(genes)-1)
        #var[index] = random()*2*pi
	var[index] += pi# + (2*random()-1)*pi*0.5
	index = randint(7,len(genes)-1)
        #var[index] = random()*2*pi
	var[index] += pi# + (2*random()-1)*pi*0.5
    else:
	    # randomize 1 torsion
	index = randint(7,len(genes)-1)
        #var[index] = random()*2*pi
	var[index] = random()*2*pi - pi
    #for index in xrange(7,len(genes)-1):
    #    var[index] = random()*2*pi-pi
    newgenes = pose.genomePy.getGenesForValuesPy(var)
    pose.setGenes(newgenes)

## def pose_mutate(pose, amplitude):
##     genes = pose.getGenesPy()
##     pose.genomePy.setGenes(genes)
##     var = pose.genomePy.getVariablesPy()
##     which = randint(0, len(genes)-5)
    
##     if which == 0:
##         dt = amplitude * random_direction()
##         for i in xrange(3):
##             var[i] += dt[i]
##     which -= 1

##     if which == 0:
##         gr = gyration_radius(pose) 
##         if gr>sys.float_info.epsilon:
##             amp = amplitude / gr
##             v = random_direction()
##             rotation = [amp*v[0], amp*v[1], amp*v[2]]
##             newquat = quaternion_increment(var[3:7], rotation)
##             var[3:7] = newquat
##     which -= 1

##     if which < len(genes)-7:
##         var[7+which] = random()*2*pi-pi

##     newgenes = pose.genomePy.getGenesForValuesPy(var)
##     pose.setGenes(newgenes)

def shelfPose(rmsdCalc, poseBuffer, pose, solutions):
    # check RMSD with solution and avoid them since we got stuck
    if pose._score < 0: return
    rmsdCalc.setRefCoords(pose.phenotype[1])
    #for sol in solutions:
    #	  rmsdCalc.setRefCoords(sol.phenotype[1])
    #	  rmsd = rmsdCalc.computeRMSD(pose.phenotype[1])
    #	  if rmsd < 2.0:
    #	      print 'Not Shelving pose as already seen', -sol._score, -pose._score
    #	      return
    #ind = 0
    scoreDiff = [0]*len(poseBuffer)
    for i in range(len(poseBuffer)):
        rmsd = rmsdCalc.computeRMSD(poseBuffer[i].phenotype[1])
	if rmsd < 2.:
	    if pose._score > poseBuffer[i]._score:
	        poseBuffer[i] = pose
	    return
        scoreDiff[i] = pose._score - poseBuffer[i]._score
        #if poseBuffer[i]._score > pose._score:
	#    ind += 1

    #import pdb; pdb.set_trace()
    poseBuffer[np.argmax(scoreDiff)] = pose
    #if  pose._score > poseBuffer[ind]._score:
    #    poseBuffer[ind] = pose
    print 'Shelf add', -pose._score, [(i, -x._score) for i,x in enumerate(poseBuffer) if x._score>0]

def MCsearch(pose, maxEvals=200000, maxTry=2500, temperature=1.2, amplitude=2, flip=True, noImprove=450, debug=False):

    from ADFRcc.adfr import BFGSQuat
    minikw = {"aMaxSteps" : 20,
              "xtol_abs" : -1,
              "xtol_rel" : -1,
              "ftol_abs" : -1,
              "ftol_rel" : -1,
              "trajFilename" : None
              }
    miniHardkw = {"aMaxSteps" : 60,
              "xtol_abs" : -1,
              "xtol_rel" : -1,
              "ftol_abs" : -1,
              "ftol_rel" : -1,
              "trajFilename" : None
              }

    poseBuffer = [pose.clone() for i in range(10)]

    # minimize th pose
    bfgs = BFGSQuat(pose)
    score, steps = minimizeBFGS(bfgs, pose, **minikw)

    # best score sceen in this search
    bestScore = pose._score      # best in local search
    bestPose = pose.clone()

    # gbestPose is only used for printing progress and check that
    # we do not lose the solution
    gbestPose = pose.clone()     # global best
    gbestScore = gbestPose._score # globally best
    # list of ssolutions found during this search
    solutions = []
    
    rmsdCalc = RMSDCalculator()
    done = False
    countNoImprove = 0
    countTry = 0

    while True:
        # create pose to try MC moves
        candidate = pose.clone()
        #assert candidate._score==pose._score
        
        # apply MC move to candidate
        nn = random()
        if nn < 0.3:
            pose_mutate_orientation(candidate, amplitude, flip)
        elif nn < 0.6:
            pose_mutate_position(candidate, amplitude)
        else:
            pose_mutate_torsion(candidate)
        #print 'step', nn, pose._score, candidate._score

        # minimize candidate
        score, steps = minimizeBFGS(bfgs, candidate, **miniHardkw)

        if metropolis_oracle(-bestPose._score, -candidate._score, temperature):
            # move is accepted
            sdiff = candidate._score - bestPose._score
            if sdiff > 0.0:
                # move is improving best seen 
                bestPose = candidate.clone()
                #if bestPose._score > bestScore:
                bestScore = bestPose._score
                if sdiff > 0.3:
                    # shelf the pose that lead to this improvement
                    shelfPose(rmsdCalc, poseBuffer, pose, solutions)
                if bestScore > gbestScore:
                    gbestScore = bestScore
                    gbestPose = bestPose.clone()
                print 'LOG', -bestPose._score, pose.genomePy.scorer.getNumEvals(), countTry, countNoImprove, maxEvals, maxTry, bestScore, gbestScore
                sys.stdout.flush()
                #import pdb; pdb.set_trace()
                countTry = 0
            pose = candidate
        else:
            countTry += 1
            if countTry > maxTry:
                countNoImprove += 1
                if countNoImprove > noImprove:
                    status = 'maxNoImprovement %d reached'%noImprove
                    break # end this line of search
                poses = [x for x in poseBuffer if x._score > 0]
                #l = len([x for x in poseBuffer if x._score > 0])
                solutions.append(bestPose)
                if debug:
                    print 'add solution', bestPose._score, len(solutions), countTry
                if len(poses)-1 > 0:
                    scores = [-x._score for x in poses]
                    proba = scores / np.sum(scores)
                    #solutions.append(bestPose)
                    bestPose = pose = poseBuffer[choice(range(len(scores)), 1, p=proba)]
                    #bestPose = pose = poseBuffer[(randint(0, l-1))].clone()
                    if debug:
                        print 'UNSHELF', countTry, -pose._score
                    countTry = 0
                else:
                    #if pose._score < 0.0:
                    pose.randomize()
                    score, steps = minimizeBFGS(bfgs, pose, **minikw)
                    bestPose = pose.clone()
                    if debug:
                        print 'RANODOMIZE', score
                    countTry = 0

        if pose.genomePy.scorer.getNumEvals() > maxEvals:
            status = 'maxEvals %d reached'%maxEvals
            break

    solutions.append(bestPose)

    best = solutions[0]
    for i, s in enumerate(solutions[1:]):
        if s._score > best. _score: best = s
	print s._score
    assert best._score==gbestScore

    print 'SEARCH_END', -best._score, countTry, countNoImprove, pose.genomePy.scorer.getNumEvals()
    return best, status

if __name__=="__main__":
    import sys, os
    from ADFR.utils.optParser import ArgParser
    from ADFR.GA import Population
    from ADFR.utils.cluster import clusterPoses, oneCluster
    from mglutil.math.rmsd import RMSDCalculator

    parser = ArgParser('OneGA')
    args = parser.parse_args()

    logFile = args['logfile']
    if logFile:
        f = open(logFile, 'w')
        path = os.path.split(logFile)[0]
    else:
        f = None
        path = None

    debug = False

    pdbid = os.path.basename(args['receptorMapsFile']).split('_')[0][:4]

    adfr = _mkADFR(args['receptorMapsFile'],
                   args['ligandFile'],
                   args['seedValue'])

    print 'seed', adfr.initialSeed

    if args['maxEvals']:
        maxEvals = args['maxEvals']
    if maxEvals==2500000:
        maxEvals = 200000 + 60000*adfr.ligandFT.torTree.torsdof
    #maxTry = 1000 + 50*adfr.ligandFT.torTree.torsdof
    maxTry = 500
    pose,status = search(adfr, maxEvals=maxEvals, maxTry= maxTry, debug=debug)

    if args['reference']:
        ref = Read(args['reference'])
        rmsdCalc0 = RMSDCalculator(ref._ag.getCoords())
    
        from mglutil.math.rmsd import MorpicRMSD
        _autos = adfr.ligand.getAutomorphisms()
        assert len(_autos)>0

        rmsdCalcMorph = MorpicRMSD(_autos)
        rmsdCalcMorph.setRefCoords(ref._ag.getCoords())
    else:
        rmsdCalcMorph = rmsdCalc0 = ref = None
        
    if pose._score >0:
        if f:
	    f.write("END %s %9.5f %5.2f %8d %8d %d %d %s\n"%(
		    pdbid, -pose._score, rmsdCalcMorph.computeRMSD(pose.phenotype[1]),
		    pose.genomePy.scorer.getNumEvals(), adfr.initialSeed, maxEvals, maxTry, status))
	print 'END %s score: %9.5f RMSD: %5.2f evals: %8d seed: %8d maxEvals: %d maxTry: %d status: "%s"\n'%(
		pdbid, -pose._score, rmsdCalcMorph.computeRMSD(pose.phenotype[1]),
		pose.genomePy.scorer.getNumEvals(), adfr.initialSeed, maxEvals, maxTry,status)

    if f:
        f.close()

    if path is not None:
        adfr.writeSolution(pose, adfr.ligandFT.torTree, '%s_%d_out.pdbqt'%(
            os.path.join(path,pdbid), adfr.initialSeed), 0)
    else:
        adfr.writeSolution(pose, adfr.ligandFT.torTree, '%s_%d_out.pdbqt'%(pdbid, adfr.initialSeed), 0)
#adpy fw3.py -t ../6ggd_rec.trg -l ../6ggd_lig.pdbqt --seed 3
#adpy fw3.py -t ../1xd0_protein.trg -l ../1xd0_ligand.pdbqt --seed 3
#time adpy fw3.py -t astex/targets/1kzk_rec.trg -l 1kzk.pdb_lig_001.pdbqt
#export f="1g9v"; time adpy fw3.py -t astex/targets/${f}_recGrad.trg -l  astex/ligands_randomized/${f}_random.pdbqt -r astex/ligands/${f}_lig.pdbqt
#export f="1gkc"; time adpy fw3.py -t astex/targets/${f}_recGrad.trg -l  astex/ligands_randomized/${f}_random.pdbqt -r astex/ligands/${f}_lig.pdbqt
#export f="1w1p"; time adpy fw3.py -t astex/targets/${f}_recGrad.trg -l  astex/ligands_randomized/${f}_random.pdbqt -r astex/ligands/${f}_lig.pdbqt
#export f="1kzk"; time adpy fw3.py -t astex/targets/${f}_recGrad.trg -l  astex/ligands_randomized/${f}_random.pdbqt -r astex/ligands/${f}_lig.pdbqt
#export f="1sj0"; time adpy fw3.py -t astex/targets/${f}_recGrad.trg -l  astex/ligands_randomized/${f}_random.pdbqt -r astex/ligands/${f}_lig.pdbqt

#for f in `cat astex.list`; do adpy fw3.py -t astex/targets/${f}_recGrad.trg -l  astex/ligands_randomized/${f}_random.pdbqt -r astex/ligands/${f}_lig.pdbqt | grep END0; done

#print 'seed', adfr.initialSeed, time()-t0
## END0 1sj0   246 -19.15093  0.81  51   985459   3 980351
## END1 1sj0   246 -17.74063  1.83  51   985459   3 980351
## END2 1sj0   246 -17.26760  5.32  51   985459   3 980351
## END3 1sj0   246 -15.48212  5.49  51   985459   3 980351
## END4 1sj0   246 -14.25362  6.98  51   985459   3 980351
## END5 1sj0   246 -14.09393  6.86  51   985459   3 980351
## END6 1sj0   246 -13.85789  7.11  51   985459   3 980351
## END7 1sj0   246 -13.47605  5.67  51   985459   3 980351
## END8 1sj0   246 -13.13968  7.80  51   985459   3 980351

#END0 1kzk   346 -19.64613  0.72 1387383 316403

#export f="1w1p"; time adpy dock.py -t astex/targets/${f}_recGrad.trg -l  astex/ligands_randomized/${f}_random.pdbqt -r astex/ligands/${f}_lig.pdbqt

#export f="1kzk"; time adpy search2.py -t astex/targets/${f}_recGrad.trg -l  astex/ligands_randomized/${f}_random.pdbqt -r astex/ligands/${f}_lig.pdbqt
#adpy search2.py -t ../1xd0_protein.trg -l ../1xd0_ligand.pdbqt -r ../1xd0_ligand.pdbqt
#adpy search2.py -t ../1xd0_protein.trg -l ../1xd0_ligand.pdbqt -r ../1xd0_ligand.pdbqt
