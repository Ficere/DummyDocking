################################################################################
##
## This library is free software; you can redistribute it and/or
## modify it under the terms of the GNU Lesser General Public
## License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## 
## This library is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## Lesser General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public
## License along with this library; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
##
## (C) Copyrights Dr. Michel F. Sanner and TSRI 2016
##
################################################################################

############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner and TSRI 2020
#
#########################################################################

import csv
from collections import OrderedDict

import numpy as np

from ADFR import molToAtomSetStatic
from ADFRcc.adfr import AtomSet, PairwiseScorer, Parameters

parameters = Parameters.getParameters()
feCoeffVdw = parameters.feCoeffVdw
feCoeffHbond = parameters.feCoeffHbond
feCoeffEstat = parameters.feCoeffEstat
feCoeffDesolv = parameters.feCoeffDesolv
feCoeffTors = parameters.feCoeffTors

class ListInteractionEnergyTerms(object):
    """Object for obtaining lists of interaction energies terms between 2
    prody atomSets.
    An ADFRcc PairWiser scorer object is built and used to score the
    interaction between the prody 2 atom sets. Energy terms for atom pairs
    can be retrieved based on: distance or using prody selection strings.

    the list of energetic terms can be saved as a csv file using
    self.saveCSV(filename, pairwiseTerms)

    examples:

    # 1 - get all pairs withing 6.0 Angstroms and report energetic terms
    # for these pairs at the atomic level
    pairs = self.getPairs(distance=6.0)
    pairwiseTerms = self.getInteractions(pairs)
    self.saveCSV('atomicPWT.csv', pairwiseTerms)
    
    # 2 - get pairs between ligand oxygen and nitrogen atoms and all
    # receptor atoms  and report energetic terms for these pairs per
    # residues in atomset1
    pairs = self.getPairs(sel2='element O N')
    pairwiseTermsRes = self.getInteractions(pairs, mode='perSet1Residue')
    self.saveCSV('recResPWT.csv', pairwiseTermsRes)
    """

    def __init__(self, atoms1, atoms2):
        """
        atoms1 and atoms2 are prody atoms set.
        """
        super(ListInteractionEnergyTerms, self).__init__()
        self.scorer = None
        self._atoms1 = atoms1
        self._atoms2 = atoms2
        recAS = AtomSet(molToAtomSetStatic(rec))
        ligAS = AtomSet(molToAtomSetStatic(lig))

        scorer = PairwiseScorer()
        scorer.initialize(recAS, ligAS)
        score = scorer.calculateScores()
        self.scorer = scorer

        self._distance = scorer.getDistanceMatrix()
        self._vdwArray = scorer.getVdwEnergyMatrix()*feCoeffVdw
        self._eArray = scorer.getEstatEnergyMatrix()*feCoeffEstat
        self._hArray = scorer.getHbEnergyMatrix()*feCoeffHbond
        self._dsArray = scorer.getSolvEnergyMatrix()*feCoeffDesolv
        self._names1 = self._atoms1.getNames()
        self._names2 = self._atoms2.getNames()
        self._resnames1 = self._atoms1.getResnames()
        self._resnames2 = self._atoms2.getResnames()
        self._resnums1 = self._atoms1.getResnums()
        self._resnums2 = self._atoms2.getResnums()
        self._icodes1 = self._atoms1.getIcodes()
        self._icodes2 = self._atoms2.getIcodes()
        self._chids1 = self._atoms1.getChids()
        self._chids2 = self._atoms2.getChids()

    def getPairs(self, distance=None, sel1='all', sel2='all'):
        # return list of i,j pairs of receptor-ligand atoms within distance
        # i and j are 0-based atom indices in the prody atom sets
        # use 'all' to refer to all atoms in a given set
        self._sel1 = sel1
        self._sel2 = sel2
        if distance is not None:
            if sel1 == 'all':
                self._sel1 = 'd<%g'%distance
            else:
                self._sel1 += ' and d<%g'%distance
            if sel2 == 'all':
                self._sel2 = 'd<%g'%distance
            else:
                self._sel2 += ' and d<%g'%distance
            
        if sel1=='all':
            sel1Ind = range(len(self._atoms1))
        else:
            sel1Ind = self._atoms1.select(sel1).getIndices()

        if sel2=='all':
            sel2Ind = range(len(self._atoms2))
        else:
            sel2Ind = self._atoms2.select(sel2).getIndices()
        pairs = []
        for i in sel1Ind:
            for j in sel2Ind:
                if distance is not None and self._distance[i,j]<=distance:
                    pairs.append((i,j))
        return pairs

    def getInteractions(self, pairs, mode='perAtom'):
        # pairs or 0-based indices into atoms1 and atoms2
        # mode can be:
        #   "perAtom" to report terms between each pair of atoms
        #   "perSet1Residue" or "perSet2Residue" to report terms summed up over
        #     residues in either set1 or set2.
        # return a list with each entry being  a list containing:
        #   - (chid1, resname1, resnum1, icode1, fullname1)
        #   - (chid2, resname2, resnum2, icode2, fullname2)
        #   - (vdw, elec, hb, ds, total)
        values = []
        if mode == 'perAtom':
            for i,j in pairs:
                vdw  = self._vdwArray[i,j]
                elec = self._eArray[i,j]
                hb   = self._hArray[i,j]
                ds   = self._dsArray[i,j]
                tot  = vdw+elec+hb+ds
                key1 = '%s:%s`%d%s:%s'%(
                    self._chids1[i], self._resnames1[i], self._resnums1[i],
                    self._icodes1[i], self._names1[i])
                key2 = '%s:%s`%d%s:%s'%(
                    self._chids2[j], self._resnames2[j], self._resnums2[j],
                    self._icodes2[j], self._names2[j])
                values.append(
                    [self._chids1[i], self._resnames1[i], self._resnums1[i],
                      self._icodes1[i], self._names1[i], key1,
                     self._chids2[j], self._resnames1[j], self._resnums2[j],
                      self._icodes2[j], self._names1[j], key2,
                     vdw, elec, hb, ds, tot])
            return values

        else: # we need to sum over residues
            d = OrderedDict() # used to sum over residues
            for i,j in pairs:
                vdw  = self._vdwArray[i,j]
                elec = self._eArray[i,j]
                hb   = self._hArray[i,j]
                ds   = self._dsArray[i,j]
                tot  = vdw+elec+hb+ds
                if mode == 'perSet1Residue':
                    key = '%s:%s`%d%s'%(self._chids1[i], self._resnames1[i],
                                        self._resnums1[i], self._icodes1[i])

                elif mode == 'perSet2Residue':
                    key = '%s:%s`%d%s'%(self._chids2[j], self._resnames2[j],
                                        self._resnums2[j], self._icodes2[j])

                else:
                    raise ValueError('Invalid mode "%s" for ListInteractionEnergyTerms.getInteractions().\nUse "perAtom", "perSet1Residue" or "perSet2Residue'%mode)
                try:
                    d[key][12] += vdw
                    d[key][13] += elec
                    d[key][14] += hb
                    d[key][15] += ds
                    d[key][16] += tot
                except KeyError:
                    if mode == 'perSet1Residue':
                        d[key] = [
                            self._chids1[i], self._resnames1[i],
                            self._resnums1[i], self._icodes1[i],
                            self._names1[i], key,
                            '*', '*', '*', '*', '*', self._sel1,
                            vdw, elec, hb, ds, tot]
                    elif mode == 'perSet2Residue':
                        d[key] = [
                            '*', '*', '*', '*', '*', self._sel2,
                            self._chids2[j], self._resnames2[j],
                            self._resnums2[j], self._icodes2[j],
                            self._names2[j], key,
                            vdw, elec, hb, ds, tot]

            values = []
            for k,v in d.items():
                values.append(v)
            return values

    def saveCSV(self, filename, pairwiseTerms):
        """
        save pariwise terms as a csv file
        """
        # field names  
        fields = ['chid1', 'resname1', 'resnum1', 'icode1', 'name1', 'fullname1',
                  'chid2', 'resname2', 'resnum2', 'icode2', 'name2', 'fullname2',
                  'vdw', 'elec', 'hb', 'ds', 'total']

        # writing to csv file  
        with open(filename, 'w') as csvfile:  
            # creating a csv writer object  
            csvwriter = csv.writer(csvfile)  
        
            # writing the fields  
            csvwriter.writerow(fields)  
        
            # writing the data rows  
            csvwriter.writerows(pairwiseTerms) 

if __name__=='__main__':
    from MolKit2 import Read
    rec = Read('/home/sanner/local/python/AGFRPaper/4ek3_rec/4EK3_rec.pdbqt')
    lig = Read('/home/sanner/local/python/AGFRPaper/4EK4_lig.pdbqt')

    IElister = ListInteractionEnergyTerms(rec._ag, lig._ag)
    pairs = IElister.getPairs(distance=100.0)
    pairwiseTermsAtoms = IElister.getInteractions(pairs)
    IElister.saveCSV('atomicPWT.csv', pairwiseTermsAtoms)
    for v in pairwiseTermsAtoms[:10]:
        print(v)
    print()
    
    pairwiseTermsRes = IElister.getInteractions(pairs, mode='perSet1Residue')
    IElister.saveCSV('recResPWT.csv', pairwiseTermsRes)
    for v in pairwiseTermsRes[:10]:
        print(v)

    pairs = IElister.getPairs(distance=8.0, sel2='element O N')
    pairwiseTermsRes = IElister.getInteractions(pairs, mode='perSet1Residue')
    IElister.saveCSV('recResONPWT.csv', pairwiseTermsRes)
    for v in pairwiseTermsRes[:10]:
        print(v)
