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

#############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner and TSRI 2015
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/ADFR/bin/runAGFR.py,v 1.5 2017/05/05 00:21:02 annao Exp $
#
# $Id: runAGFR.py,v 1.5 2017/05/05 00:21:02 annao Exp $
#
#
# Usage: pythonsh runADFR.py lig_random.pdbqt mapFolder -o file.log -ref lig_xtal.pdbqt
#
import sys, os

from ADFR.utils.runAGFR import runAGFR
from ADFR.utils.optParser import AGFRParser
from glob import glob

parser = AGFRParser()
cnfFile = None
# check if argument list contains configuration file
args = sys.argv[1:]
if len(args):
    for i, arg in enumerate(args):
        if arg in ['--config', '-F']:
            cnfFile = args[i+1]
            break
if cnfFile:
    # read arguments from the file, pass them to the parser
    from ADFR.utils.optParser import readConfigFile
    cnfArgs = readConfigFile(cnfFile)
    kw =  vars(parser.parse_args(cnfArgs))
else:
    # process command line arguments
    kw =  vars(parser.parse_args())

if kw['testBatch']:
    class DontPrint(object):
        def write(*args): pass

    dp = DontPrint()
    out = sys.stdout
    err = sys.stderr
    sys.stdout = dp
    sys.stderr = dp
    try:
        runner = runAGFR()
        runner(**kw)
    except:
        sys.stdout = out
        sys.stderr = err
        print 'failed for', kw['receptorFile']
        #sys.exit(1)
else:
    folder = None
    onException = kw.pop('onException')
    receptor = kw.get('receptorFile', None)
    if receptor:
        # check if folder:
        if os.path.isdir(receptor):
            folder = receptor
        if folder:
            filenames = glob(os.path.join(folder, "*.pdb"))
            filenames.extend(glob(os.path.join(folder, "*.ent.gz")))
            filenames.extend(glob(os.path.join(folder, "*.pdbqt")))
            if len(filenames)==0:
                print 'ERROR: no suitable files found in %s'% folder
            filenames.sort()
        else:
            filenames = [receptor]
        for filename in filenames:
            try:
                print 'processing:',os.path.split(filename)[1] 
                kw['receptorFile'] = filename
                runner = runAGFR()
                runner(**kw)
            except Exception, e:
                if onException=='raise':
                    raise
                elif onException=='stop':
                    print e
                    import sys
                    sys.exit(1)
                elif onException=='report':
                    print e
