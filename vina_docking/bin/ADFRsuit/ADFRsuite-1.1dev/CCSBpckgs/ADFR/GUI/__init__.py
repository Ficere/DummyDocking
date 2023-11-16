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
# $Header: /mnt/raid/services/cvs/ADFR/GUI/__init__.py,v 1.2 2016/12/07 00:38:32 sanner Exp $
#
# $Id: __init__.py,v 1.2 2016/12/07 00:38:32 sanner Exp $
#
import os
from glob import glob
from PySide import QtCore, QtGui

import ADFR
import ADFR.GUI
ICONPATH = os.path.join(ADFR.__path__[0], 'GUI', 'Icons')

_AGFRGUI_debug = False # set to True to avoid @waiting to mask exceptions

def waiting_effects(function):
    def new_function(*args, **kwargs):
        #QtGui.QApplication.setOverrideCursor(QtGui.QCursor(QtCore.Qt.IBeamCursor))
        if _AGFRGUI_debug:
            return function(*args, **kwargs)
        else:
            QtGui.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
            try:
                return function(*args, **kwargs)
            except Exception as e:
                print "error in function" , function
                print("Error {}".format(e.args[0]))
                raise e
            finally:
                QtGui.QApplication.restoreOverrideCursor()
    return new_function

def makeGridGUI(app):
    import sys, PySide
    from PySide import QtGui

    from ADFR.GUI.gridGUI import GridGUI
    from ADFR.utils.optParser import AGFRParser
    qtpath = None
    if sys.platform == "win32":
        qtpath = os.path.join(os.path.dirname(PySide.__file__), 'plugins')
    else:
        root = os.getenv('ADS_ROOT')
        if root:
            qtpath = os.path.join(root, 'plugins')
        if sys.platform.startswith('linux'):
            app.setStyle('Cleanlooks')
            app.setFont (QtGui.QFont ("Times", 10, QtGui.QFont.Normal))
        if sys.platform.startswith('darwin') :
            app.setStyle('Cleanlooks')
            w = app 
            palette = QtGui.QPalette(w.palette())
            palette.setColor(palette.Window,  "#EFF0F1") #"#ECF0F1")
            palette.setColor(palette.Base, QtCore.Qt.white)
            palette.setColor(palette.Button, "#F2F3F4")
            palette.setColor(palette.AlternateBase, QtCore.Qt.white)
            w.setPalette(palette)

    if qtpath:
        app.addLibraryPath(qtpath)
    else:
        print ("Warning AGFRgui.py: Failed to load PySide plugins library")
    parser = AGFRParser(prog='agfrgui')
    options = vars(parser.parse_args())
    _debug = options.pop('debug')
    if _debug>0:
        ADFR.GUI._AGFRGUI_debug = True
        #options['pyShell'] = True
    folder = None
    receptor = options.get('receptorFile', None)
    if receptor:
        if os.path.isdir(receptor):
            folder = receptor
    if folder:
        filenames = glob(os.path.join(folder, "*.pdb"))
        filenames.extend(glob(os.path.join(folder, "*.ent.gz")))
        filenames.extend(glob(os.path.join(folder, "*.pdbqt")))
        if len(filenames)==0:
            print 'ERROR: no suitable files found in %s'% folder
        options['receptorFile'] = filenames[0]
        #print [os.path.split(name)[1] for name in filenames]
        widget = GridGUI(files=filenames)
    else:
        widget = GridGUI()
        
    

    widget.show()
    widget.configure(options)
    return widget
