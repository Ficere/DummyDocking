##############################################################################
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
# $Header: /mnt/raid/services/cvs/ADFR/GUI/gridGUI.py,v 1.110.2.17 2018/03/29 21:27:38 annao Exp $
#
# $Id: gridGUI.py,v 1.110.2.17 2018/03/29 21:27:38 annao Exp $
#
import os, sys, weakref, tempfile, shutil, numpy, string, threading
from math import sqrt
from random import random
from time import time, sleep
from PySide import QtCore, QtGui

from MolKit2 import getCrystalMatesMatrices, Read
from MolKit2.selection import Selection, SelectionSet
from MolKit2.molecule import MultiMolecule

from ADFR import checkLigandFile
from ADFR.utils.runAGFR import errorCodes
from ADFR.utils.maps import flexResStr2flexRes
import ADFR.GUI
from ADFR.GUI import ICONPATH, waiting_effects
from ADFR.GUI.pairwiseScore import PairwiseScoreTable, molToAtomSetStatic
from ADFRcc.adfr import AtomSet, PairwiseScorer

from mglutil.util.packageFilePath import getResourceFolderWithVersion,findFilePath
from mglutil.util.callback import CallbackFunction
from mglutil.preferences import UserPreference, PreferenceItem

from PmvApp.Pmv import MolApp
from PmvApp.Pmv import RefreshDisplayEvent
from prody import Atom

from AutoSite.fillBuriedness import Buriedness
from AutoSite.utils.clusterTPoints import DensityClustering
from AutoSite.scoreClusters import scoreClusters
#from AutoSite.ASfeaturepoints import featurePts

from DejaVu2.Cylinders import Cylinders
from DejaVu2.Spheres import Spheres

LeftMargin = rightMargin = spacing = 8
topMargin = bottomMargin = 0


use_ipython_shell=False

try:
    from IPython.qt.console.rich_ipython_widget import RichIPythonWidget
    from IPython.qt.inprocess import QtInProcessKernelManager
    from IPython.lib import guisupport
    use_ipython_shell=True
except ImportError:
    print "Warning: failed to import IPython.qt module. Cannot use IPython shell widget,\nusing PyShell instead."

def ribbonArrowGeom(arc, normal, width, headLength=2):
    headLength = 2
    # make verices for body of arrow
    w = width*0.5
    normal = numpy.array(normal)
    arc =  numpy.array(arc)
    left = arc[:-headLength] - w*normal
    right = arc[:-headLength] + w*normal

    # add 2 vertices for arrow head base
    p1 = arc[-headLength-1] - width*normal
    p2 = arc[-headLength-1] + width*normal

    # calculate arrow tip
    p3 = arc[-1] # tip of the arrow
    vertices = numpy.array(left.tolist() + right.tolist() + [p1,p2,p3])

    # build the faces
    faces = []
    n = len(left)
    for i in range(n-1):
        faces.append( (i, i+n, i+1) )
        faces.append( (i+1, i+n, i+n+1) )

    # add arrow head with 3 triangles ((1,2,5) (2,3,5) and (3,4,5)
    #               5
    #           1 2   3 4
    nn = 2*len(left)
    faces.append((n-1, nn+2, nn))
    faces.append((n-1, 2*n-1, nn+2))
    faces.append((2*n-1, nn+1, nn+2))

    return vertices, faces

class MyQLineEdit(QtGui.QLineEdit):
    recall = QtCore.Signal(str)

    def keyPressEvent(self, e):
        key = e.key()
        if key == QtCore.Qt.Key_Up:
            self.recall.emit('up')
        elif key == QtCore.Qt.Key_Down:
            self.recall.emit('down')
        else:
           QtGui.QLineEdit.keyPressEvent(self, e)

class LogBox(QtGui.QWidget):

    def __init__(self, title, widget, parent=None):
        super(LogBox, self).__init__(parent)
        self.setWindowTitle(title)
        layout = QtGui.QVBoxLayout()
        self.text = QtGui.QTextEdit()
        self.text.setReadOnly(True)
        layout.addWidget(self.text)
        self.setLayout(layout)
        rel_pos = widget.pos()
        pos = widget.mapToGlobal(rel_pos)
        self.move(pos.x()+20,pos.y()+20)
        self.resize(400,600)

class DockingTestParams( QtGui.QWidget):

    def __init__(self, parent=None):
        super(DockingTestParams, self).__init__(parent=parent)
        l = QtGui.QFormLayout(self)

        w = self.maxEvalsWidget = QtGui.QSpinBox()
        w.setMinimum(1)
        w.setMaximum(99999999)
        w.setValue(1000)
        l.addRow('max evaluations (thousands):', w)

        w = self.maxGensWidget = QtGui.QSpinBox()
        w.setMinimum(1)
        w.setMaximum(99999999)
        w.setValue(50)
        l.addRow('max generations:', w)

        w = self.maxNoImprovWidget = QtGui.QSpinBox()
        w.setMinimum(1)
        w.setMaximum(100)
        w.setValue(5)
        l.addRow('max no improvment:', w)

        w = self.seedWidget = QtGui.QSpinBox()
        w.setMaximum(999999999)
        w.setValue(random()*10000)
        l.addRow('random seed:', w)

        w = self.popSizeWidget = QtGui.QSpinBox()
        w.setMinimum(1)
        w.setMaximum(1000)
        w.setValue(150)
        l.addRow('population size:', w)


class DockingTestParamsDialog(QtGui.QDialog):

    def __init__(self, title='No Name', parent=None):
        super(DockingTestParamsDialog, self).__init__(parent)
        self.setWindowTitle(title)

        # place the widget close to where we clicked
        rel_pos = QtGui.QCursor.pos()
        pos = self.mapToGlobal(rel_pos)
        self.move(pos.x()+20, pos.y()+15)
        self.paramsWidget = DockingTestParams()
        self.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)

        self.button = QtGui.QPushButton("Run")
        self.button.clicked.connect(self.accept)

        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.paramsWidget)
        layout.addWidget(self.button)
        # Set dialog layout
        self.setLayout(layout)

    
class BoxParameters(QtGui.QGroupBox):

    def __init__(self, app, parent=None):
        super(BoxParameters, self).__init__(parent)

        self.setTitle('box parameters')
        self.app = weakref.ref(app)
        self.cmdHistory = []
        self.historyPointer = 0
        #from DejaVu2.Box import NiceBox
        #self.box = NiceBox('psudoBox')
        #self.box.addToViewer(app.viewer)
        box = app.boxGeom
        self.data = {}

        gridLayout = QtGui.QGridLayout()
        gridLayout.setContentsMargins(LeftMargin,topMargin,rightMargin,bottomMargin)
        gridLayout.setSpacing(spacing)

        gridLayout.addWidget(QtGui.QLabel("x:"), 0, 1, 1, 1, QtCore.Qt.AlignCenter)
        gridLayout.addWidget(QtGui.QLabel("y:"), 0, 2, 1, 1, QtCore.Qt.AlignCenter)
        gridLayout.addWidget(QtGui.QLabel("z:"), 0, 3, 1, 1, QtCore.Qt.AlignCenter)
        gridLayout.addWidget(QtGui.QLabel("center:"), 1, 0, 1, 1, QtCore.Qt.AlignRight)
        gridLayout.addWidget(QtGui.QLabel("size:"), 2, 0, 1, 1, QtCore.Qt.AlignRight)
        gridLayout.addWidget(QtGui.QLabel("spacing:"), 3, 0, 1, 1, QtCore.Qt.AlignRight)
        gridLayout.addWidget(QtGui.QLabel("smoothing:"), 3, 2, 1, 1, QtCore.Qt.AlignRight)
        gridLayout.addWidget(QtGui.QLabel("cmd:"), 4, 0, 1, 1, QtCore.Qt.AlignRight)

        self._initialized = False
        self.centerSpinBoxes = []
        values = box.center
        for i in range(3):
            widget = QtGui.QDoubleSpinBox()
            widget.valueChanged.connect(self.centerChanged)
            widget.setDecimals(3)
            widget.setMinimum(-9999.)
            widget.setMaximum(9999.)
            widget.setValue(values[i])
            self.centerSpinBoxes.append(widget)
            gridLayout.addWidget(widget, 1, i+1, 1, 1)

        self.sizeSpinBoxes = []
        values = box.sides
        for i in range(3):
            widget = QtGui.QDoubleSpinBox()
            widget.valueChanged.connect(self.sizeChanged)
            widget.setDecimals(3)
            widget.setMinimum(app._spacing)
            widget.setMaximum(9999.)
            widget.setValue(values[i])
            self.sizeSpinBoxes.append(widget)
            gridLayout.addWidget(widget, 2, i+1, 1, 1)

        self.spacingWidget = widget = QtGui.QDoubleSpinBox()
        widget.setDecimals(3)
        widget.setValue(app._spacing)
        widget.setMinimum(0.00001)
        widget.setSingleStep(0.025)
        widget.valueChanged.connect(self.spacingChanged)
        gridLayout.addWidget(widget, 3, 1, 1, 1)
 
        self.smoothWidget = widget = QtGui.QDoubleSpinBox()
        widget.setDecimals(3)
        widget.setValue(app._smooth)
        widget.setMinimum(0.00001)
        widget.valueChanged.connect(self.smoothChanged)
        widget.setSingleStep(0.1)
        gridLayout.addWidget(widget, 3, 3, 1, 1)
        
        self.cmdWidget = MyQLineEdit()
        self.cmdWidget.recall.connect(self.handleRecall)
        gridLayout.addWidget(self.cmdWidget, 4, 1, 1, 3)
        self.cmdWidget.returnPressed.connect(self.executeCmd)

        self.setLayout(gridLayout)
        self._initialized = True

    def sizeChanged(self, value=None):
        if not self._initialized: return
        app = self.app()
        #padding = app.paramsWidget.gridPaddingWidget.value()
        self.data["lengths"] = [x.value() for x in self.sizeSpinBoxes]
        l1, l2, l3 = self.data['lengths']
        #app.agfr.setBox(["user", x, y, z, l1, l2, l3], padding, app._spacing)
        #padding = app.agfr.padding
        #app._baseSize[:] = app.agfr.boxLengths - 2*padding
        #app._baseSize[:] = [x-2*padding for x in lengths]
        #app.boxGeom.setSides(*app.agfr.boxLengths)
        #app.paramsWidget.gridPaddingWidget.setValue(padding)
        #app.onBoxChange()
        app.boxGeom.setSides(l1,l2,l3)

    def centerChanged(self, value=None):
        if not self._initialized: return
        app = self.app()
        center = self.data['center'] = [x.value() for x in self.centerSpinBoxes]
        #app.agfr.setBoxCenter(center)
        #app.boxGeom.setCenter( *center)
        #app.onBoxChange()
        app.boxGeom.setCenter(*center)

        
    def spacingChanged(self, value):
        #app = self.app()
        #app._spacing = value
        #app.agfr.setSpacing(value)
        #app.onBoxChange()
        self.data['spacing'] = value

    def smoothChanged(self, value):
        #app = self.app()
        #app._smooth = value
        self.data['smooth'] = value

    def applyChange(self):
        if not len(self.data): return
        app = self.app()
        center = self.data.get("center", None)
        if center is not None:
            app.agfr.setBoxCenter(center)
            app.boxGeom.setCenter( *center)
        spacing = self.data.get("spacing", None)
        if spacing is not None:
            app._spacing = value
            app.agfr.setSpacing(value)
        smooth = self.data.get("smooth", None)
        if smooth is not None:
            app._smooth = value
        lengths = self.data.get("lengths", None)
        if lengths is not None:
            x, y, z = app.agfr.boxCenter
            l1, l2, l3 = lengths
            padding = app.agfr.padding
            app.agfr.setBox(["user", x, y, z, l1, l2, l3], padding, app._spacing)
            app._baseSize[:] = app.agfr.boxLengths - 2*padding
            #app._baseSize[:] = [x-2*padding for x in lengths]
            app.boxGeom.setSides(*app.agfr.boxLengths)
        padding = app.agfr.padding
        app.paramsWidget.gridPaddingWidget.setValue(padding)
        app.onBoxChange()
        self.data = {}

    def _cancel(self):
        if len(self.data):
            app = self.app()
            # set the box to the original data
            app.boxGeom.setCenter(*app.agfr.boxCenter)
            app.boxGeom.setSides(*app.agfr.boxLengths)
    ##
    ## methods for cmd entry
    def handleRecall(self, what):
        if what=='up':
            if self.historyPointer > 0:
                self.historyPointer -= 1
                self.cmdWidget.clear()
                self.cmdWidget.insert(self.cmdHistory[self.historyPointer])
        elif what=='down':
            if self.historyPointer < len(self.cmdHistory):
                self.historyPointer += 1
                self.cmdWidget.clear()
                if self.historyPointer < len(self.cmdHistory):
                    self.cmdWidget.insert(self.cmdHistory[self.historyPointer])

    def _getIntOtFloat(self, s):
        try:
            value = int(s)
            return 'int', value
        except ValueError:
            try:
                value = float(s) 
                return 'float', value
            except ValueError:
                return None, None

    def executeCmd(self):
        cmd = self.cmdWidget.text()
        self.cmdHistory.append(cmd)
        self.historyPointer += 1
        self.cmdWidget.clear()
        cmds = cmd.split(';')
        # move size command to the front for performance reasons
        otherCmds = []
        splitCmd = None
        for cmd in cmds:
            c = cmd.split()[0].lower()
            if c=='size'[:len(c)] or c=='dimensions'[:len(c)]:
                splitCmd = cmd
            else:
               otherCmds.append(cmd)
        if splitCmd:
            otherCmds.insert(0, splitCmd)
        cmds = otherCmds
        
        app = self.app()
        for cmd in cmds:
            words = cmd.split()
            if words[0].lower()=='center'[:len(words[0])]:
                if len(words)==4:
                    size = [float(x) for x in words[1:4]]
                    self._initialized = False # to avoid 2 callbacks
                    for i in range(3):
                        if i==2:
                            self._initialized = True
                        self.centerSpinBoxes[i].setValue(size[i])
                else:
                    msgBox = QtGui.QMessageBox(self)
                    msgBox.setText("ERROR: Bad syntax. expected center x y z")
                    msgBox.exec_()

            elif words[0].lower()=='size'[:len(words[0])]:
                size = None
                if len(words)==2:
                    vtype, size = self._getIntOtFloat(words[1])
                    size = [size, size, size]                
                elif len(words)==4:
                    s = self.spacingWidget.value()
                    sx = float(words[1])
                    sy = float(words[2])
                    sz = float(words[3])
                    nx = round(sx/s)
                    ny = round(sy/s)
                    nz = round(sz/s)
                    size = [nx*s, ny*s, nz*s]
                else:
                    msgBox = QtGui.QMessageBox(self)
                    msgBox.setText("ERROR: bad syntax. expected size x or size x y z")
                    msgBox.exec_()

                if size is not None:
                    
                    self._initialized = False
                    for i in range(3):
                        if i==2:
                            self._initialized = True
                        self.sizeSpinBoxes[i].setValue(size[i])

            elif words[0].lower()=='dimensions'[:len(words[0])]:
                size = None
                if len(words)==2:
                    vtype, size = self._getIntOtFloat(words[1])
                    size = size*self.spacingWidget.value() 
                    size = [size, size, size]
                elif len(words)==4:
                    vtype, size = self._getIntOtFloat(words[1])
                    nx = size
                    ny = int(words[2])
                    nz = int(words[3])
                    if 2*(nx/2)-nx==0: #event number
                        nx += 1
                    if 2*(ny/2)-ny==0: #event number
                        ny += 1
                    if 2*(nz/2)-nz==0: #event number
                        nz += 1
                    s = self.spacingWidget.value()
                    size = [nx*s, ny*s, nz*s]
                else:
                    msgBox = QtGui.QMessageBox(self)
                    msgBox.setText("ERROR: bad syntax. expected nx, ny nz as int for nb grid points")
                    msgBox.exec_()

                if size is not None:
                    self._initialized = False
                    for i in range(3):
                        if i==2:
                            self._initialized = True
                        self.sizeSpinBoxes[i].setValue(size[i])
            else:
                msgBox = QtGui.QMessageBox(self)
                msgBox.setText("ERROR: %s is not a valid keyword. Use center, size, or dimensions"%words[0])
                msgBox.exec_()

class BoxParametersDialog(QtGui.QDialog):
    closedSignal = QtCore.Signal()

    def __init__(self, app,  title='No Name',  parent=None):
        super(BoxParametersDialog, self).__init__(parent)
        self.setWindowTitle(title)

        # place the widget close to where we clicked
        rel_pos = QtGui.QCursor.pos()
        pos = self.mapToGlobal(rel_pos)
        self.move(pos.x()+20, pos.y()+15)

        self.boxParamsWidget = BoxParameters(app, parent)
        self.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok,
            QtCore.Qt.Horizontal, self)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        #self.buttonBox.button(self.buttons.Cancel).clicked.connect(self.cancel_cb)
        #self.buttonBox.button(self.buttons.Ok).clicked.connect(self.ok_cb)
        #self.buttonBox.button(self.buttons.Ok).setDefault(False)

        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.boxParamsWidget)
        layout.addWidget(self.buttonBox)
        # Set dialog layout
        self.setLayout(layout)
        #self.setWindowFlags( self.windowFlags() & ~QtCore.Qt.WindowCloseButtonHint)
        # create a Box

    def reject(self):
        QtGui.QDialog.reject(self)
        self.boxParamsWidget._cancel()
        self.closeEvent()
        
    def closeEvent(self, evnt=None):
        self.closedSignal.emit()

    def accept(self):
        QtGui.QDialog.accept(self)
        self.boxParamsWidget.applyChange()

    def keyPressEvent(self, evt):
        #this is to prevent the dialog's closing when 'return' is pressed in the LineEdit widget.
        if evt.key() == QtCore.Qt.Key_Enter or evt.key() == QtCore.Qt.Key_Return:
            return
        QtGui.QDialog.keyPressEvent(self, evt)

class CheckParamsInBoxDialog(QtGui.QDialog):
    finished = QtCore.Signal(int, bool)
    
    def __init__(self, msgtxt, buttons=["OK"], auto="do not ask me again", title="", parent=None):
        super(CheckParamsInBoxDialog, self).__init__(parent)
        self.setWindowTitle(title)
        self.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
        layout = QtGui.QVBoxLayout()
        label = QtGui.QLabel("WARNING: %s" % msgtxt)
        buttonBox = QtGui.QDialogButtonBox(QtCore.Qt.Horizontal)
        self.buttons = []
        for i , btxt in enumerate(buttons):
            bn = QtGui.QPushButton(self.tr(btxt))
            bn.clicked.connect(CallbackFunction(self.button_cb, i))
            buttonBox.addButton(bn, QtGui.QDialogButtonBox.ActionRole)
            self.buttons.append(bn)
        self.checkbox = QtGui.QCheckBox()
        self.checkbox.setText(auto)
        layout.addWidget(label)
        layout.addWidget(buttonBox)
        layout.addWidget(self.checkbox)
        self.setLayout(layout)

    def button_cb(self, val):
        self.finished.emit(val, self.checkbox.isChecked())
        QtGui.QDialog.accept(self)
        
    def reject(self):
        self.finished.emit(-1, self.checkbox.isChecked())
        QtGui.QDialog.reject(self)

    def exec_(self, *args, **kwargs):
        """
        Override the exec_ method to return the value of the checkbox
        """
        #print "EXEC", args, kwargs
        return QtGui.QMessageBox.exec_(self, *args, **kwargs), self.checkbox.isChecked()
    
    
class SettingsDialog(QtGui.QDialog):
    def __init__(self, app, title="User settings", parent=None):
        super(SettingsDialog, self).__init__(parent)
        self.setWindowTitle(title)
        # place the widget close to where we clicked
        rel_pos = QtGui.QCursor.pos()
        pos = self.mapToGlobal(rel_pos)
        self.move(pos.x()+20, pos.y()+15)
        self.buttonBox = QtGui.QDialogButtonBox()
        self.buttonBox.accepted.connect(self.accept)

        layout = QtGui.QVBoxLayout()
        self.groupBox1 = QtGui.QGroupBox(self)
        self.groupBox1.setTitle(self.tr("Flexible Residues Out Of Box"))
        layout1 = QtGui.QVBoxLayout(self.groupBox1)
        
        layout.addWidget(self.buttonBox)
        # Set dialog layout
        self.setLayout(layout)

class readLogThread(QtCore.QThread):
    progress = QtCore.Signal(str) # create a custom signal we can subscribe
                                  # to to emit update commands
    progressBar = QtCore.Signal(float)

    errorSignal = QtCore.Signal(str)

    def __init__(self, process, logFile, app, postProcess=None, parent=None):
        super(readLogThread,self).__init__(parent)
        self.app = app
        self.postProcess = postProcess
        self.exiting = False
        self.process = process
        self.logFile = logFile
        self.nbLogLines = 0

    def run(self):
        emitThreadDone = True
        emitMsg = 'readLogThread THREAD DONE'
        while True:
            #print 'IN RUN', self.process.poll() is None
            if self.process.poll() == None:
                self.msleep(500)
                #print 'reading STDOUT'
                #for line in iter(self.process.stdout.readline, ""):
                #    self.progress.emit('out: '+line[:-1])
                #print 'reading STDERR'
                #for line in iter(self.process.stderr.readline, ""):
                #    self.progress.emit('err: '+line[:-1])
                #self.msleep(1)
                #print 'reading LOG'
                try:
                    f = open(self.logFile)
                    #self.progressBar.emit(1)
                except IOError:
                    continue
                lines = f.readlines()
                f.close()
                if len(lines)>self.nbLogLines:
                    #print 'NBL', self.nbLogLines
                    text = ''.join(lines[self.nbLogLines:])
                    line = lines[-1]
                    #if len(line)>83 and line[26]=='%':
                    #    print 'FUGU', float(line.split()[2][:-1])/100.
                    if  len(line)>26 and line[26]=='%':
                        self.progressBar.emit( float(line.split()[2][:-1]) )
                    #else:
                    #    self.progress.emit( 0. )
                    #self.progress.emit(text)
                    self.nbLogLines = len(lines)
            else:
                output = self.process.stderr.readlines()
                _errors = []
                _warnings = []
                if len(output):
                    for err in output:
                        if len(err) > 1:
                            if err.find("WARNING") > 0:
                                _warnings.append(err)
                            else:
                                _errors.append(err)
                    if len(_warnings):
                        msg = "AutoGrid folder %s\nWARNING message:\n%s"%(
                            self.app._gc.folder, ''.join(_warnings))
                        print msg
                    if len(_errors):
                        errmsg = ''.join(_errors)
                        if errmsg.find("no closestH atom was found") >= 0:
                            errmsg = "It seems that hydrogen atoms are missing in the receptor."
                        msg = "AutoGrid failed to run in folder %s\nError message:\n%s"%(
                            self.app._gc.folder, errmsg)
                        self.errorSignal.emit(msg)
                        emitThreadDone = True
                        emitMsg = "AutoGrid failed"
                        print msg
                if not len(_errors):
                    self.progressBar.emit(100.0)
                    if self.postProcess:
                        self.progressBar.emit(0.0)
                        emitThreadDone = False
                        self.runPostProcess()
                break
        if emitThreadDone:
            #self.progress.emit('readLogThread THREAD DONE')
            self.progress.emit(emitMsg)
            
    def runPostProcess(self):
        try:
            self.postProcess()
        except Exception as e:
            print("Error {}".format(e.args[0]))
            raise e
        #print "postprocess thread done"
        self.progressBar.emit(100.0)
        self.progress.emit('readLogThread THREAD DONE')

class progressLogFile(QtCore.QThread):
    progressBar = QtCore.Signal(float)
    progress = QtCore.Signal(str)
    def __init__(self, filename, maxlines, app, parent=None, doneMsg=None):
        super(progressLogFile ,self).__init__(parent)
        self.logFile = filename
        self.maxlines = maxlines
        self._stop = False
        self.app = app
        self.doneMsg=doneMsg
        #print "Log THREAD maxlines:", maxlines

    def run(self):
        fileobj = None
        while True:
            if self._stop:
                #print "progressLogFile stop."
                #self.progressBar.emit(100.)
                if fileobj:
                    try:
                        fileobj.close()
                    except:
                        pass
                #self.progressBar.emit('0')
                return
            try:
                with open(self.logFile)as fileobj:
                    lines = fileobj.readlines()
            except IOError:
                self.msleep(300)
                continue
            #lines = f.readlines()
            #f.close()
            nlines = len(lines)
            #print "AAAA lines:", nlines, lines
            if nlines < self.maxlines:
                val = nlines*100./self.maxlines
                if val == 0.0 : val = 5. # just starting the progress line
                #print "progress:" , val
                self.progressBar.emit(val)
            else:
                self.progressBar.emit(100.)
                if self.doneMsg:
                    self.progress.emit(self.doneMsg)
                return
            self.msleep(300)

    def stop(self):
        self._stop = True


class MapsFolderDialog(QtGui.QDialog):

    def __init__(self, name, folder=None, parent=None):
        super(MapsFolderDialog, self).__init__(parent)
        self.setWindowTitle("target file name")

        labelName = QtGui.QLabel("name:")
        self.nameEntry = QtGui.QLineEdit()
        self.nameEntry.setText(name)
        
        label = QtGui.QLabel('Enter file name for saving maps')
        browseButton = QtGui.QPushButton('browse')
        browseButton.clicked.connect(self.browseForFolder)
        
        self.folderPath = QtGui.QLineEdit()
        self.folderPath.returnPressed.connect(self.setPath)
        self.folderPath.setMinimumWidth(400)
        if folder is None:
            self.folderPath.setText(os.getcwd())
        else:
            self.folderPath.setText(folder)

        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok |
                                                QtGui.QDialogButtonBox.Cancel)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        #self.buttonBox.buttons()[0].setDisabled(True)
        
        layouth1 = QtGui.QHBoxLayout()
        layouth1.setContentsMargins(LeftMargin,topMargin,rightMargin,bottomMargin)
        layouth1.setSpacing(spacing)
        layouth1.addWidget(labelName)
        layouth1.addWidget(self.nameEntry)
        layouth = QtGui.QHBoxLayout()
        layouth.setContentsMargins(LeftMargin,topMargin,rightMargin,bottomMargin)
        layouth.setSpacing(spacing)
        layouth.addWidget(browseButton)
        layouth.addWidget(self.folderPath)
        
        layout = QtGui.QVBoxLayout()
        layout.setContentsMargins(LeftMargin,topMargin,rightMargin,bottomMargin)
        layout.setSpacing(spacing)
        layout.addWidget(label)
        layout.addLayout(layouth1)
        layout.addLayout(layouth)
        layout.addWidget(self.buttonBox)
        self.setLayout(layout)

    def browseForFolder(self):
        dialog = QtGui.QFileDialog(self)
        dialog.setFileMode(QtGui.QFileDialog.Directory)
        dialog.setOption(QtGui.QFileDialog.ShowDirsOnly)
        if dialog.exec_():
            fileNames = dialog.selectedFiles()
            self.folderPath.setText(fileNames[0])
            self.setPath()

    def setPath(self):
        folder = self.folderPath.text()
        if not os.path.exists(folder):
            msgBox = QtGui.QMessageBox(self)
            msgBox.setText("ERROR: Folder %s does not exist"%folder)
            msgBox.exec_()
            self.folderPath.setText('')
            return
        if not os.path.isdir(folder):
            msgBox = QtGui.QMessageBox(self)
            msgBox.setText("ERROR: %s is not a folder"%folder)
            msgBox.exec_()
            self.folderPath.setText('')
            return
        self.buttonBox.buttons()[0].setDisabled(False)
        
class AtomTypeSelector(QtGui.QDialog):

    def __init__(self, atomList, parent=None):
        super(AtomTypeSelector, self).__init__(parent)
        self.setWindowTitle("select atom types")

        self.allButtonOn = QtGui.QCheckBox("select all")
        self.allButtonOn.clicked.connect(self.selectAll)

        self.allButtonOff = QtGui.QCheckBox("deselect all")
        self.allButtonOff.clicked.connect(self.deselectAll)

        self.atypeList = QtGui.QListWidget()
        self.items = []
        self.atomList = atomList
        for atype, checked in atomList.items():
            item = QtGui.QListWidgetItem(atype)
            self.items.append(item)
            if checked:
                item.setCheckState(QtCore.Qt.Checked)
            else:
                item.setCheckState(QtCore.Qt.Unchecked)
            self.atypeList.addItem(item)
                
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok|QtGui.QDialogButtonBox.Cancel)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.cancel_cb)

        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.allButtonOn)
        layout.addWidget(self.allButtonOff)
        layout.addWidget(self.atypeList)
        layout.addWidget(self.buttonBox)
        # Set dialog layout
        self.setLayout(layout)
        self.atypeList.itemClicked.connect(self.itemClicked_cb)
        self.okbutton = self.buttonBox.button(QtGui.QDialogButtonBox.Ok)

    def selectAll(self):
        for item in self.items:
            item.setCheckState(QtCore.Qt.Checked)
        self.okbutton.setDisabled(False)
         
    def deselectAll(self):
        for item in self.items:
            item.setCheckState(QtCore.Qt.Unchecked)
        self.okbutton.setDisabled(True)

    def itemClicked_cb(self, item):
        nitems = len(self.items)
        count = 0
        for item in self.items:
            if item.checkState()==QtCore.Qt.Unchecked:
                count+=1
        if count == nitems: #all unchecked, disable "OK" button
            self.okbutton.setDisabled(True)
        else:
            self.okbutton.setDisabled(False)

    def getTypesString(self):
        s = ''
        for item in self.items:
            if item.checkState()==QtCore.Qt.Checked:
                s += '%s '%item.text()
        return s[:-1]
    
    def cancel_cb(self):
        n = 0
        for atype, checked in self.atomList.items():
            item = self.items[n]
            if checked:
                item.setCheckState(QtCore.Qt.Checked)
            else:
                item.setCheckState(QtCore.Qt.Unchecked)
            n+=1
        self.reject()

## attempt to change sort of table to keep ALL at the top
## and have small clusters at the bottom, but it does not seem to
## work.
class TableWidgetItem(QtGui.QTableWidgetItem):
    def __lt__(self, other):
        #import pdb; pdb.set_trace()
        str1 = self.data(QtCore.Qt.DisplayRole)
        if str1.startswith('NA1') or str1.startswith('all'):
            n1 = 9999999999.
        elif str1.startswith('NA2') or str1.startswith('small'):
            n1 = -1.
        else:
            n1 = float(str1)
                
        str2 = other.data(QtCore.Qt.DisplayRole)
        if str2.startswith('NA1') or str2.startswith('all'):
            n2 = 9999999999.
        elif str2.startswith('NA2') or str2.startswith('small'):
            n2 = -1.
        else:
            n2 = float(str2)
        #print str1, str2, n1, n2
        return (n1 < n2)

        
class ADFRGridMapParametersWidget(QtGui.QWidget):
    
    def __init__(self, app, parent=None):

        super(ADFRGridMapParametersWidget, self).__init__(parent)
        self._tpok = False
        self._gridok = True
        self._boxok = False
        self._covbondok = False
        self.app = weakref.ref(app)
        
        ## create group for input
        ##
        ## self.groupBox1 = QtGui.QGroupBox(self)
        ## self.groupBox1.setTitle(self.tr("load"))
        ## layout = QtGui.QHBoxLayout()
        ## self.groupBox1.setLayout(layout)
        ## w = QtGui.QPushButton(QtGui.QIcon(os.path.join(ICONPATH, 'rec1.png')), '')
        ## w.setIconSize(QtCore.QSize(48, 48))
        ## layout.addWidget(w)
        ## self.loadRecButton = w

        ## w = QtGui.QPushButton(QtGui.QIcon(os.path.join(ICONPATH, 'target.png')), '')
        ## w.setIconSize(QtCore.QSize(48, 48))
        ## layout.addWidget(w)
        ## self.loadMapsButton = w
        
        ## w = QtGui.QPushButton(QtGui.QIcon(os.path.join(ICONPATH, 'ligandNoArrow.png')), '')
        ## w.setIconSize(QtCore.QSize(48, 48))
        ## layout.addWidget(w)
        ## self.loadLigButton = w
        
        ## w = QtGui.QPushButton(QtGui.QIcon(os.path.join(ICONPATH, 'folder.png')), '')
        ## w.setIconSize(QtCore.QSize(48, 48))
        ## layout.addWidget(w)
        ## self.setFolderButton = w
        
        ## w = QtGui.QPushButton(QtGui.QIcon(os.path.join(ICONPATH, 'www.png')), '')
        ## w.setIconSize(QtCore.QSize(48, 48))
        ## layout.addWidget(w)
        ## self.fetchPDBButton = w
        
        ## self.groupBox1 = QtGui.QGroupBox(self)
        ## self.groupBox1.setTitle(self.tr("receptor"))
        ## #self.groupBox1.setContentsMargins(LeftMargin,6,rightMargin,6) # causes button to overlap with top of the box
        ## layout = QtGui.QHBoxLayout()
        ## layout.setContentsMargins(LeftMargin,topMargin,rightMargin,bottomMargin)
        ## #layout.setSpacing(spacing)
        ## self.groupBox1.setLayout(layout)
        ## w = QtGui.QPushButton("Open ...")
        ## layout.addWidget(w)
        ## self.loadRecButton = w
        
        ## w = QtGui.QPushButton("target file(.trg) ...")
        ## layout.addWidget(w)
        ## self.loadMapsButton = w
        
        ## ## create group for ligand
        ## ##
        ## self.groupBox2 = QtGui.QGroupBox(self)
        ## self.groupBox2.setTitle(self.tr("[ligand]"))
        ## #self.groupBox2.setContentsMargins(LeftMargin,6,rightMargin,6)
        ## layout = QtGui.QHBoxLayout()
        ## layout.setContentsMargins(LeftMargin,topMargin,rightMargin,bottomMargin)
        ## #layout.setSpacing(spacing)
        ## self.groupBox2.setLayout(layout)
        ## w = QtGui.QPushButton("Open ...")
        ## layout.addWidget(w)
        ## self.loadLigButton = w

        ## create group for placing the box
        ##
        self.groupBox3 = QtGui.QGroupBox(self)
        self.groupBox3.setTitle(self.tr("docking box"))
        layout1 = QtGui.QVBoxLayout(self.groupBox3) # VBox because we might add manual grid widget below later
        layout1.setContentsMargins(LeftMargin,topMargin,rightMargin,bottomMargin)
        layout1.setSpacing(spacing)
        self.groupBox3Layout = layout1
        #layout = QtGui.QGridLayout() # grid layout to hold bar and padding
        layout = QtGui.QHBoxLayout() # grid layout to hold bar and padding
        layout1.addLayout(layout)
        
        self.toolBar1 = QtGui.QToolBar()
        self.toolBar1.setOrientation(QtCore.Qt.Horizontal)
        layout.addWidget(self.toolBar1)#, 0, 0, 1, 1)
        
        #self.ASGridAct = act = QtGui.QAction(
        #    QtGui.QIcon(os.path.join(ICONPATH, 'center_grid_AutoSite.png')),
        #    'AutoSite', self)
        #act.setStatusTip('place grid based on AutoSite prediction')
        #act.setDisabled(True)
        #self.toolBar1.addAction(act)
        
        self.recGridAct = act = QtGui.QAction(
            QtGui.QIcon(os.path.join(ICONPATH, 'center_grid_REC.png')),
            'box entire receptor', self)
        act.setStatusTip('make the box cover the receptor')
        act.setDisabled(True)
        self.toolBar1.addAction(act)

        self.ligGridAct = act = QtGui.QAction(
            QtGui.QIcon(os.path.join(ICONPATH, 'center_grid_LIG.png')),
            'box entire ligand', self)
        act.setStatusTip('make the box cover the ligand')
        #act.setCheckable(True)
        act.setDisabled(True)
        self.toolBar1.addAction(act)

        self.TPGridAct = act = QtGui.QAction(
            QtGui.QIcon(os.path.join(ICONPATH, 'center_grid_TRANS.png')),
            'box all displayed ligand binding pockets', self)
        act.setStatusTip('make the box cover the visible ligand binding pockets')
        #act.setCheckable(True)
        act.setDisabled(True)
        self.toolBar1.addAction(act)

        self.residuesGrdAct = act = QtGui.QAction(
            QtGui.QIcon(os.path.join(ICONPATH, 'center_grid_RES.png')),
            'box around selected residues', self)
        act.setStatusTip('set box around selected residues')
        act.setCheckable(True)
        act.setDisabled(True)
        self.toolBar1.addAction(act)

        self.manualGrdAct = act = QtGui.QAction(
            QtGui.QIcon(os.path.join(ICONPATH, 'center_grid_HAND.png')),
            'set grid manually', self)
        act.setStatusTip('size and place box manually')
        act.setCheckable(True)
        act.setDisabled(True)
        self.toolBar1.addAction(act)

        # add grid padding widget
        #w = QtGui.QLabel('padding:')
        #w.setAlignment(QtCore.Qt.AlignRight)
        #layout.addWidget(w)#, 0, 1, 1, 1)

        self.gridPaddingWidget = w = QtGui.QDoubleSpinBox()
        w.setDecimals(1)
        w.setValue(4.0)
        w.setMaximumSize(60,20)
        #layout.addWidget(w)#, 0, 2, 1, 1)

        flayout = QtGui.QFormLayout()
        flayout.addRow('padding:', self.gridPaddingWidget)
        layout.addLayout(flayout)
        
        #act.triggered.connect(self.gridGUI.setGridFullLigand)
        act.setDisabled(True)
        self.toolBar1.addAction(act)
        self.groupBox3.setDisabled(True)
        
        ## create group for flexible residues
        ##
        self.groupBox4 = QtGui.QGroupBox(self)
        self.groupBox4.setTitle(self.tr("[0 flexible residue(s)]"))
        layout = QtGui.QGridLayout(self.groupBox4)
        layout.setContentsMargins(LeftMargin,topMargin,rightMargin,bottomMargin)
        layout.setSpacing(spacing)

        w = self.flexResWidget = MyQLineEdit()
        layout.addWidget(self.flexResWidget, 0, 0, 1, 1)
        #self.flexResWidget.setStyleSheet("QLineEdit{background: grey;}")

        self.setFlexResButton = QtGui.QPushButton(QtGui.QIcon(
            os.path.join(ICONPATH, 'flex_res_in_grid.png')),'')
        layout.addWidget(self.setFlexResButton, 0, 1, 1, 1)
        self.setFlexResButton.setCheckable(True)
        self.groupBox4.setDisabled(True)

        self.setFlexResGrdButton = QtGui.QPushButton(QtGui.QIcon(
            os.path.join(ICONPATH, 'resize_grid_per_FlexRes.png')),'')
        layout.addWidget(self.setFlexResGrdButton, 0, 2, 1, 1)
        self.setFlexResGrdButton.setDisabled(True)
        self.groupBox4.setDisabled(True)

        ## create group for covalent docking residues
        ##
        self.covDockingGroupBox = QtGui.QGroupBox(self)
        self.covDockingGroupBox.setTitle(self.tr("[covalent docking]"))
        self.covDockingGroupBox.setCheckable(True)
        self.covDockingGroupBox.setChecked(False)
        self.covDockingGroupBox.setDisabled(True)
        
        layout = QtGui.QVBoxLayout(self.covDockingGroupBox)
        layout.setContentsMargins(LeftMargin,topMargin,rightMargin,bottomMargin)
        layout.setSpacing(spacing)

        layout1 = QtGui.QHBoxLayout()
        lab1 =  QtGui.QLabel('atom1:')
        lab2 =  QtGui.QLabel('atom2:')
        self.covAt1Lab = QtGui.QLabel('       ')
        self.covAt2Lab = QtGui.QLabel('       ')
        layout1.addWidget(lab1)
        layout1.addWidget(self.covAt1Lab)
        layout1.addWidget(lab2)
        layout1.addWidget(self.covAt2Lab)
        layout.addLayout(layout1)
        
        ## create group for translational points
        ##
        self.groupBox5 = QtGui.QGroupBox(self)
        self.groupBox5.setTitle(self.tr("ligand binding pockets"))
        #self.groupBox5.setStyleSheet("QGroupBox {  margin-top: 11px}")

        layout = QtGui.QVBoxLayout(self.groupBox5)
        layout.setContentsMargins(LeftMargin,0,rightMargin,0)
        layout.setSpacing(0)

        layout2 = QtGui.QHBoxLayout()
        layout2.setContentsMargins(LeftMargin,0,rightMargin,0)
        layout2.setSpacing(spacing)
        asVer = "1.0"
        if self.app()._autoSite2: asVer = "1.1"
        self.computeTPointsButton = w = QtGui.QPushButton("compute pockets [AutoSite %s]" % asVer)
        w.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        layout2.addWidget(w)
        self.autoSiteVButton = w = AutoSiteVersionSwitch(version=asVer)
        w.setSizePolicy(QtGui.QSizePolicy.Maximum, QtGui.QSizePolicy.Maximum)#QtGui.QSizePolicy.Expanding)
        layout2.addWidget(w)
        layout.addLayout(layout2)

        self.clustersWidget = w = QtGui.QTableWidget(self.groupBox5)
        w. setSelectionMode(QtGui.QAbstractItemView.SelectionMode.NoSelection) 
        w.setColumnCount(5)
        w.verticalHeader().hide()
        w.setMinimumSize(100, 100)
        w.setHorizontalHeaderLabels(
            ['fills', 'AS Score', '#Points', 'RadGyr','Buriedness'])
        w.horizontalHeader().setResizeMode(0, QtGui.QHeaderView.Stretch)
        #self.clustersWidget.resizeColumnsToContents()
        #self.clustersWidget.horizontalHeader().sortIndicatorChanged.connect(self.sortClusters)
        layout.addWidget(w)

        layout1 = QtGui.QHBoxLayout()
        layout1.setContentsMargins(LeftMargin,0,rightMargin,0)
        layout1.setSpacing(spacing)
        self.selectAll = QtGui.QPushButton("select all")
        self.deselectAll = QtGui.QPushButton("deselect all")
        self.invertSelection = QtGui.QPushButton("invert")
        layout1.addWidget(self.selectAll)
        layout1.addWidget(self.deselectAll)
        layout1.addWidget(self.invertSelection)
        layout.addLayout(layout1)

        # create set of buttons for feature points
        #fplayout = QtGui.QHBoxLayout()
        #fplayout.setContentsMargins(LeftMargin,topMargin,rightMargin,bottomMargin)
        #fplayout.setSpacing(spacing)
        #layout.addLayout(fplayout)

        # compute FP button
        #self.featurePtsButton = w = QtGui.QPushButton("compute feature points")
        #w.setDisabled(True)
        #fplayout.addWidget(w)

        # pick residues from list contact FPs
        #self.setFlexResButtonFromFPButton = QtGui.QPushButton(QtGui.QIcon(
        #    os.path.join(ICONPATH, 'cube.png')),'')
        #fplayout.addWidget(self.setFlexResButtonFromFPButton)
        #self.setFlexResButtonFromFPButton.setCheckable(True)

        self.groupBox5.setDisabled(True)

        ## create group for computing affinity maps
        ##
        self.groupBox6 = QtGui.QGroupBox(self)
        self.groupBox6.setTitle(self.tr("affinity maps"))
        layout = QtGui.QGridLayout(self.groupBox6)
        layout.setContentsMargins(0,0,rightMargin,0)
        layout.setSpacing(0)

        layout1 = QtGui.QHBoxLayout()
        layout1.setContentsMargins(LeftMargin,0,rightMargin,0)
        layout1.setSpacing(spacing)
        
        self.compAllButton = button1 = QtGui.QRadioButton("for all atom types")
        button1.setChecked(True)
        button1.clicked.connect(self.configureEditButton)
        layout1.addWidget(button1)

        self.wMapParamB = wMapB = WaterMapParamsButton("water map setting")
        layout1.addWidget(wMapB)
        layout.addLayout(layout1,0,0,1,3)

        self.compLigandAtypes = button2 = QtGui.QRadioButton(" for types:")
        button2.clicked.connect(self.configureEditButton)
        layout.addWidget(button2, 1, 0, 1, 1)

        self.atypesLabel = QtGui.QLabel('C OA HD')
        layout.addWidget(self.atypesLabel, 1, 1, 1, 1, QtCore.Qt.AlignLeft)

        self.editAtypesButton = w = QtGui.QPushButton("edit ...")
        w.setDisabled(True)
        layout.addWidget(w, 1, 2, 1, 1)
        
        # group for adding gradient:
        self.addGradGroupBox = QtGui.QGroupBox(self)
        self.addGradGroupBox.setContentsMargins(LeftMargin,topMargin,rightMargin,bottomMargin)
        self.addGradGroupBox.setTitle(self.tr("create gradients inside receptor"))
        self.addGradGroupBox.setCheckable(True)
        self.addGradGroupBox.setChecked(False)

        self.spacerLabel = QtGui.QLabel(' ')
        layout.addWidget(self.spacerLabel, 2, 0, 1, 3)

        gradLayout = QtGui.QGridLayout(self.addGradGroupBox)
        gradLayout.setContentsMargins(LeftMargin,4,rightMargin,0)
        gradLayout.setSpacing(spacing)
        w =  QtGui.QLabel("define exterior as:")
        gradLayout.addWidget(w,0,0,1,1)
        self.gradLargeCl = b1 = QtGui.QRadioButton("largest cluster of negative values" )
        b1.setChecked(True)
        b1.clicked.connect(self.configureGradCutOffBox)
        gradLayout.addWidget(b1, 1, 0, 1,2)
        self.useGradCutOff = b2 = QtGui.QRadioButton("all clusters larger than:")
        b2.clicked.connect(self.configureGradCutOffBox)
        b2.setChecked(False)
        gradLayout.addWidget(b2, 2, 0,1,1)
        self.gradCutOffBox = w = QtGui.QSpinBox()
        w.setDisabled(True)
        gradLayout.addWidget(w, 2, 1, 1, 1, QtCore.Qt.AlignLeft)
        layout.addWidget(self.addGradGroupBox, 3, 0, 1, 3)
        self.groupBox6.setDisabled(True)

        #hlayout = QtGui.QHBoxLayout()
        #hlayout.setContentsMargins(LeftMargin,0,rightMargin,0)
        #hlayout.setSpacing(spacing)
        #self.toolBar2 = QtGui.QToolBar()

        #self.TPOKAct = act = QtGui.QAction(
        #    QtGui.QIcon(os.path.join(ICONPATH, 'cancel24.png')),
        #    'select at least one binding pockets overlaping with the box', self)
        #act.setStatusTip('select at least one binding pockets overlaping with the box')
        #self.toolBar2.addAction(act)

        #self.gridOKAct = act = QtGui.QAction(
        #    QtGui.QIcon(os.path.join(ICONPATH, 'ok24.png')),
        #    'Docking box covers moving atoms', self)
        #act.setStatusTip('docking box covers moving receptor atoms')
        #act.setDisabled(True)
        #self.toolBar2.addAction(act)
        
        #if sys.platform=='darwin':
        #    self.toolBar2.setMaximumWidth(100)
        #elif sys.platform=='linux2':
        #    self.toolBar2.setMaximumWidth(70)
        #elif sys.platform=='win32':
        #    self.toolBar2.setMaximumWidth(80)
        #else:
        #    self.toolBar2.setMaximumWidth(100)

        #hlayout.addWidget(self.toolBar2)

        w = self.computeGridsButton = QtGui.QPushButton("generate target file ... ")
        w.setDisabled(True)
        #hlayout.addWidget(w)
        #layout.addLayout(hlayout, 4, 0, 1, 3)
        layout.addWidget(w, 4, 0, 1, 3)

        self.progressBar = QtGui.QProgressBar()
        layout.addWidget(self.progressBar, 5, 0, 1, 3)
        self.progressBar.setMinimum(0)
        #self.progressBar.setMaximum(100)
        self.progressBar.setValue(0)
        self.progressBar.setContentsMargins(LeftMargin,topMargin,rightMargin,bottomMargin)
        
        ## w = self.dockButton = QtGui.QPushButton('dock ...')
        ## w.clicked.connect(self.displayDockingPanel)
        ## w.setDisabled(True)
        ## layout.addWidget(w, 4, 0, 1, 3)
        
        mainLayout = QtGui.QGridLayout()
        #mainLayout.setContentsMargins(2,2,2,2)
        #mainLayout.setSpacing(2)
        #mainLayout.setContentsMargins(LeftMargin,10,rightMargin,10)
        #mainLayout.setSpacing(spacing)
        #mainLayout.addWidget(self.groupBox1, 0, 0, 1, 1)
        #mainLayout.addWidget(self.groupBox2, 0, 1, 1, 1)
        mainLayout.addWidget(self.groupBox3, 1, 0, 1, 2)
        mainLayout.addWidget(self.groupBox4, 2, 0, 1, 2)
        mainLayout.addWidget(self.covDockingGroupBox, 3, 0, 1, 2)
        mainLayout.addWidget(self.groupBox5, 4, 0, 1, 2)
        mainLayout.addWidget(self.groupBox6, 5, 0, 1, 2)
        self.setLayout(mainLayout)

    ## def displayDockingPanel(self):
    ##     from ADFR.GUI.ADFRgui import SinglLocalDockingWidget
    ##     widget = self.dockingWidget = SinglLocalDockingWidget(self.app().viewer)
    ##     widget.inputWidget.ligandEntryWidget.setText(self.app().ligand.filename)
    ##     widget.inputWidget.mapsEntryWidget.setText(self.app().destinationFolder+'.zip')
    ##     widget.inputWidget.outputNameWidget.setText('NoName')
    ##     widget.inputWidget.gaNbWidget.setValue(10)
    ##     widget.inputWidget.refLigWidget.setText(self.app().ligand.filename)


            
    ## attempt to hijack sorting of clusters but
    ## sorting still happens before this function is called
    def sortClusters(self, logicalIndex, order):
        pass
        #print ';SORT', logicalIndex, order
        
    def configureGenerateMaps(self):
        if self._tpok and self._gridok and self._boxok:
            self.computeGridsButton.setDisabled(False)
            
        else:
            self.computeGridsButton.setDisabled(True)

    def disableGridCheck(self, value):
        self.gridOKAct.setDisabled(value)
        
    def handleTPointsSignal(self, value, msg='translation points needed'):
        if value:
            self.TPGridAct.setDisabled(False)
            #self.TPOKAct.setDisabled(False)
            #self.featurePtsButton.setDisabled(False)
            self.TPOKAct.setStatusTip('ligand binding pocket is OKAY')
            self.TPOKAct.setText('ligand binding pocket is OKAY')
            self.TPOKAct.setIcon(QtGui.QIcon(os.path.join(ICONPATH, 'ok24.png')))
            self._tpok = True
        else:
            self.TPGridAct.setDisabled(True)
            #self.featurePtsButton.setDisabled(True)
            #self.TPOKAct.setDisabled(True)
            #print "handleTPointsSignal, msg", msg
            self.TPOKAct.setStatusTip(msg)
            self.TPOKAct.setText(msg)
            self.TPOKAct.setIcon(QtGui.QIcon(os.path.join(ICONPATH, 'cancel24.png')))
            self._tpok = False
        self.configureGenerateMaps()
        
    def handleGridOKSignal(self, value):
        if value:
            self.gridOKAct.setIcon(QtGui.QIcon(os.path.join(ICONPATH, 'ok24.png')))
            self.gridOKAct.setText('docking box covers moving receptor atoms')
            self.gridOKAct.setStatusTip('docking box covers moving receptor atoms')
            self._gridok = True
        else:
            self._gridok = False
            #self.showFPAction.setDisabled(True)
            #self.TSpheresC.Set(vertices=[])
            #self.TSpheresO.Set(vertices=[])
            #self.TSpheresH.Set(vertices=[])
            self.gridOKAct.setIcon(QtGui.QIcon(os.path.join(ICONPATH, 'cancel24.png')))
            self.gridOKAct.setText('docking box does NOT cover moving receptor atoms')
            self.gridOKAct.setStatusTip('docking box does NOT cover moving receptor atoms')
        self.configureGenerateMaps()

    def configureEditButton(self):
        if self.compAllButton.isChecked():
            self.editAtypesButton.setDisabled(True)
        elif self.compLigandAtypes.isChecked():
            self.editAtypesButton.setDisabled(False)
            self.atypesLabel.setDisabled(False)

    def configureGradCutOffBox(self):
        if self.gradLargeCl.isChecked():
            self.gradCutOffBox.setDisabled(True)
        elif self.useGradCutOff.isChecked():
            self.gradCutOffBox.setDisabled(False)

class AutoSiteVersionSwitch(QtGui.QPushButton):
    # A tool button with drop-down menu to select AutoSite version and
    # some AutoSite parameters
    autosite2 = QtCore.Signal(bool)
    ligSizeChanged = QtCore.Signal(int)
    usePepScore = QtCore.Signal(bool)
    
    def __init__(self, version="1.0", parent=None):
        super(AutoSiteVersionSwitch, self).__init__(parent)
        self.setIcon(QtGui.QIcon(os.path.join(ICONPATH, "model.png")))
        self.setMenu(QtGui.QMenu(self))
        action = QtGui.QWidgetAction(self)
        # Add a groupBox with widgets to for the dropdown menu
        menuGrBox = QtGui.QGroupBox()
        gbLayout = QtGui.QVBoxLayout()
        lout1 = QtGui.QHBoxLayout()
        # radiobutton group to select the version
        rbGroup = QtGui.QButtonGroup(self)
        self.rb1 = rb1 = QtGui.QRadioButton("1.0")
        self.rb2 = rb2 = QtGui.QRadioButton("1.1")
       
        rbGroup.addButton(rb1)
        rbGroup.addButton(rb2)
        lout1.addWidget(QtGui.QLabel("version: "))
        lout1.addWidget(rb1)
        lout1.addWidget(rb2)
        lout1.setContentsMargins(0,0,0,0)
        gbLayout.addLayout(lout1)
        rb1.toggled.connect(self.setAutoSiteVersion)
        # add widgets to set AutoSite2 (version 1.1) parameters
        lout2 = QtGui.QHBoxLayout()
        self.pepCkeckB = pb = QtGui.QCheckBox("peptide scoring func")
        pb.stateChanged.connect(self.useAutoSitePepScore)
        lout2.addWidget(pb)
        lout2.setContentsMargins(0,0,0,0)
        lout2.addStretch()
        pb.setChecked(True)
        #pb.setLayoutDirection(QtCore.Qt.RightToLeft) #puts a space in front of the text label.
        pb.setContentsMargins(0,0,0,0)
        gbLayout.addLayout(lout2)
        
        lout3 =  QtGui.QHBoxLayout()
        lab = QtGui.QLabel("ligand size: ")
        self.sizeSB = sb = QtGui.QSpinBox()
        sb.setMinimum(1)
        sb.setMaximum(999999)
        sb.setValue(500)
        lout3.addWidget(lab)
        lout3.addWidget(sb)
        lout3.setContentsMargins(0,0,0,0)
        sb.valueChanged.connect(self.sizeChanged)
        gbLayout.addLayout(lout3)

        if version == "1.0":
            rb1.setChecked(True)
            self.pepCkeckB.setEnabled(False)
            self.sizeSB.setEnabled(False)
        else:
            rb2.setChecked(True)
            self.pepCkeckB.setEnabled(True)
            self.sizeSB.setEnabled(True)
        menuGrBox.setLayout(gbLayout)
        action.setDefaultWidget(menuGrBox)
        self.menu().addAction(action)
        
    def setAutoSiteVersion(self):
        autoSite1 = self.rb1.isChecked()
        autoSite2 = self.rb2.isChecked()
        self.autosite2.emit(autoSite2)
        # enable/disable  widgets to set AutoSite2 parameters 
        self.pepCkeckB.setEnabled(autoSite2)
        self.sizeSB.setEnabled(autoSite2)
        
    def sizeChanged(self, val):
        #print "Lig size: ", val
        self.ligSizeChanged.emit(val)

    def useAutoSitePepScore(self, val):
        # returns True or False
        #print "useAutoSitePepScore val:", val
        self.usePepScore.emit(self.pepCkeckB.isChecked())
    
class WaterMapParamsButton(QtGui.QPushButton):
    entropyChanged = QtCore.Signal(float)
    weightChanged = QtCore.Signal(float)
    
    def __init__(self, text, parent=None):
        super(WaterMapParamsButton, self).__init__(text, parent)
        self.setMenu(QtGui.QMenu(self))
        action = QtGui.QWidgetAction(self)
        menuGrBox = QtGui.QGroupBox()
        gbLayout = QtGui.QGridLayout()
        self.wWeight = w = QtGui.QDoubleSpinBox()
        w.setSingleStep(.1)
        w.valueChanged.connect(self.weightChanged.emit)
        w.setValue(0.6)
        lab1 =  QtGui.QLabel("weight: ")
        gbLayout.addWidget(lab1, 0, 0, 1, 1)
        gbLayout.addWidget(w, 0, 1, 1, 1)
        self.wEntropy = w = QtGui.QDoubleSpinBox()
        w.valueChanged.connect(self.entropyChanged.emit)
        w.setMinimum(-1000.)
        w.setSingleStep(.1)
        w.setValue(-0.2)
        lab2 =  QtGui.QLabel("entropy: ")
        gbLayout.addWidget(lab2, 1, 0, 1, 1)
        gbLayout.addWidget(w, 1, 1, 1, 1)
        menuGrBox.setLayout(gbLayout)
        action.setDefaultWidget(menuGrBox)
        self.menu().addAction(action)

        
class CovalentLigAtomsList(QtGui.QDialog):
    """A dialog to select/add atoms to ignore in covalent docking  """
    resSelected = QtCore.Signal(str, bool)
    oncloseSignal = QtCore.Signal(bool)
    
    def __init__(self, parent=None):
        super(CovalentLigAtomsList, self).__init__(parent)
        self.setWindowTitle("define covalent bond")

        # place the widget close to where we clicked
        rel_pos = QtGui.QCursor.pos()
        pos = self.mapToGlobal(rel_pos)
        self.move(pos.x()+20, pos.y()+15)

        self.app = None
        self.immediate = True
        self.resList = []
        self.items = {}
        self.buildForm()

    def addItems(self, atomList):
        hv = atomList.getHierView()
        for ch in hv.iterChains():
            for res in ch.iterResidues():
                chid, resname, resnum, icode = res.getChid(), res.getResname(), res.getResnum(), res.getIcode()
                if icode:
                    resName = "%s:%s%d%s:"%(chid, resname, resnum, icode)
                else:
                    resName = "%s:%s%d:"%(chid, resname, resnum)
                icode = res.getIcode()
                if icode:
                    resName = "%s:%s"%(resName, icode)
                atomsStr = ""
                atinds = []
                for a in res.iterAtoms():
                    atname = a.getName()
                    atomsStr += "%s,"% atname
                    atinds.append(a.getIndex())
                resStr = resName + atomsStr
                self.resList.append(resStr)
                self.items[resStr] = atinds
                item = QtGui.QListWidgetItem(resStr)
                #item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
                self.resListWidget.addItem(item)
        
    def fill(self, atomList):
        self.resList = []
        self.resListWidget.clear()
        self.items = {}
        self.addItems(atomList)

    def buildForm(self):
        layout = QtGui.QVBoxLayout()
        layout1 = QtGui.QHBoxLayout()
        layout1.setContentsMargins(LeftMargin,topMargin,rightMargin,bottomMargin)
        layout1.setSpacing(spacing)
        
        self.covAt1Lab = QtGui.QLabel('atom1:')
        self.covAt1w = QtGui.QLineEdit()
        self.pickCovAt1w =  QtGui.QPushButton(QtGui.QIcon(
            os.path.join(ICONPATH, 'pickAtom.png')),'')
        self.pickCovAt1w.setCheckable(True)

        self.covAt2Lab = QtGui.QLabel('atom2:')
        self.covAt2w = QtGui.QLineEdit()
        self.pickCovAt2w =  QtGui.QPushButton(QtGui.QIcon(
            os.path.join(ICONPATH, 'pickAtom.png')),'')
        self.pickCovAt2w.setCheckable(True)
        
        layout1.addWidget(self.covAt1Lab)
        layout1.addWidget(self.covAt1w)
        layout1.addWidget(self.pickCovAt1w)
        layout1.addWidget(self.covAt2Lab)
        layout1.addWidget(self.covAt2w)
        layout1.addWidget(self.pickCovAt2w)
        layout.addLayout(layout1)
        
        layout2 = QtGui.QHBoxLayout()
        layout2.setContentsMargins(LeftMargin,topMargin,rightMargin,bottomMargin)
        layout2.setSpacing(spacing)
        covResLabel = QtGui.QLabel('limit cov. lig.:')
        self.covalentResEdit = QtGui.QLineEdit()
        layout2.addWidget(covResLabel)
        layout2.addWidget(self.covalentResEdit)
        layout.addLayout(layout2)
        
        layout3 =  QtGui.QHBoxLayout() # contains the list widget with checkbuttons
        layout3.setContentsMargins(LeftMargin,topMargin,rightMargin,bottomMargin)
        layout3.setSpacing(spacing)
        self.resListWidget = QtGui.QListWidget()
        layout3.addWidget(self.resListWidget)
        layout.addLayout(layout3)
        
        #layout4 = QtGui.QHBoxLayout() # contains All & None buttons
        #self.allButton = QtGui.QPushButton("All")
        #self.allButton.clicked.connect(self.selectAll)
        #self.noneButton = QtGui.QPushButton("None")
        #self.noneButton.clicked.connect(self.deselectAll)
        #layout4.addWidget(self.allButton)
        #layout4.addWidget(self.noneButton)

        #layout.addLayout(layout4)
        self.startOverButton =  QtGui.QPushButton("Start Over")
        self.startOverButton.clicked.connect(self.clearForm)
        self.buttonBox = QtGui.QDialogButtonBox(QtCore.Qt.Horizontal, self)
        self.buttonBox.addButton(self.startOverButton, QtGui.QDialogButtonBox.ActionRole)
        self.buttonBox.addButton(QtGui.QDialogButtonBox.Ok)
        self.buttonBox.addButton(QtGui.QDialogButtonBox.Cancel)
        self.buttonBox.button(self.buttonBox.Ok).setDefault(False)
        self.startOverButton.setDefault(False)
        self.buttonBox.button(self.buttonBox.Cancel).setDefault(False)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox)
        self.setLayout(layout)

    def keyPressEvent(self, evt):
        # This prevents the closing of the dialog when "return" key is pressed in the QLineEdit widget
        if evt.key() == QtCore.Qt.Key_Enter or evt.key() == QtCore.Qt.Key_Return:
            return
        QtGui.QDialog.keyPressEvent(self, evt)

    def getAtomList(self):
        resDict = {}
        resAtomList = []
        atinds = []
        for i in range(self.resListWidget.count()):
            item = self.resListWidget.item(i)
            txt = str(item.text())
            resAtomList.append(txt)
            atinds.extend(self.items[txt])
            txtList = txt.split(":")
            chid = txtList[0]
            if not resDict.has_key(chid):
                resDict[chid] = txtList[1]
            else:
                resDict[chid] += ","+txtList[1]
        resStr = ""
        for chid, res in resDict.items():
            if len(resStr): resStr+=";"
            resStr += "%s:%s"%(chid,res)
        return resStr, resAtomList, atinds

    def accept(self):
        self.oncloseSignal.emit(True)
        QtGui.QDialog.accept(self)
        self.app().freezeUI(False)
        #super(CovalentLigAtomsList, self).accept()
        
    def reject(self):
        self.oncloseSignal.emit(False)
        super(CovalentLigAtomsList, self).reject()
        self.app().freezeUI(False)

    def clearForm(self):
        self.covAt1w.setText("")
        self.covAt2w.setText("")
        self.covalentResEdit.setText("")
        self.resListWidget.clear()
        self.items = {}
        self.pickCovAt2w.setChecked(False)
        #self.pickCovAt1w.click()
        self.pickCovAt1w.setChecked(True)
        self.app().viewer.processPicking = self.app().pickFirstCovalentBondAtom



from ADFR.utils.runAGFR import runAGFR
from .statusBar import MyStatusBar

class Drawer(QtGui.QWidget):
    def __init__(self, parent=None):
        super(Drawer, self).__init__(parent)
        self.setWindowFlags(QtCore.Qt.Tool | QtCore.Qt.FramelessWindowHint |
                            QtCore.Qt.WindowStaysOnTopHint) # need to be on top on mac
        
        self._animation = QtCore.QPropertyAnimation(self, "size")
        self._animation.setDuration(200)
        self._animation.valueChanged.connect(self.setFixedSize)
        self.hide()

    def setToolBar(self, tb):
        self.tb = tb

    def open_(self):
        h = self.tb.getHeight()
        w = self.tb.getWidth()
        if self.tb._orient=='vertical':
            self._animation.setStartValue(QtCore.QSize(0, h))
            self._animation.setEndValue(QtCore.QSize(w, h))
        else:
            self._animation.setStartValue(QtCore.QSize(h, 0))
            self._animation.setEndValue(QtCore.QSize(w, h))
        self._animation.setDirection(QtCore.QAbstractAnimation.Forward)
        self._animation.start()
        self.show()

    def close(self):
        self._animation.setDirection(QtCore.QAbstractAnimation.Backward)
        self._animation.start()

class SmallToolBar(QtGui.QFrame):
#class SmallToolBar(QtGui.QWidget):

    def __init__(self, size, orientation, parent=None):
        super(SmallToolBar, self).__init__(parent)
        assert orientation in ['vertical', 'horizontal']

        #self.setWindowFlags(QtCore.Qt.FramelessWindowHint)
        self.setFrameStyle(QtGui.QFrame.StyledPanel | QtGui.QFrame.Plain)
        self._orient = orientation
        self._size = size
        self._widgets = []
        if orientation=='vertical':
            self.layout = QtGui.QVBoxLayout()
        else:
            self.layout = QtGui.QHBoxLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(self.layout)
        self.setStyleSheet("background-color: white")

    def addWidget(self, widget):
        widget.setIconSize(QtCore.QSize(self._size, self._size))
        widget.setFixedSize(QtCore.QSize(self._size, self._size))
        self._widgets.append(widget)
        self.layout.addWidget(widget)

    def reset(self):
        for w in self._widgets:
            if w.isCheckable():
                w.setChecked(False)
        
    def getWidth(self):
        if self._orient=='vertical':
            return self._size + 4
        else:
            return (4+self._size)*len(self._widgets)
        
    def getHeight(self):
        if self._orient=='vertical':
            return (4+self._size)*len(self._widgets)
        else:
            return self._size + 4

class LigandToolBar(SmallToolBar):

    def __init__(self, size, orientation, app, parent=None):
        super(LigandToolBar, self).__init__(size, orientation, parent)
        self.app = weakref.ref(app)

        # show hide ligand
        button = QtGui.QToolButton(self)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(os.path.join(ICONPATH, 'hidden.png')),
                       QtGui.QIcon.Normal,  QtGui.QIcon.On )
        icon.addPixmap(QtGui.QPixmap(os.path.join(ICONPATH, 'visible.png')),
                       QtGui.QIcon.Normal,  QtGui.QIcon.Off)
        button.setIcon(icon)
        button.setCheckable(True)
        button.setStyleSheet("QToolButton { border: 2px solid #8f8f91; border-radius: 6px;}")
        button.toggled.connect(self.app().showHideLigand)
        self.addWidget(button)
        
        # focus on ligand
        button = QtGui.QPushButton(self)
        button.setIcon(QtGui.QIcon(os.path.join(ICONPATH, 'focus.png')))
        button.setStyleSheet("QToolButton { border: 2px solid #8f8f91; border-radius: 6px;}")
        button.released.connect(self.app().focusLigand)
        self.addWidget(button)

        # label atoms buttons
        for name, prop in [('labName.png', 'name'),
                           ('labType.png', 'AD_element'),
                           ('labq.png', 'charge')]:
            button = QtGui.QToolButton(self)
            button.setIcon(QtGui.QIcon(os.path.join(ICONPATH, name)))
            button._label = prop
            #button.setStyleSheet("QToolButton { border: 2px solid #8f8f91; border-radius: 6px;\nQToolButton:checked {background-color: yellow;}")
            button.toggled.connect(self.label)
            button.setCheckable(True)
            self.addWidget(button)

        # view/edit torsion tree
        button = QtGui.QToolButton(self)
        button.setIcon(QtGui.QIcon(os.path.join(ICONPATH, 'torTree.png')))
        button.setCheckable(True)
        #button.setStyleSheet("QToolButton { border: 2px solid #8f8f91; border-radius: 6px;\nQToolButton:checked {background-color: yellow;}")
        button.toggled.connect(self.app().showHideTorsionTree)
        button.setToolTip("view/edit torsion")
        self.addWidget(button)

        # save ligand
        button = QtGui.QPushButton(self)
        button.setIcon(QtGui.QIcon(os.path.join(ICONPATH, 'folder.png')))
        button.setStyleSheet("QToolButton { border: 2px solid #8f8f91; border-radius: 6px;}")
        button.released.connect(self.saveLigand)
        self.addWidget(button)

        # toggle rotatable bonds
        #button = QtGui.QToolButton(self)
        #button.setIcon(QtGui.QIcon(os.path.join(ICONPATH, 'toggleBonds.png')))
        #button.setCheckable(True)
        #button.setStyleSheet("QToolButton { border: 2px solid #8f8f91; border-radius: 6px;\nQToolButton:checked {background-color: yellow;}")
        #button.toggled.connect(self.app().toggleRotatableBonds)
        #self.addWidget(button)

    def saveLigand(self):
        ADlig = self.app().agfr._adLigs[0]
        fdialog = QtGui.QFileDialog()
        fdialog.setDefaultSuffix("pdbqt")
        name = ADlig.molPrep.name+'.pdbqt'
        filepath = fdialog.getSaveFileName(self, "ligand filename", name, 
                                           filter="PDBQT files (*.pdbqt);; Any files (*)")[0]

        if filepath is None: # cancelled dialog
            return
        elif filepath[-6:] != '.pdbqt':
            filepath += '.pdbqt'
        lines = ADlig.getPDBQTlines()
        f = open(filepath, 'w')
        [f.write('%s\n'%l) for l in lines]
        f.close()
         
    def label(self, checked):
        button = self.sender()
        if not checked:
            if button._label=='charge':
                self.app().pmv.unlabelAtoms(
                    self.app().ligand.select(), propName=button._label, formatString='%.2f')
            else:
                self.app().pmv.unlabelAtoms(
                    self.app().ligand.select(), propName=button._label, formatString=None)

            return

        for bt in self._widgets[2:5]:
            if not bt._label == button._label:
                bt.setChecked(False)

        if button._label=='charge':
            self.app().pmv.labelAtoms(self.app().ligand.select(), propName=button._label,
                                      formatString='%.2f', fontSize=(0.25, 0.25, 0.1))
        else:
            self.app().pmv.labelAtoms(self.app().ligand.select(), propName=button._label,
                                      formatString=None, fontSize=(0.35, 0.35, 0.1))
        

class TopBar(QtGui.QWidget):

    def __init__(self, parent=None):
        super(TopBar, self).__init__(parent)
        self.app = weakref.ref(parent)
        vlayout = QtGui.QVBoxLayout()
        hlayout = QtGui.QHBoxLayout()
        vlayout.addLayout(hlayout)
        ## input group defining receptor from file folder or fetch from web
        ##
        size = 40
        self.groupBox1 = QtGui.QGroupBox(self)
        self.groupBox1.setTitle(self.tr("load receptor"))
        layout = QtGui.QHBoxLayout()
        self.groupBox1.setLayout(layout)
        w = QtGui.QPushButton(QtGui.QIcon(os.path.join(ICONPATH, 'rec1.png')), '')
        w.setIconSize(QtCore.QSize(size, size))
        layout.addWidget(w)
        self.loadRecButton = w
        
        w = QtGui.QPushButton(QtGui.QIcon(os.path.join(ICONPATH, 'folder.png')), '')
        w.setIconSize(QtCore.QSize(size, size))
        layout.addWidget(w)
        self.setFolderButton = w
        
        w = QtGui.QPushButton(QtGui.QIcon(os.path.join(ICONPATH, 'www.png')), '')
        w.setIconSize(QtCore.QSize(size, size))
        layout.addWidget(w)
        self.fetchPDBButton = w
        hlayout.addWidget(self.groupBox1)

        self.groupBox2 = QtGui.QGroupBox(self)
        self.groupBox2.setTitle(self.tr("target"))
        layout = QtGui.QHBoxLayout()
        self.groupBox2.setLayout(layout)
        w = QtGui.QPushButton(QtGui.QIcon(os.path.join(ICONPATH, 'target.png')), '')
        w.setIconSize(QtCore.QSize(size, size))
        layout.addWidget(w)
        hlayout.addWidget(self.groupBox2)
        self.loadTargetButton = w
        
        self.groupBox3 = QtGui.QGroupBox(self)
        self.groupBox3.setTitle(self.tr("ligand"))
        layout = QtGui.QHBoxLayout()
        self.groupBox3.setLayout(layout)
        w = QtGui.QPushButton(QtGui.QIcon(os.path.join(ICONPATH, 'ligandNoArrow.png')), '')
        w.setIconSize(QtCore.QSize(size, size))
        layout.addWidget(w)
        hlayout.addWidget(self.groupBox3)
        self.loadLigButton = w

        size = 24
        hlayout2 = QtGui.QHBoxLayout()
        vlayout.addLayout(hlayout2)
        self.groupBox4 = QtGui.QGroupBox(self)
        layout = QtGui.QHBoxLayout()
        self.groupBox4.setLayout(layout)
        w = QtGui.QPushButton(QtGui.QIcon(os.path.join(ICONPATH, 'previous.png')), '')
        w.setIconSize(QtCore.QSize(size, size))
        w.setDisabled(True)
        self.previousWidget = w
        w.pressed.connect(self.app().gotoPreviousFile)
        layout.addWidget(w)

        w = QtGui.QPushButton(QtGui.QIcon(os.path.join(ICONPATH, 'next.png')), '')
        w.setIconSize(QtCore.QSize(size, size))
        w.setDisabled(True)
        self.nextWidget = w
        w.pressed.connect(self.app().gotoNextFile)
        layout.addWidget(w)
        hlayout2.addWidget(self.groupBox4)
        
        size = 26
        w = QtGui.QPushButton(QtGui.QIcon(os.path.join(ICONPATH, 'autofocus.png')), '')
        w.setIconSize(QtCore.QSize(size, size))
        w.setCheckable(True)
        w.setChecked(True)
        self.autofocusButton = w
        hlayout2.addWidget(w)
        
        hlayout2.addStretch()

        PMVICONPATH = findFilePath('Icons', 'PmvApp.GUI')
        w = QtGui.QPushButton(QtGui.QIcon(os.path.join(PMVICONPATH, 'PyShell.png')), '')
        w.setIconSize(QtCore.QSize(size, size))
        w.setCheckable(True)
        w.setChecked(False)
        w.toggled.connect(self.app().showPythonShell)
        self.PyShellButton = w
        hlayout2.addWidget(w)

        w = QtGui.QPushButton(QtGui.QIcon(os.path.join(ICONPATH, 'cogwheel.png')), '')
        w.setIconSize(QtCore.QSize(size, size))
        hlayout2.addWidget(w)
        w.setCheckable(True)
        w.setChecked(False)
        self.userPreferencesButton = w
        self.setLayout(vlayout)

class GridGUI(QtGui.QWidget):

    jobDone = QtCore.Signal()
    flexResChanged = QtCore.Signal()
    TPointsOK = QtCore.Signal(bool, str)
    gridOK = QtCore.Signal(bool)
    AutoSiteDone = QtCore.Signal()
    
    def __init__(self, parent=None, pmv=None, pmvViewer=None, eventHandler=None, files=None):
        super(GridGUI, self).__init__(parent)
        try:
            import mslib
            self._hasMSMS = True
        except RuntimeError:
            self._hasMSMS = False
        self._labDisplayMode = 0 # 0: hidden, 1: residues in box, 2: flexRes
        self._MSMSDisplayMode = 0 # 0: hidden, 1: surface in box, 3: complete surface
        self._hasReceptorSurface = False

        self._filesToLoopOver = None # use to support inputFolder cmdline argument
        self._currentFileIndex = 0
        self._haveTrg = False # set to True if target file is loaded or computed
        self._ligandMultiMol = None # is set to multimol when multi mol lig is read
        self._openDialogs = [] # use to close open dialogues. When a dialog is posted
        # it should add itself to that list and remove itself when dismissed 
        
        self._LL = [0, 0, 0] # lower left box corner
        self._UR = [1, 1, 1] # upper right box corner
        self._baseSize = [0,0,0]
        self._smooth = 0.5
        self._spacing = 0.375
        self.carbon_cutoff = -0.3
        self._clusterItems = []
        self.smallClusterCutOff = 10
        
        self._noFlexResParse = False # set to True to set entry without parsing string
        self._flexRes = [] # list of prody residues for flexible residues
        self._flexResSCAtoms = None # MolKit2 selection with all moving side chain atoms in flexible residues
        self._flexResAtoms = []
        self._covalentBondAtoms = [None, None] # will hold 2 atoms
        self._torsionAtom = None
        self._covalentDocking = False
        self._covalentLigAtoms = None
        self.agfr = runAGFR()
        self.jobDone.connect(self.finishUp)
        self.flexResChanged.connect(self.onFlexResChanged)
        
        self.agfr.getAllADatomTypes()
        self.PDBprofile = None
        self.receptor = None # will be he prody molecule for the receptor
        self.ligand = None # will be he prody molecule for the ligand
        self._mapTypes = ['C', 'OA', 'HD'] # this Gui attribute will hold a list of
        # map types selected by the user (when "for types" radio button is checked)
        # self.agfr.setMaptypes() function should be called before we compute the grids.
        self.agfr.setMapTypes(self._mapTypes)
        
        self._autoSite2 =  self.agfr.autoSite2 #False # run AutoSite1.0 (faster version) by default
        # AutoSite2 parameters (can be set by the user):
        self._ligandSize = 500
        self._pepScore = True 
        if pmv is None:
            self.createPmvApp()
        else:
            assert isinstance(Pmv, MolApp)
            self.pmv = pmv
        self.pmv.gui = weakref.ref(self)
        self.buildUI(pmvViewer=pmvViewer, parent=parent)

        self.setFiles(files)
        
        self.viewer.pickedAtomsFuncs.append(self.printPickedAtoms)
        self.viewer.pickedBondsFuncs.append(self.printPickedBonds)
        self.defaultProcessPicking = self.viewer.processPicking
        
        self.AutoSiteDone.connect(self.fillFillsTable)

        self.tmpFolder = os.path.join(getResourceFolderWithVersion(), 'tmp')
        if not os.path.exists(self.tmpFolder):
            os.mkdir(self.tmpFolder)
        elif not os.path.isdir(self.tmpFolder):
            msgBox = QtGui.QMessageBox(self)
            msgBox.setText("ERROR: %s is not a folder, please remove this file"%self.tmpFolder)
            msgBox.setStandardButtons(QtGui.QMessageBox.Ok)
            msgBox.exec_()

        #def _flexResChanged(frstr):
        #    self.paramsWidget.flexResWidget.setText(frstr)
        #    self.setFlexRes()
        self.setAppTitle()
        self.covLigAtomsDialog = None
        #outside = numpy.load('outside.npy')
        #from DejaVu2.Spheres import Spheres
        #sph = Spheres('outside', centers=outside, radii=(0.15,),
        #              inheritMaterial=False, materials=[[1,0,0]])
        #self.viewer.AddObject(sph)
        self.trgFiles = []
        self.trgMapGui = None
        #self.setStyleSheet("QLAYOUT { margin-top: 0; margin-bottom: 0}")
        self._pdbSites = {}
        self._altlocResidues = {}
        self._mutatedResidues = {}
        self._gridGuiOptions = {}
        self._xmMats = None # list of instance matrices for crystal mates
        self.preferences = UserPreference()
        self.preferences.addItem(PreferenceItem('Flexible Residues Out Of Box', "Always ask", validValues=["Always ask", "Adjust box", "Make sidechains rigid", "None"], category="General", doc="Specify action when selected flexible residues are found outside of the grid box."))
        self.preferences.addItem(PreferenceItem('Number of Flexible Residues', "Always ask", validValues=["Always ask", "None"], category="General" , doc="Show a warning when there are more flexible residues selected then the AutoDock can handle."))
        self.preferences.loadSettings()
        self.missingResLines = None
        self.missingResSpheres = None
        self.missingResLabels = None

        self.analyserWidget = None
        self._resOutsideBoxDialog = None
        self._flexResDialog = None
        self._dockedLigand = None

    def autofocus(self):
        """returns True is autofocus is on"""
        return self.topBar.autofocusButton.isChecked()
    
    def setFiles(self, files):
        print 'FILES', files
        self._filesToLoopOver = files
        if not self._filesToLoopOver:
            self.topBar.nextWidget.setDisabled(True)
            self.topBar.previousWidget.setDisabled(True)
        else:
            self.topBar.nextWidget.setDisabled(False)
        
    def closeEventNoAsk(self, event): # use for testing the gui
        self.close()

    def closeEvent(self, event):
        if ADFR.GUI._AGFRGUI_debug:
            #close all the windows:
            for w in self._openDialogs:
                w.close()
            self.close()
            return
        
        self.setWindowState(QtCore.Qt.WindowActive)
        self.activateWindow()
        self.raise_()
        reply = QtGui.QMessageBox.question(self, 'Message',
            "Are you sure to quit?", QtGui.QMessageBox.Yes |
            QtGui.QMessageBox.No, QtGui.QMessageBox.No)

        if reply == QtGui.QMessageBox.Yes:
            # We need to loop over a copy of the list because onClose with remove
            # elements from the list causing the for loop to not complete
            for w in self._openDialogs:
                #print 'closing', w
                w.close()
            event.accept()
        else:
            event.ignore()
        
    def setAppTitle(self):
        rname = lname = 'None'
        if self.PDBprofile:
            app = 'PDBprofiler'
            rname = self.PDBprofile._pdbid
        else:
            app = 'AGFR'
            if self.receptor:
                rname = self.receptor.name
            if self.ligand:
                lname = self.ligand.name
        import mglutil
        from ADFR import __version__, __subversion__, __revision__
        title = "%s %d.%d%s (built: %s) receptor: %s   ligand: %s"%(
            app, __version__,__subversion__,__revision__,
            mglutil.__revision__,rname, lname)

        try:
            from ADFRcc.adfrcc import get_ompthread_count
            nthreads = get_ompthread_count()
            if nthreads > 1:
                title += "   OpenMP: ON %d threads"% nthreads
            else:
                title += " OpenMP OFF"
        except ImportError:
                print "Failed to import get_ompthread_count() from ADFRcc.adfrcc"
        self.setWindowTitle(title)

    ## def resizeEvent(self, event):
    ##     print "RESIZE EVENT", event
    ##     if self.isVisible():
    ##         p1 = self.splitter.widget(0).sizePolicy()
    ##         p2 = self.splitter.widget(1).sizePolicy()
    ##         print "stretch:", p1.horizontalStretch(), p2.horizontalStretch(), self.splitter.size()
    ##         if p1.horizontalStretch() == 0 and p2.horizontalStretch() == 0:
    ##             self.splitter.setStretchFactor(0, 2)
    ##             self.splitter.setStretchFactor(1, 10)
    ##             #self.splitter.setSizes([1, 100])
    ##     super(GridGUI, self).resizeEvent(event)

    def printPickedAtoms(self, selections, operation):
        msg = None
        if selections.nbAtoms():
            msg = 'you picked atom(s): '
            for sel in selections:
                mol = sel.getAtomGroup().getMolecule()
                for atom in sel:
                    msg += mol.atomFullName(atom)+ ' '
        if msg:
            self.displayMsg(msg)

    def printPickedBonds(self, bonds, operation):
        msg = None
        if len(bonds):
            msg = 'you picked bond(s): '
            for mol, at1, at2 in bonds:
                msg += '%s - %s '%(mol.atomFullName(at1), mol.atomFullName(at2))
        if msg:
            self.displayMsg(msg)
            
    def displayMsg(self, msg):
        # Make sure that the message fits inside the self.statusbar._message label and does not stretch the widget. Too long messages would cause Seg fault.
        fm = QtGui.QFontMetrics(self.statusbar._message.font())
        msgW = fm.width(msg)
        labelW = self.statusbar._message.width()
        if msgW > labelW:
            #truncate the message, put "... " at the end
            ff = float(msgW)/labelW
            nchar4label = int(len(msg)/ff)-5
            msg = msg[:nchar4label]+"... "
        self.statusbar.showMessage(msg, timeout=5000)

    def pickFirstCovalentBondAtom_cb(self):
        #print "pickFirstCovalentBondAtom_cb"
        w = self.covLigAtomsDialog
        w.pickCovAt2w.setChecked(False)
        if w.pickCovAt1w.isChecked():
            #print 'pick cov1'
            self.viewer.processPicking = self.pickFirstCovalentBondAtom
        else:
            #print 'pick default'
            self.viewer.processPicking = self.defaultProcessPicking

    def pickSecondCovalentBondAtom_cb(self):
        #print "pickSecondCovalentBondAtom_cb"
        w = self.covLigAtomsDialog
        w.pickCovAt1w.setChecked(False)
        if w.pickCovAt2w.isChecked():
            self.viewer.processPicking = self.pickSecondCovalentBondAtom
        else:
            self.viewer.processPicking = self.defaultProcessPicking
            #print 'pick default'

    def pickFirstCovalentBondAtom(self, pick):
        selections = self.viewer.findPickedAtoms(pick)
        if selections.nbAtoms()==1:
            mol = selections[0].getAtomGroup().getMolecule()
            atname = [mol.atomFullName(x) for x in selections[0]][0]
            ag = selections[0]._ag
            self._covalentBondAtoms[0] = Atom(ag, selections[0].getIndices()[0], 0)
            # check that the 2 atoms share a bond
            # if not do not set 1 and display error message and return
            if self._covalentBondAtoms[1] is not None:
                # check that they are bonded
                bonded = False
                for neighborAtom in self._covalentBondAtoms[1].iterBonded():
                    if neighborAtom == self._covalentBondAtoms[0]:
                        bonded = True
                        break
                if not bonded:
                    reply = self.showNoBondDialog(mol, self._covalentBondAtoms[0], self._covalentBondAtoms[1])
                    if reply == QtGui.QMessageBox.No:
                        self.statusbar.showMessage(
                            'atom %s does not share a bond with atom %s'%(
                                atname, mol.atomFullName(self._covalentBondAtoms[1])),
                            level='error', timeout=3000)
                        self._covalentBondAtoms[0] = None
                        return
            self.covBondAtom1Sphere.Set(vertices=[self._covalentBondAtoms[0].getCoords()], visible=True)
            w = self.covLigAtomsDialog
            w.covAt1w.setText(atname)
            self.statusbar.showMessage(
                'atom %s set at first covalent bond atom'%atname)#, timeout=3000)

            if self._covalentBondAtoms[1] is None:
                w.pickCovAt2w.click()
            self.getCovalentBondAtoms()
            self.updateStatusBar()
            if self._covalentBondAtoms[0] is not None and self._covalentBondAtoms[1] is not None:# and self.paramsWidget._covbondok:
                self.checkExcludedCovalenLigAtoms()
                self.viewer.processPicking = self.addRemoveCovLigAtoms 

    def pickSecondCovalentBondAtom(self, pick):
        selections = self.viewer.findPickedAtoms(pick)
        if selections.nbAtoms()==1:
            mol = selections[0].getAtomGroup().getMolecule()
            atname = [mol.atomFullName(x) for x in selections[0]][0]
            ag = selections[0]._ag
            self._covalentBondAtoms[1] = Atom(ag, selections[0].getIndices()[0], 0)
            self.covBondAtom2Sphere.Set(vertices=[self._covalentBondAtoms[1].getCoords()], visible=True)
            # check that the 2 atoms share a bond
            # if not do not set 1 and display error message and return
            if self._covalentBondAtoms[0] is not None:
                # check that they are bonded
                bonded = False
                for neighborAtom in self._covalentBondAtoms[1].iterBonded():
                    if neighborAtom == self._covalentBondAtoms[0]:
                        bonded = True
                        break
                if not bonded:
                    # show a dialog to ask the user if the second atom should be used as covalent bond:
                    reply = self.showNoBondDialog(mol, self._covalentBondAtoms[0], self._covalentBondAtoms[1])
                    if reply == QtGui.QMessageBox.No:
                        # would be nice to put text in RED with warning icon
                        # in status bar
                        self.statusbar.showMessage(
                            'atom %s does not share a bond with atom %s'%(
                                atname, mol.atomFullName(self._covalentBondAtoms[1])), level='error', timeout=3000)
                        self._covalentBondAtoms[1] = None
                        self.covBondAtom2Sphere.Set(visible = False)
                        return
            w = self.covLigAtomsDialog
            w.covAt2w.setText(atname)
            self.statusbar.showMessage(
                'atom %s set at second covalent bond atom'%atname)#, timeout=3000)

            if self._covalentBondAtoms[0] is None:
                w.pickCovAt1w.click()
            else:
                w.pickCovAt2w.setChecked(False)
                self.viewer.processPicking = self.defaultProcessPicking
            self.getCovalentBondAtoms()
            self.updateStatusBar()
            if self._covalentBondAtoms[0] is not None and self._covalentBondAtoms[1] is not None :#and self.paramsWidget._covbondok:
                self.checkExcludedCovalenLigAtoms()
                self.viewer.processPicking = self.addRemoveCovLigAtoms
            #print "DONE picking second atom"

    def showNoBondDialog(self, mol, at1, at2):
        # show a dialog to ask the user if the second atom should be used as covalent bond:
        # find the distance between he atoms(to display it in the dialog)
        p1 = at1.getCoords()
        p2 = at2.getCoords()
        import math
        distance = math.sqrt(pow(p2[0]-p1[0], 2)+pow(p2[1]-p1[1], 2)+pow(p2[2]-p1[2], 2))
        msg =  "<p>atom %s does not share a bond with atom %s.</p>" \
            "<p>The distance between atoms is %.2f. Do you want to use this atom in covalent docking?</p>" %(mol.atomFullName(at1), mol.atomFullName(at2), distance)
        return QtGui.QMessageBox.question(self, "Covalent docking", msg, QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
        

    def getFirstCovalentBondAtom_cb(self):
        pass
    
    def firstCovalentBondAtomChanged_cb(self):
        #print "firstCovalentBondAtomChanged_cb"
        w = self.covLigAtomsDialog
        at1Str = w.covAt1w.text().encode('ascii', 'replace')
        if not len(at1Str):
           self._covalentBondAtoms[0] = None
           w.pickCovAt1w.setChecked(True)
           if w.pickCovAt2w.isChecked():
               w.pickCovAt2w.setChecked(False)
           self.viewer.processPicking = self.pickFirstCovalentBondAtom
           self.covBondAtom1Sphere.Set(visible=False)
           self.clearCovDockingExcludedAtoms()
           self.updateStatusBar() 

    def getSecondCovalentBondAtom_cb(self):
        pass

    def secondCovalentBondAtomChanged_cb(self):
        w = self.covLigAtomsDialog
        at2Str = w.covAt2w.text().encode('ascii', 'replace')
        #print "secondCovalentBondAtomChanged_cb"
        if not len(at2Str):
           self._covalentBondAtoms[1] = None
           w.pickCovAt2w.setChecked(True)
           if w.pickCovAt1w.isChecked():
               w.pickCovAt1w.setChecked(False)
           self.viewer.processPicking = self.pickSecondCovalentBondAtom
           self.covBondAtom2Sphere.Set(visible=False)
           self.clearCovDockingExcludedAtoms()
           self.updateStatusBar()
           
    def getLimitResName_cb(self):
        w = self.covLigAtomsDialog
        resName = w.covalentResEdit.text().encode('ascii', 'replace')
        _resName = resName.upper()
        if  resName != _resName : w.covalentResEdit.setText(_resName)
        if not len(resName): 
            self.clearCovDockingExcludedAtoms()
            return
        #try:
        #    residues = flexResStr2flexRes(_resName)
        #except:
        #    reply = QtGui.QMessageBox.critical(self, "Invalid residue string",
        #             'Residue string should be in the form:  chid1:res1,res2..;chid2:res1,..' ,
        #                                       QtGui.QMessageBox.Ok)
        #    w.covalentResEdit.setText("")
        #    return

        selStr, flexResAtoms  = self.parseResString(_resName)
        if not len(selStr):
            w.covalentResEdit.setText("")
            return
        self.hideCovLigandAtoms()
        self._covalentLigAtoms = None
        self._newExcluded = []
        self.covLigAtomsDialog.resListWidget.clear()
        self.getCovalentBondAtoms()
        self.updateStatusBar()
        if self._covalentBondAtoms[0] is not None and self._covalentBondAtoms[1] is not None :#and self.paramsWidget._covbondok:
            self.checkExcludedCovalenLigAtoms()
            self.viewer.processPicking = self.addRemoveCovLigAtoms

    def setCovalentDockingMode(self, val):
        # callback of the [covalent docking] group check button
        #import pdb; pdb.set_trace()
        self._covalentDocking = val
        self.agfr.covalentBond = val
        self.showCovBondAtoms.setDisabled(not val)
        if  self._covalentBondAtoms[0] is not None and self._covalentBondAtoms[1] is not None:
            self.covBondAtom1Sphere.Set(visible=val)
            self.covBondAtom2Sphere.Set(visible=val)
            if self._covalentLigAtoms is not None:
                if val:
                    self.showCovLigandAtoms()
                else: self.hideCovLigandAtoms()
        if val: # show dialog
            self.freezeUI(True)
            if not self.covLigAtomsDialog:
                self.covLigAtomsDialog = w = CovalentLigAtomsList()
                w.resSelected.connect(self.showCovLigandAtoms_cb)
                w.oncloseSignal.connect(self.covLigAtomsDialogClosed)
                w.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
                w.app = weakref.ref(self)
                w.pickCovAt1w.clicked.connect(self.pickFirstCovalentBondAtom_cb)
                #w.covAt1w  LineEdit widget
                w.covAt1w.returnPressed.connect(self.getFirstCovalentBondAtom_cb)
                w.covAt1w.textChanged.connect(self.firstCovalentBondAtomChanged_cb)
                w.pickCovAt2w.clicked.connect(self.pickSecondCovalentBondAtom_cb)
                #w.covAt1w. LineEdit widget
                w.covAt2w.returnPressed.connect(self.getSecondCovalentBondAtom_cb)
                w.covAt2w.textChanged.connect(self.secondCovalentBondAtomChanged_cb)

                w.covalentResEdit.returnPressed.connect(self.getLimitResName_cb)
                self._openDialogs.append(w)                
            else:
                w = self.covLigAtomsDialog
            if self._covalentBondAtoms[0] is None:
                w.pickCovAt1w.setChecked(True)
                if w.pickCovAt2w.isChecked():
                    w.pickCovAt2w.setChecked(False)
                self.viewer.processPicking = self.pickFirstCovalentBondAtom
            elif self._covalentBondAtoms[1] is None:
                w.pickCovAt2w.setChecked(True)
                if w.pickCovAt1w.isChecked():
                    w.pickCovAt1w.setChecked(False)
                self.viewer.processPicking = self.pickSecondCovalentBondAtom
            self.covLigAtomsDialog.show()
            self.paramsWidget.computeGridsButton.setDisabled(True)
        else: # remove dialog
            self.freezeUI(False)
            if self.covLigAtomsDialog and self.covLigAtomsDialog.isVisible():
                self.covLigAtomsDialog.reject()
            else:
                self.covLigAtomsDialog = None
            self._covalentBondAtoms[0] = None
            self._covalentBondAtoms[1] = None
            self._torsionAtom = None
            self.hideCovLigandAtoms()
            self._covalentLigAtoms = None
            self._newExcluded = []
            self.viewer.processPicking = self.defaultProcessPicking
            self.paramsWidget.covAt1Lab.setText('       ')
            self.paramsWidget.covAt2Lab.setText('      ')
            self.paramsWidget.computeGridsButton.setDisabled(False)
            self._openDialogs.remove(self.covLigAtomsDialog)
        self.updateStatusBar()

    def checkExcludedCovalenLigAtoms(self):
        #print "checkExcludedCovalenLigAtoms"
        w = self.covLigAtomsDialog
        self.hideCovLigandAtoms()
        w.resListWidget.clear()
        self._covalentLigAtoms = None
        self._newExcluded = []
        #self.covLigAtomsDialog.covalentResEdit.setText("")
        self.covLigAtomsDialog.resListWidget.clear()
        if len(self.agfr.covalentBondToExclude):
            #import pdb; pdb.set_trace()
            # do not add self._covalentBondAtoms to the list of excluded atoms
            sel = self.receptor.select("index %d %d" %(self._covalentBondAtoms[0].getIndex(), self._covalentBondAtoms[1].getIndex()))
            self.agfr.covalentBondToExclude = self.agfr.covalentBondToExclude - sel
            if len(self.agfr.covalentBondToExclude):
                w.fill(self.agfr.covalentBondToExclude)
                resStr, atmsList, atinds = self.covLigAtomsDialog.getAtomList()
                self._covalentLigAtoms = Selection(self.receptor._ag, atinds, '')
                self._newExcluded = [set(atinds)]
                self.showCovLigandAtoms()
                w.pickCovAt1w.setChecked(False)
                w.pickCovAt2w.setChecked(False)
                self.viewer.processPicking = self.addRemoveCovLigAtoms
            
    def covLigAtomsDialogClosed(self, res):
        self.paramsWidget.covDockingGroupBox.setDisabled(False)
        if res: #OK
            #import pdb; pdb.set_trace()
            if self._covalentBondAtoms[0] and self._covalentBondAtoms[1]:
                resStr, resList, atinds = self.covLigAtomsDialog.getAtomList()
                if len(atinds):
                    self._covalentLigAtoms = Selection(self.receptor._ag, atinds, '')
                    #toRemoveAtoms = self._covalentLigAtoms[0]
                else:
                    self._covalentLigAtoms = None
                    #toRemoveAtoms = []
                #self._covalentLigAtoms = self.selectionFromResString(resList)
                #for sel in self._covalentLigAtoms:
                #    self.receptor._ag._flags['stippled'][sel.getIndices()] = True
                #self.receptor.app().displayLines.refreshDisplay(self.receptor)
                #print "COV LIG ATOMS DIALOG:", resList
                w = self.paramsWidget
                w.covAt1Lab.setText(self.covLigAtomsDialog.covAt1w.text())
                w.covAt2Lab.setText(self.covLigAtomsDialog.covAt2w.text())
                #at1 = self._covalentBondAtoms[0].getSerial()
                #at2 = self._covalentBondAtoms[1].getSerial()
                #torAt = self._torsionAtom.getSerial()
                
                #self.agfr.setCovalentDocking(torAt, [at1, at2], toRemoveAtoms=toRemoveAtoms)
            else:
                self.agfr.covalentBond = False
                self.covLigAtomsDialog.clearForm()
                self.paramsWidget.computeGridsButton.setDisabled(False)
                self.paramsWidget.covDockingGroupBox.setChecked(False)
                #self.clearCovDockingExcludedAtoms()
                self.paramsWidget.covAt1Lab.setText('       ')
                self.paramsWidget.covAt1Lab.setText('       ')
        else: #cancel
            # clear the text edit areas for atom 1 and 2
            self.agfr.covalentBond = False
            self.agfr.data['covalentBondAtom1'] = None
            self.agfr.data['covalentBondAtom2'] = None
            self.agfr.data['covalentBondTorsionAtom'] = None
            self._covalentBondAtoms[0] = None
            self._covalentBondAtoms[1] = None
            self.paramsWidget.covAt1Lab.setText('       ')
            self.paramsWidget.covAt1Lab.setText('       ')
            self.covLigAtomsDialog.clearForm()
            self._torsionAtom = None
            self.hideCovLigandAtoms()
            self._covalentLigAtoms = None
            self._newExcluded = []
            self.paramsWidget.computeGridsButton.setDisabled(False)
            self.paramsWidget.covDockingGroupBox.setChecked(False)
            #self.covBondAtom1Sphere.Set(visible=False)
            #self.covBondAtom2Sphere.Set(visible=False)
            #self.updateStatusBar()
        self.viewer.processPicking = self.defaultProcessPicking

    def clearCovDockingExcludedAtoms(self):
        #self.agfr.data['covalentBond'] = []
        self.agfr.data['covalentRes'] = []
        self.agfr.covalentBondToExclude = []
        self.hideCovLigandAtoms()
        self._covalentLigAtoms = None
        self._newExcluded = []
        self.covLigAtomsDialog.resListWidget.clear()

    @waiting_effects
    def getCovalentBondAtoms(self):
        #print "getCovalentBondAtoms"
        #import pdb; pdb.set_trace()
        if self._covalentDocking:
            self.agfr.covalentBond = True
            if self._covalentBondAtoms[0] and self._covalentBondAtoms[1]:
                self.agfr.covalentBond = True
                at1 = self._covalentBondAtoms[0].getSerial()
                at2 = self._covalentBondAtoms[1].getSerial()
                torAt = None
                for at in self._covalentBondAtoms[0].iterBonded():
                    if at != self._covalentBondAtoms[1] and at.getElement() != "H":
                        self._torsionAtom = at
                        torAt = at.getSerial()
                        break
                if torAt:
                    resStr = str(self.covLigAtomsDialog.covalentResEdit.text())
                    if not len(resStr): resStr=None
                    self.agfr.setCovalentDocking(torAt, [at1, at2], resStr)
                else:
                    # FIX ME: invalid atom1, clear the form , issue Error in the status bar
                    at1Name = self.receptor.atomFullName(self._covalentBondAtoms[0])
                    reply = QtGui.QMessageBox.critical(self, "Invalid first atom",
                'invalid atom1 %s. Can not find the receptor atom used to compute the torsion angle for the covalent bond.'% at1Name ,
                QtGui.QMessageBox.Cancel | QtGui.QMessageBox.Retry)
                    if reply == QtGui.QMessageBox.Cancel:
                        self.setCovalentDockingMode(False)
                    else:
                        self.covLigAtomsDialog.covAt1w.setText("")
                        self.clearCovDockingExcludedAtoms()
            else:
                self.clearCovDockingExcludedAtoms()
        else:
            self.agfr.covalentBond = False
            self.clearCovDockingExcludedAtoms()
            
    def selectionFromResString(self, resStrList):
        #resStrList looks like this: [A:CYS164:CA,C,O;A:ASP165,...,..."
        #import pdb;pdb.set_trace()
        _covalentLigAtoms = []
        chDict = {}
        for resStr in resStrList:
            chid, resstr, atoms = resStr.split(":")
            if not chDict.has_key(chid):
                chDict[chid] = {}
            if not chDict[chid].has_key(resstr):
                chDict[chid][resstr] = ""
            atnames = ""
            for at in atoms.replace(',', ' ').split():
                atnames += at+ " "
            chDict[chid][resstr]+=atnames
        for chid in  chDict.keys():
            for resstr in chDict[chid].keys():
                atnames = chDict[chid][resstr].strip()
                #resnum = "".join([ss for ss in resstr if ss.isdigit()])
                resnum = []
                for i in range(len(resstr)-1, -1, -1):
                    if resstr[i].isdigit(): resnum.insert(0, resstr[i])
                    else: break
                selstr = "chain %s resnum `%s` name %s" %(chid, "".join(resnum), atnames)
                _covalentLigAtoms.append(self.receptor.select(selstr))
        return _covalentLigAtoms

    def addRemoveCovLigAtoms(self, pick):
        #import pdb; pdb.set_trace()
        selections = self.viewer.findPickedAtoms(pick)
        if selections.nbAtoms()!=1: return
        mol = selections[0].getAtomGroup().getMolecule()
        atname = [mol.atomFullName(x) for x in selections[0]][0]
        ag = selections[0]._ag
        at = Atom(ag, selections[0].getIndices()[0], 0)
        atInd = at.getIndex()
        if atInd == self._covalentBondAtoms[0].getIndex() or atInd == self._covalentBondAtoms[1].getIndex():return 
        excluded = []
        if self._covalentLigAtoms is not None:
            if len(self._covalentLigAtoms):
                excluded.extend(self._covalentLigAtoms.getIndices().tolist())
        excluded = set(excluded)
        if atInd in excluded:
            self.removeLigAtom(atInd, excluded)
        else:
            self.addLigAtom(at, excluded)

    @waiting_effects        
    def addLigAtom(self, at, excluded):
        newList = set([])
        aInd = at.getIndex()
        #import pdb; pdb.set_trace()
        if at.numBonds()==0:
            newList = set([aInd])
        else:
            for natom in at.iterBonded():
                naInd = natom.getIndex()
                #print "2 atoms:", aInd, naInd, at.getName(), natom.getName()
                if naInd not in excluded and naInd not in newList:
                    inSubTree = self.receptor._traverse(
                        at, natom, [aInd, naInd])
                    #print "inSubTree:", inSubTree
                    newList = newList | set(inSubTree)
                #print "newList:", newList
        if len(newList):
            at1, at2 = self._covalentBondAtoms
            inds = set([at1.getIndex(), at2.getIndex()])
            newList = newList-inds-excluded
            newExcluded = Selection(self.receptor._ag,
                                                    list(newList), "")
            self._newExcluded.append(newList)
            self.covLigAtomsDialog.addItems(newExcluded)
            resStr, atmsList, atinds = self.covLigAtomsDialog.getAtomList()
            self._covalentLigAtoms = Selection(self.receptor._ag, atinds, '')
            self.showCovLigandAtoms()

    @waiting_effects
    def removeLigAtom(self, atInd, excluded):
        if not self._covalentLigAtoms:
            return
        atInSet = None
        for i, atset in enumerate(self._newExcluded):
            if atInd in atset:
                atInSet = i
                break
        if atInSet == None: return
        setToRemove = self._newExcluded.pop(atInSet)
        #print "removeLigAtom:", "setToRemove:", setToRemove,
        #print "_newExcluded:", self._newExcluded
        newExcluded = []
        self.covLigAtomsDialog.resListWidget.clear()
        self.covLigAtomsDialog.items = {}
        map(newExcluded.extend, self._newExcluded)
        if len(newExcluded):
            self.covLigAtomsDialog.addItems(Selection(self.receptor._ag, list(newExcluded), ""))
            #unlabel all covligatoms
            if len(self._covalentLigAtoms):
                self.receptor.app().unlabelAtoms.updateModel(self._covalentLigAtoms)
                
            resStr, atmsList, atinds = self.covLigAtomsDialog.getAtomList()
            #self._covalentLigAtoms = self.selectionFromResString(atmsList)
            self._covalentLigAtoms = Selection(self.receptor._ag, atinds, '')
            # show new atoms
            self.showCovLigandAtoms()
        else:
            self.hideCovLigandAtoms()
            self._covalentLigAtoms = None
            
    def showCovLigandAtoms_cb(self, selstr, display):
        # callback of the covLigAtomsDialog list widget and "All" and "None"
        # buttons
        #import pdb; pdb.set_trace()
        #print "showCovLigandAtoms_cb",
        resStr, atmsList, atinds = self.covLigAtomsDialog.getAtomList()
        if not len(atmsList):
            self.hideCovLigandAtoms()
            self._covalentLigAtoms = None
        else:
            self._covalentLigAtoms = Selection(self.receptor._ag, atinds, '')
            if display == False:
                selList = self.selectionFromResString(selstr.split(";"))
                for sel in selList:
                    self.receptor.app().unlabelAtoms.updateModel(sel)
            self.showCovLigandAtoms()

    def showCovLigandAtoms(self):
        #import pdb; pdb.set_trace()
        if not self._covalentLigAtoms: return
        coords = []
        self.receptor._ag._flags['stippled'][:] = False
        if len(self._covalentLigAtoms):
            self.receptor._ag._flags['stippled'][self._covalentLigAtoms.getIndices()] = True
            coords.extend(self._covalentLigAtoms.getCoords().tolist())
            self.receptor.app().labelAtoms.updateModel(self._covalentLigAtoms, propName='name')
        self.receptor.app().displayLines.refreshDisplay(self.receptor)
        gc = self.receptor.geomContainer
        geom = gc.geoms[self.receptor.app().labelAtoms._atomSetName]
        self.atomsToIgnoreSpheres.Set(vertices=coords, visible=True)
        geom.Set(translation=[[0.3, 0.3, 0.3]])
        self.receptor.app().labelAtoms.refreshDisplay(self.receptor)

    def hideCovLigandAtoms(self):
        if not self._covalentLigAtoms: return
        if len(self._covalentLigAtoms): 
            self.receptor._ag._flags['stippled'][self._covalentLigAtoms.getIndices()] = False
            self.receptor.app().unlabelAtoms.updateModel(self._covalentLigAtoms)
        self.receptor.app().displayLines.refreshDisplay(self.receptor)
        gc = self.receptor.geomContainer
        geom = gc.geoms[self.receptor.app().labelAtoms._atomSetName]
        self.atomsToIgnoreSpheres.Set(visible=False)
        geom.Set(translation=[[0., 0., 0.]])
        self.receptor.app().unlabelAtoms.refreshDisplay(self.receptor)
        
    def createPmvApp(self, eventHandler=None):
        pmv = MolApp(eventHandler=eventHandler)
        pmv.trapExceptions = False
        #pmv.lazyLoad('bondsCmds', package='PmvApp')
        #pmv.lazyLoad('fileCmds', package='PmvApp')
        #pmv.lazyLoad('displayCmds', package='PmvApp')
        #pmv.lazyLoad('editCmds', package='PmvApp')
        #pmv.lazyLoad("colorCmds", package="PmvApp")
        #pmv.lazyLoad("selectionCmds", package="PmvApp")
        #pmv.lazyLoad('deleteCmds', package='PmvApp')
        #pmv.lazyLoad('labelCmds', package='PmvApp')

        #pmv.lazyLoad('msmsCmds', package='PmvApp')
        #pmv.lazyLoad('displayHyperBallsCmds', package='PmvApp')
        #pmv.lazyLoad('interactionsCmds', package='PmvApp')
        #pmv.lazyLoad('coarseMolecularSurfaceCmds', package='PmvApp')

        pmv.setOnAddObjectCmd('Molecule', [pmv.displayLines, pmv.colorByAtomType])
        self.pmv = pmv

    def setGridVisible(self, value):
        if not self.boxGeom: return
        # value is 0 for unchecked and 2 for checked for checkbox
        # not(value==0) make it work for 0, 1, 2, False, True
        self.boxGeom.master.Set(visible = not(value==0))
        for c in self.boxGeom.master.children:
            if c.name=='faces':
                c.Set(visible = 0)
            else:
                c.Set(visible = not(value==0))

    def freezeUI(self, freeze):
        if freeze:
            self.stack.setDisabled(True)
        else:
            self.stack.setDisabled(False)

    def buildUI(self, pmvViewer=None, parent=None):

        self.topBar = TopBar(self)

        self.stack = QtGui.QStackedWidget(self)
        w = self.paramsWidget = ADFRGridMapParametersWidget(self, parent)
        self.stack.addWidget(w)
        self.topBar.loadRecButton.clicked.connect(self.getReceptorFilename)
        self.topBar.loadTargetButton.clicked.connect(self.getReceptorMapsFilename)
        self.topBar.loadLigButton.clicked.connect(self.getLigandFilename)
        self.topBar.setFolderButton.clicked.connect(self.setFolder)
        self.topBar.fetchPDBButton.clicked.connect(self.fetchPDB)
        self.topBar.userPreferencesButton.clicked.connect(self.showUserPreferences)
        w.editAtypesButton.clicked.connect(self.pickATypes)
        w.computeGridsButton.clicked.connect(self.computeGrids)
        w.computeTPointsButton.clicked.connect(self.computeTPoints)
        w.clustersWidget.itemClicked.connect(self.showHideCluster)
        #w.featurePtsButton.clicked.connect(self.computeFPoints)
        w.flexResWidget.returnPressed.connect(self.setFlexRes)
        w.setFlexResButton.clicked.connect(self.showFlexResChooser)
        #w.setFlexResButton.clicked.connect(self.showFlexResChooserFromFP)
        w.setFlexResGrdButton.clicked.connect(self.setGridFlexRes)
        
        w.covDockingGroupBox.toggled.connect(self.setCovalentDockingMode)

        w.gridPaddingWidget.valueChanged.connect(self.paddingChanged)
        w.compAllButton.clicked.connect(self.toggleComputeAll)
        w.recGridAct.triggered.connect(self.setGridFullReceptor)
        w.ligGridAct.triggered.connect(self.setGridFullLigand)
        w.TPGridAct.triggered.connect(self.setGridClusterPoints)
        w.residuesGrdAct.triggered.connect(self.showResiduesGridControls)
        w.manualGrdAct.triggered.connect(self.showManualGridControls)
        w.selectAll.clicked.connect(self.selectAllPockets)
        w.deselectAll.clicked.connect(self.deselectAllPockets)
        w.invertSelection.clicked.connect(self.invertPocketSelection)

        # AutoSite params button:
        w.autoSiteVButton.autosite2.connect(self.setAutoSiteVersion)
        w.autoSiteVButton.ligSizeChanged.connect(self.setAutoSiteLigSize)
        w.autoSiteVButton.usePepScore.connect(self.setAutoSitePepScoreFunc)

        # water map settings button:
        w.wMapParamB.entropyChanged.connect(self.setWaterMapEntropy)
        w.wMapParamB.weightChanged.connect(self.setWaterMapWeight)
        ## if sys.platform=='darwin':
        ##     self.paramsWidget.setMaximumWidth(450)
        ## elif sys.platform=='linux2':
        ##     self.paramsWidget.setMaximumWidth(350)
        ## elif sys.platform=='win32':
        ##     self.paramsWidget.setMaximumWidth(450)
        ## else:
        ##     self.paramsWidget.setMaximumWidth(400)
        #self.TPointsOK.connect(w.handleTPointsSignal)
        #self.gridOK.connect(w.handleGridOKSignal)
        
        ## create the vertical tool bar
        ##
        self.toolBarV = QtGui.QToolBar(self)
        self.toolBarV.setOrientation(QtCore.Qt.Vertical)

        act = self.focusBoxAction = QtGui.QAction(
            QtGui.QIcon(os.path.join(ICONPATH, 'focus.png')),
            'center view on docking box', self)
        act.setStatusTip('focus on docking box')
        act.triggered.connect(self.focusOnGrid)
        self.toolBarV.addAction(act)

        self.toolBarV.addSeparator()

        act = self.showReceptorAction = QtGui.QAction(
            QtGui.QIcon(os.path.join(ICONPATH, 'receptor2_NO.png')),
            'show receptor', self)
        act.setStatusTip('show/hide receptor')
        act.triggered.connect(self.showHideReceptor)
        act.setCheckable(True)
        act.setDisabled(True)
        self.toolBarV.addAction(act)
        
        act = self.showReceptorSurfaceAction = QtGui.QAction(
            QtGui.QIcon(os.path.join(ICONPATH, 'ms.png')),
            'show receptor surface', self)
        act.setStatusTip('show/hide receptor surface')
        act.triggered.connect(self.showHideReceptorSurface)
        act.setCheckable(True)
        act.setDisabled(True)
        self.toolBarV.addAction(act)
        
        act = self.showResLabAction = QtGui.QAction(
            QtGui.QIcon(os.path.join(ICONPATH, 'labels.png')),
            'show receptor residue labels', self)
        act.setStatusTip('show/hide labels ')
        act.triggered.connect(self.showHideResidueLabels)
        act.setCheckable(True)
        act.setDisabled(True)
        self.toolBarV.addAction(act)

        self.toolBarV.addSeparator()
        
        act = self.editLigandAct = QtGui.QAction(
            QtGui.QIcon(os.path.join(ICONPATH, 'ligand.png')),
            '', self)
        act.setStatusTip('show ligand tools')
        act.toggled.connect(self.toggleLigandEdit)
        act.setCheckable(True)
        act.setDisabled(True)
        self.toolBarV.addAction(act)
                
        # create drawer for ligand editing
        self.drawer = Drawer(self)
        self.drawer.raise_()
        #tb_size= 36
        self.ligand_tb = LigandToolBar(36, 'vertical', self)
        #self.ligand_tb = SmallToolBar(tb_size, 'vertical')
        content_lay = QtGui.QVBoxLayout()
        content_lay.setContentsMargins(0, 0, 0, 0)
        content_lay.addWidget(self.ligand_tb)
        self.drawer.setToolBar(self.ligand_tb)
        self.drawer.setLayout(content_lay)

        act = self.showCovBondAtoms = QtGui.QAction(
            QtGui.QIcon(os.path.join(ICONPATH, 'covalentBond.png')),
            'show/hide atoms forming covalent bond', self)
        act.setStatusTip('show/hide covalent bond atoms')
        act.triggered.connect(self.showHideCovBond)
        #act.setCheckable(True)
        #act.setChecked(False)
        act.setDisabled(True)
        self.toolBarV.addAction(act)

        self.toolBarV.addSeparator()

        act = self.showGridBoxAction = QtGui.QAction(
            QtGui.QIcon(os.path.join(ICONPATH, 'cube.png')),
            'show docking box and binding pockets fill points', self)
        act.setStatusTip('show/hide labels docking box')
        act.triggered.connect(self.showHideGridBox)
        act.setCheckable(True)
        act.setChecked(True)
        act.setEnabled(False)
        self.toolBarV.addAction(act)

        act = self.showTrgGuiAction = QtGui.QAction(
            QtGui.QIcon(os.path.join(ICONPATH,'iso.png')),
            'trg map gui', self)
        act.triggered.connect(self.openTrgMapGui)
        act.setDisabled(True)
        act.setCheckable(True)
        act.setChecked(False)
        self.toolBarV.addAction(act)

        act = self.showXtalMates = QtGui.QAction(
            QtGui.QIcon(os.path.join(ICONPATH,'xtal.png')),
            'show/hide Xtal mates', self)
        act.triggered.connect(self.toggleXtalMates)
        act.setDisabled(True)
        act.setCheckable(True)
        act.setChecked(False)
        self.toolBarV.addAction(act)
                
        ## act = self.editLigandAct = QtGui.QAction(
        ##     QtGui.QIcon(os.path.join(ICONPATH,'ligandTools.png')),
        ##     'editLigand', self)
        ## act.toggled.connect(self.toggleLigandEdit)
        ## act.setCheckable(True)
        ## act.setDisabled(True)
        ## act.setChecked(False)
        ## self.toolBarV.addAction(act)
                
        act = self.testDockAct = QtGui.QAction(
            QtGui.QIcon(os.path.join(ICONPATH,'ad4.png')),
            'test run', self)
        act.triggered.connect(self.testDocking)
        act.setDisabled(True)
        self.toolBarV.addAction(act)

        act = self.viewInteractionsAct = QtGui.QAction(
            QtGui.QIcon(os.path.join(ICONPATH,'score.png')),
            'intreactions panel', self)
        act.triggered.connect(self.viewInteractions)
        act.setCheckable(True)
        act.setChecked(False)
        act.setDisabled(True)
        self.toolBarV.addAction(act)

        ## act = self.showFPAction = QtGui.QAction(
        ##     "FP",
        ##     #QtGui.QIcon(os.path.join(ICONPATH, 'cube.png')),
        ##     #'show/hide feature points',
        ##     self)
        ## act.setStatusTip('show/hide feature points')
        ## act.triggered.connect(self.showHideFP)
        ## act.setCheckable(True)
        ## act.setDisabled(True)
        ## self.toolBarV.addAction(act)
        spacer = QtGui.QWidget(self)
        spacer.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        self.toolBarV.addWidget(spacer)

        self.toolBarV.addSeparator()

        act = self.showReport = QtGui.QAction(
            QtGui.QIcon(os.path.join(ICONPATH,'report.png')),
            'show report', self)
        act.triggered.connect(self.showHideReport)
        self.toolBarV.addAction(act)

        #act = self.nextFileAct = QtGui.QAction(
        #    QtGui.QIcon(os.path.join(ICONPATH,'next.png')),
        #    'show/hide setting', self)
        #act.triggered.connect(self.gotoNextFile)
        #act.setDisabled(True)
        #self.toolBarV.addAction(act)

        #act = self.showSettings = QtGui.QAction(
        #    QtGui.QIcon(os.path.join(ICONPATH,'cogwheel.png')),
        #    'show/hide setting', self)
        #act.triggered.connect(self.displaySettings)
        #act.setDisabled(True)
        #self.toolBarV.addAction(act)
        
        from PmvApp.GUI.Qt.PmvGUI import PmvViewer
        if pmvViewer is None:
            ## create the Pmv Viewer
            self.viewer = PmvViewer(self.pmv, master=self)
            # overwrite dragSelectFunction in viewer to selecting atoms
            def noSelect(*args, **kw):
                #print 'caught selection attempt'
                pass
            self.viewer.dragSelectFunc = noSelect
        else:
            assert isisntance(pmvViewer, PmvViewer)
            self.viewer = viewer
        self.viewer.cameras[0].setMinimumWidth(450)

        from DejaVu2.Points import Points
        self.TPoints = Points('tpoints', visible=0, inheritMaterial=False,
                              materials=[(1,0,0)], inheritPointWidth=False,
                              pointWidth=4.)
        self.viewer.AddObject(self.TPoints)
        self.boxGeom = None
        self.anchorAtomGeom = None
        self.rotBondsArcs = None
        
        # Spheres to display covalent bond atoms
        self.covBondAtom1Sphere = None
        self.covBondAtom2Sphere = None
        self.atomsToIgnoreSpheres = None

        # Lines and spheres to display gaps (missing residues)
        #self.missingResLines = None
        #self.missingResSpheres = None
        #self.missingResLabels = None

        # feature points spheres
        ## self.TSpheresC = Spheres('tsphC', visible=0, inheritMaterial=False,\
        ##                       materials=[(0,1,0)], inheritPointWidth=False,\
        ##                       pointWidth=4.)
        ## self.viewer.AddObject(self.TSpheresC)
        ## self.TSpheresO = Spheres('tsphO', visible=0, inheritMaterial=False,\
        ##                       materials=[(0,0,1)], inheritPointWidth=False,\
        ##                       pointWidth=4.)
        ## self.viewer.AddObject(self.TSpheresO)
        ## self.TSpheresH = Spheres('tsphH', visible=0, inheritMaterial=False,\
        ##                       materials=[(1,1,1)], inheritPointWidth=False,\
        ##                       pointWidth=4.)
        ## self.viewer.AddObject(self.TSpheresH)

        # place widgets
        topLayout = QtGui.QVBoxLayout()
        topLayout.setContentsMargins(LeftMargin,topMargin,rightMargin,bottomMargin)
        topLayout.setSpacing(spacing)
        
        mainLayout = QtGui.QHBoxLayout()
        mainLayout.setContentsMargins(LeftMargin,topMargin,rightMargin,bottomMargin)
        mainLayout.setSpacing(spacing)
        topLayout.addLayout(mainLayout)
        self.splitter = splitter = QtGui.QSplitter()
        mainLayout.addWidget(splitter)
        
        #splitter.addWidget(self.paramsWidget)
        frame = QtGui.QFrame()
        frame.setContentsMargins(0,0,0,0)
        vlayout = QtGui.QVBoxLayout()
        frame.setLayout(vlayout)
        vlayout.setContentsMargins(0,0,0,0)
        splitter.addWidget(frame)
        vlayout.addWidget(self.topBar)
        vlayout.addWidget(self.stack)
        #self.stack.setCurrentWidget(self.paramsWidget)
        
        self.viewWidget = QtGui.QWidget()
        viewLayout = QtGui.QHBoxLayout()

        # view layout contains the vertical tool bar and the Camera
        viewLayout.setContentsMargins(LeftMargin,topMargin,rightMargin,bottomMargin)
        viewLayout.setSpacing(spacing)
        
        self.viewWidget.setLayout(viewLayout)
        viewLayout.addWidget(self.toolBarV)

        # vlayout contains the camera and optionally a table below
        self.vlayout = QtGui.QVBoxLayout()
        self.vlayout.setContentsMargins(0,0,0,0)
        viewLayout.addLayout(self.vlayout)
        self.vsplitter = vsplitter = QtGui.QSplitter(QtCore.Qt.Vertical)
        vsplitter.setContentsMargins(0,0,0,0)
        self.vlayout.addWidget(vsplitter)
        vsplitter.addWidget(self.viewer.cameras[0])
        
        self.pairsTable = PairwiseScoreTable()
        self.pairsTable.setMaximumHeight(0)
        self.pairsTable.table.itemSelectionChanged.connect(self.onSelectionChange)
        self.pairsTable.tableUpdated.connect(self.onTablechange)
        self._interCyl = Cylinders('interCyl', visible=0, inheritMaterial=False)
        self._interSph = Spheres('interSph', visible=0, inheritMaterial=False)
        self.viewer.AddObject(self._interCyl)
        self.viewer.AddObject(self._interSph)
        vsplitter.addWidget(self.pairsTable)
        
        splitter.addWidget(self.viewWidget)
        
        #splitter.setStretchFactor(0, 0.1)
        #splitter.setStretchFactor(1, 0.9)# needs an integer
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 5)

        vsplitter.setStretchFactor(0, 5)
        vsplitter.setStretchFactor(1, 3)
        
        # MS 8/25/16 commented out because it triggered a warning from QT
        # QLayout::addChildLayout: layout "" already has a parent NOT SURE WHY:(
        #mainLayout.addLayout(viewLayout)

        #self.statusbar = bar = QtGui.QStatusBar(self)
        self.statusbar = bar = MyStatusBar(self)
        bar.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Fixed)
        topLayout.addWidget(bar)
        #self.statusbar.showMessage("load a receptor, define a docking box, select translational points to enable maps calculations")
        self.statusbar.showMessage("Welcome to AGFR gui, Please load a receptor")
        self.boxOkLabel = label = QtGui.QLabel("Box")
        label.setFrameStyle(QtGui.QFrame.StyledPanel | QtGui.QFrame.Sunken)
        sample_palette = QtGui.QPalette() 
        sample_palette.setColor(QtGui.QPalette.Window, QtCore.Qt.red)
        label.setAutoFillBackground(True)
        label.setPalette(sample_palette)
        #label.show()
        #layout.addWidget(label)
        bar.addPermanentWidget (label)

        self.tpOkLabel = label1 = QtGui.QLabel("Pocket")
        label1.setFrameStyle(QtGui.QFrame.StyledPanel | QtGui.QFrame.Sunken)
        sample_palette = QtGui.QPalette() 
        sample_palette.setColor(QtGui.QPalette.Window, QtCore.Qt.gray)
        label1.setAutoFillBackground(True)
        label1.setPalette(sample_palette)
        bar.addPermanentWidget (label1)

        self.frOkLabel = label2 = QtGui.QLabel("FR")
        label2.setFrameStyle(QtGui.QFrame.StyledPanel | QtGui.QFrame.Sunken)
        sample_palette = QtGui.QPalette() 
        sample_palette.setColor(QtGui.QPalette.Window, QtCore.Qt.gray);
        label2.setAutoFillBackground(True)
        label2.setPalette(sample_palette)
        bar.addPermanentWidget(label2)

        self.covLabel = label3 = QtGui.QLabel("Cov")
        label3.setFrameStyle(QtGui.QFrame.StyledPanel | QtGui.QFrame.Sunken)
        sample_palette = QtGui.QPalette() 
        sample_palette.setColor(QtGui.QPalette.Window, QtCore.Qt.gray);
        label3.setAutoFillBackground(True)
        label3.setPalette(sample_palette)
        bar.addPermanentWidget(label3)
        self.setLayout(topLayout)

    def placeDrawer(self):
        rel_pos = self.toolBarV.actionGeometry(self.editLigandAct).topLeft()
        pos = self.toolBarV.mapToGlobal(rel_pos)
        self.drawer.move(pos.x()+ self.toolBarV.sizeHint().width(), pos.y())
        
    def resizeEvent(self, event):
        self.placeDrawer()
        super(GridGUI, self).resizeEvent(event)

    def moveEvent(self, event):
        self.placeDrawer()
        super(GridGUI, self).moveEvent(event)

    def toggleLigandEdit(self, checked):
        if checked:
            #self.drawer.raise_()
            self.placeDrawer()
            self.drawer.open_()
        else:
            self.drawer.close()

    def showHideLigand(self, checked):
        if checked:
            self.pmv.showMolecules(self.ligand, show=True)
        else:
            self.pmv.showMolecules(self.ligand, show=False)

    def focusLigand(self, *args, **kw):
        self.pmv.focusScene(obj=self.ligand, padding=3)

    def setRoot(self, selections, operation):
        msg = None
        if selections.nbAtoms()==1:
            adlig = self.agfr._adLigs[0]
            atom = adlig.molPrep._ag[selections[0].getIndices()[0]]
            self.agfr._adLigs[0].TTroot = self.agfr._adLigs[0].ttbuilder.buildTree(atom)
            msg = 'you picked atom: %s:%s:%s:%s'%(
                atom.getChid(), atom.getResname(), atom.getResnum(), atom.getIcode())
            self.statusbar.showMessage(msg, level='info', timeout=5000)
            self._colorLigByRigidBody()
            self.anchorAtomGeom.Set(vertices=[atom.getCoords()])
            self._ligandRoot = atom
            self._ligandChanged = True
            
    def _colorLigByRigidBody(self):
        adlig = self.agfr._adLigs[0]
        # build colors for per-rigid-body coloring
        from PmvApp.pmvPalettes import Rainbow, RainbowSortedKey
        palette = [Rainbow[x] for x in RainbowSortedKey]
        colors = numpy.zeros( (len(adlig.molPrep._ag),3), 'f')
        for node in adlig.ttbuilder._allRigidBodies:
            colors[node.atoms] = palette[node.num%len(palette)]
        self.pmv.customColor(adlig.molPrep.select(), colors)
        
    def showHideTorsionTree(self, checked):
        self.anchorAtomGeom.Set(visible=checked)
        if checked:
            self._oldPickMode = {'vertices':'atoms', 'parts':'bonds', 'both':'both'}[self.viewer.pickLevel]
            self.viewer.setPickingMode('both')
            self.statusbar.showMessage('Left click on an atom to define a new root or a rotatable bonds to freeze/free them',
                                       level='info', timeout=5000)
            # assign left pick to pickRoot
            self._savedPickingAtomsFuncs = self.viewer.pickedAtomsFuncs[:]
            self.viewer.pickedAtomsFuncs = [self.setRoot]
            self._savedPickingBondsFuncs = self.viewer.pickedBondsFuncs[:]
            self.viewer.pickedBondsFuncs = [self.toggleRotatableBond]
            self._colorLigByRigidBody()
            self._restoreColors = self.pmv.undo.cmdStack[-1][0][0]
            self.showRotatableBonds(checked)
        else:
            self.viewer.setPickingMode(self._oldPickMode)
            # restore ligand coloring
            cmd, args, kw = self._restoreColors
            cmd(*args, **kw)
            # restore picking functions
            self.viewer.pickedAtomsFuncs = self._savedPickingAtomsFuncs
            self.viewer.pickedBondsFuncs = self._savedPickingBondsFuncs
            del self._savedPickingAtomsFuncs
            del self._savedPickingBondsFuncs
            del self._restoreColors
            self.showRotatableBonds(checked)
            
    def showRotatableBonds(self, checked):
        self.rotBondsArcs.Set(visible=checked)
        ## if checked:
        ##     self.statusbar.showMessage('Left click on bond or arrow to freeze/free rotatable bond',
        ##                                level='info', timeout=5000)
        ##     self._savedPickingFuncs = self.viewer.pickedBondsFuncs[:]
        ##     self.viewer.pickedBondsFuncs = [self.toggleRotatableBond]
        ##     self.viewer.setPickingMode('bonds')
        ## else:
        ##     self.viewer.pickedBondsFuncs = self._savedPickingFuncs
        ##     self.viewer.setPickingMode('atoms')
    
    def toggleRotatableBond(self, selections, operation):
        mol, at1, at2 = selections[0]
        key = '%d %d'%(at1.getIndex(),at2.getIndex())
        if not self._frozenBonds.has_key(key): return
        isFrozen = self._frozenBonds[key]
        if isFrozen:
            self._frozenBonds[key] = False
            self.agfr._adLigs[0].ttbuilder._rotatableBonds[key] = False
        else:
            self._frozenBonds[key] = True
            self.agfr._adLigs[0].ttbuilder._rotatableBonds[key] = True
        # recolor
        colors = []
        nbtri = self.rotBondsArcs._trianglesPerArrow
        frozenBonds = []
        for i, j in self._allRotBonds:
            if i>j:
                tmp = i
                i = j
                j = tmp
            if self._frozenBonds['%d %d'%(i,j)]:
                colors.extend( [(1.,0.,0.)]*nbtri )
                frozenBonds.append((i,j))
            else:
                colors.extend( [(0.,1.,0.)]*nbtri )
        self.rotBondsArcs.Set(materials=colors)
        self.agfr._adLigs[0].TTroot = self.agfr._adLigs[0].ttbuilder.buildTree(
            self._ligandRoot, frozenBonds=frozenBonds)
        self._colorLigByRigidBody()
        self._ligandChanged = True
        
    def createPyShellWidget(self):
        #from mglutil.util.packageFilePath import findFilePath
        #PMVICONPATH = findFilePath('Icons', 'PmvApp.GUI')
        #pyShell = self.togglePyShellAct =  QtGui.QAction(
        #QtGui.QIcon(os.path.join(PMVICONPATH, 'PyShell.png')),
        #    'Python shell', self)
        #pyShell.setStatusTip('Show/hide Python shell')
        #pyShell.connect(QtCore.SIGNAL('triggered()'), self.showPythonShell)
        #self.topBar.PyShellButton.pressed.connect(pyShell)
        #pyShell.setCheckable(True)
        #pyShell.setChecked(False)
        if not use_ipython_shell:        
            # create the PyShell widget
            from pyshell import PyShell
            self.pyShellWidget = PyShell(self)
        else :
            #create ipython embeded shell
            kernel_manager = QtInProcessKernelManager()
            kernel_manager.start_kernel()
            kernel = kernel_manager.kernel
            kernel.gui = 'qt4'
            #other variable to pass to the context ?
            kernel.shell.push({'self':self})
            #do you want pylab ?
            #kernel.shell.run_cell("import matplotlib")
            #kernel.shell.run_cell("matplotlib.rcParams['backend.qt4']='PySide'")
            #kernel.shell.run_cell("%pylab")
#            kernel.shell
            #matplotlib.rcParams['backend.qt4']='PySide'
            kernel_client = kernel_manager.client()
            kernel_client.start_channels()
            app = guisupport.get_app_qt4()
            def stop():
                kernel_client.stop_channels()
                kernel_manager.shutdown_kernel()
                app.exit()
    
            self.pyShellWidget = RichIPythonWidget()
            self.pyShellWidget.kernel_manager = kernel_manager
            self.pyShellWidget.kernel_client = kernel_client
            self.pyShellWidget.exit_requested.connect(stop)

        self.pyShellDialog = QtGui.QDialog(self.parent())
        self.pyShellDialog.setModal(False)
        # remove close button, for some reason on linux we have to set minimize
        # to get the close to disappear, try windowflags.py
        self.pyShellDialog.setWindowFlags(QtCore.Qt.Dialog |
                                          QtCore.Qt.WindowMinimizeButtonHint)
        hl = QtGui.QHBoxLayout()
        hl.addWidget(self.pyShellWidget)
        hl.setContentsMargins(0,0,0,0)
        self.pyShellDialog.setLayout(hl)
        self.pyShellDialog.resize(800, 480)
        
    def showPythonShell(self, checked):
        """toggles show/hide of the Python Shell. """
        if not checked:
            self.pyShellDialog.hide()
            if self.pyShellDialog in self._openDialogs:
                self._openDialogs.remove(self.pyShellDialog)
        else:
            self.pyShellDialog.show()
            self._openDialogs.append(self.pyShellDialog)
            

    def getReceptorMapsFilename(self):
        filename, selfilter = QtGui.QFileDialog().getOpenFileName(
            self, self.tr("Read Receptor Maps"), '',
            self.tr("trg Files (*.trg);; All files (*)"))
        if filename:
            if self.receptor:
                self.deleteReceptor()
                self.agfr.reset()
                self.reset()
            self.getReceptorMaps(filename.encode('ascii', 'replace'))

    def createGridBoxGeom(self):
        if not self.boxGeom:
            from DejaVu2.Box import NiceBox
            self.boxGeom = NiceBox('gridOutline')
            self.boxGeom.addToViewer(self.viewer)
            self.setGridVisible(True)
            self.showGridBoxAction.setEnabled(True)

    def createAnchorAtomGeom(self):
        from DejaVu2.Spheres import Spheres
        from DejaVu2.IndexedPolygons import IndexedPolygons
        from DejaVu2.Arcs3D import Arcs3D, arcVertices
        if not self.anchorAtomGeom:
            self.anchorAtomGeom = Spheres('rootAtom', visible=0, inheritMaterial=False,
                                          materials=[(1,0,1)], inheritFrontPolyMode=False,
                                          frontPolyMode='line', quality=2,
                                          inheritLineWidth=0, lineWidth=1)
        self.viewer.AddObject(self.anchorAtomGeom, parent=self.ligand.geomContainer.masterGeom)
        #if not self.rotBondsArcs:
        lig = self.ligand
        coords = lig._ag.getCoords()
        #import pdb; pdb.set_trace()
        colors = []
        verts = []
        faces = []
        off = 0
        #import pdb; pdb.set_trace()
        #for bondstr, rot in self.agfr._adLigs[0].ttbuilder._rotatableBonds.items():
        self._frozenBonds = {}
        for i, j in self.agfr._adLigs[0].ttbuilder._allRotatableBonds:
            if i>j:
                tmp = i
                i = j
                j = tmp
            if lig._ag._bondOrder['%d %d'%(i,j)] !=1: continue
            self._frozenBonds['%d %d'%(i,j)] = False
            c1 = coords[i]
            c2 = coords[j]
            center = ((c1[0]+c2[0])*.5, (c1[1]+c2[1])*.5, (c1[2]+c2[2])*.5)
            norm = 1./sqrt((c2[0]-c1[0])*(c2[0]-c1[0]) +
                           (c2[1]-c1[1])*(c2[1]-c1[1]) +
                           (c2[2]-c1[2])*(c2[2]-c1[2]))
            normal = ((c2[0]-c1[0])*norm, (c2[1]-c1[1])*norm, (c2[2]-c1[2])*norm)
            v,f = ribbonArrowGeom(arcVertices(center, normal, .5, 270, 15), normal, 0.3)
            verts.extend(v)
            faces.extend((numpy.array(f, 'i')+off).tolist())
            off += len(v)
            colors.extend( [(0.,1.,0.)]*len(f))

        self.rotBondsArcs = IndexedPolygons('rotBondsRibbonArrows', vertices=verts, faces=faces,
                                            inheritMaterial=0, 
                                            inheritCulling=0, culling='none', visible=0)
        if len(colors):
            self.rotBondsArcs.Set(materials=colors)
            self.rotBondsArcs._trianglesPerArrow = len(f)

        self._ligandRoot = self.agfr._adLigs[0]._root
        # makea copy of all rotatble bonds so we can keep drawing the arrows correctly
        # after freezing some
        self._allRotBonds = self.agfr._adLigs[0].ttbuilder._allRotatableBonds[:]
        self.viewer.AddObject(self.rotBondsArcs, parent=self.ligand.geomContainer.masterGeom)
            
    def createCovAtmsSpheresGeometries(self):
        from DejaVu2.Spheres import Spheres
        if not self.covBondAtom1Sphere:
            # Spheres to display covalent bond atoms
            self.covBondAtom1Sphere = Spheres('covBondAtom1', visible=0,
                                              inheritMaterial=False, 
                                              materials=[(1, 0.6, 0.5)], radii=0.2)
            self.viewer.AddObject(self.covBondAtom1Sphere,
                                  parent=self.receptor.geomContainer.geoms['master'])
            self.covBondAtom2Sphere = Spheres('covBondAtom2', visible=0,
                                              inheritMaterial=False, 
                                              materials=[(1, 0.6, 0.5)], radii=0.2)
            self.viewer.AddObject(self.covBondAtom2Sphere,
                                  parent=self.receptor.geomContainer.geoms['master'])

        if not self.atomsToIgnoreSpheres:
            self.atomsToIgnoreSpheres = Spheres('atomsToIgnore', visible=0,
                                                inheritMaterial=False, 
                                                materials=[(1, 0.8, 0.1)], radii=0.1)
            self.viewer.AddObject(self.atomsToIgnoreSpheres,
                                  parent=self.receptor.geomContainer.geoms['master'])
            
        
    @waiting_effects
    def getReceptorMaps(self, filename):
        from zipfile import ZipFile
        from ADFRcc.adfr import GridMap
        zf = ZipFile(filename)
        self._targetFile = filename
        if not filename in self.trgFiles:
            self.trgFiles.append(os.path.abspath(filename))
        zf.extractall(self.tmpFolder)
        folder = os.path.join(self.tmpFolder, os.path.splitext(os.path.basename(filename))[0])
        self.agfr.data = {}
        # find receptor file name i.e PDBQT file which is NOT receptor.pdbqt
        f = open(os.path.join(folder, 'data.pkl'))
        import pickle
        _data = pickle.load(f)
        f.close()
        #import pdb; pdb.set_trace()
        receptorFile = _data['inputReceptor']
        flexResStr = _data.get('flexResStr', None)
        ## for filename in zf.namelist():
        ##     if filename.endswith('.pdbqt'):
        ##         filename = os.path.split(filename)[1]
        ##         recname = os.path.splitext(filename)[0]
        ##         if recname != 'receptor':
        ##             break
        self._clusterItems = []
        #self.loadReceptor(os.path.join(folder, '%s.pdbqt'%recname))
        if self._covalentLigAtoms:
            self.hideCovLigandAtoms()
            self._covalentLigAtoms = None
            self.paramsWidget.covDockingGroupBox.setChecked(False)
            self.agfr.covalentBondToExclude = []
        self.loadReceptor(os.path.join(folder, receptorFile))
        # Find out all map types and load 'e' map for map info 
        atomTypes = ''
        nb = 0
        _mapLoaded = False
        for mapFileName in zf.namelist():
            w = mapFileName.split('.')
            if w[-1]=='map':
                if self.agfr.ADAtomTypes.has_key(w[-2]):
                    atomTypes += '%s '%w[-2]
                    nb += 1
                if not _mapLoaded and not (w[-2]=='e' or w[-2]=='d'): # read 'e' map to get center, box size, spacing, and flexres
                    _map = GridMap()
                    _f, name = os.path.split(mapFileName)
                    _map.loadFromMapFile(w[-2], folder, name)
                    center = _map.getCenterPy()
                    nx, ny, nz = _map.getNumGridPointsPy()
                    self._spacing = spacing = _map.getDistBetweenGridPoints()
                    sx = spacing*(nx-1)
                    sy = spacing*(ny-1)
                    sz = spacing*(nz-1)
                    #print 'ASDASda', mapFileName, center, sx, sy, sz
                    _mapLoaded = True

        #print 'SADASDA', nb, atomTypes, 'F', flexResStr, 'F', center, sx, sy, sz, spacing
        # set box
        self.boxGeom.setCenter( *center)
        self.boxGeom.setSides( sx, sy, sz)
        self._baseSize[:] = (sx, sy, sz)
        padding = 0.0
        self.agfr.setBox(["user", center[0], center[1], center[2], sx, sy, sz], padding, spacing)
        self.paramsWidget.gridPaddingWidget.setValue(padding)
        npts = self.agfr.boxSize
        self.paramsWidget.gradCutOffBox.setMaximum((npts[0]+1)*(npts[1]+1)*(npts[2]+1) )
        # set flexRes
        if flexResStr:
            self.paramsWidget.flexResWidget.setText(flexResStr)
        else:
            self.paramsWidget.flexResWidget.setText("")
        self.setFlexRes()
        # set atom types
        if nb==len(self.agfr.ADAtomTypes): # all types
            self.paramsWidget.atypesLabel.setText('')
            self.paramsWidget.compAllButton.setChecked(True)
        else:
            self.paramsWidget.compLigandAtypes.setChecked(True)
            self._mapTypes = atomTypes.split()
            self.setLigAtomTypesLabel(atomTypes)
        # set AutoSiteVersion
        aSiteVer = _data.get('AutoSiteVersion', "1.0")
        if aSiteVer == "1.1":
            val = _data.get('ligandSize', None)
            if val is not None:
                self.paramsWidget.autoSiteVButton.sizeSB.setValue(int(val))
                self._ligandSize = int(val)
            val = _data.get('pepScore', True)
            self._pepScore = val # True or False
            self.paramsWidget.autoSiteVButton.pepCkeckB.setChecked(val)
            self._autoSite2=True
            self.agfr.setAutoSiteVersion("1.1", ligandSize=self._ligandSize, pepScore=self._pepScore)
            self.paramsWidget.autoSiteVButton.rb2.setChecked(True)   
        else:
            self._autoSite2=False
            self.agfr.setAutoSiteVersion("1.0")
            self.paramsWidget.autoSiteVButton.rb1.setChecked(True)   
        # set water map parameters:
        val = _data.get("wMapEntropy", None)
        if val is not None:
            self.paramsWidget.wMapParamB.wEntropy.setValue(float(val))
            self.agfr.setWmapParams(ENTROPY=float(val))
        val = _data.get("wMapWeight", None)
        if val is not None:
            self.paramsWidget.wMapParamB.wWeight.setValue(float(val))
            self.agfr.setWmapParams(weight=float(val))
        # do the maps have gradients:
        val = _data.get("mapGradients", False)
        self.paramsWidget.addGradGroupBox.setChecked(val)
        # set TPoints
        tpFile = os.path.join(folder, 'translationPoints.npy')
        if os.path.exists(tpFile):
            tp = numpy.load(tpFile)
            self._clusters = [[], [[], tp, [], []]]
            self._allClusterColors = [] #[[0, 1, 0]]
            self.colorMap = [[0,1,0]]
            clustersWidget = self.paramsWidget.clustersWidget
            self.clearClustersWidget()
            clustersWidget.setHorizontalHeaderLabels(
                ['fills', 'AS Score', '#Points', 'RadGyr','Buriedness'])

            #newItem = QtGui.QListWidgetItem()
            clustersWidget.insertRow(0)
            newItem = TableWidgetItem('from file')
            newItem.setCheckState(QtCore.Qt.Checked)
            clustersWidget.setItem(0, 0, newItem)
            self._clusterItems.append(newItem)
            newItem = TableWidgetItem(str(len(tp)))
            clustersWidget.setItem(0, 2, newItem)
            self.TPoints.Set(visible=True, vertices=tp, materials=[[0, 1, 0]])
            self.agfr.setFillPoints(tp)
            #self.TPointsOK.emit(True, "")
        else:
            if _data.has_key('covalentBond'):
                inds = _data['covalentBond']
                at1Str = _data['covalentBondAtom1']
                at2Str = _data['covalentBondAtom2']
                tStr = _data['covalentBondTorsionAtom']
                tind = int(tStr.split("(")[1].split(")")[0])
                self._covalentBondAtoms =[None, None]
                w = self.paramsWidget
                w.covDockingGroupBox.setChecked(True)
                if self.covLigAtomsDialog:
                    w.covDockingGroupBox.setChecked(True)
                    QtGui.QDialog.accept(self.covLigAtomsDialog)
                w.covAt1Lab.setText(at1Str)
                w.covAt2Lab.setText(at2Str)
                id1 = self.receptor.select('serial %d'%inds[0]).getIndices()[0]
                self._covalentBondAtoms[0] = Atom(self.receptor._ag, id1, 0)
                id2 = self.receptor.select('serial %d'%inds[1]).getIndices()[0]
                self._covalentBondAtoms[1] = Atom(self.receptor._ag, id2, 0)
                id3 = self.receptor.select('serial %d'%tind).getIndices()[0]
                self._torsionAtom = Atom(self.receptor._ag, id3, 0)
                self.covBondAtom1Sphere.Set(vertices=[self._covalentBondAtoms[0].getCoords()], visible=True)
                self.covBondAtom2Sphere.Set(vertices=[self._covalentBondAtoms[1].getCoords()], visible=True)

                if _data.has_key('covalentRes'):
                    covRes = str(_data['covalentRes'])

                ignAtInds = _data.get('covalentLigandAtomIndices', [])
                toRemoveAtoms = None
                if len(ignAtInds):
                    selstr= "index " + "%d "*len(ignAtInds) % tuple(ignAtInds)
                    toRemoveAtoms = self.receptor.select(selstr)
                    self._covalentLigAtoms = toRemoveAtoms
                    self.showCovLigandAtoms()
                self.agfr.setCovalentDocking(tind, inds, toRemoveAtoms=toRemoveAtoms)
                self.TPoints.Set(visible=False)
                self.agfr.setFillPoints([])
        self.updateStatusBar()
        self.showTrgGuiAction.setDisabled(False)
        if os.path.exists(folder):
            shutil.rmtree(folder)
        self._haveTrg = True
        if self.ligand:
            self.testDockAct.setDisabled(False)
            self.viewInteractionsAct.setDisabled(False)
            
    def addGeomsForGaps(self, gapsPerCh):
        # show missing residues as dotted lines with spheres and labels for each missing residue,
        # gapsPerCh is a dictionary; {Chid: [["Nterm", [res1obj, missingResChar1,.., missingResChar#, res2obj], "CTerm"], ....], }
        if not len(gapsPerCh):
            if self.missingResLines:
                self.missingResLines.Set(visible=False)
            if self.missingResSpheres:
                self.missingResSpheres.Set(visible=False)
            if self.missingResLabels:
                self.missingResLabels.Set(visible=False)
            return
        if not self.missingResLines:
            from DejaVu2.IndexedPolylines import IndexedPolylines
            from DejaVu2.Spheres import Spheres
            self.missingResLines = IndexedPolylines('gaps', visible=0,
                                                    materials=[[1., 0.65, 0]], inheritMaterial=False,
                                                    stippleLines=True, inheritStippleLines=False)
            self.viewer.AddObject(self.missingResLines)
            
            self.missingResSpheres = Spheres('gapRes', visible=0, inheritMaterial=False,
                                             materials=[(1., 0.65, 0)], quality=2, radii=0.2)
            self.viewer.AddObject(self.missingResSpheres)
            from DejaVu2.glfLabels import GlfLabels
            self.missingResLabels = GlfLabels('gapLabels', fontStyle='solid3d',
                                     inheritMaterial=False, materials=[[1. ,0.65, 0]], 
                                     fontTranslation=(0,0,.1), fontScales=(0.3, 0.3, 0.1),
                                     pickable=0, visible = False)
            
            self.viewer.AddObject(self.missingResLabels)
        lCoords = []
        lFaces = []
        nface = 0
        spCoords = []
        labels = []
        mol = None
        #import pdb; pdb.set_trace()
        if self.receptor: # this is the case  when a pdbqt file(created by pdb analyzer) is loaded
            mol = self.receptor
        elif self.PDBprofile and self.PDBprofile._currentBM:
            mol = self.PDBprofile._currentBM
        if not mol :return
        for chid, chgaps in gapsPerCh.items():
            #print " GAPS in chain:", chid, chgaps
            for gapData in chgaps[1]:
                r1Name, r1Num, r1Icode = gapData[0]
                r2Name, r2Num, r2Icode = gapData[-1]
                # check if theese residues exists:
                if r1Icode:
                    res1 = mol._ag.select("chid %s resnum `%d` icode %s" %(chid, r1Num, r1Icode))
                else:
                    res1 = mol._ag.select("chid %s resnum `%d`" %(chid, r1Num))
                if r2Icode:
                    res2 = mol._ag.select("chid %s resnum `%d` icode %s" %(chid, r2Num, r2Icode))
                else:
                    res2 = mol._ag.select("chid %s resnum `%d`" %(chid, r2Num))
                        
                if not res1 or not res2:
                    continue
                #import pdb; pdb.set_trace()
                misRes = gapData[1:-1]
                try:
                    ccoords = res1.select('name C').getCoords()[0]
                except:
                    ccoords = res1.select('name CA').getCoords()[0]
                ncoords = res2.select('name N').getCoords()[0]
                #print "coords in _showGaps:", ccoords, ncoords
                lCoords.extend([ccoords, ncoords])
                lFaces.append([nface, nface+1])
                nface += 2
                # add spheres for missing residues
                coords = []
                x1, y1, z1 = ccoords
                x2, y2, z2 = ncoords
                nres = len(misRes)
                # this computes coords for equally spaced spheres placed over the dotted lines. One sphere per each
                # missing residue
                # for i in range(nres):
                #     t = (i+1)*1./(nres+1)
                #     d1, d2, d3 = t*(x2-x1), t*(y2-y1), t*(z2-z1)
                #     coords.append([x1+d1, y1+d2, z1+d3])
                #labels.extend(misRes)
                # place one sphere in in middle of the dotted line and a label containing letters for all missing residues
                d1, d2, d3 = (x2-x1)/2., (y2-y1)/2., (z2-z1)/2.
                coords.append([x1+d1, y1+d2, z1+d3])
                # create a label for the sphere (if it is longer than 10 characters make shorten it:ABC...EFG"
                misResStr="".join(misRes)
                if len(misResStr) > 10:
                    misResStr = misResStr[:3]+"..."+misResStr[-3:]
                spCoords.extend(coords)
                labels.append(misResStr)
            # Add speres and labels for N and C temini
            terms = []
            if len(chgaps[0]):
                terms.append([chgaps[0], "N"])# N term
            if len(chgaps[-1]): # C term
                terms.append([chgaps[-1], "C"])
            #import pdb; pdb.set_trace()
            for term, lastAtName in terms:
                trName, trNum, trIcode = term[0] # last residue
                if trIcode:
                    mtres = mol._ag.select("chid %s resnum `%d` icode %s" %(chid, trNum, trIcode))
                else:
                    mtres = mol._ag.select("chid %s resnum `%d`" %(chid, trNum))
                if not mtres: continue
                mtresStr = term[1] # missing residues string (1letter name per each missing res)
                # find N , CA or C atom in the last residue:
                atomSet = None
                for atName in [ lastAtName, "CA"]:
                    atomSet = mtres.select("name %s"% atName)
                    if atomSet: break
                if not atomSet: continue
                nbCoords = []
                for atom in atomSet.iterAtoms():
                    atCoords = atom.getCoords()
                    # find neighbors:
                    for nb in atom.iterBonded():
                        nbCoords.append(numpy.subtract(nb.getCoords(), atCoords))
                    break # should be one atom in atomSet
                vsum = numpy.add.reduce(nbCoords)
                norm = numpy.linalg.norm(vsum)
                tspCoords = atCoords-3*vsum/norm
                spCoords.append(tspCoords)
                lCoords.extend([tspCoords, atCoords])
                lFaces.append([nface, nface+1])
                nface += 2
                #vsum
                if len(mtresStr) > 10:
                    mtresStr = mtresStr[:3]+"..."+mtresStr[-3:]
                labels.append(mtresStr)
        #import pdb; pdb.set_trace()
        if len(lCoords):
            self.missingResLines.Set(vertices=lCoords, faces=lFaces, visible=True)
            self.missingResSpheres.Set(vertices=spCoords, visible=True)
            tr = numpy.ones( (len(spCoords), 3), 'f') * 0.5
            self.missingResLabels.Set(vertices=spCoords, labels=labels, labelTranslation=tr,
                                      visible=True)
        else:
            self.missingResLines.Set(vertices=[], faces=[], visible=False)
            self.missingResSpheres.Set(vertices=[], visible=False)
            self.missingResLabels.Set(vertices=[], labels=[], 
                                      visible=False)
        
    def setFolder(self):
        folder = QtGui.QFileDialog.getExistingDirectory(
            self, self.tr("select folder"), '')
        if folder:
            from glob import glob           
            filenames = glob(os.path.join(folder, "*.pdb"))
            filenames.extend(glob(os.path.join(folder, "*.ent.gz")))
            filenames.extend(glob(os.path.join(folder, "*.pdbqt")))
            if len(filenames)==0:
                msgBox = QtGui.QMessageBox(self)
                msgBox.setWindowTitle("no molecules")
                msgBox.setText("The folder %s contains no suitable molecules"%folder)
                msgBox.setStandardButtons(QtGui.QMessageBox.Ok)
                msgBox.setDefaultButton(QtGui.QMessageBox.Ok)
                ret =  msgBox.exec_()
                return
            self.setFiles(filenames)
            self._currentFileIndex = -1
            self.gotoNextFile() # go to files[0]
            
    def fetchPDB(self, name=None):
        if name is None:
            from mglutil.gui.Qt.fetchMolGUI import FetchGUI
            name, ans = FetchGUI.getName(parent=self, formats=['pdb'])
            if not ans: return
            names = [n for n in name[0].split(" ") if len(n)]
            ext = name[1]
            #print "fetchFile:", names, ext
            mol = None
            #import pdb; pdb.set_trace()
            from prody import fetchPDB
            molFiles = []
            for molName in names:
                molfile = fetchPDB(molName)
                if molfile:
                    molFiles.append(molfile)
                else:
                    msg = "Could not fetch %s molecule" % (molName +".pdb")
                    errmsg = ""
                    if len(msg):
                        errmsg = errmsg+"\n" + msg
                    if len(errmsg):        
                        # pop up a warning dialog, ask the user to retry
                        msgBox = QtGui.QMessageBox(QtGui.QMessageBox.Warning,
                                                   "Fetch Failed", errmsg,
                                                   QtGui.QMessageBox.NoButton, self)
                        msgBox.addButton("OK", QtGui.QMessageBox.AcceptRole)
                        #msgBox.addButton("Close", QtGui.QMessageBox.RejectRole)
                        #if msgBox.exec_() == QtGui.QMessageBox.AcceptRole:
                        msgBox.exec_()
            if len(molFiles)== 1:
                self.loadReceptor(molFiles[0])
            elif len(molFiles) > 1:
                self.setFiles(molFiles)
                self._currentFileIndex = -1
                self.gotoNextFile() # go to molFiles[0]
            
    def getReceptorFilename(self):
        filename, selfilter = QtGui.QFileDialog.getOpenFileName(
            self, self.tr("Read Receptor"), '',
            self.tr("receptor files (*.pdbqt *.pdb *.pdb.gz *.ent *.ent.gz);; All files (*)"))
        if filename:
            self.loadReceptor(filename.encode('ascii', 'replace'))

    def isLigand(self, filename):
        if checkLigandFile(filename):
            msgBox = QtGui.QMessageBox(self)
            msgBox.setWindowTitle("invalid receptor file")
            msgBox.setText("The file %s contains a torsion tree indicating this is a ligand."%filename)
            #msgBox.setInformativeText("No torsion tree defined in this file, still use it as a ligand?")
            msgBox.setStandardButtons(QtGui.QMessageBox.Ok)
            msgBox.setDefaultButton(QtGui.QMessageBox.Ok)
            ret =  msgBox.exec_()
            return True
        return False
    
    @waiting_effects
    def loadReceptor(self, filename):
        #logbox = LogBox('processing %s'%filename, self)
        #logbox.show()
        if self._filesToLoopOver and filename not in self._filesToLoopOver:
            self._filesToLoopOver = None
        end = self.agfr.loadReceptor(filename, processReceptor=False)
        #disable or enable save receptor button here
        if self.receptor:
            self.deleteReceptor()
        else:
            if len (self.pmv.Mols):
                [self.pmv.deleteMolecule(m) for m in self.pmv.Mols[:]]
                self.reset()
        if self.agfr._pdbInput:
            #import pdb; pdb.set_trace()
            self.setGridVisible(False)
            self.PDBprofile = self.agfr.PDBprofile
            self.setAppTitle()
            for n, name in enumerate(self.PDBprofile._assemblies.keys()):
                chset, mol = self.PDBprofile._assemblies[name]
                opac = numpy.ones( len(mol._ag) )*0.5
                segA = mol._ag.select('segment A')
                if segA:
                    opac[segA.getIndices()] = 1.
                mol._ag.setData('opacity_lines', opac)
                self.pmv.addMolecule(mol)

            # create analyser widget 
            self._pdbMol = self.agfr.receptor
            from pdbAnalyser import PDBAnalyserWidget
            self.analyserWidget = PDBAnalyserWidget(self, self.PDBprofile, None)

            # show Biomolecule 1 by default
            tw = self.analyserWidget.treeWidgets
            if tw.has_key('Biomol 1'):
                #tw['Biomol 1'].buildUI()
                self.analyserWidget.treesNoteBook.setCurrentWidget(tw['Biomol 1'])
            else:
                #tw.values()[0].buildUI()
                self.analyserWidget.treesNoteBook.setCurrentWidget(tw.values()[0])
            self.analyserWidget.switchBM()
            self.addGeomsForGaps(self.PDBprofile._gapsPerChain)
            self.stack.addWidget(self.analyserWidget)
            self.stack.setCurrentWidget(self.analyserWidget)
            # disable widgets in the vertical toolbar :
            self.showReceptorAction.setDisabled(True)
            self.showReceptorSurfaceAction.setDisabled(True)
            self.showResLabAction.setDisabled(True)
            self.editLigandAct.setChecked(False)
            self.editLigandAct.setDisabled(True)
            self.showCovBondAtoms.setDisabled(True)
            self.showGridBoxAction.setDisabled(True)
            self.showTrgGuiAction.setDisabled(True)
            self.showXtalMates.setDisabled(True)
            self.testDockAct.setDisabled(True)
            self.viewInteractionsAct.setDisabled(True)
        else:
            mols = self.pmv.readMolecules(filename)
            self.setReceptor(mols[0])
            self.PDBprofile = None
            self.setAppTitle()
            # hide missing residues geometry
            if self.missingResLines:
                self.missingResLines.Set(visible=False)
            if self.missingResSpheres:
                self.missingResSpheres.Set(visible=False)
            if self.missingResLabels:
                self.missingResLabels.Set(visible=False)
            self.stack.setCurrentWidget(self.paramsWidget)
            
        #logbox.close()

    @waiting_effects
    def finalizeMolecules(self):
        # called when the PDB analyzer is done and we want to proceed AFGRgui
        mol = self.agfr._BMProcessed
        if mol is None:
            mol = self.agfr.processBM()

        if self.agfr._processedReceptor is None:
            atoms = mol.select('receptor')
            if atoms is None:
                print 'WARNING: no atoms selected for the receptor'
                return
            receptor = self.agfr.processReceptor(atoms, charges=self.analyserWidget.getCharges())
            self.agfr._processedReceptor = receptor
            if self.agfr._saveFilename is None:
                self.agfr._saveFilename = self.PDBprofile._pdbid
            filename = self.agfr.saveReceptorAsPDBQT(receptor, self.agfr._saveFilename+'_rec.pdbqt')
        else:
            filename = self.agfr._processedReceptor.filename

        if len(self.agfr._adLigs) == 0: # we do not have processed ligands
            ligs = self.agfr.saveLigands()
            if len(ligs) > 1:
                msgBox = QtGui.QMessageBox(self)
                msgBox.setText("INFO: ligand atoms form %d connected fragments.\nA separate ligand file was created for each fragment. The first ligand will be loaded in th maps GUI"%len(ligs))
                msgBox.exec_()
                
        else:
            ligs = self.agfr._processedLigands
            
        # this works but we do not have to the PDBQT source needed to run AutoGrid
        # so we write the PDBQT file and read it
        # self.pmv.addMolecule(receptor)
        # self.setReceptor(rreceptor)
        recmol = self.pmv.readMolecules(filename)[0]
        recmol.pdbHeader = self.agfr.PDBprofile._receptor.pdbHeader.copy()
        self.setReceptor(recmol)
        
        # cleanup the Viewer
        for n, name in enumerate(self.PDBprofile._assemblies.keys()):
            chset, mol = self.PDBprofile._assemblies[name]
            self.pmv.deleteMolecule(mol)
            # delete rotamol molecules
            for rotamol in self.analyserWidget.treeWidgets[name]._rotamols:
                self.pmv.deleteMolecule(rotamol._mol)
            del self.analyserWidget.treeWidgets[name]._rotamols
            del self.analyserWidget.treeWidgets[name]._mol

        # switch to AGFRgui
        self.stack.setCurrentWidget(self.paramsWidget)

        # load a ligand if none was specified on cmdline and we have some from pdbanalyzer
        if self.agfr._ligFilename is None: # no ligand specified on cmdline
            if ligs: # we have processed ligands we use the first one
                self.pmv.addMolecule(ligs[0])
                self.setLigand(ligs[0])

        # add gap geometries for context
        self.addGeomsForGaps(self.PDBprofile._gapsPerChain)

        chids = set(recmol._ag.getChids())
        if len(self.PDBprofile._sites):
            for k,v in self.PDBprofile._sites.items():
                #self._pdbSites[name] = [(res.getChid(), res.getResnum(), res.getIcode()) for res in v[0]]
                chid, resname, resnum = k.split(':')
                if chid in chids:
                    self._pdbSites[k] = v[0]
        
        # create structures for residues which sidechains that should
        # potentially made flexible
        chids = set(self.agfr.receptor._ag.getChids())
        self._altlocResidues = self.analyserWidget.getAltlocResidues()
        #print "ALTLOCK RESIDUES:", self._altlocResidues
        self._mutatedResidues = {}
        for key, resnames in self.PDBprofile._resToMutate.items():
            chid, resname, resnum, icode,  = key.split(':')
            if chid not in chids: continue
            self._mutatedResidues[key] = resnames
        del self.analyserWidget

        # enable the button to show crystal mates if info avaiable
        if len(recmol.pdbHeader['SCALE'])>0:
            self.showXtalMates.setDisabled(False)
        
    def getLigandFilename(self):
        """Callback of the 'load ligand' pushbutton"""
        filename, selfilter = QtGui.QFileDialog.getOpenFileName(
            self, self.tr("Read Ligand"), '',
            self.tr("ligand files (*.pdbqt *.pdb *.ent *.ent.gz *.mol2 *.sdf *.cif);; All files (*)"))
        if filename:
            self.loadLigand(filename.encode('ascii', 'replace'))

    def loadLigand(self, filename):
        #if not checkLigandFile(filename):
        #    msgBox = QtGui.QMessageBox(self)
        #    msgBox.setWindowTitle("invalid ligand file")
        #    msgBox.setText("The file %s does not contain a torsion tree indicating this is a receptor."%filename)
            #msgBox.setInformativeText("No torsion tree defined in this file, still use it as a ligand?")
        #    msgBox.setStandardButtons(QtGui.QMessageBox.Ok)
        #    msgBox.setDefaultButton(QtGui.QMessageBox.Ok)
        #    ret =  msgBox.exec_()
        #    return
        if self.ligand:
            self.pmv.deleteMolecule(self.ligand)
        _mols = self.pmv.readMolecules(filename)
        #self.setLigand(mols[0])
        #_mol = Read(filename)

        self.agfr._processedLigands = self.agfr.processLigands(_mols[0]._ag)
        mol = self.agfr._processedLigands[0]
        self.pmv.addMolecule(mol)
        self.pmv.undisplayLines(_mols[0])
        self.setLigand(mol)
        #self.viewer.Reset_cb()
        #self.viewer.Normalize_cb()
        #self.viewer.Center_cb()
        if isinstance(_mols[0], MultiMolecule) or _mols[0]._multi=='conformations':
            self._ligandMultiMol = _mols[0]
            #self.nextFileAct.setDisabled(False)
            self.topBar.nextWidget.setDisabled(False)
            
    def setReceptor(self, mol):
        self.receptor = mol
        self.agfr.setReceptor(mol)
        self.createGridBoxGeom()
        self.createCovAtmsSpheresGeometries()
        self._hasReceptorSurface = False
        self.setAppTitle()
        self.showReceptorAction.setDisabled(False)
        self.showReceptorSurfaceAction.setDisabled(False)
        self.showReceptorAction.setChecked(True)
        self.showResLabAction.setDisabled(False)
        #self.defaultReceptorBox()

        self.setGridVisible(True)
        self.setGridFullReceptor()
        self.showGridBoxAction.setEnabled(True)
        self.showGridBoxAction.setChecked(True)

        self.paramsWidget.groupBox3.setDisabled(False)
        #self.paramsWidget.ASGridAct.setDisabled(False)
        self.paramsWidget.recGridAct.setDisabled(False)
        self.paramsWidget.manualGrdAct.setDisabled(False)
        self.paramsWidget.residuesGrdAct.setDisabled(False)
        self.paramsWidget.groupBox4.setDisabled(False)
        self.paramsWidget.groupBox5.setDisabled(False)
        self.paramsWidget.groupBox6.setDisabled(False)
        self.paramsWidget.covDockingGroupBox.setDisabled(False)
        self.viewer.Reset_cb()
        self.viewer.Normalize_cb()
        self.viewer.Center_cb()
        
    def setLigand(self, mol):
        self._ligandRoot = None
        self._frozenBonds = {}
        self._ligandChanged = False
        self.ligand = mol
        self.agfr.setLigand(mol)
        self.setAppTitle()
        atypes = list(set(mol._ag.getData('AD_element')))
        if 'C' not in atypes:
            atypes.insert(0, 'C')
        self.pmv.displaySB(mol)
        self.pmv.undisplayLines(mol)
        self.pmv.customColor(mol.select('element C'), [(1.,1.,0.)])
        self.paramsWidget.groupBox3.setDisabled(False)
        self.createAnchorAtomGeom()
        self.anchorAtomGeom.Set(vertices=[self.agfr._adLigs[0]._root.getCoords()], radii=0.5)
        self.paramsWidget.ligGridAct.setDisabled(False)
        self.paramsWidget.editAtypesButton.setDisabled(False)
        self._mapTypes = atypes
        self.setLigAtomTypesLabel(' '.join(atypes))
        self.paramsWidget.compLigandAtypes.setChecked(True)
        self.editLigandAct.setDisabled(False)
        self.ligand_tb._widgets[0].setChecked(True)
        if self._haveTrg:
            self.testDockAct.setDisabled(False)
            self.viewInteractionsAct.setDisabled(False)
        self.viewInteractionsAct.setChecked(False)
        self.viewInteractions()

    def deleteReceptor(self):
        if self.receptor:
            self.pmv.deleteMolecule(self.receptor)
            self.receptor = None
            #remove/hide geoms associated with this receptor:
            # hide missing residues geometry
            if self.missingResLines:
                self.missingResLines.Set(visible=False)
            if self.missingResSpheres:
                self.missingResSpheres.Set(visible=False)
            if self.missingResLabels:
                self.missingResLabels.Set(visible=False)
            if self.TPoints:
                self.TPoints.Set(visible=False, vertices=[])
            if self.agfr.data.has_key('flexResStr') and self.agfr.data['flexResStr']:
                self.agfr.data.pop('flexResStr')
                self._flexRes = [] 
                self._flexResSCAtoms = None
                self._flexResAtoms = []
                self.paramsWidget.flexResWidget.setText("")
                self.paramsWidget.setFlexResGrdButton.setDisabled(True)
            self._hasReceptorSurface = False
            self._MSMSDisplayMode = 0

        
    def setLigAtomTypesLabel(self, atomTypesStr):
        #atomTypesStr = atomTypesStr.encode('ascii', 'replace')
        if len(atomTypesStr) > 30:
            lab = atomTypesStr[:10] + ' ... ' + atomTypesStr[-11:-1]
            self.paramsWidget.atypesLabel.setText(lab)
        else:
            self.paramsWidget.atypesLabel.setText(atomTypesStr)
        
    def updateDisplay(self):
        # handle lines of receptor
        if self.showReceptorAction.isChecked():
            self.receptor.geomContainer.geoms['master'].Set(visible=True)
            resWithAtomsInBox = self.getResWithAtomsInBox()
            gc = self.receptor.geomContainer
            opac = numpy.ones( len(gc.allCoords) )*0.3
            if resWithAtomsInBox:
                opac[resWithAtomsInBox.getIndices()] = 1.
            #gc.geoms['singleBonds'].Set(
            #    opacity=opac, transparent=True, inheritMaterial=False, polyFace='front')
            self.receptor._ag.setData('opacity_lines', opac)
            self.pmv.displayLines(self.receptor)
        else:
            self.receptor.geomContainer.geoms['master'].Set(visible=False)
            ## if self.receptor:
            ##     self.pmv.undisplayLines(self.receptor)
            ## if self.covBondAtom1Sphere:
            ##     self.covBondAtom1Sphere.Set(visible=False)
            ##     self.covBondAtom2Sphere.Set(visible=False)
            ##     self.atomsToIgnoreSpheres.Set(visible=False)

        # handle residue labels
        if self._labDisplayMode==0 and self.receptor: # hide all residue labels
            self.pmv.unlabelResidues(self.receptor)
        elif self._labDisplayMode==1: # label residues in box
            # hide existing labels
            self.pmv.unlabelResidues(self.receptor)
            # put labels that are now inside the box
            atoms =  self.selectCAInBox()
            if atoms:
                self.pmv.labelResidues(SelectionSet([atoms]))
        elif self._labDisplayMode==2: # label flexres
            #import pdb; pdb.set_trace()
            if self.receptor:
                self.pmv.unlabelResidues(self.receptor)
            if self._flexResAtoms:
                self.pmv.labelResidues(SelectionSet([self._flexResAtoms]))

        # handle receptor surface
        #print self._MSMSDisplayMode, self._hasReceptorSurface
        if self._hasMSMS:
            if self._MSMSDisplayMode==0 and self.receptor: # hide all residue labels
                if self._hasReceptorSurface:
                    self.pmv.undisplayMSMS(self.receptor)
            else:
                if not self._hasReceptorSurface:
                    self.pmv.computeMSMS(self.receptor, surfName="%s_surface"%self.receptor.name, density=6.0)
                    self._hasReceptorSurface = True
                if self._MSMSDisplayMode==1: # surface for atoms in the box
                    # hide existing surface
                    self.pmv.undisplayMSMS(self.receptor)
                    self.pmv.displayMSMS(self.selectAtomsInBox())
                elif self._MSMSDisplayMode==2: # entire surface
                    self.pmv.displayMSMS(self.receptor)
            
    def updateStatusBar(self):
        from time import time
        t0 = time()
        errCodes = self.agfr.checkComputeGrids()
        statusMsg = ""
        #print "UPDATE STATUS BAR:", errCodes
        self.paramsWidget._tpok = True
        self.paramsWidget._gridok = True
        self.paramsWidget._boxok = True
        self.paramsWidget._covbondok = True
        boxOK =  QtCore.Qt.green
        boxOkTip = "docking box overlaps receptor"
        tpOK = QtCore.Qt.green
        tpOkTip = "ligand binding pocket points overlap with box"
        frOK = None # QtCore.Qt.green
        frOkTip = "Number of flexible residue(s)"
        ready = False
        numFR = 0
        if not self._flexResSCAtoms:
            frOK = QtCore.Qt.gray
            frOkTip = "no flexible residues selected"
        if self._covalentDocking:
            tpOK = QtCore.Qt.gray
            tpOkTip = "covalent docking mode. No binding pocket points needed"
            covDock = QtCore.Qt.green
            covDockTip = "covalent docking is enabled"
        else:
            covDock = QtCore.Qt.gray
            covDockTip = "covalent docking is disabled"
        for err, msg in errCodes:
            if err in [104, 109]: # covalent bond atoms are missing or not in the box
                covDock = QtCore.Qt.red
                covDockTip = msg
                self.paramsWidget._covbondok = False
            elif err in [105, 106]: # no translational points or the points are not in the box
                tpOK = QtCore.Qt.red
                self.paramsWidget._tpok = False
                tpOkTip = msg
            elif err == 101: #no box , TP are not OK either
                boxOK =  QtCore.Qt.red
                self.paramsWidget._boxok = False
                boxOkTip = msg
                tpOK = QtCore.Qt.red
                self.paramsWidget._tpok = False
                tpOkTip = msg
            elif err == 102: #receptor is not in the box
                boxOK =  QtCore.Qt.red
                self.paramsWidget._boxok = False
                boxOkTip = msg
            elif err in [103, 108]:
                if err == 103: # flex residue(s) are outside the box.
                    self.paramsWidget.setFlexResGrdButton.setDisabled(False)
                    prefVal = self.preferences['Flexible Residues Out Of Box']['value']
                    if prefVal == "Always ask":
                        # open a dialog asking the user whether to adjust the box or make residues that are outside the box rigid.
                        self._resOutsideBoxDialog = dial = CheckParamsInBoxDialog(msg, title="Selected residues outside the box",
                                    buttons=["Adjust the box to include\n all flexible sidechains",
                                             "Make sidechains outside\n the box rigid",
                                             "I'll fix it later"])
                        dial.finished.connect(self.flexResOutOfBox_cb)
                        self._openDialogs.append(dial)
                        # if flexible residues chooser dialog is open, then disable it so that it does not allow to check any
                        # more residues and does not pop up this dialog.
                        if self._flexResDialog and self._flexResDialog.isVisible():
                            self._flexResDialog.setEnabled(False)
                        res = dial.show()
                        msg = ""
                    elif prefVal == "Adjust box":
                        self.adjustGridBoxFlexRes()
                        frOK = QtCore.Qt.green
                        msg = ""
                    elif prefVal == "Make sidechains rigid":
                        self.makeResOutsideBoxRigid()
                        frOK = QtCore.Qt.green
                        msg = ""
                    else: #"None"
                        frOK = QtCore.Qt.red
                        self.paramsWidget._gridok = False
                        frOkTip = msg

            elif err == 110: # specified residue cannot be made flexible
                res = msg[1] # this should be a tuple (chainID, resName, resNum)
                dg = QtGui.QMessageBox(QtGui.QMessageBox.Warning, "Flexible residues error", msg[0])
                dg.exec_()
                chid, rs, rsnum = res
                #Uncheck the incorrect residue in the Flex Res dialog (if it is displayed) 
                if hasattr(self, "_flexResDialog") and self._flexResDialog.isVisible():
                    resstr = "%s%d"%(rs, rsnum)
                    resitems = self._flexResDialog.resTreeWidget.treeWidget.findItems(resstr, QtCore.Qt.MatchContains | QtCore.Qt.MatchRecursive)
                    for resitem in resitems:
                        if resitem.text(0) == resstr:
                            resitem.setCheckState(0, QtCore.Qt.Unchecked)
                else:
                    # remove the residue from the list of flex resudues in the text entry of the "flexible residue(s)" group
                    flexResStr= self.agfr.data['flexResStr'] # this string contains the residue that needs to be removed
                    newFlexResStr = ""
                    for _chid, residues in flexResStr2flexRes(flexResStr):
                        resstr = ""
                        if _chid == chid:
                            for _res in residues:
                                if _res[0] == rs and _res[1] == rsnum: 
                                    continue
                                resstr +="%s%d,"%(_res[0], _res[1])
                            if len(resstr):
                               newFlexResStr += _chid + ":" + resstr[:-1]+ ";" 
                        else:
                            newFlexResStr += _chid + ":"
                            for _res in residues:
                                resstr +="%s%d,"%(_res[0], _res[1])
                            newFlexResStr += resstr[:-1]+ ";"
                    if len(newFlexResStr):
                        newFlexResStr = newFlexResStr[:-1] #remove trailing ";"
                    # set the "flexible residue(s)" text line widget with the new string and update.
                    self.paramsWidget.flexResWidget.setText(newFlexResStr)
                    self.setFlexRes()
                #this should fix the problem, so do not change the FlexRes status button to red
                msg=""
            elif err == 0:
                ready = True
            if len(msg):
                statusMsg += msg +"; "
        if self._flexResAtoms:
            numFR = self._flexResAtoms.getHierView().numResidues()
            if frOK != QtCore.Qt.red: # flex residues are inside the box - test for the number of the FR
                if numFR > 0 and numFR<9:
                    frOK = QtCore.Qt.green
                elif numFR >8 and numFR < 16:
                    #frOK = QtCore.Qt.yellow
                    frOK = QtGui.QColor("#FFA500").lighter() # light orange
                elif numFR > 15:
                    #frOK = QtCore.Qt.red
                    frOK = QtGui.QColor("#FF8C00") # darker orange
            frOkTip = "%s: %d" %("Number of flexible side chains:" , numFR) 
        
        statusMsg = statusMsg[:-2]# remove last "; "
        #print "status Message:", statusMsg
        palette = QtGui.QPalette() 
        palette.setColor(QtGui.QPalette.Window, boxOK)
        self.boxOkLabel.setPalette(palette)
        self.boxOkLabel.setToolTip(boxOkTip)
        
        palette = QtGui.QPalette() 
        palette.setColor(QtGui.QPalette.Window, tpOK)
        self.tpOkLabel.setPalette(palette)
        self.tpOkLabel.setToolTip(tpOkTip)
        
        palette = QtGui.QPalette() 
        palette.setColor(QtGui.QPalette.Window, frOK)
        self.frOkLabel.setPalette(palette)
        self.frOkLabel.setToolTip(frOkTip)
        if numFR > 0:
            self.frOkLabel.setText("FR-%d"%numFR)
        else:
            self.frOkLabel.setText("FR")
        self.statusbar.showMessage(statusMsg)
        
        palette = QtGui.QPalette() 
        palette.setColor(QtGui.QPalette.Window, covDock)
        self.covLabel.setPalette(palette)
        self.covLabel.setToolTip(covDockTip)
        if ready:
            self.paramsWidget.computeGridsButton.setDisabled(False)
        else:
            self.paramsWidget.computeGridsButton.setDisabled(True)
        #print "UPDATE STATUS BAR, time" , time()-t0
            
    def onBoxChange(self, invalidateClusterPoints=True):
        cx, cy, cz = self.boxGeom.center
        sx, sy, sz = self.boxGeom.sides
        self._LL = [cx - (sx*.5), cy - (sy*.5), cz - (sz*.5)]
        self._UR = [cx + (sx*.5), cy + (sy*.5), cz + (sz*.5)]
        coords = []
        if len(self._clusterItems):
            coords = self.getClusterPoints()
        self.agfr.setFillPoints(coords)
        npts = self.agfr.boxSize
        self.paramsWidget.gradCutOffBox.setMaximum((npts[0]+1)*(npts[1]+1)*(npts[2]+1) )
        self.updateStatusBar()
        #if invalidateClusterPoints:
        #    self.clearClustersWidget()
        #    self.TPoints.Set(visible=0)
        #    self.paramsWidget.computeGridsButton.setDisabled(True)
        #    self.paramsWidget.TPGridAct.setDisabled(True)
        #    self._clusterItems = []
        self.updateDisplay()

    def showHideCovBond(self):
        if self._covalentBondAtoms[0] is not None and self._covalentBondAtoms[1] is not None:
            #val = self.showCovBondAtoms.isChecked()
            val = False
            if self.covBondAtom1Sphere.visible and self.covBondAtom2Sphere.visible:
                val = True
            self.covBondAtom1Sphere.Set(visible=not val)
            self.covBondAtom2Sphere.Set(visible=not val)
        
    def showHideReceptor(self):
        self.updateDisplay()

    def showHideReceptorSurface(self):
        if self._hasMSMS:
            self._MSMSDisplayMode = (self._MSMSDisplayMode+1)%3
            self.updateDisplay()
        else:
            msg = """MSMS/MSLIB is not enabled in your installation of MGLTools2.\nHence MSMS-based molecular surfaces cannot be calculated.\nMSMS/MSLIB is free for academic but requires a license for commercial use.\nPlease contact Dr. Michel F. Sanner (sanner@scripps.edu) for a license dor commercial use of MSMS/MSLIB.\nOnce you obtain a license key, you can enable the molecular surface features through the help menu in the Pmv application"""

            msgBox = QtGui.QMessageBox()
            msgBox.setText(msg)
            msgBox.exec_()
            
    def showHideGridBox(self):
        self.setGridVisible(self.showGridBoxAction.isChecked())
        self.TPoints.Set(visible=self.showGridBoxAction.isChecked())
        
    def showHideResidueLabels(self):
        self._labDisplayMode = (self._labDisplayMode+1)%3
        if self._labDisplayMode==2 and self._flexResSCAtoms is None:
            self._labDisplayMode=0
        self.updateDisplay()
        
    def showHideFP(self):
        if self.showFPAction.isChecked():
            self.TSpheresC.Set(visible=True)
            self.TSpheresO.Set(visible=True)
            self.TSpheresH.Set(visible=True)
        else:
            self.TSpheresC.Set(visible=False)
            self.TSpheresO.Set(visible=False)
            self.TSpheresH.Set(visible=False)
        
    def focusOnGrid(self):
        if self.boxGeom:
            self.pmv.focusScene(obj=self.boxGeom.corners)
        else:
            self.pmv.focusScene(obj=self.PDBprofile._currentBM)
        ## oldRecVis = oldLigVis = None
        ## if self.receptor:
        ##     oldRecVis = self.receptor.geomContainer.geoms['master'].visible
        ##     self.receptor.geomContainer.geoms['master'].Set(visible=False)
        ## if self.ligand:
        ##     oldLigVis = self.ligand.geomContainer.geoms['master'].visible
        ##     self.ligand.geomContainer.geoms['master'].Set(visible=False)
        ## oldGridVis = self.showGridBoxAction.isChecked()
        ## if not oldGridVis:
        ##     self.setGridVisible(True)
        ## rot = self.viewer.rootObject.rotation
        ## self.viewer.Reset_cb()
        ## self.viewer.Normalize_cb()
        ## self.viewer.Center_cb()
        ## self.viewer.rootObject.Set(rotation=rot)
        ## self.setGridVisible(oldGridVis)
        ## if oldRecVis is not None:
        ##     self.receptor.geomContainer.geoms['master'].Set(visible=oldRecVis)
        ## if oldLigVis is not None:
        ##     self.ligand.geomContainer.geoms['master'].Set(visible=oldLigVis)

    def openTrgMapGui(self):
        #if not len(self.trgFiles): return
        if self.showTrgGuiAction.isChecked():
            #show trg gui 
            if not self.trgMapGui:
                # create popup GUI for displaying maps in a tree widget and
                # showing isocontours
                from trgMapsGUI import TrgMapsGui
                self.trgMapGui = TrgMapsGui(parent=self, app=self)
                self.trgMapGui.viewer = self.viewer
                for filename in self.trgFiles:
                    self.trgMapGui.addTrg(filename)#, boxGeom=self.boxGeom)
                self.trgMapGui.addButton.setDisabled(True)
                self.trgMapGui.removeButton.setDisabled(True)
                def onClose():
                    self.showTrgGuiAction.setChecked(False)
                self.trgMapGui.closedSignal.connect(onClose)
            else:
                # the GUI has been created. Update its tree widget with the
                # trg names from self.trgFiles
                trgInGui = [val['filename'] for val in self.trgMapGui.trg_data.values()]
                #print "TRG in mapGUI:", trgInGui
                for filename in self.trgFiles:
                    if not filename in trgInGui:
                        self.trgMapGui.addTrg(filename)
            self.trgMapGui.show()
            self._openDialogs.append(self.trgMapGui)
        else:
            # hide trg gui
            if self.trgMapGui and self.trgMapGui.isVisible():
                self.trgMapGui.hide()
                if self.trgMapGui in self._openDialogs:
                    self._openDialogs.remove(self.trgMapGui)
    def showHideReport(self):
        print 'Not Yet'

    def reset(self):
        self.receptor = None
        self._haveTrg = False # set to True if target file is loaded
        #self.setGridVisible(False)

    def manageButtons(self):
        if self._currentFileIndex > 0:
            self.topBar.previousWidget.setDisabled(False)
        if self._currentFileIndex < len(self._filesToLoopOver)-1:
            self.topBar.nextWidget.setDisabled(False)

    @waiting_effects
    def gotoNextFile(self):
        self.topBar.previousWidget.setDisabled(True)
        self.topBar.nextWidget.setDisabled(True)
        if self._filesToLoopOver:
            self._currentFileIndex += 1
            self.manageButtons() 
            if self._currentFileIndex<len(self._filesToLoopOver):
                if self.receptor:
                        self.deleteReceptor()
                currentMolecules = self.pmv.Mols[:]
                try:
                    #print 'loading', self._filesToLoopOver[self._currentFileIndex]
                    [self.pmv.deleteMolecule(m) for m in currentMolecules]
                    self.agfr.reset()
                    self.reset()
                    #import pdb; pdb.set_trace()
                    self.loadReceptor(self._filesToLoopOver[self._currentFileIndex])
                except Exception, e:
                    print e

        if self._ligandMultiMol:
            if isinstance(self._ligandMultiMol, MultiMolecule):
                self.ligand_tb.reset()
                if self.ligand:
                    self.pmv.deleteMolecule(self.ligand)
                _mol = self._ligandMultiMol
                while True:
                    ok = _mol.gotoNext()
                    try:
                        self.agfr._processedLigands = self.agfr.processLigands(_mol._ag)
                        break
                    except:
                        #FIXME .. usually happend with unsupported atom types.
                        # maybe offer to fix bad atoms
                        pass
                if _mol.curMolIndex() > 0:
                    self.topBar.previousWidget.setDisabled(False)
                if _mol.curMolIndex() < _mol.numMols():
                    self.topBar.nextWidget.setDisabled(False)
                self.pmv.addMolecule(self.agfr._processedLigands[0])
                self.setLigand(self.agfr._processedLigands[0])
            else: # multi conf
                ind = self._ligandMultiMol.curMolIndex()
                ok = self.pmv.multiGoto(self._ligandMultiMol, ind+1)
                if ok:
                    # self._ligandMultiMol is not display. It holds the
                    # multimol. We have to update self.ligand which is the
                    # visible ligand molecule
                    c = self._ligandMultiMol._ag.getCoords()
                    self.ligand._ag.setCoords(c)
                    self.ligand.geomContainer.allCoords[:] = c
                    self.ligand.setBondorder(self.ligand._ag._bondOrder)
                    if self.ligand.geomContainer.geoms.has_key('aromaticSpheres'):
                        bonds = self.ligand.select().getBonds()
                        if len(bonds[4]):
                            self.pmv.displaySB.updateAromaticArcs(self.ligand, bonds[4])
                    event = RefreshDisplayEvent(molecule=self.ligand)
                    self.pmv.eventHandler.dispatchEvent(event)
                    self.viewInteractions()

                if self._ligandMultiMol.curMolIndex() > 0:
                    self.topBar.previousWidget.setDisabled(False)
                if self._ligandMultiMol.curMolIndex() < self._ligandMultiMol.numMols()-1:
                    self.topBar.nextWidget.setDisabled(False)
        else:
            pass

    @waiting_effects
    def gotoPreviousFile(self):
        self.topBar.previousWidget.setDisabled(True)
        self.topBar.nextWidget.setDisabled(True)
        if self._filesToLoopOver:
            self._currentFileIndex -= 1
            self.manageButtons() 
            if self._currentFileIndex>=0:
                self.topBar.nextWidget.setDisabled(False)
                currentMolecules = self.pmv.Mols[:]
                try:
                    #print 'loading', self._filesToLoopOver[self._currentFileIndex]
                    [self.pmv.deleteMolecule(m) for m in currentMolecules]
                    self.agfr.reset()
                    self.reset()
                    #import pdb; pdb.set_trace()
                    self.loadReceptor(self._filesToLoopOver[self._currentFileIndex])
                except Exception, e:
                    print e

        if self._ligandMultiMol:
            if isinstance(self._ligandMultiMol, MultiMolecule):
                self.ligand_tb.reset()
                if self.ligand:
                    self.pmv.deleteMolecule(self.ligand)
                _mol = self._ligandMultiMol
                while True:
                    ok = _mol.gotoPrevious()
                    try:
                        self.agfr._processedLigands = self.agfr.processLigands(_mol._ag)
                        break
                    except:
                        #FIXME .. usually happend with unsupported atom types.
                        # maybe offer to fix bad atoms
                        pass
                if _mol.curMolIndex() > 0:
                    self.topBar.previousWidget.setDisabled(False)
                if _mol.curMolIndex() < _mol.numMols():
                    self.topBar.nextWidget.setDisabled(False)
                self.pmv.addMolecule(self.agfr._processedLigands[0])
                self.setLigand(self.agfr._processedLigands[0])
            else: # multi conf
                ind = self._ligandMultiMol.curMolIndex()
                ok = self.pmv.multiGoto(self._ligandMultiMol, ind-1)
                if ok:
                    c = self._ligandMultiMol._ag.getCoords()
                    self.ligand._ag.setCoords(c)
                    self.ligand.geomContainer.allCoords[:] = c
                    self.ligand.setBondorder(self.ligand._ag._bondOrder)
                    if self.ligand.geomContainer.geoms.has_key('aromaticSpheres'):
                        bonds = self.ligand.select().getBonds()
                        if len(bonds[4]):
                            self.pmv.displaySB.updateAromaticArcs(self.ligand, bonds[4])
                    event = RefreshDisplayEvent(molecule=self.ligand)
                    self.pmv.eventHandler.dispatchEvent(event)
                    self.viewInteractions()

                if self._ligandMultiMol.curMolIndex() > 0:
                    self.topBar.previousWidget.setDisabled(False)
                if self._ligandMultiMol.curMolIndex() < self._ligandMultiMol.numMols()-1:
                    self.topBar.nextWidget.setDisabled(False)
                    
        else:
            pass
        
    def displaySettings(self):
        print 'NOT YET'
        
    def _mkADFR(self, ligandFile, seed):
        from time import time
        t0 = time()
        if ligandFile is None:
            ligandFile = '%s_lig.pdbqt'%self.ligand.name
            self.agfr.saveLigandAsPDBQT(self.ligand.name)
            self.ligand.filename = ligandFile

        from ADFR import getLigandFromFile
        ligand, error = getLigandFromFile(ligandFile)

        from ADFR.utils.maps import MapsFile, flexResStr2flexRes
        mf = MapsFile(self._targetFile)
        mf.unzipMaps()
        mapsFolder = mf.getMapsFolder()
        receptorFilename = os.path.join(mf.getMapsFolder(),
                                        mf.getReceptorFilename())
        receptor = Read(receptorFilename)
        flexResStr = mf.getFlexResStr()
        covalentRec = mf.getCovalentBond()

        if covalentRec is not None:
            covalentRec += [d.get('covalentBondTorsionAtoms')[-1]]
        else:
            covalentRec = None    
        tPtsFile = os.path.join(mapsFolder, "translationPoints.npy")
        if os.path.exists(tPtsFile):
            tpointsFilename = os.path.join(mapsFolder, "translationPoints.npy")
        else:
            tpointsFilename = None
        mapFilesRoot = 'rigidReceptor'

        from ADFR import ADFR
        adfr = ADFR(ligand, mapsFolder, logFile=None,
                    receptor=receptor,
                    flexibleResidues=flexResStr,
                    mapFilesRoot=mapFilesRoot,
                    tpointsFilename=tpointsFilename,
                    seedValue=seed,
                    ## covalentIndices=covalentIndices,
                    ## FTRecsrc=args.get('FTRecsrc', None),
                    ## fixedRoot=fixedRoot,
                    ## neighborSearchCutoff=neighborSearchCutoff
                    )
        print time()-t0
        return adfr

    def viewInteractions(self):
        if not self.viewInteractionsAct.isChecked():
            self.pairsTable.setMaximumHeight(0)
            self._interCyl.Set(visible=0)
            self._interSph.Set(visible=0)
            return

        rec = self.receptor
        lig = self.ligand
        recAS = AtomSet(molToAtomSetStatic(rec))
        ligAS = AtomSet(molToAtomSetStatic(lig))

        scorer = PairwiseScorer()
        scorer.initialize(recAS, ligAS)
        score = scorer.calculateScores()

        self.pairsTable.setScorer(scorer, rec._ag, lig._ag)
        self.pairsTable.setMaximumHeight(self.viewer.cameras[0].height()/2)

        self.showGridBoxAction.setChecked(False)
        self.showHideGridBox()
        if self.autofocus():
            self.pmv.focusScene(obj=self.ligand, padding=1)

    def onTablechange(self):
        self._interCyl.Set(visible=0)
        self._interSph.Set(visible=0)

    def onSelectionChange(self):
        from DejaVu2.stipples import stippleLines
        rows = set([x.row() for x in self.pairsTable.table.selectedItems()])
        if len(rows)==0:
            self.onTablechange()
            return
        c1 = self.receptor._ag.getCoords()
        c2 = self.ligand._ag.getCoords()
        c = []
        f = []
        col = []
        n = 0
        for row in rows:
            i,j =  self.pairsTable._pairs[int(self.pairsTable.table.item(row,0).text())]
            c.append(c1[i])
            c.append(c2[j])
            f.append((n, n+1))
            text = self.pairsTable.table.item(row,3).text()
            if text=='.' or float(text)>0.:
                col.append([1,0,0])
            else:
                col.append([0,1,0])
            n += 2

        v1, f1, fcol, vcol = stippleLines(numpy.array(c), f, col, segLen=0.3, spaceLen=0.2)

        self._interCyl.Set(vertices=v1, faces=f1, radii=(0.03,), visible=1, 
                           materials=vcol, tagModified=True) 
        self._interSph.Set(vertices=v1, radii=(0.03,), materials=vcol,
                           visible=1, tagModified=True)
        if self.autofocus():
            self.pmv.focusScene(obj=self._interCyl, padding=5)
        
        #import pdb; pdb.set_trace()
        
        ## adfr = self._mkADFR(self.agfr.ligand.filename, 0)
        ## ind = adfr.createPopulation(1)[0]
        ## ind.setGenes(ind.genomePy.getIdentityGenesPy())
        ## print 'SCORE', -ind.score()
        ## print "RR-L ---------------------------------------------------------------"
        ## print "name elem ene       x         y          z"
        ## ADelem = adfr.ligandFT.mol._ag._data['AD_element']
        ## genome = ind.genomePy#GenomePy(adfr.FT, adfr.scorer, scaleRE=scaleRE)
        ## RRL = genome.scorer.getLrrGridScorer().getScoreArray()
        ## for i, a1 in enumerate(adfr.ligandFT.mol.select()):
        ##     x, y ,z = a1.getCoords()
        ##     print "%4s %2s %9.3f %9.3f %9.3f %9.3f"%(
        ##         a1.getName(), ADelem[a1.getIndex()], RRL[i],x, y, z)
        ## print "RR-L END------------------------------------------------------------"
        ## return ind._score

    def testDocking(self):
        params = DockingTestParamsDialog(title="docking test parameters")
        params.exec_()
        maxEvals = params.paramsWidget.maxEvalsWidget.value()
        maxGens = params.paramsWidget.maxGensWidget.value()
        noImprov = params.paramsWidget.maxNoImprovWidget.value()
        seed = params.paramsWidget.seedWidget.value()
        popSize = params.paramsWidget.popSizeWidget.value()

        reference = ligandFile = self.agfr.ligand.filename
        adfr = self._mkADFR(ligandFile, seed)
        pop = adfr.createPopulation(popSize)
        adfr.createGA(pop, reference, RMSDMatching='1to1', savebestLig=True)
        if self._dockedLigand:
            self.pmv.deleteMolecule(self._dockedLigand)
        self._dockedLigand = lig = adfr.ligand
        self.pmv.addMolecule(lig)
        self.pmv.displaySB(lig)
        self.pmv.customColor(lig.select('element C'), [(.8,0.3,0.)],
                             geomsToColor=['sb',])
        
        from threading import Thread
        thread = Thread(target=adfr.ga.evolve, kwargs={'maxGens':maxGens, 'verbose':1, 'maxEvals':maxEvals*1000, 'noImproveStop':noImprov})
        thread.start()
        #adfr.ga.evolve(maxGens=maxGens, maxEvals=maxEvals*1000,
        #               noImproveStop=noImprov, verbose=3)
        
        self._posePlayingTimer = QtCore.QTimer(self)
        cb = CallbackFunction(self.playNextPose, adfr, lig)
        self._posePlayingTimer.timeout.connect(cb)
        self._posePlayingTimer.start(500)

    def playNextPose(self, adfr, lig):
        coords = adfr.ga.getBestLigCoords()
        lig._ag.setCoords(coords[1])
        lig.geomContainer.allCoords[:] = coords[1]
        event = RefreshDisplayEvent(molecule=lig)
        self.pmv.eventHandler.dispatchEvent(event)
        self.pmv.gui().viewer.selectionEventHandler()
        if adfr.ga.status == 'ended':
            self._posePlayingTimer.stop()
            msgBox = QtGui.QMessageBox(self)
            msgBox.setText('Docking search finished with status:\n "%s"\n'%adfr.ga._endingStatus)
            msgBox.setStandardButtons(QtGui.QMessageBox.Ok)
            ret =  msgBox.exec_()
            
    def toggleXtalMates(self):
        if self.showXtalMates.isChecked():
            if self._xmMats is None:
                self._xmMats = getCrystalMatesMatrices(self.receptor, cutoff=5.0)
            self.receptor.geomContainer.geoms['master'].Set(instanceMatrices=self._xmMats)
        else:
            self.receptor.geomContainer.geoms['master'].Set(instanceMatrices=numpy.identity(4))

    def selectAtomsInBox(self):
        # select receptor atoms in the docking box, sets inbox__ flag
        if self.receptor is not None:
            ag = self.receptor._ag
            ag.setFlags('inbox__', [False]*len(ag))
            selstr = 'x>%f and x<%f and y>%f and y<%f and z>%f and z<%f'%(
                self._LL[0], self._UR[0], self._LL[1],
                self._UR[1], self._LL[2], self._UR[2])
            inbox = self.receptor.select(selstr)
            ag[inbox.getIndices()].setFlags('inbox__', True)
            return self.receptor.select(selstr)

    def selectCAInBox(self):
        # get list of CA atoms for residues with moving side chain atoms in box
        if self.receptor is not None:
            ag = self.receptor._ag
            atoms = self.selectAtomsInBox()
            if atoms:
                ag.setFlags('inboxsel__', [False]*len(ag))
                for res in atoms.getHierView().iterResidues():
                    if res.select('ca inbox__'):
                        res.setFlags('inboxsel__',True)
                return self.receptor.select('inboxsel__')

                # this is much slower
                ## serials = set(atoms.getResindices())
                ## selstr = 'resindex %s'%' '.join([str(x) for x in atoms.getSerials()])
                ## #selstr = 'same residue as (serial %s) and ca'%serials
                ## inBox = self.receptor.select(selstr)
                ## caInBox = inBox.select('ca')
                ## if caInBox is not None:
                ##     return caInBox

    def getResWithAtomsInBox(self):
        # select all moving receptor side chain atoms inside the box
        # then expand to full residues

        #inBoxAtoms = self.selectAtomsInBoxNoNH()
        inBoxAtoms = self.selectAtomsInBox()
        if inBoxAtoms:
            ag = self.receptor._ag
            ag.setFlags('inboxsel__', [False]*len(ag))
            for res in inBoxAtoms.getHierView().iterResidues():
                res.setFlags('inboxsel__',True)
            return self.receptor.select('inboxsel__')
            # much slower version
            ## serials = set(inBoxAtoms.getResindices())
            ## selStr = 'resindex %s'%" ".join([str(x) for x in serials])
            ## resWithAtomsInBox = self.receptor.select(selStr)
            ## print 'S2', time()-t0
            ## return resWithAtomsInBox

    ## def selectMovingSideChainAtomsInBox(self):
    ##     # get list of atoms beyond CB for residues with moving side chain atoms in box
    ##     if self.receptor is not None:
    ##         atoms = self.selectAtomsInBoxNoNH()
    ##         if atoms:
    ##             serials = set(atoms.getResindices())
    ##             selstr = 'resindex %s and not bb and not name CB'%' '.join([str(x) for x in atoms.getSerials()])
    ##             inBox = self.receptor.select(selstr)
    ##             if inBox is not None:
    ##                 return inBox

    def selectAtomsInBoxNoNH(self):
        # select moving side chains atoms in the docking box
        if self.receptor is not None:
            selstr = '(not bb) and (not name CB) and (not resname PRO ALA GLY) and '
            selstr += 'x>%f and x<%f and y>%f and y<%f and z>%f and z<%f'%(
                self._LL[0], self._UR[0], self._LL[1],
                self._UR[1], self._LL[2], self._UR[2])
            inBox = self.receptor.select(selstr)
            if inBox is not None:
                return self.removeNH(inBox)
        
    def removeNH(self, atoms):
        # remove H atoms attached to N atoms, WHY ???? (MS)
        toRemove = []
        for atom in atoms:
            if atom.getElement()=='H':
                for natom in atom.iterBonded():
                    if natom.isbackbone:
                        toRemove.append(atom.getIndex())
        atoms = atoms - Selection(self.receptor._ag, toRemove, '')
        return atoms

    def pickATypes(self):
        self.paramsWidget.compLigandAtypes.setChecked(True)
        self.paramsWidget.atypesLabel.setDisabled(False)
        d = {}
        atstr = self._mapTypes
        for at in self.agfr.ADAtomTypes.keys():
            if at in atstr:
                d[at] = True
            else:
                d[at] = False
        dialog = AtomTypeSelector(d, self)
        result = dialog.exec_()
        atypesstr = dialog.getTypesString()
        self._mapTypes = atypesstr.split()
        for tt in ["OA", "HD"]: # these map types should always be in the map types list, since
            # they are going to be used to create W map file.
            if tt not in self._mapTypes:
                self._mapTypes.append(tt)
                #atypesstr += " %s"%tt # not sure if we need to display this map type in the
                # label since the user did not select it .
        self.setLigAtomTypesLabel(atypesstr)

    def getFlexResDialogInput(self):
        # convert selection string into dict with key 'Chid:ResnameResnum' and
        # value True is this residue in 
        resNames = []
        df = {}
        for res in self._flexRes:
            df['%c:%s%d'%(res.getChid(), res.getResname(), res.getResnum())] = True

        d = {}
        for atom in self.selectCAInBox():
            key = '%s:%s%d'%(atom.getChid(), atom.getResname(), atom.getResnum())
            if df.has_key(key):
                d[key] = True
            else:
                d[key] = False
        return d

    def flexResSelStrFromAtoms(self, atoms):
        # build flexRes selection string from an atom selection
        d = {}
        for atom in atoms:
            chid = atom.getChid()
            resnum = atom.getResnum()
            resname = atom.getResname()
            icode = atom.getIcode()
            if icode is None: icode=""
            if not d.has_key(chid):
                #d[chid] = ['%s%d'%(atom.getResname(),atom.getResnum())]
                d[chid] = [(resname, resnum, icode)]
            else:
                #d[chid].append('%s%d'%(atom.getResname(), atom.getResnum()))
                d[chid].append((resname, resnum, icode))
        chids = d.keys()
        chids.sort()
        s = ''
        descr = []
        for chid in chids:
            s += '%c:'%chid
            resl = {}
            for res in d[chid]:
                #s += '%s,'%res
                resstr = '%s%d%s'%res
                if not resl.has_key(resstr):
                    s += '%s%d%s,'%res
                    #resl.append((res[:3], res[3:]))
                    resl[resstr] = (res[0], str(res[1]), res[2])
            descr.append( (chid, resl.values()) )
            s = s[:-1]
            s += ';'
        s = s[:-1]
        return s, descr
    
    #def showFlexResChooserFromFP(self):
    #    print 'ASDAS'

    def showFlexResChooser(self):

        def setFR():
            #caAtoms = self._flexResDialog.resTreeWidget.getAtoms()
            #if caAtoms:
            #    caAtoms = caAtoms.select('ca')
            #    allResAtoms = self.receptor.select('same residue as index %s'%' '.join(
            #    [str(x) for x in caAtoms.getIndices()]))
            from time import time
            allResAtoms = self._flexResDialog.resTreeWidget.getAtoms()
            if allResAtoms:
                self._flexResAtoms = allResAtoms
                self._flexResSCAtoms = self.removeNH(allResAtoms.select('not bb and not ca'))
                #selStr, flexResDescr = self.flexResSelStrFromCA(caAtoms)
                selStr, flexResDescr = self.flexResSelStrFromAtoms(allResAtoms)
            else:
                self._flexResAtoms = []
                self._flexResSCAtoms = []
                selStr = ''
            if self._flexResAtoms:
                self.paramsWidget.setFlexResGrdButton.setDisabled(False)
            self._noFlexResParse = True
            self.paramsWidget.flexResWidget.setText(selStr)
            self._noFlexResParse = False
            self.agfr.setFlexResidues(selStr)
            self.flexResChanged.emit()
            
        def accept():
            self.paramsWidget.flexResWidget.setDisabled(False)
            self.paramsWidget.setFlexResGrdButton.setDisabled(False)
            self.paramsWidget.setFlexResButton.setChecked(False)
            self._openDialogs.remove(self._flexResDialog)
            self.freezeUI(False)
                
        def onClose():
            self.paramsWidget.setFlexResButton.setChecked(False)
            accept()

        if self.paramsWidget.setFlexResButton.isChecked():
            self.paramsWidget.flexResWidget.setDisabled(True)
            self.paramsWidget.setFlexResGrdButton.setDisabled(True)
            from ADFR.GUI.ResiduesTree import ResiduesTreeDialog
            caInbox = self.selectCAInBox()
            checked = {}
            for ca in caInbox:
                checked['%s:%d%s'%(ca.getChid(),ca.getResnum(), ca.getIcode())] = False
            if self._flexResSCAtoms:
                selstr, descr = self.flexResSelStrFromAtoms(self._flexResSCAtoms)
                for chid, residues in descr:
                    for resname, resnum, icode in residues:
                        checked['%s:%s%s'%(chid,resnum, icode)] = True

            self._flexResDialog = ResiduesTreeDialog(
                self, self.receptor.select(), title="flexible receptor chains",
                lig=self.ligand, checkedItems=checked, flexResiduesNames=self.agfr.flexResiduesNames,
                altloc=self._altlocResidues,
                mutated=self._mutatedResidues,
                mse=self.agfr._MSEtoMET)

            self._flexResDialog.resTreeWidget.resChanged.connect(setFR)
            self._flexResDialog.accepted.connect(accept)
            self._flexResDialog.closedSignal.connect(onClose)
            self._flexResDialog.show()
            self._openDialogs.append(self._flexResDialog)
            self.freezeUI(True)
        else:
            self.paramsWidget.flexResWidget.setDisabled(False)
            self.paramsWidget.setFlexResGrdButton.setDisabled(False)
            self._flexResDialog.accept()
            self._openDialogs.remove(self._flexResDialog)
            
    ## def defaultReceptorBox(self):
    ##     coords = self.receptor._ag.getCoords()
    ##     mini = numpy.min(coords, 0)
    ##     maxi = numpy.max(coords, 0)
    ##     lengths = (maxi-mini) + 2*self.paramsWidget.gridPaddingWidget.value()
    ##     center = mini + 0.5*lengths
    ##     self._baseSize[:] = lengths - 2*self.paramsWidget.gridPaddingWidget.value()
    ##     self.boxGeom.setCenter( *center)
    ##     self.boxGeom.setSides( *lengths)
    ##     self.onBoxChange()

    def showResiduesGridControls(self):
        def setBox():
            atoms = self._resControlWidget.resTreeWidget.getAtoms()
            if atoms:
                coords = atoms.getCoords()
                padding = self.paramsWidget.gridPaddingWidget.value()
                self.agfr.setBoxForCoords(coords, padding, self._spacing)
                self.agfr.data['boxMode'] =  "residues"
                self._baseSize[:] = self.agfr.boxLengths - 2*padding
                self.boxGeom.setCenter( *self.agfr.boxCenter)
                self.boxGeom.setSides( *self.agfr.boxLengths)
                self.onBoxChange()
            else:
                self.setGridFullReceptor()
        def accept():
            self.paramsWidget.toolBar1.setDisabled(False)
            self.paramsWidget.residuesGrdAct.setChecked(False)
            self._openDialogs.remove(self._resControlWidget)
            self.freezeUI(False)

        def onClose():
            self.paramsWidget.toolBar1.setDisabled(False)
            accept()

        if self.paramsWidget.manualGrdAct.isChecked():
            self.paramsWidget.manualGrdAct.setChecked(False)
            self.showManualGridControls()

        if self.paramsWidget.residuesGrdAct.isChecked():
            from ADFR.GUI.ResiduesTree import ResiduesTreeDialog
            
            self._resControlWidget = ResiduesTreeDialog(
                self, self.receptor.select(), self.ligand,
                title="grid receptor amino acids", sites=self._pdbSites)
            self._resControlWidget.resTreeWidget.resChanged.connect(setBox)
            self._resControlWidget.accepted.connect(accept)
            self._resControlWidget.closedSignal.connect(onClose)
            self._resControlWidget.show()
            self.paramsWidget.toolBar1.setDisabled(True)
            self._openDialogs.append(self._resControlWidget)            
            self.freezeUI(True)
        else:
            self.paramsWidget.toolBar1.setDisabled(False)
            self._resControlWidget.accept()
            
    def showManualGridControls(self):
        def accept():
            self.paramsWidget.toolBar1.setDisabled(False)
            self.paramsWidget.manualGrdAct.setChecked(False)
            self.paramsWidget.gridPaddingWidget.setDisabled(False)
            self._openDialogs.remove(self._manualControlDialog)
        def onClose():
            self.paramsWidget.toolBar1.setDisabled(False)
            self.paramsWidget.gridPaddingWidget.setDisabled(False)
            accept()
        if self.paramsWidget.manualGrdAct.isChecked():
            self.paramsWidget.gridPaddingWidget.setDisabled(True)
            self.paramsWidget.toolBar1.setDisabled(True)
            self._manualControlDialog = BoxParametersDialog(self)
            self._manualControlWidget = self._manualControlDialog.boxParamsWidget
            #self.paramsWidget.groupBox3Layout.addWidget(self._manualControlWidget)
            self._manualControlDialog.accepted.connect(accept)
            self._manualControlDialog.closedSignal.connect(onClose)
            self._manualControlDialog.show()
            self._openDialogs.append(self._manualControlDialog)
        else:
            self.paramsWidget.gridPaddingWidget.setDisabled(False)
            self.paramsWidget.toolBar1.setDisabled(False)
            self._manualControlDialog.accept()
            #self._manualControlWidget.deleteLater()
            self._openDialogs.remove(self._manualControlDialog)

    def setGridFlexRes(self):
        if not self._flexResAtoms: return
        #coords = self._flexResAtoms.getCoords()
        #mini = numpy.min(coords, 0)
        #maxi = numpy.max(coords, 0)
        if self.paramsWidget.manualGrdAct.isChecked():
            self.paramsWidget.manualGrdAct.setChecked(False)
            self.showManualGridControls()

        padding = self.paramsWidget.gridPaddingWidget.value()
        ## nc1 = numpy.min( (mini-padding, numpy.array(self._LL)), 0) # new corner1
        ## nc2 = numpy.max( (maxi+padding, numpy.array(self._UR)), 0) # new corne2
        ## nc = 0.5*(nc1 + nc2) # new center
        ## nl = (nc2 - nc1)
        ## self._baseSize[:] = nl - 2*self.paramsWidget.gridPaddingWidget.value()
        ## self.boxGeom.setCenter( *nc)
        ## self.boxGeom.setSides( *nl)
        flexResStr = self.agfr.data['flexResStr']
        self.agfr.setBox(["residues", flexResStr], padding, self._spacing)
        self._baseSize[:] = self.agfr.boxLengths - 2*padding
        self.boxGeom.setCenter( *self.agfr.boxCenter)
        self.boxGeom.setSides( *self.agfr.boxLengths)
        self.onBoxChange()

    def flexResOutOfBox_cb(self, val, doNotAsk):
        """callback of the pop up dialog when some of the flex residues are found out of the grid box """
        #print "flexResOutOfBox button pressed:", val  
        # val is in [0 -adjust box, 1 - make sidechains rigid, 2 - fix later] 
        if val == 0:
            self.adjustGridBoxFlexRes()
        elif val == 1:
            self.makeResOutsideBoxRigid()
        else: #fix later button or  "X" - close dialog 
            # set the value temporary so that the updateStatusBar does not call the dialog recursively
            self.preferences['Flexible Residues Out Of Box']['value'] = "None"
            self.updateStatusBar()

        if doNotAsk:
            if val == 0:
                self.preferences['Flexible Residues Out Of Box']['value'] = "Adjust box"
            elif val == 1:
                self.preferences['Flexible Residues Out Of Box']['value'] = "Make sidechains rigid"
            else:
                self.preferences['Flexible Residues Out Of Box']['value'] = "None"
        else:
            self.preferences['Flexible Residues Out Of Box']['value'] = "Always ask"
        if self._resOutsideBoxDialog and self._resOutsideBoxDialog in self._openDialogs:
            self._openDialogs.remove(self._resOutsideBoxDialog)
        # enable the flexible dialog if it is opened
        if self._flexResDialog and self._flexResDialog.isVisible():
            self._flexResDialog.setEnabled(True)

            
    def adjustGridBoxFlexRes(self):
        # make the box big enough to include all flex residues:
        #import pdb; pdb.set_trace()
        llx1, lly1, llz1 = self._LL
        urx1, ury1, urz1 = self._UR
        padding = self.agfr.padding
        flexResStr = self.agfr.data['flexResStr']
        self.agfr.setBox(["residues", flexResStr], padding, self._spacing)
        cx, cy, cz = self.agfr.boxCenter
        sx, sy, sz = self.agfr.boxLengths
        llx2, lly2, llz2 = [cx - (sx*.5), cy - (sy*.5), cz - (sz*.5)]
        urx2, ury2, urz2 = [cx + (sx*.5), cy + (sy*.5), cz + (sz*.5)]
        llx, lly, llz = [min(llx1, llx2), min(lly1, lly2),  min(llz1, llz2)]
        urx, ury, urz = [max(urx1, urx2),  max(ury1, ury2), max(urz1, urz2)]
        
        sx, sy, sz = [urx-llx, ury-lly, urz-llz]
        cx, cy, cz = [llx+(sx*.5), lly+(sy*.5), llz+(sz*.5)]
        self.agfr.setBoxCenter([cx, cy, cz])
        self.agfr.computeGridSize([sx, sy, sz])
        self._baseSize[:] = self.agfr.boxLengths - 2*padding
        self.boxGeom.setCenter( *self.agfr.boxCenter)
        self.boxGeom.setSides( *self.agfr.boxLengths)
        
        self.onBoxChange()

    def makeResOutsideBoxRigid(self):
        if not self._flexResAtoms:return
        #print "in makeResOutsideBoxRigid:", len(self._flexResAtoms)
        from ADFR.utils.MakeGrids import splitFlexRes
        flexResStr = self.agfr.data['flexResStr']
        flexresList =  flexResStr2flexRes(flexResStr) # fix it -- does not return Icode
        outsideRes = []
        for frchain in flexresList:
            ch = frchain[0]
            for fr in frchain[1]:
                fratoms = self.receptor.select('chid %s resnum `%d`'% (ch, fr[1]))
                #import pdb; pdb.set_trace()

                icode = fratoms.getIcodes()[0]
                inside, outside, inds = self.pointsInBox(fratoms.getCoords())
                if len(outside):
                    self._flexResAtoms = self._flexResAtoms - fratoms
                    outsideRes.append('%s:%s%d%s'%(ch, fr[0], fr[1], icode))
        #print "in makeResOutsideBoxRigid: 1", len(self._flexResAtoms)

        if len(self._flexResAtoms):
            self._flexResSCAtoms = self.removeNH(self._flexResAtoms.select('not bb and not ca'))
            selStr, flexResDescr = self.flexResSelStrFromAtoms(self._flexResAtoms)
        else:
            self._flexResAtoms = []
            self._flexResSCAtoms = []
            selStr = ''
                                      
        self._noFlexResParse = True
        self.paramsWidget.flexResWidget.setText(selStr)
        if hasattr(self, "_flexResDialog") and self._flexResDialog.isVisible():
            for resstr in outsideRes:
                resitem = self._flexResDialog.resTreeWidget._resToItem.get(resstr, None)
                self._flexResDialog.resTreeWidget._suspend = True
                resitem.setCheckState(0, QtCore.Qt.Unchecked)
                self._flexResDialog.resTreeWidget._suspend = False
        self._noFlexResParse = False
        self.agfr.setFlexResidues(selStr)
        self.flexResChanged.emit()

    def setGridFullReceptor(self):
        if self.paramsWidget.manualGrdAct.isChecked():
            self.paramsWidget.manualGrdAct.setChecked(False)
            self.showManualGridControls()
        padding = self.paramsWidget.gridPaddingWidget.value()
        self.agfr.setBox(["receptor"], padding, self._spacing)
        self._baseSize[:] = self.agfr.boxLengths - 2*padding
        self.boxGeom.setCenter( *self.agfr.boxCenter)
        self.boxGeom.setSides( *self.agfr.boxLengths)
        self.onBoxChange()
        
    def setGridFullLigand(self):
        if self.paramsWidget.manualGrdAct.isChecked():
            self.paramsWidget.manualGrdAct.setChecked(False)
            self.showManualGridControls()
        padding = self.paramsWidget.gridPaddingWidget.value()
        self.agfr.setBox(["ligand"], padding, self._spacing)
        self._baseSize[:] = self.agfr.boxLengths - 2*padding
        self.boxGeom.setCenter( *self.agfr.boxCenter)
        self.boxGeom.setSides( *self.agfr.boxLengths)
        self.onBoxChange()

    def setGridClusterPoints(self):
        if self.paramsWidget.manualGrdAct.isChecked():
            self.paramsWidget.manualGrdAct.setChecked(False)
            self.showManualGridControls()
        if self.paramsWidget.residuesGrdAct.isChecked():
            self.paramsWidget.residuesGrdAct.setChecked(False)
            self.showResiduesGridControls()
        coords = self.getClusterPoints()
        padding = self.paramsWidget.gridPaddingWidget.value()
        self.agfr.setBoxForCoords(coords, padding, self._spacing)
        self._baseSize[:] = self.agfr.boxLengths - 2*padding
        self.boxGeom.setCenter( *self.agfr.boxCenter)
        self.boxGeom.setSides( *self.agfr.boxLengths)
        self.onBoxChange(invalidateClusterPoints=False)

    def paddingChanged(self, value):
        self.agfr.setPadding(value)
        self._baseSize[:] = self.agfr.boxLengths - 2*value
        self.boxGeom.setSides(*self.agfr.boxLengths)
        self.onBoxChange()

    def setFlexRes(self):
        ## expect string of type chid1:Resnum1,Resnum2,Resnum3 ;  chid2:Resnum1,Resnum2,Resnum3
        ## if ":" is missing chid will match all chains
        ## residues can also be specified using ResnameResnum and can beseparated by space or ,
        if self._noFlexResParse:
            return
        self._flexRes = []
        #self._flexResStr = [] # list of (
        self._flexResAtoms = None # MolKit2 selection with all atoms in flexible residues
        self._flexResSCAtoms = None # MolKit2 selection with all moving side chain atoms in flexible resid

        frstr = self.paramsWidget.flexResWidget.text().encode('ascii', 'replace')
        if len(frstr)==0:
            self.agfr.setFlexResidues(None)
            #disable button that sets the box around flex residues
            self.paramsWidget.setFlexResGrdButton.setDisabled(True)
            self.flexResChanged.emit()
            return
        selStr, flexResAtoms  = self.parseResString(frstr)
        self._flexResAtoms = Selection(self.receptor._ag, flexResAtoms, '')
        if self._flexResAtoms:
            self._flexResSCAtoms = self.removeNH(self._flexResAtoms.select('not bb and not ca'))
        #print "FRA", len(self._flexResAtoms), len(self._flexResSCAtoms), self._flexResSCAtoms.getResnames()
        # check that grid covers flexRec atoms
        #flexRecInside, outsideAtoms = self.atomInBox(self._flexResSCAtoms)
        #if not flexRecInside:
        #    outSideResnames = ''
        #    resnums = outsideAtoms.getResnums()
        #    chids = outsideAtoms.getChids()
        #    done = {}
        #    outSideRes = []
        #    for chid, resnum in zip(chids, resnums):
        #        key = '%s:%d'%(chid, resnum)
        #        if done.has_key(key): continue
        #        done[key] = True
        #        res = hv.getResidue(chid, resnum)
        #        outSideRes.append(res)
        #        outSideResnames += '%s%d '%(res.getResname(), res.getResnum())
        #    msgBox = QtGui.QMessageBox(self)
        #    msgBox.setText('atoms for residues %s are outside the grid\nWould you like change the grid to cover these atoms ? If you choose NO this(ese) residue(s) will bre removed from the flexible residues list'%outSideResnames)
        #    msgBox.setStandardButtons(QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
        #    ret =  msgBox.exec_()

        #    if ret == QtGui.QMessageBox.Yes: # change grid to cover flexres atoms
        #        coords = self._flexResSCAtoms.getCoords()
        #        mini = numpy.min(coords, 0)
        #        maxi = numpy.max(coords, 0)
        #        center = [x.value() for x in self.centerSpinBoxes]
        #        size = [x.value() for x in self.sizefSpinBoxes]
        #        corner1 = center[0]-size[0]*.5, center[1]-size[1]*.5, center[2]-size[2]*.5
        #        corner2 = center[0]+size[0]*.5, center[1]+size[1]*.5, center[2]+size[2]*.5
        #        nc1 = numpy.min( (mini, corner1), 0) # new corner1
        #        nc2 = numpy.max( (maxi, corner2), 0) # new corne2
        #        nc = 0.5*(nc1 + nc2) # new center
        #        for i in range(3):
        #            self.centerSpinBoxes[i].setValue(nc[i])
        #            self.sizefSpinBoxes[i].setValue(nc2[i]-nc1[i]+2*self.gridPaddingWidget.value())
        #            self._baseSize = [x.value() for x in self.sizefSpinBoxes]
        #    else: # remove offending residues from flexRes list
        #        for res in outSideRes:
        #            self._flexRes.remove(res)
        #        selStr, descr = self.buildFlexRec(self._flexRes)
        #        #import pdb; pdb.set_trace()
        #        self.paramsWidget.flexResWidget.setText(selStr)
        #        self.flexResDialog.fill(self.getFlexResDialogInput())
        #        self.setFlexRes()

        #print self._flexRes
        #print self._flexResStr
        #print selStr
        # string might have been expanded to add residues names
        self.paramsWidget.flexResWidget.setText(selStr)
        self.agfr.setFlexResidues(selStr)
        if self._flexResAtoms:
            self.paramsWidget.setFlexResGrdButton.setDisabled(False)
        self.flexResChanged.emit()

    def parseResString(self, frstr):
        selStr = ""
        flexResStr = ""
        errors = []
        hv = self.receptor._ag.getHierView()
        #import pdb; pdb.set_trace()
        flexResAtoms = [] # list of atom indices for all atoms in FlexRes
        for expr in frstr.split(';'): #Loop over expressions delimited by ;
            dum = expr.split(':')
            if len(dum)==1:
                #chids = numpy.unique(self.receptor._ag.getChids())
                #residues = dum[0]
                #selStrTmp = ""
                errors.append('ERROR: invalid syntax\nShould be: chain1:resName1Resnum1,...,resNameNResnumN;\nchain2:resName1Resnum1, ...')
                continue
            else:
                chids = [dum[0]]
                residues = dum[1]
                selStrTmp = "%s:"%dum[0]
            for res in residues.replace(',', ' ').split():
                if res[0] in string.digits: # number only
                    resname = None
                    try:
                        resnum = int(res)
                    except ValueError:
                        
                        errors.append('ERROR: invalid syntax\n expect a residue number but got "%s"\n'%res)
                        continue
                else: # name and number
                    resname = res[:3]
                    icode = ""
                    try:
                        resnum = int(res[3:])
                    except ValueError:
                        resnumStr = ""
                        for i, st in enumerate(res[3:]):
                            if st.isdigit(): resnumStr += st
                            else: break
                        icode = res[3+i:]
                        if not len(resnumStr) or not icode.isalpha():
                            errors.append('ERROR: invalid syntax\n expect a residue number but got "%s"\n'%res[3:])
                            continue
                        resnum = int(resnumStr)
                    
                for chid in chids:
                    chainResStr = ""
                    if len(icode):
                        res = hv.getResidue(chid, resnum, icode=icode)
                    else:
                        res = hv.getResidue(chid, resnum)
                    if res is None:
                        errors.append('Residue %d not found in chain %s'%(
                            resnum, chid))
                    else:
                        flexResAtoms.extend(res.getIndices())
                        if resname:
                            if res.getResname()==resname:
                                self._flexRes.append(res)
                                #self._flexResStr.append( (res.getResname(), res.getResnum()) )
                                selStrTmp += "%s%s%s,"%(resname,resnum,icode)
                            else:
                                errors.append('Residue %d is %s which does not match specified %s'%(
                                    resnum, res.getResname(), resname))
                        else:
                            self._flexRes.append(res)
                            #self._flexResStr.append( (res.getResname(), res.getResnum()) )
                            selStrTmp += "%s%s%s,"%(res.getResname(),resnum, icode)
            if selStrTmp and selStrTmp[-1]==',': # i.e. there is at least 1 residue in that chain
                selStr += selStrTmp[:-1]+';' # remove trailing ","

        selStr = selStr[:-1] # remove trailing ";"
        if len(errors):
            msgBox = QtGui.QMessageBox(self)
            msgBox.setText('\n'.join(errors))
            msgBox.setStandardButtons(QtGui.QMessageBox.Ok)
            ret =  msgBox.exec_()
        return selStr, flexResAtoms
        
    def atomInBox(self, atoms):
        # check is atoms are in the grid box, If not identify atoms outside
        center = self.boxGeom.center
        size = self.boxGeom.sides
        p1x, p1y, p1z = center[0]-size[0]*.5, center[1]-size[1]*.5, center[2]-size[2]*.5
        p2x, p2y, p2z = center[0]+size[0]*.5, center[1]+size[1]*.5, center[2]+size[2]*.5
        coords = atoms.getCoords()
        mini = numpy.min(coords, 0)
        maxi = numpy.max(coords, 0)
        # MS ??? not sure this is right.
        if mini[0]<p1x or mini[1]<p1y or mini[2]<p1z or \
            maxi[0]>p2x or maxi[1]>p2y or maxi[2]>p2z:
            outside = atoms.select('x<%f or x>%f or y<%f or y>%f or z<%f or z>%f'%(
                p1x, p2x, p1y, p2y, p1z, p2z))
            return False, outside
        return True, None

    def toggleComputeAll(self):
        if self.paramsWidget.compAllButton.isChecked():
            #self.atypesEntry.setDisabled(True)
            self.paramsWidget.atypesLabel.setDisabled(True)
        else:
            #self.atypesEntry.setDisabled(False)
            self.paramsWidget.atypesLabel.setDisabled(True)

#######################################################################END
    
            
    ## def setAtypes(self):
    ##     #atstring = self.atypesEntry.text()
    ##     atstring = self.paramsWidget.atypesLabel.text()
    ##     atstring = atstring.replace(',', ' ')
    ##     atypes = []
    ##     for atype in atstring.split():
    ##         if not self.agfr.ADAtomTypes.has_key(atype):
    ##             msgBox = QtGui.QMessageBox(self)
    ##             msgBox.setText("ERROR: %s is not a valid AutoDock5 atom type"%atype)
    ##             msgBox.setStandardButtons(QtGui.QMessageBox.Ok)
    ##             ret =  msgBox.exec_()
    ##         else:
    ##             atypes.append(atype)
    ##     #self.atypesEntry.setText(' '.join(atypes))
    ##     self.paramsWidget.atypesLabel.setText(' '.join(atypes))
    ##     self.compLigandAtypes.setChecked(True)

    def onFlexResChanged(self):
        if self._flexResSCAtoms:
            fratypes = numpy.unique(self._flexResSCAtoms.getData('AD_element')).tolist()
            self._mapTypes = numpy.unique(self._mapTypes + fratypes).tolist()
            self.setLigAtomTypesLabel(' '.join(self._mapTypes))
            self.paramsWidget.groupBox4.setTitle(self.tr("[%d flexible residue(s)]"%len(self._flexResAtoms.select('ca'))))
        else:
            self.paramsWidget.groupBox4.setTitle(self.tr("[0 flexible residue(s)]"))
        self.updateStatusBar()        
        self.pmv.undisplaySB(self.receptor)
        if self._flexResSCAtoms:
            self.pmv.displaySB(self._flexResAtoms)
            self.pmv.customColor(self._flexResAtoms.select('element C'), [(1.,0.5,0.)],
                                 geomsToColor=['sb',])
            gc = self.receptor.geomContainer
            # FIXME for some reason the sticks and balls get the transparency of the lines
            gc.geoms['sticks'].Set(transparent=False)
            gc.geoms['balls'].Set(transparent=False)
            #gc.geoms['sticks'].Set(opacity=0.8)
            #gc.geoms['balls'].Set(opacity=0.8)
        if self._flexResAtoms:
            numFR = self._flexResAtoms.getHierView().numResidues()
            if numFR > 8:
                if self.preferences['Number of Flexible Residues']['value'] == "Always ask":
                    msgDial = CheckParamsInBoxDialog("There are currently, %d flexible sidechains. AutoDock can realistically handle ~4 and ADFR up to 15."%numFR, auto="Do not show this again", title="Number of Flexible Residues")
                    res = msgDial.exec_()
                    if res[1] == True:
                        self.preferences['Number of Flexible Residues']['value'] = "None"
        self.updateDisplay()

    def showUserPreferences(self):
        def onClose():
            self.topBar.userPreferencesButton.setChecked(False)

        if self.topBar.userPreferencesButton.isChecked():
            from mglutil.gui.Qt.preferencesGui import PreferencesTreeDialog
            self._userPreferencesDialog =  PreferencesTreeDialog(self.preferences, self, title="User Preferences")
            self._userPreferencesDialog.closedSignal.connect(onClose)
            self._userPreferencesDialog.show()
            self._openDialogs.append(self._userPreferencesDialog)
        else:
            self._userPreferencesDialog.reject()
            if self._userPreferencesDialog in self._openDialogs:
                self._openDialogs.remove(self._userPreferencesDialog)

                
    def clearClustersWidget(self):
        clustersWidget = self.paramsWidget.clustersWidget
        clustersWidget.clear()
        clustersWidget.setHorizontalHeaderLabels(
            ['fills', 'AS Score', '#Points', 'RadGyr','Buriedness'])
        for i in range(clustersWidget.rowCount()):
            clustersWidget.removeRow(0)
        clustersWidget.resizeColumnsToContents()
            
##    @waiting_effects
    def AutoSiteFill(self, spacing=None, folder=None, background=False):

        if spacing is None:
            spacing = self._spacing
        self.__spacing = spacing
        self.__folder = folder
        #print "AutoSiteFill", "center", self.agfr.boxCenter[:], "boxGUI center",self.boxGeom.center, "length", length, "boxGUI sides", self.boxGeom.sides
        flexResStr = None
        if self._flexResAtoms:
            flexResStr = self.agfr.data['flexResStr']
        #background=False
        #import pdb; pdb.set_trace()
        gc, status = self.agfr.runAutoSite(flexResStr=flexResStr, spacing=spacing,
                                            background=background, outlev=2)
        process = gc.process
        self._gc = gc
        if gc is None: return
        if background:
            logFile = os.path.join(self.agfr.tmpFolder, gc.logFile)
            logReader = readLogThread(process, logFile, self, postProcess=self._afterAutoSite)
            logReader.progress.connect(self.updateASText, QtCore.Qt.QueuedConnection)
            logReader.progressBar.connect(self.updateASProgress, QtCore.Qt.QueuedConnection)
            logReader.errorSignal.connect(self.handleAutoGridError, QtCore.Qt.QueuedConnection)
            logReader.start()
            self.logReader = logReader # to avoid garbage collection
            self.paramsWidget.computeTPointsButton.setDisabled(True)
        else:
            self._afterAutoSite()

    def updateASText(self,text):
        if text=='readLogThread THREAD DONE':
            self.paramsWidget.progressBar.setValue(0)
        elif text == 'AutoGrid failed':
            self.paramsWidget.progressBar.setValue(0)
            self.paramsWidget.computeTPointsButton.setDisabled(False)
        else:
            self.statusbar.showMessage(text, level='busy')

    def updateASProgress(self, value):
        self.paramsWidget.progressBar.setValue(value)
            
    def _afterAutoSite(self):
        gc = self._gc
        self.logReader.progress.emit('getting high affinity points ...')
        spacing = self.__spacing
        if self._autoSite2 == True:
            self.agfr.setAutoSiteVersion("1.1", ligandSize=self._ligandSize, pepScore=self._pepScore)
            #run AutoSite2 for pocket detection
            #print "AUTOSITE 2", self._ligandSize, self._pepScore
            clustersorted, clPropsorted, dcl = self.agfr.afterAutoSite2(gc, ligandSize=self._ligandSize, pepScore=self._pepScore, verbose=True)
        else:
            self.agfr.setAutoSiteVersion("1.0")
            clustersorted, clPropsorted, dcl = self.agfr.afterAutoSite(gc, spacing=spacing, cutoff=10, verbose=True)
        #gc.getASPoints()
        #dcl = DensityClustering([spacing,spacing,spacing])
        #dcl.findClustersD(gc._indices, cVolcut=0)
        del self.__spacing    
        #self.logReader.progress.emit('found %d TPoints in %d clusters'%(len(gc._coords),
        #                                                                len(dcl._clusters)))
        npoints = sum([len(cl[1]) for cl in clustersorted])
        
        self.logReader.progress.emit('found %d TPoints in %d clusters'%( npoints, len(clustersorted)))
        self.colorMap = [(0,1,0), (1,0.5,0), (1,1,0), (0,1,1), (1,0,1), ]
        #colors = [(1,1,1)]*len(gc._indices)
        ncol = len(gc._indices)
        colors = numpy.ones(3*ncol).reshape(ncol,3) 
        nbLargeCl = 0
        #print "len clustersorted:", len(clustersorted)
        #print "num coords :", len(gc._coords)
        #print "len dcl._clen, dcl._clusters:", len(dcl._clen) , dcl._clen, len(dcl._clusters)
        
        #self._clusters = [[gc._indices, gc._coords, gc._potential,gc._atype]]
        self._clusters = [[]]
        #lclusters, clProp = scoreClusters(self.receptor, dcl, gc)
        #clPropsort = []
        #if len(lclusters)!=0 :
        #    tmpclsort = sorted(lclusters,key=lambda x:x[4],reverse = True)
        #    clPropsort = sorted(clProp,key=lambda x:x[5],reverse = True)

        if not self._autoSite2:
            for x in clustersorted:
                if len(x[0])>=self.smallClusterCutOff:
                    insidePts, outPoints, inds = self.pointsInBox(x[1])
                    
                    #self._clusters.append([x[0], x[1], x[2], x[3]])
                    self._clusters.append([x[0][inds], insidePts , x[2][inds], x[3][inds]])
            smallClustersCoords = []
            smallClustersInds = []
            smallClustersPots = []
            smallClustersAtypes = []
            for length, clinds in zip(dcl._clen, dcl._clusters):
                if length<self.smallClusterCutOff:
                    color = [1,0,0]
                    smallClustersInds.extend( gc._indices[clinds] )
                    smallClustersCoords.extend( gc._coords[clinds] )
                    smallClustersPots.extend( gc._potential[clinds] )
                    smallClustersAtypes.extend( gc._atype[clinds] )
                else:
                    color = self.colorMap[nbLargeCl%len(self.colorMap)]
                    nbLargeCl += 1
                for ind in clinds:
                    colors[ind] = color
            if len(smallClustersCoords):
                insideSmPts, outPoints, inds = self.pointsInBox(smallClustersCoords)
                self._clusters.append( [smallClustersInds[inds], insideSmPts, smallClustersPots[inds], smallClustersAtypes][inds])
            else:
                self._clusters.append([[],[],[],[]])
            # remove the the points that are outside the box
            insidePts, outPoints, insideInds = self.pointsInBox(gc._coords)
            self._clusters[0] = [gc._indices[insideInds], insidePts, gc._potential[insideInds], gc._atype[insideInds]]
            self._allClusterColors = colors[insideInds]
            
        else: # AutoSite2
            from AutoSite.shrink import shrinkPocket
            
            for x in clustersorted:
                shrinkedPoints = shrinkPocket(x[1],finalsize=len(x[1])/5)
                
                self._clusters.append([x[0], self.pointsInBox(shrinkedPoints)[0], x[2] ,x[3]])
            self._clusters.append([[],[],[],[]])
            self._allClusterColors = []
        
        self._ASclProp =  clPropsorted
        self._ASfolder = self.__folder
        del self.__folder

        self.paramsWidget.computeTPointsButton.setDisabled(False)
        self.AutoSiteDone.emit()
        if hasattr(self.agfr, "tmpFolder"):
            if self.agfr.tmpFolder and os.path.exists(self.agfr.tmpFolder):
                shutil.rmtree(self.agfr.tmpFolder)

    def fillFillsTable(self):
        clustersWidget = self.paramsWidget.clustersWidget
        self.clearClustersWidget()

        clProp = self._ASclProp
        folder = self._ASfolder
        del self._ASclProp
        del self._ASfolder
        clustersWidget.setHorizontalHeaderLabels(
            ['fills', 'AS Score', '#Points', 'RadGyr','Buriedness'])

        self._clusterItems = []
        #print 'time 6', time()-t0
        #t0 = time()
        # select clusters with more than self.smallClusterCutOff points
        row=0
        for inds, coords, pots, atypes in self._clusters[1:]:
            if len(coords)==0:
                continue
            clustersWidget.insertRow(row)
            #if row == 0:
            #    data = ['all', '-', len(inds), '-', '-']
                
            if row < len(self._clusters)-2:
                #data = [row+1, clProp[row][5], len(inds), clProp[row][3], clProp[row][4]]
                data = [row+1, clProp[row][5], len(coords), clProp[row][3], clProp[row][4]]
            else:
                #data = ['small clusters (<%d)'%self.smallClusterCutOff, '-', len(inds), '-','-']
                data = ['small clusters (<%d)'%self.smallClusterCutOff, '-', len(coords), '-','-']

            for col, val in enumerate(data):
                #if col !=0 and row ==1:
                #    clustersWidget.insertColumn(j)
                if isinstance(val, float):
                    newItem = TableWidgetItem("%.2f"%val)
                elif isinstance(val, int):
                    newItem = TableWidgetItem(str(val))
                elif isinstance(val, str):
                    newItem = QtGui.QTableWidgetItem(val)
                clustersWidget.setItem(row, col, newItem)
                if col == 0:
                    self._clusterItems.append(newItem)
                    newItem.setCheckState(QtCore.Qt.Unchecked)
                    #if row==0:
                    #    newItem.setCheckState(QtCore.Qt.Checked)
                    #    self._allClustersItem = newItem
            row += 1

        #vlabs = ['All']+[str(i) for i in range(row-1)]
        #print 'VLABS', vlabs
        #clustersWidget.setHorizontalHeaderLabels(vlabs)
        clustersWidget.resizeColumnsToContents()
        #clustersWidget.setSortingEnabled(True)

        #self.TPoints.Set(visible=True, vertices=gc._coords, materials=colors)
        #self.TPointsOK.emit(True)

        # show best scoring cluster by default
        self._clusterItems[0].setCheckState(QtCore.Qt.Checked)
        self.showHideCluster(self._clusterItems[0])
        #print 'removing', folder
        shutil.rmtree(folder)

    def selectAllPockets(self):
        clustersWidget = self.paramsWidget.clustersWidget
        if clustersWidget.rowCount() == 0: return
        colors = []
        coords = []
        for n in range(clustersWidget.rowCount()):
            item = clustersWidget.item(n,0)
            name = item.text()
            item.setCheckState(QtCore.Qt.Checked)
            try: cn = int(name)
            except:
                #if name.startswith('small') or name.startswith('from file'):
                cn = len(self._clusterItems)
            if not len(self._allClusterColors):
                c = self._clusters[cn][1]
                colors.extend( [self.colorMap[(cn-1)%len(self.colorMap)]]*len(c) )
                coords.extend(c)
        if len(self._allClusterColors):
            coords = self._clusters[0][1]
            colors = self._allClusterColors
        #print "IN SELECT ALL POCKETS, len(coords):", len(coords)
        self.TPoints.Set(visible=True, vertices=coords,
                         materials=colors)
        self.agfr.setFillPoints(coords)
        self.updateStatusBar()

    def deselectAllPockets(self):
        clustersWidget = self.paramsWidget.clustersWidget
        if clustersWidget.rowCount() == 0: return
        for n in range(clustersWidget.rowCount()):
            item = clustersWidget.item(n,0)
            name = item.text()
            item.setCheckState(QtCore.Qt.Unchecked)
        self.TPoints.Set(visible=False)
        self.agfr.setFillPoints([])
        self.updateStatusBar()

    def invertPocketSelection(self):
        clustersWidget = self.paramsWidget.clustersWidget
        if clustersWidget.rowCount() == 0: return
        colors = []
        coords = []
        for n in range(clustersWidget.rowCount()):
            item = clustersWidget.item(n,0)
            name = item.text()
            if item.checkState() == QtCore.Qt.CheckState.Checked:
                item.setCheckState(QtCore.Qt.Unchecked)
            else:
                item.setCheckState(QtCore.Qt.Checked) 
                if name.startswith('small'):
                    cn = len(self._clusterItems)
                else:
                    cn = int(name)
                c = self._clusters[cn][1]
                
                colors.extend( [self.colorMap[(cn-1)%len(self.colorMap)]]*len(c) )
                coords.extend(c)
        if not len(colors): colors = [[1,1,1]]
        self.TPoints.Set(visible=True, vertices=coords,
                         #materials=self._allClusterColors)
                         materials=colors)
        self.agfr.setFillPoints(coords)
        self.updateStatusBar()
        
    def computeTPoints(self):
        # create scratch folder name
        #self.clearClustersWidget()
        #self._clusterItems = []
        temp_name = next(tempfile._get_candidate_names())
        folder = os.path.join(self.tmpFolder, temp_name)
        os.mkdir(folder)
        self.statusbar.showMessage('running AutoSite', 'busy')
        self.AutoSiteFill(spacing=1.0, folder=folder, background=True)
        ## self.statusbar.showMessage('AutoSite is Done')

    ## def computeFPoints(self):
    ##     #self.TPoints.Set(visible=0)
    ##     #self.TPointsOK.emit(False)
    ##     clusterPoints = self.getClusterPoints()
    ##     clusterPotentials = self.getClusterPotentials()
    ##     clusterAtomtypes = self.getClusterAtypes()
    ##     if len(clusterPoints)==0:
    ##         msgBox = QtGui.QMessageBox(self)
    ##         msgBox.setText("ERROR: no translation points selected")
    ##         msgBox.exec_()
    ##         return
    ##     ccoords = [x for i,x in enumerate(clusterPoints) if clusterAtomtypes[i] == 'C']
    ##     cpot = [x for i,x in enumerate(clusterPotentials) if clusterAtomtypes[i] == 'C']
    ##     ocoords = [x for i,x in enumerate(clusterPoints) if clusterAtomtypes[i] == 'O']
    ##     opot = [x for i,x in enumerate(clusterPotentials) if clusterAtomtypes[i] == 'O']
    ##     hcoords = [x for i,x in enumerate(clusterPoints) if clusterAtomtypes[i] == 'H']
    ##     hpot = [x for i,x in enumerate(clusterPotentials) if clusterAtomtypes[i] == 'H']

    ##     fp = featurePts(self.receptor,ccoords, cpot, ocoords, opot, hcoords, hpot)

    ##     self.TSpheresC.Set(visible=True, vertices=fp.cfp, materials=[(0,1,0)], radii=[0.3])
    ##     self.TSpheresO.Set(visible=True, vertices=fp.ofp, materials=[(1,0,0)], radii=[0.3])
    ##     self.TSpheresH.Set(visible=True, vertices=fp.hfp, materials=[(1,1,1)], radii=[0.3])
    ##     self.showFPAction.setDisabled(False)
    ##     self.showFPAction.setChecked(True)
        
        ## fresdict = {}
        ## for x in fp.bsres:
        ##     xsplit = x.split(':')
        ##     fresdict.setdefault(xsplit[0], [])
        ##     fresdict[xsplit[0]].append(xsplit[1])

        ## flexResStr = ''
        ## for k,v in fresdict.iteritems():
        ##     resstr =  "".join(str(x)+',' for x in v)
        ##     flexResStr = flexResStr+k+':'+resstr+';'
            

        ## #print flexResStr
        ## self.paramsWidget.flexResWidget.setText(flexResStr)
        ## self.setFlexRes()
        ## #self.flexResChanged.emit() 
        
    def pointsInBox(self, pts):
        inside = []
        outside = []
        insideInds = []
        llx, lly, llz = self._LL
        urx, ury, urz = self._UR        
        #for x,y,z in pts:
        for i, pt  in enumerate(pts):
            x, y, z = pt
            if x > llx and x < urx and y > lly and y < ury and z > llz and z < urz:
                inside.append( (x,y,z) )
                insideInds.append(i)
            else:
                outside.append( (x,y,z) )
        return inside, outside, insideInds

    def checkTPInsideBox(self, gridPoints):
        llx, lly, llz = self._LL
        urx, ury, urz = self._UR
        for x,y,z in gridPoints:
            if x > llx and x < urx and y > lly and y < ury and z > llz and z < urz:
                return True
        return False

    def showHideCluster(self, item):
        clustersWidget = self.paramsWidget.clustersWidget
        coords = []
        ## if id(item) == id(self._allClustersItem): # all
        ##     if item.checkState() == QtCore.Qt.CheckState.Unchecked: # we clicked to uncheck
        ##         self.TPoints.Set(visible=False)
        ##     else: # we clicked to check ALL clusters
        ##         for n in range(clustersWidget.rowCount()):
        ##             if id(clustersWidget.item(n,0))== id(self._allClustersItem):
        ##                 continue
        ##             clustersWidget.item(n,0).setCheckState(QtCore.Qt.Unchecked)
        ##         coords = self._clusters[0][1]
        ##         self.TPoints.Set(visible=True, vertices=coords,
        ##                          materials=self._allClusterColors)
                
        #else:
            #self._allClustersItem.setCheckState(QtCore.Qt.Unchecked)
        coords = []
        colors = []
        for n in range(clustersWidget.rowCount()):
            item = clustersWidget.item(n,0)
            if item.checkState()==QtCore.Qt.Checked:
                #if id(clustersWidget.item(n,0))== id(self._allClustersItem):
                #    continue
                name = item.text()
                if name.startswith('small') or name.startswith('from'):
                    cn = len(self._clusterItems)
                else:
                    cn = int(name)
                c = self._clusters[cn][1]
                coords.extend( c )
                colors.extend( [self.colorMap[(cn-1)%len(self.colorMap)]]*len(c) )

        if coords:
            self.TPoints.Set(visible=True, vertices=coords, materials=colors)
            self.paramsWidget.TPGridAct.setDisabled(False)
        else:
            self.TPoints.Set(visible=False)
            self.paramsWidget.flexResWidget.setText('')
            self.setFlexRes()
        self.agfr.setFillPoints(coords)
        self.updateStatusBar()
            
    def isBoxOnReceptor(self):
        if not self.gridGUI.isBoxOnReceptor():
            msgBox = QtGui.QMessageBox(self)
            msgBox.setText("WARNING: Grid box does not overlap with receptor")
            msgBox.exec_()
            return False
        else:
            return True

    def getClusterPoints(self):
        clusterPoints = [] # cluster points coordinates
        
        for i, item in enumerate(self._clusterItems):
            if item.checkState() == QtCore.Qt.Checked:
                #clusterPoints.extend( self._clusters[i][1].tolist() )
                points = self._clusters[i+1][1]
                clusterPoints.extend( self._clusters[i+1][1] )
        return clusterPoints

    def getClusterPotentials(self):
        item0 = self._clusterItems[0]
        if item0.checkState() == QtCore.Qt.Checked: # all
            clusterPotentials = self._clusters[0][2] # cluster points potentials
        else:
            clusterPotentials = [] # cluster points coordinates
            i = 1
            #import pdb; pdb.set_trace()
            for item in self._clusterItems[1:]:
                if item.checkState() == QtCore.Qt.Checked:
                    #clusterPoints.extend( self._clusters[i][1].tolist() )
                    clusterPotentials.extend( self._clusters[i][2] )
                i += 1
        return clusterPotentials

    def getClusterAtypes(self):
        item0 = self._clusterItems[0]
        if item0.checkState() == QtCore.Qt.Checked: # all
            clusterAtomtypes = self._clusters[0][3] # cluster points atom types
        else:
            clusterAtomtypes = [] # cluster points coordinates
            i = 1
            #import pdb; pdb.set_trace()
            for item in self._clusterItems[1:]:
                if item.checkState() == QtCore.Qt.Checked:
                    #clusterPoints.extend( self._clusters[i][1].tolist() )
                    clusterAtomtypes.extend( self._clusters[i][3] )
                i += 1
        return clusterAtomtypes

    def setAutoSiteVersion(self, val):
        # val is True or False
        self._autoSite2=val
        self.paramsWidget.computeTPointsButton.setText("compute pockets [AutoSite 1.%d]"%val)
        
    def setAutoSiteLigSize(self, val):
        #val is an integer
        self._ligandSize = val
        self.paramsWidget.computeTPointsButton.setText("compute pockets [AutoSite 1.%d]"%val)
        
    def setAutoSitePepScoreFunc(self, val):
        # val is True or False
        self._pepScore = val

    def setWaterMapEntropy(self, val):
        self.agfr.setWmapParams(ENTROPY=val)

    def setWaterMapWeight(self, val):
        self.agfr.setWmapParams(weight=val)

    def computeGrids(self):
        self.__to = time()
        recName = os.path.splitext(os.path.split(self.receptor.filename)[1])[0]
        if os.path.isdir(recName):
            recName = os.path.abspath(os.path.dirname(recName))
        # ask for file name
        filename, filtr = QtGui.QFileDialog.getSaveFileName(
            self, "Target file name", recName, "Target file (*.trg)")
        if not filename: return
        fi = QtCore.QFileInfo(filename)
        self.destinationFolder = fi.completeBaseName().encode('ascii', 'replace')       
        #self.destinationFolderPath = fi.path().encode('ascii', 'replace')
        qd = QtCore.QDir()
        self.destinationFolderPath = qd.toNativeSeparators(fi.path()).encode('ascii', 'replace')

        folder = os.path.join(self.destinationFolderPath, self.destinationFolder)
        if not self._covalentDocking: 
            fillPoints = self.pointsInBox(self.getClusterPoints())[0]
            if len(fillPoints)==0:
                msgBox = QtGui.QMessageBox(self)
                msgBox.setText("ERROR: no binding pocket points selected")
                msgBox.exec_()
                return
            
            self.agfr.setFillPoints(fillPoints)
        else: #covalent docking
            at1, at2 = self._covalentBondAtoms
            at3 = self._torsionAtom
            torsionAtIndex = at3.getSerial()
            self.agfr.data['covalentBondAtom1'] = self.receptor.atomFullName(at1)
            self.agfr.data['covalentBondAtom2'] = self.receptor.atomFullName(at2)
            self.agfr.data['covalentBondTorsionAtom'] = '%s (%d)'%( self.receptor.atomFullName(at3), torsionAtIndex)
        
            self.agfr.data['covalentAtomsCoords'] = at1.getCoords().tolist()+ at2.getCoords().tolist()+ \
                                           at3.getCoords().tolist()
            self.agfr.data['covalentBond'] = [at1.getSerial(), at2.getSerial()]
            self.agfr.covalentBond = True
            self.agfr.data['covalentRes'] = None
            if self._covalentLigAtoms:
                ## if len(self._covalentLigAtoms):
                ##     newRecSel = self.receptor.select() - self._covalentLigAtoms
                ##     self.agfr.receptor = newRecSel.toMolecule(self.receptor.name)
                ##     print "AGFR cov docking:", "old rec", len(self.receptor._ag), "new rec", len(self.agfr.receptor._ag), "atms to exclude", len(self._covalentLigAtoms)
                self.agfr.covalentBondToExclude = self._covalentLigAtoms
            
        if self.paramsWidget.compAllButton.isChecked():
                self.agfr.setMapTypes("all")
        else:
            self.agfr.setMapTypes(self._mapTypes)
        flexResStr = None
        if self._flexResAtoms:
            flexResStr = self.agfr.data['flexResStr']
        self.paramsWidget.computeGridsButton.setDisabled(True)
        #self.paramsWidget.computeGridsButton.setText('computing affinity maps')
        self.statusbar.showMessage('computing affinity maps', level='busy')
        gc, status = self.agfr.computeGrids(folder, flexResStr,
                                            self._spacing, background=True)
        if status is not 0:
            self.statusbar.showMessage('%s' %(status,) , level='error')
            self.paramsWidget.computeGridsButton.setDisabled(False)
            return
        process  = gc.process
        self._gc = gc
        self._selStr=flexResStr
        logFile = os.path.join(folder, gc.logFile)
        logReader = readLogThread(process, logFile, self, postProcess=self.mkTrg)
        logReader.progress.connect(self.updateText, QtCore.Qt.QueuedConnection)
        logReader.progressBar.connect(self.updateProgress, QtCore.Qt.QueuedConnection)
        logReader.errorSignal.connect(self.handleAutoGridError, QtCore.Qt.QueuedConnection)

        if not logReader.isRunning():
            logReader.start()
        self.logReader = logReader
        #self.paramsWidget.computeGridsButton.setText('computing maps ...')
        self.statusbar.showMessage('computing affinity maps', level='busy')
        

    def handleAutoGridError(self, msg):
        msgBox = QtGui.QMessageBox()
        msgBox.setText(msg)
        msgBox.exec_()
        self.statusbar.showMessage(msg, level='error')


    def mkTrg(self):
        self.__to = time()
        folder = os.path.join(self.destinationFolderPath , self.destinationFolder)
        logthread = None
        logFileName = None
        if self.paramsWidget.addGradGroupBox.isChecked():
            if self.paramsWidget.gradLargeCl.isChecked():
                cutoffval = -1
            else:
                cutoffval = int(self.paramsWidget.gradCutOffBox.value())
            self.agfr.cutOffValue = cutoffval
            #print "Adding gradients to maps, cutoff:", cutoffval
            #logFileName = next(tempfile._get_candidate_names())
            # create a temporary file that will be used by the C++ function computing receptor gradient
            # to output the maptype names as it processes them.
            # This file will be read by progressLogFile.run() method which will update the progress bar.
    
            logFileName = tempfile.mkstemp(text=True)[1]
            ntypes = 0
            # find out how many map types will be processed:
            for atype in self.agfr.data['mapTypes']:
                if atype not in ['e', 'd', 'sd']:
                    ntypes += 1
            # this is a new QThread object
            logthread = progressLogFile(logFileName, ntypes, self, doneMsg="generating .trg file ...")
            #self.logReader.progress.emit("adding gradients to maps and generating .trg file ...")
            logthread.progress.connect(self.updateText, QtCore.Qt.QueuedConnection)
            logthread.progressBar.connect(self.updateProgress, QtCore.Qt.QueuedConnection)
            logthread.progress.emit("adding gradient to maps ...")
            logthread.start()
            self.logthread = logthread
            self.agfr.generateTrgFile(self._gc, folder, self._selStr, addGradients=True, logFileName=logFileName)#, callback=self.updateProgress)
            
        else:
            self.logReader.progress.emit('generating .trg file ...')
            self.agfr.generateTrgFile(self._gc, folder, self._selStr)
        self.paramsWidget.computeGridsButton.setDisabled(False)
        trgFile = os.path.abspath(os.path.join(self.agfr.destinationFolderPath,  self.agfr.destinationFolder+'.trg'))
        if os.path.exists(trgFile):
            if not trgFile in self.trgFiles:
                self.trgFiles.append(trgFile)
                if self.trgMapGui and self.trgMapGui.isVisible():
                    self.trgMapGui.addTrg(trgFile)
        self.showTrgGuiAction.setDisabled(False)
        if logthread and logthread.isRunning():
            logthread.stop()
        if logFileName and os.path.exists(logFileName):
            #print "trying to remove tmp log file:", logFileName
            try:
                os.remove(logFileName)
            except:
                #On windows could not remove the file: Get an error
                # that the file is in use by other process
                #print "could not remove omp logfile file:" , logFileName
                pass


    def addFlexRecHeader(self, gc, selStr):
        gc.addFlexRecHeader('FLEXRES "%s"'%selStr)

    @waiting_effects
    def finishUp(self):
        self.paramsWidget.computeGridsButton.setText('generate target file ... ')
        filename = os.path.join(self.destinationFolderPath,
                                self.destinationFolder)+'.trg'
        self._targetFile = filename
        filesize = os.path.getsize(filename)/1048576. # in Mb
        self.statusbar.showMessage('Generated target file %s (%.2fMb) in %.2f(s)'%(
            filename, filesize, time()-self.__to))
        del self.__to
        self._haveTrg = True
        if self.ligand:
            self.testDockAct.setDisabled(False)
            self.viewInteractionsAct.setDisabled(False)
        
    def updateText(self,text):
        #print '++++++++++++++++++++++++++++++++++++++++++'
        #print text
        #self.text.append(text)
        if text=='readLogThread THREAD DONE':
            self.paramsWidget.progressBar.setValue(0)
            self.jobDone.emit()
        else:
            #self.paramsWidget.computeGridsButton.setText(text)
            self.statusbar.showMessage(text, level='busy')
            
    def updateProgress(self, value):
        self.paramsWidget.progressBar.setValue(value)

    def configure(self, options):
        """Update the gui with the parameters from options dictionary
        DEFAULT OPTIONS: {'receptorFile': None,  'ligandFile': None, 'receptorMaps': None,
                   'flexres': None, 'boxMode': None,
                   'padding': 4.0,  'smooth': 0.5, 'spacing': 0.375,
                   'waterWeight': 0.6, 'waterEntropy': -0.2,
                   'receptorGradient': True, 'recGradVolCut': None,  
                   'mapTypes': 'all',   
                   'autoSiteVersion': 1.0, 'pepScore': False, 'ligandSize': 500}
        """
        self.agfr.proccessPreLoadCmdLineOptions(options)

        self._cmdlineOption = options

        #pyShell = options.get('pyShell', None)
        #if pyShell:
        self.createPyShellWidget()

        trgFile = options.get('receptorMaps', None)
        receptorFile = options.get('receptorFile', None)
        if trgFile:
            self.getReceptorMaps(trgFile)
        elif receptorFile:
            options.pop("receptorFile")
            self._gridGuiOptions = options
            self.loadReceptor(receptorFile)
        else:
            self.setGridGuiOptions(options)

        if self.agfr._ligFilename:
            self.loadLigand(self.agfr._ligFilename)

            #the trg file contains receptor, ligand and all the parameters data

    def setGridGuiOptions(self, options):
        flexResStr = options.get('flexres', None)
        if flexResStr is not None and  self.receptor is not None:
            self._noFlexResParse = False
            self.paramsWidget.flexResWidget.setText(flexResStr)
            self.setFlexRes()
            self._noFlexResParse = True
            
        padding = options.get('padding', self.agfr.padding)
        spacing = options.get('spacing', self.agfr.spacing)
        smooth =  options.get('smooth', 0.5)
        if padding != self.paramsWidget.gridPaddingWidget.value():
            self.paramsWidget.gridPaddingWidget.setValue(padding)
        self._spacing = spacing
        self._smooth = smooth
        boxMode = options.get('boxMode', None)
        """
          -b boxmode [modes ...], --boxMode boxmode [modes ...]
                        docking box definition mode. Mode can be:
                         receptor: smallest box encompassing the entire receptor (default).
                         ligand  : smallest box encompassing a specified ligand (see -l/--ligand)
                         fill    : smallest box encompassing fill points (see fills)
                         residues chid resnamesResnums: smallest box encompassing the specified residues
                         user cx cy cz sx sy sz: box centered at (cx, cy, cz) of size sx, sy, sz (units Angstroms)
                         user centerMode sx sy sz:  centerMode can be: receptor, ligand, fill, residues
        """
        
        if boxMode is not None:
            if boxMode[0] == "ligand" and ligandFile:
                self.setGridFullLigand()
            elif boxMode[0] == "receptor" and self.receptor:
                self.setGridFullReceptor()
            elif boxMode[0] == "residues" and self.receptor:
                if len (boxMode) == 1:
                    self.agfr.myprint( "\nWARNING: can't set boxMode 'user %s', no residue selection string is specified" % boxMode[0])
                else:
                    try:
                        self.agfr.setBox(boxMode, padding, self._spacing)
                        self._baseSize[:] = self.agfr.boxLengths - 2*padding
                        
                        self.boxGeom.setCenter( *self.agfr.boxCenter)
                        self.boxGeom.setSides( *self.agfr.boxLengths)
                        self.onBoxChange()
                    except:
                        self.agfr.myprint( "\nWARNING: can't set boxMode 'user %s', invalid residue selection string  %s " % (boxMode[0], boxMode[1]))
            elif boxMode[0] == "user":
                setmode = True
                mode1 = boxMode[1]
                if mode1 in ["receptor","ligand","fill","residues"]:
                    if mode1 in ["receptor", "residues"] and self.receptor is None:
                        self.agfr.myprint( "\nWARNING: can't set boxMode 'user %s', no receptor file is specified." % mode1)
                        setmode = False
                    elif mode1 == "ligand" and ligandFile is None:
                        self.agfr.myprint( "\nWARNING: can't set boxMode 'user ligand', no ligand file is specified.")
                        setmode = False
                    elif mode1 == "fill":
                        self.agfr.myprint( "\nWARNING: can't set boxMode 'user fill', no fill points have been computed yet")
                        setmode = False
                else:
                    try:
                        boxVals = [float(x) for x in boxMode[1:]]
                        assert len(boxVals) == 6
                    except:
                        self.agfr.myprint( "WARNING: Invalid boxmode %s.\nExpected receptor, ligand, fill, residues or 6 floats (cx,cy,cz,sx,sy,sz)\nfor box center and side lengths."% (boxMode,))
                        setmode = False
                if setmode:
                    try:
                        self.agfr.setBox(boxMode, padding, spacing)
                        center = self.agfr.boxCenter
                        sx, sy, sz = self.agfr.boxLengths
                        self.boxGeom.setCenter( *center)
                        self.boxGeom.setSides( sx, sy, sz)
                        self._baseSize[:] = (sx, sy, sz)
                        self.paramsWidget.gridPaddingWidget.setValue(self.agfr.padding)
                    except:
                        self.agfr.myprint( "WARNING: Failed to set boxmode 'user' %s" % (mode1,))
            elif boxMode[0] is not None:
                self.agfr.myprint( "WARNING: Failed to set boxmode %s -not supported for AGFR GUI"%boxMode[0])

        ### -ng
        recGrad = options.get('receptorGradient') #True default
        self.paramsWidget.addGradGroupBox.setChecked(recGrad)
        if recGrad:
            volCut = options.get('recGradVolCut', None)
            if volCut is not None:
                self.paramsWidget.gradLargeCl.setChecked(False)
                self.paramsWidget.useGradCutOff.setChecked(True)
                self.paramsWidget.gradCutOffBox.setValue(float(volCut))
                self.paramsWidget.gradCutOffBox.setDisabled(False)
            else:
                self.paramsWidget.gradLargeCl.setChecked(True)
                self.paramsWidget.useGradCutOff.setChecked(False)
                self.paramsWidget.gradCutOffBox.setDisabled(True)
                
        ### autosite version (1.0, 1.1)
        asversion = options.get('autoSiteVersion', None)
        if asversion is not None:
            if asversion == 1.1:
                self.paramsWidget.autoSiteVButton.rb2.setChecked(True)
                # -ls , --ligandSize
                # 'ligandSize': 500
                ligSize = int(options.get('ligandSize', 500))
                if ligSize != self._ligandSize:
                    self.paramsWidget.autoSiteVButton.sizeSB.setValue(ligSize)
                    self._ligandSize = ligSize                    
                # -ps, --pepScore
                #'pepScore': False
                pepSc = options.get('pepScore', True)
                if pepSc != self._pepScore:
                    self._pepScore = pepSc # True or False
                    self.paramsWidget.autoSiteVButton.pepCkeckB.setChecked(pepSc)

            else:
                self.paramsWidget.autoSiteVButton.rb2.setChecked(False)
        ## maptypes -m {all, ligand}
        maptypes = options.get('mapTypes', 'all')
        if maptypes == "ligand" and self.ligand:
            self.paramsWidget.compLigandAtypes.setChecked(True)
            self.paramsWidget.compAllButton.setChecked(False)
        else:
            self.paramsWidget.compLigandAtypes.setChecked(False)
            self.paramsWidget.compAllButton.setChecked(True)
        
        #--waterWeight
        #'waterWeight': 0.6,
        wW = options.get('waterWeight', None)
        if wW is not None and self.agfr._wMapWeight != wW:
            self.paramsWidget.wMapParamB.wWeight.setValue(wW)
        #--waterEntropy
        # 'waterEntropy': -0.2,
        wE = options.get('waterEntropy', None)
        if wE is not None and self.agfr._wMapEntropy != wE:
            self.paramsWidget.wMapParamB.wEntropy.setValue(wE)
        self._gridGuiOptions= {}



        
if __name__=='__main__':
    import sys

    from optparse import OptionParser
    parser = OptionParser(usage="usage: %prog ligand mapFolder [options] filename",
                      version="%prog 0.1")
    parser.add_option("-r", "--receptor",
                      action="store", # optional because action defaults to "store"
                      dest="receptorFile",
                      help="receptor PDBQT file",)
    parser.add_option("-l", "--ligand",
                      action="store", # optional because action defaults to "store"
                      dest="ligandFile",
                      help="ligand PDBQT file",)
    parser.add_option("-t", "--trg",
                      action="store", # optional because action defaults to "store"
                      dest="receptorMaps",
                      help=".trg file containing receptor maps and receptor structure",)

    (options, args) = parser.parse_args()
    app = QtGui.QApplication(sys.argv)
    widget = GridGUI()
    if len(sys.argv)==2:
        widget.pmv.readMolecules(sys.argv[1])
    #widget.resize(1000,400)
    widget.show()
    if options.receptorFile:
        widget.loadReceptor(options.receptorFile)
    if options.ligandFile:
        widget.loadLigand(options.ligandFile)
    if options.receptorMaps:
        widget.getReceptorMaps(options.receptorMaps)
    
    #widget.viewer.OneRedraw()
    widget.raise_()
    timer = QtCore.QTimer()
    timer.singleShot(200, widget.viewer.OneRedraw)
    timer.start(1)

    sys.exit(app.exec_())
