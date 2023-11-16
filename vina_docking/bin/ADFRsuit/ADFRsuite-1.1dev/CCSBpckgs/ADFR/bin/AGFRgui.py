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
#
# Author: Michel Sanner
#
# $Header: /mnt/raid/services/cvs/ADFR/bin/AGFRgui.py,v 1.3.2.1 2017/08/03 01:34:11 mgltools Exp $
#
# $Id: AGFRgui.py,v 1.3.2.1 2017/08/03 01:34:11 mgltools Exp $
#

import sys, os

from ADFR.GUI import makeGridGUI
from ADFR.GUI.gridGUI import GridGUI

from PySide import QtCore, QtGui

# use 4.2 by default
#from ADFRcc import setForceFieldVersion
#parameters = setForceFieldVersion('4.2')

sys.setrecursionlimit(10000)
app = QtGui.QApplication(sys.argv)

widget = makeGridGUI(app)

#widget.show()
widget.raise_()
timer = QtCore.QTimer()
timer.singleShot(200, widget.viewer.OneRedraw)
timer.start(1)

sys.exit(app.exec_())
