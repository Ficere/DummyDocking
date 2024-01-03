import os
import numpy as np

from ADFR import molToAtomSetStatic
from ADFR.GUI import ICONPATH
from ADFRcc.adfr import AtomSet, PairwiseScorer, Parameters

parameters = Parameters.getParameters()
feCoeffVdw = parameters.feCoeffVdw
feCoeffHbond = parameters.feCoeffHbond
feCoeffEstat = parameters.feCoeffEstat
feCoeffDesolv = parameters.feCoeffDesolv
feCoeffTors = parameters.feCoeffTors

from PySide import QtCore, QtGui
from ADFRcc.adfr import Parameters


class MyTableItem(QtGui.QTableWidgetItem):
    """Subclass QTableWidgetItem to allow numericla sortng of values"""

    def __lt__(self, other):
        try:
            return float(self.data(QtCore.Qt.EditRole)) < float(
                other.data(QtCore.Qt.EditRole)
            )
        except ValueError:
            # float will fail for cells below the precision (containing ".")
            return self.data(QtCore.Qt.EditRole) < other.data(QtCore.Qt.EditRole)


class PairwiseScoreTable(QtGui.QWidget):
    tableUpdated = QtCore.Signal()

    def __init__(self, precision=2, parent=None):
        super(PairwiseScoreTable, self).__init__(parent)
        self.scorer = None
        self.setPrecision(precision)
        self.buildUI()

    def setScorer(self, scorer, atoms1, atoms2):
        self.scorer = scorer
        self._distance = scorer.getDistanceMatrix()
        self._vdwArray = scorer.getVdwEnergyMatrix() * feCoeffVdw
        self._eArray = scorer.getEstatEnergyMatrix() * feCoeffEstat
        self._hArray = scorer.getHbEnergyMatrix() * feCoeffHbond
        self._dsArray = scorer.getSolvEnergyMatrix() * feCoeffDesolv
        self._atomSet1 = scorer.getAtomSet1().atomSetStatic  ## ADFRcc AtomSetStatic
        self._atomSet2 = scorer.getAtomSet2().atomSetStatic  ##
        self._atoms1 = atoms1  ## prody AtomSet
        self._atoms2 = atoms2
        vdw = np.sum(self._vdwArray)
        elec = np.sum(self._eArray)
        hb = np.sum(self._hArray)
        ds = np.sum(self._dsArray)
        self.setEnergyLabel(vdw, elec, hb, ds)
        self.updateTable()

    def setPrecision(self, precision):
        self._precision = precision
        self._format = "%." + str(precision) + "f"
        self._toosmall = "0."
        for n in range(precision):
            self._toosmall += "0"
        # print 'FOGO', self._precision, self._format, self._toosmall
        self.updateTable()

    def setEnergyLabel(self, vdw, elec, hb, ds):
        self.energyLabel.setText(
            "ene.: %.2f vdw: %.2f el: %.2f hb: %.2f ds: %.2f "
            % (vdw + elec + hb + ds, vdw, elec, hb, ds)
        )

    def buildUI(self):
        layout = QtGui.QVBoxLayout()
        hlayout = QtGui.QHBoxLayout()

        self.distanceWidget = w = QtGui.QDoubleSpinBox()
        w.setMaximum(99999.0)
        w.setMinimum(0.0)
        w.setValue(5.0)
        w.setMaximumSize(60, 20)
        w.valueChanged.connect(self.updateTable)
        label = QtGui.QLabel("dist. cut.:")
        hlayout.addWidget(label)
        hlayout.addWidget(self.distanceWidget)

        self.precisionWidget = w = QtGui.QSpinBox()
        label = QtGui.QLabel("precision:")
        w.setMaximum(8)
        w.setMinimum(1)
        w.setValue(2)
        w.setMaximumSize(40, 20)
        w.valueChanged.connect(self.setPrecision)
        hlayout.addWidget(label)
        hlayout.addWidget(self.precisionWidget)

        ## from mglutil.util.packageFilePath import findFilePath
        ## PMVICONPATH = findFilePath('Icons', 'ADFR.GUI')
        ## w = QtGui.QPushButton(QtGui.QIcon(os.path.join(PMVICONPATH, 'table-icon-png-2.png')), '')
        ## w.setIconSize(QtCore.QSize(24, 24))
        ## w.setCheckable(True)
        ## w.setChecked(False)

        ## w = QtGui.QComboBox(self)
        ## w.addItems(['per atom', 'per ligand residue', 'per receptor residue'])
        ## self.tableModeWidget = w

        ## hlayout.addWidget(w)
        # horizontalSpacer = QtGui.QSpacerItem(1, 1, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        # hlayout.insertSpacerItem(4,horizontalSpacer)
        hlayout.addStretch()

        self.energyLabel = QtGui.QLabel("ene.:")
        hlayout.addWidget(self.energyLabel)

        # hlayout.addStretch()

        # self.focusButton = button = QtGui.QPushButton(self, checkable=True)
        # button.setIcon(QtGui.QIcon(os.path.join(ICONPATH, 'focus.png')))
        # button.setStyleSheet("QToolButton { border: 2px solid #8f8f91; border-radius: 6px;}")
        # button.setChecked(True)
        # hlayout.addWidget(button)

        # build table of pairwise interactions
        self.table = QtGui.QTableWidget(1, 9)
        self.table.verticalHeader().hide()
        self.table.horizontalHeader().setResizeMode(QtGui.QHeaderView.Stretch)
        self.columnLabels = [
            "num",
            "at1",
            "at2",
            "total",
            "vdw",
            "elec",
            "hb",
            "ds",
            "dist",
        ]
        self.table.setHorizontalHeaderLabels(self.columnLabels)
        # layout.addLayout(formLay)
        layout.addLayout(hlayout)
        layout.addWidget(self.table)
        self.setLayout(layout)

    def getPairs(self):
        # return list of i,j pairs based on current cutoff values
        dcut = self.distanceWidget.value()
        pairs = np.nonzero(self._distance < dcut)
        return zip(pairs[0], pairs[1])

    def getValueString(self, value):
        valstr = self._format % value
        if (
            valstr == self._toosmall
            or valstr == "+" + self._toosmall
            or valstr == "-" + self._toosmall
        ):
            valstr = "."
        return valstr

    def updateTable(self):
        if self.scorer is None:
            return
        # from time import time
        # t0 = time()

        # forget existing rows
        self.table.setRowCount(0)

        pairs = self.getPairs()
        if len(pairs) == 0:
            return

        self.table.setRowCount(len(pairs))
        names1 = self._atoms1.getNames()
        names2 = self._atoms2.getNames()
        resnames1 = self._atoms1.getResnames()
        resnames2 = self._atoms2.getResnames()
        resnums1 = self._atoms1.getResnums()
        resnums2 = self._atoms2.getResnums()
        from DejaVu2.colorTool import GreenWhiteRed, Map

        ramp = GreenWhiteRed(upper=1.0)

        maxi = max(
            max(self._vdwArray.flat),
            max(self._eArray.flat),
            max(self._hArray.flat),
            max(self._dsArray.flat),
        )
        mini = min(
            min(self._vdwArray.flat),
            min(self._eArray.flat),
            min(self._hArray.flat),
            min(self._dsArray.flat),
        )
        limit = max(abs(maxi), abs(mini))
        vdw = [self._vdwArray[x] for x in pairs]
        elec = [self._eArray[x] for x in pairs]
        hb = [self._hArray[x] for x in pairs]
        ds = [self._dsArray[x] for x in pairs]
        tot = vdw + elec + hb + ds

        vdwc = Map(vdw, ramp, -limit, limit) * 255
        elecc = Map(elec, ramp, -limit, limit) * 255
        hbc = Map(hb, ramp, -limit, limit) * 255
        dsc = Map(ds, ramp, -limit, limit) * 255
        totc = Map(tot, ramp, -limit, limit) * 255

        # self.table.setSortingEnabled(False)
        n = 0  # index in pairs
        k = 0  # index in table
        keptPairs = []
        for i, j in pairs:
            strs = [
                self.getValueString(tot[n]),
                self.getValueString(vdw[n]),
                self.getValueString(elec[n]),
                self.getValueString(hb[n]),
                self.getValueString(ds[n]),
            ]
            if set(strs) == set(["."]):
                n += 1
                continue  # all entries are below precision (i.e. '.'), we skip
            # print k, i, j, strs
            self.table.setItem(
                k, 0, MyTableItem(str(k))
            )  # hidden colum to find pairs after sort
            self.table.setItem(
                k,
                1,
                QtGui.QTableWidgetItem(
                    "%s`%s:%s" % (resnames1[i], resnums1[i], names1[i])
                ),
            )
            self.table.setItem(
                k,
                2,
                QtGui.QTableWidgetItem(
                    "%s`%s:%s" % (resnames2[j], resnums2[j], names2[j])
                ),
            )

            w = MyTableItem(strs[0])
            w.setBackground(QtGui.QColor(*totc[n]))
            self.table.setItem(k, 3, w)
            if vdw[n]:
                w = MyTableItem(strs[1])
                w.setBackground(QtGui.QColor(*vdwc[n]))
                self.table.setItem(k, 4, w)
            if elec:
                w = MyTableItem(strs[2])
                w.setBackground(QtGui.QColor(*elecc[n]))
                self.table.setItem(k, 5, w)
            if hb:
                w = MyTableItem(strs[3])
                w.setBackground(QtGui.QColor(*hbc[n]))
                self.table.setItem(k, 6, w)
            if ds:
                w = MyTableItem(strs[4])
                w.setBackground(QtGui.QColor(*dsc[n]))
                self.table.setItem(k, 7, w)
            self.table.setItem(k, 8, MyTableItem("%.2f" % self._distance[i, j]))
            keptPairs.append((i, j))
            n += 1
            k += 1
            # print n, i, j, names1[i], names2[j], vdw+elec+hb+ds

        self._pairs = keptPairs
        self.table.setRowCount(len(keptPairs))
        self.table.resizeColumnsToContents()
        self.table.setColumnHidden(0, True)
        self.table.sortItems(3, QtCore.Qt.AscendingOrder)
        self.table.setSortingEnabled(True)
        self.tableUpdated.emit()
        # print 'UPDATE TIME', time()-t0


if __name__ == "__main__":
    import sys
    from MolKit2 import Read

    rec = Read("4ek4_rec.pdbqt")
    lig = Read("4ek4_lig.pdbqt")

    recAS = AtomSet(molToAtomSetStatic(rec))
    ligAS = AtomSet(molToAtomSetStatic(lig))

    scorer = PairwiseScorer()
    scorer.initialize(recAS, ligAS)
    score = scorer.calculateScores()
    scorer.printDebugDescription()

    app = QtGui.QApplication(sys.argv)

    w = PairwiseScoreTable(scorer, rec._ag, lig._ag)
    w.resize(400, 800)
    w.show()
    sys.exit(app.exec_())
