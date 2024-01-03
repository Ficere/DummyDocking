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

########################################################################
#
# Date: 2014 Authors: Michel Sanner
#
#    sanner@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Michel Sanner and TSRI
#
#########################################################################

from PySide import QtGui, QtCore
import os, sys
from mglutil.util.packageFilePath import findFilePath

PMVICONPATH = findFilePath("Icons", "PmvApp.GUI")
import urllib2


class FetchGUI(QtGui.QDialog):
    """This class provides GUI for fetching a molecule from PDB"""

    def __init__(self, parent=None, formats=["pdb", "mmtf"]):
        super(FetchGUI, self).__init__(parent)

        self.pdbidEditW = QtGui.QLineEdit("")
        self.pdbidEditW.returnPressed.connect(self.accept)

        hLayout = QtGui.QHBoxLayout()
        # self.kwEditW = QtGui.QLineEdit("")
        # self.kwEditW.returnPressed.connect(self.showTable)
        self.searchButton = QtGui.QPushButton()
        self.searchButton.setIcon(
            QtGui.QIcon(os.path.join(PMVICONPATH, "search-icon.png"))
        )

        self.searchButton.clicked.connect(self.showTable)
        self.searchButton.setToolTip("Open a keyword search widget")
        self.searchButton.setCheckable(True)
        self.searchButton.setChecked(False)
        label0 = QtGui.QLabel("Pdb ID:")
        hLayout.addWidget(label0)
        hLayout.addWidget(self.pdbidEditW)
        self.pdbidEditW.setSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        self.searchButton.setSizePolicy(
            QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed
        )
        # label1 = QtGui.QLabel("keyword:")
        # hLayout.addWidget(label1)
        hLayout.addWidget(self.searchButton)
        spacer = QtGui.QWidget()
        spacer.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        hLayout.addWidget(spacer)
        self.layout_ = layout = QtGui.QFormLayout()
        layout.setFieldGrowthPolicy(
            QtGui.QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow
        )  # (2|1))
        layout.addRow(hLayout)
        self.formats = formats
        self.formatBox = None
        if len(formats) > 1:
            self.formatBox = QtGui.QComboBox()
            self.formatBox.addItems(self.formats)
            layout.addRow("Format:", self.formatBox)
        self.buttons = QtGui.QDialogButtonBox(QtCore.Qt.Horizontal, self)
        fetchButton = QtGui.QPushButton(self.tr("&Fetch"))
        self.buttons.addButton(fetchButton, QtGui.QDialogButtonBox.AcceptRole)
        self.buttons.addButton(QtGui.QDialogButtonBox.Close)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        fetchButton.setDefault(False)
        layout.addRow("", self.buttons)
        self.setLayout(layout)
        self.pdbidEditW.setFocus()
        self.tableW = None

    def accept(self):
        pdbs = self.getSelectedPdbs()
        if not len(pdbs):
            flags = QtGui.QMessageBox.StandardButton.Ok
            msg = "Pdb Id is not specified"
            response = QtGui.QMessageBox.warning(self, "Warning!", msg, flags)
        else:
            pdbId = " ".join(pdbs)
            self.pdbidEditW.setText(pdbId)
            QtGui.QDialog.accept(self)

    def showTable(self):
        # opens a second dialog with an entry for keyword search and a
        # QTableWidget that will contain the results of the search.
        if self.searchButton.isChecked():
            if not self.tableW:
                self.tableW = tw = TableWidget("Search PDB", parent=self, keyword=None)
                nrows = self.layout_.rowCount()
                self.layout_.insertRow(nrows - 1, tw)
                tw.searchEnd.connect(self.setSize)
            else:
                self.tableW.show()
        else:
            self.tableW.hide()
            self.adjustSize()

    def setSize(self):
        self.adjustSize()

    def getSelectedPdbs(self):
        pdbIds = []
        if self.tableW:
            pdbIds.extend(self.tableW.getCheckedItems())
        id = self.pdbidEditW.text()
        if len(id) and id not in pdbIds:
            pdbIds.insert(0, id)

        return pdbIds

    def getSelected(self):
        selected = self.tableW.table.selectedIndexes()
        # update the entry with selected PdbIds
        pdbId = ""
        for item in selected:
            column = item.column()
            if column == 0:
                pdbId += str(item.data()) + " "
        self.pdbidEditW.setText(pdbId)
        return pdbId

    @staticmethod
    def getName(parent=None, formats=["pdb", "mmtf"]):
        dialog = FetchGUI(parent=parent, formats=formats)
        dialog.setWindowTitle(
            QtGui.QApplication.translate(
                "Dialog", "Fetch molecule", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        # dialog.setWindowModality(QtCore.Qt.WindowModal)
        result = dialog.exec_()
        # print "FetchGUI result:", result
        molname = dialog.pdbidEditW.text()
        if dialog.formatBox:
            molformat = dialog.formatBox.currentText()
        else:
            molformat = dialog.formats[0]
        return ([str(molname), str(molformat)], result == QtGui.QDialog.Accepted)

    def keyPressEvent(self, evt):
        # this is to prevent closing of the dialog when "Return" is pressed in LineEdit widget
        if evt.key() == QtCore.Qt.Key_Enter or evt.key() == QtCore.Qt.Key_Return:
            return
        QtGui.QDialog.keyPressEvent(self, evt)


import re


class TableWidget(QtGui.QWidget):
    searchResult = QtCore.Signal(str)
    searchEnd = QtCore.Signal()

    uniprotACIDpattern = re.compile(
        "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
    )

    url = "http://www.rcsb.org/pdb/rest/search"
    txtSearch = """<orgPdbQuery>
    
<queryType>org.pdb.query.simple.AdvancedKeywordQuery</queryType>

<description>Text Search for: %s</description>

<keywords>%s</keywords>

</orgPdbQuery>"""

    uniprotACSearch = """<orgPdbQuery>

<queryType>org.pdb.query.simple.UpAccessionIdQuery</queryType>

<description>Simple query for a list of Uniprot Accession IDs: %s</description>

<accessionIdList>%s</accessionIdList>

</orgPdbQuery>"""

    def __init__(self, title, keyword=None, parent=None):
        super(TableWidget, self).__init__(parent)
        self.dataList = []
        self.header = header = ["ID", "#Res", "Method", "Reso", "Date", "Title"]
        ncolumns = len(header)

        # self.setMinimumSize(700, 300)
        self.setMinimumWidth(700)
        self.setWindowTitle(title)

        hLayout = QtGui.QHBoxLayout()
        self.keywordEditW = QtGui.QLineEdit()
        if keyword is not None:
            self.keywordEditW.setText(keyword)
        # self.keywordEditW.setMaximumWidth(150)
        searchLabel = QtGui.QLabel("Keyword:")
        self.searchButton = QtGui.QPushButton("Search")
        self.searchButton.clicked.connect(self.fillTableDialog)
        self.keywordEditW.returnPressed.connect(self.fillTableDialog)
        hLayout.addWidget(searchLabel)
        hLayout.addWidget(self.keywordEditW)
        hLayout.addWidget(self.searchButton)

        self.table = table = QtGui.QTableWidget(1, ncolumns, self)
        table.setHorizontalHeaderLabels(header)
        table.horizontalHeader().setResizeMode(ncolumns - 1, QtGui.QHeaderView.Stretch)
        # table.resizeColumnsToContents()
        table.setColumnWidth(0, 50)
        table.setColumnWidth(1, 45)
        table.setColumnWidth(2, 80)
        table.setColumnWidth(3, 50)
        table.setColumnWidth(4, 80)
        table.cellClicked.connect(self.cellClicked_cb)
        lfont = QtGui.QFont()
        lfont.setPointSize(8)
        table.setFont(lfont)
        layout = QtGui.QVBoxLayout()
        layout.addLayout(hLayout)
        layout.addWidget(table)

        # enable sorting
        table.horizontalHeader().setSortIndicator(0, QtCore.Qt.AscendingOrder)
        # setSortingEnabled() calls the sort() of the table_model in the default Discending order.
        # setSortIndicator() is used to change the initial sort order to Ascending.
        # table.setSortingEnabled(True)

        # add progress bar , "close", "accept" buttons
        buttonsLayout = QtGui.QHBoxLayout()
        # self.progress = QtGui.QProgressBar(self)
        self.numIdFoundLabel = QtGui.QLabel("")

        self.toggleSelectAllB = QtGui.QPushButton("Check All")
        self.toggleSelectAllB.clicked.connect(self.toggleSelectAll)
        self.toggleSelectAllB.setSizePolicy(
            QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed
        )
        self.toggleSelectAllB.setEnabled(False)
        buttonsLayout.addWidget(self.toggleSelectAllB)
        # buttonsLayout.addWidget(self.progress)
        buttonsLayout.addWidget(self.numIdFoundLabel)
        layout.addLayout(buttonsLayout)
        self.setLayout(layout)
        self.searchResult.connect(self.fillTable)
        self.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        policy = self.sizePolicy()
        policy.setVerticalStretch(1)
        self.setSizePolicy(policy)

    def toggleSelectAll(self):
        # check or uncheck all items in the first column of the table
        text = str(self.toggleSelectAllB.text())
        if text == "Check All":
            checked = QtCore.Qt.Checked
            newtext = "Uncheck All"
        else:
            checked = QtCore.Qt.Unchecked
            newtext = "Check All"
        for i in range(self.table.rowCount()):
            item = self.table.item(i, 0)
            if item:
                item.setCheckState(checked)
        self.toggleSelectAllB.setText(newtext)

    def fillTableDialog(self):
        # display search results in the table
        keyword = str(self.keywordEditW.text())
        self.table.clearContents()
        if self.table.rowCount() > 1:
            self.table.setRowCount(1)
        if len(keyword):
            import thread

            thread.start_new_thread(self.keywordSearch, (keyword,))

    # @busy_cursor
    def keywordSearch(self, keyword):
        self.searchButton.setText("Searching....")
        self.searchButton.setDisabled(True)
        match = self.uniprotACIDpattern.match(keyword.strip())
        if match is None:
            req = urllib2.Request(self.url, data=self.txtSearch % (keyword, keyword))
        else:
            req = urllib2.Request(
                self.url, data=self.uniprotACSearch % (keyword, keyword)
            )

        f = urllib2.urlopen(req)
        self._result = f.read()
        # print 'time search', time()-t0
        pdbIds = self._result.replace("\n", ",")[:-1]
        url = (
            "http://www.rcsb.org/pdb/rest/customReport.xml?pdbids=%s&customReportColumns=structureId,residueCount,experimentalTechnique,resolution,releaseDate,structureTitle&service=wsfile&format=csv "
            % pdbIds
        )
        try:
            response = urllib2.urlopen(url)
            result = response.read()
            # self.searchResult.emit( result.replace('"','') ) # this calls self.fillTable
            self.searchResult.emit(result)  # this calls self.fillTable
        except urllib2.HTTPError:
            self.searchResult.emit('No match for: "%s"' % keyword)
        self.searchButton.setDisabled(False)
        self.searchButton.setText("Search")

    def fillTable(self, result):
        if result.startswith("No match for:"):
            return
        lines = result.split("\n")[1:-1]
        nlines = len(lines)
        self.numIdFoundLabel.setText("Found %d entities" % nlines)
        for n, line in enumerate(lines):
            ncolumns = len(self.header)
            for j, txt in enumerate(line.split('",')):
                newitem = QtGui.QTableWidgetItem(txt.replace('"', ""))
                if j == 0:
                    newitem.setCheckState(QtCore.Qt.Unchecked)
                self.table.setItem(self.table.rowCount() - 1, j, newitem)
                if j == ncolumns - 1:
                    font = QtGui.QFont()
                    font.setUnderline(True)
                    # font.setPointSize(8)
                    newitem.setFont(font)
                    newitem.setForeground(QtGui.QBrush(QtGui.QColor("blue")))
            nrows = self.table.rowCount()
            # if n < nlines-1:
            if nrows < nlines:
                self.table.insertRow(nrows)
            if nrows <= 15:  # ????
                # print "getting Size, n=", n
                h = self.getTableHeight()
                self.table.setFixedHeight(h)
            # self.progress.setValue(int(100*n/nlines))
            # self.progress.update()
            # print "Progress val:", 100*n/nlines
        # self.progress.setValue(100)
        if not self.table.isSortingEnabled():
            self.table.setSortingEnabled(True)
        self.toggleSelectAllB.setEnabled(True)
        self.adjustSize()
        self.searchEnd.emit()

    def getTableWidth(self):
        w = self.table.verticalHeader().width() + 4  # +4 seems to be needed
        for i in range(self.table.columnCount()):
            w += self.table.columnWidth(i)  # seems to include gridline
        return w

    def getTableHeight(self):
        nrows = self.table.rowCount()
        h = self.table.horizontalHeader().height()  # + 4
        for i in range(nrows):
            h += self.table.rowHeight(i)
        return h

    def cellClicked_cb(self, row, column):
        # enable "Use selected PDB IDs " button if any pdbid is selected, disable otherwise
        selecdedPdb = 0
        selected = self.table.selectedIndexes()
        for item in selected:
            if item.column() == 0:
                if item.data() is not None:
                    selecdedPdb += 1
        # open link in the browser if a cell in the last column is clicked
        if column == len(self.header) - 1:
            import webbrowser

            pdbid = str(self.table.item(row, 0).text())
            url = "http://www.rcsb.org/pdb/explore/explore.do?structureId=%s" % pdbid
            webbrowser.open(url, new=0, autoraise=True)

    def getCheckedItems(self):
        # return a list of pdbId in the table that are checked
        pdbIds = []
        for i in range(self.table.rowCount()):
            item = self.table.item(i, 0)
            if item and item.checkState() == QtCore.Qt.Checked:
                pdbIds.append(str(item.data(0)))
        return pdbIds
