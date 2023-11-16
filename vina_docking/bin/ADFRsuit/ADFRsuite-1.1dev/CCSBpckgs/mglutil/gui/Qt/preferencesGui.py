from PySide import QtCore, QtGui
from mglutil.util.callback import CallbackFunction

class PreferencesTreeWidget(QtGui.QWidget):
    def __init__(self, preferences, app, parent=None):
        super(PreferencesTreeWidget, self).__init__(parent)
        self.preferences = preferences
        layout =  QtGui.QVBoxLayout()
        self.treeWidget = treeWidget = QtGui.QTreeWidget()
        layout.addWidget(self.treeWidget)
        self.setLayout(layout)
        
        treeWidget.setColumnCount(2)
        treeWidget.headerItem().setText(0, "")
        #treeWidget.headerItem().setText(1, "current value")
        treeWidget.headerItem().setText(1, "set value")
        treeWidget.setColumnWidth(0, 300)
        treeWidget.setColumnWidth(1, 200)
        #treeWidget.setColumnWidth(2, 30)
        treeWidget.itemClicked.connect(self.onClick)
        categories = {"General": QtGui.QTreeWidgetItem(treeWidget.invisibleRootItem(), ["General"])}
        for name, prefs in preferences.items():
            category = prefs['category']
            if categories.has_key(category):
                catItem = categories[category]
            else:
                catItem = QtGui.QTreeWidgetItem(treeWidget.invisibleRootItem(),category)
            catItem.setExpanded(True)
            nameItem = QtGui.QTreeWidgetItem(catItem)
            nameItem.setText(0, name)
            nameItem.setExpanded(True)
            curVal = prefs['value']
            #print "current val for name:", name, curVal
            #nameItem.label = QtGui.QLabel(str(curVal))
            validVals = prefs['validValues']
            if len(validVals):
                nameItem.comboBox = QtGui.QComboBox()
                nameItem.comboBox.addItems([str(x) for x in validVals])
                self.treeWidget.setItemWidget(nameItem, 1, nameItem.comboBox)
                cb = CallbackFunction(self.setPrefValue, nameItem, name)
                nameItem.comboBox.currentIndexChanged.connect(cb)
                ind = validVals.index(curVal)
                nameItem.comboBox.setCurrentIndex(ind)
        self.text = QtGui.QPlainTextEdit()
        self.text.setReadOnly(True)
        self.text.setFixedHeight(100)
        self.text.setLineWrapMode(self.text.WidgetWidth)
        treeWidget.setMinimumHeight(150)
        layout.addWidget(self.text)
        
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        self.setSizePolicy(sizePolicy)
        header = self.treeWidget.header()
        header.setResizeMode(QtGui.QHeaderView.ResizeToContents)
        header.setStretchLastSection(False)

        self.saveB = QtGui.QPushButton("save preferenced")
        self.saveB.clicked.connect(self.savePreferences_cb)
        layout.addWidget(self.saveB)
        
    def getTreeWidth(self):
        header = self.treeWidget.header()
        w = self.treeWidget.frameWidth() * 2 + self.treeWidget.style().pixelMetric(QtGui.QStyle.PM_ScrollBarExtent)
        for i in range(self.treeWidget.columnCount()):
            w += header.sectionSize(i)
        return w

    def onClick(self, item, column):
        print "item clicked:", item.text(0), column
        pref = self.preferences.get(str(item.text(0)), None)
        if pref:
            self.text.setPlainText(pref['doc'])

    def setPrefValue(self, nameItem, name, index):
        self.treeWidget.setCurrentItem(nameItem)
        val = str(nameItem.comboBox.currentText())
        #val = self.preferences[str(nameItem.text(0))]['validValues'][index]
        self.preferences.set(name, val)
        #print "%s: %s"% (name,val)
        cb = self.preferences[name].get('callbackFunc')
        if cb: cb()

    def savePreferences_cb(self):
        self.preferences.saveAllSettings()

    
class PreferencesTreeDialog(QtGui.QDialog):
    closedSignal = QtCore.Signal()
    
    def __init__(self, preferences, app, title="", parent=None):
        
        super(PreferencesTreeDialog, self).__init__(parent)
        self.setWindowTitle(title)

        # place the widget close to where we clicked
        rel_pos = QtGui.QCursor.pos()
        pos = self.mapToGlobal(rel_pos)
        self.move(pos.x()+20, pos.y()+15)

        tree = self.prefTreeWidget = PreferencesTreeWidget(preferences, app , parent)
        self.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
        self.buttonBox = QtGui.QDialogButtonBox()
        self.buttonBox.accepted.connect(self.accept)

        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.prefTreeWidget)
        layout.addWidget(self.buttonBox)
        self.setLayout(layout)
        margins = self.layout().contentsMargins()
        w = tree.getTreeWidth()
        self.resize(margins.left() + margins.right()+ w,  self.minimumHeight()) # so that all columns of the tree are visible in the dialog.
            
    def reject(self):
        QtGui.QDialog.reject(self)
        self.closeEvent()
        
    def closeEvent(self, evnt=None):
        self.closedSignal.emit()
