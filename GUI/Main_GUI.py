
from PyQt5.QtWidgets import QApplication, QMainWindow, QTabWidget, QWidget
import sys
from GUI_LibraryReformat import LibraryReformatterGUI
from GUI_identification import Identification_GUI
from GUI_Quantification import Meta_Quantification_GUI 

# Assuming Identification_GUI and LibraryReformatterGUI classes are defined

class TableWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("DI_metabolome")
        self.setGeometry(150, 150, 850, 750)
        
        # Create the tab widget
        self.tabs = QTabWidget()
        self.setCentralWidget(self.tabs)
        
        # Add tabs
        
        self.Libtab = LibraryReformatterGUI()
        self.IDtab = Identification_GUI()         
        self.Quantab = Meta_Quantification_GUI()
        
        
        self.tabs.addTab(self.Libtab, "Library Reformatter")
        self.tabs.addTab(self.IDtab, "Metabolite Identification")
        self.tabs.addTab(self.Quantab, "Metabolite Quantification")
        
        
def main():
    app = QApplication(sys.argv)
    mainWindow = TableWindow()
    mainWindow.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()

