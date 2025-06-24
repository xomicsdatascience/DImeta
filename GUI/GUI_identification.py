
import sys
import os
import uuid
from PyQt5.QtWidgets import (QApplication, QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QDialog, QSpinBox, QPushButton, QCheckBox,
                             QLabel, QTextEdit, QListWidget, QLineEdit, QFileDialog, QMessageBox,QMainWindow, QTabWidget,QFormLayout,
                             QInputDialog)
from PyQt5.QtCore import Qt
from LibraryHandling import LibraryLoadingStrategy, LibraryReformat, LibrarySaveStrategy, read_path
from IdentificationMeta import QueryTargetedSpectrum, match_spectrum, group_tuples_by_same_value, filter_tuples, cosine_similarity, custom_sort, macc_score, normalize_to_100

from querylibrarymatch import get_spectra, match_and_calculate_cosine_similarity, main_processing_function, generate_plot, save_results
import pandas as pd

import os
import glob
import numpy as np
import heapq
import logging
from collections import defaultdict

class QTextEditLogger(logging.Handler, object):
    def __init__(self, widget):
        super(QTextEditLogger, self).__init__()
        self.widget = widget
        self.widget.setReadOnly(True)
    
    def emit(self, record):
        msg = self.format(record)
        self.widget.append(msg)

class StreamToLogger(object):
    """
    Fake file-like stream object that redirects writes to a logger instance.
    """
    def __init__(self, logger, log_level=logging.INFO):
        self.logger = logger
        self.log_level = log_level
        self.linebuf = ''
    
    def write(self, buf):
        for line in buf.rstrip().splitlines():
            self.logger.log(self.log_level, line.rstrip())
    
    def flush(self):
        pass



class Identification_GUI(QWidget):
    def __init__(self):
        super().__init__()
        #self.libraryFilePath = ""
        #self.figPath = ""
        #self.InputFilePaths = []
        self.initUI()
        
    def initUI(self):
        self.setWindowTitle('Metabolite Identification')
        layout = QVBoxLayout()
        self.formLayout = QFormLayout()
        
        # Create the buttons for file and folder selection
        self.libraryPathBtn = QPushButton('Browse')
        self.libraryPathBtn.clicked.connect(self.selectLibraryFile)
        self.figPathBtn = QPushButton('Browse')
        self.figPathBtn.clicked.connect(self.selectFigFolder)
        self.mzmlFilesBtn = QPushButton('Browse')
        self.mzmlFilesBtn.clicked.connect(self.selectInputFiles)
        # Create the Delete and Clear buttons
        self.deleteSelectedFilesBtn = QPushButton('Delete')
        self.deleteSelectedFilesBtn.clicked.connect(self.deleteSelectedFiles)
        self.clearBtn = QPushButton('Clear')
        self.clearBtn.clicked.connect(self.clearAllFiles)
        
        # Create labels to display the selected paths
        self.libraryPathLabel = QLabel('No file selected')
        self.figPathLabel = QLabel('No folder selected')
        self.mzmlFilesList = QListWidget()
       
        # Create input fields for analysis parameters
        self.ppmToleranceEdit = QLineEdit()
        self.minMatchedPeaksEdit = QLineEdit()
        self.precursorIonMassToleranceEdit = QLineEdit()
        self.intensityEdit = QLineEdit()  # New line edit for intensity
        self.intensityEdit.setText(f'{3e3}')  # Set default value to 1000 (10^3)
        self.cosineEdit = QLineEdit()  # New line edit for intensity
#         self.cosineEdit.setText(f'{3e3}')
        
        # self.scansEdit = QLineEdit()
        
        # Add widgets to the form layout
        self.formLayout.addRow(' Input Files (.mzML/.mzXML):', self.mzmlFilesBtn)
        self.formLayout.addRow(self.mzmlFilesList) 
        
        # Create a horizontal layout for Delete and Clear buttons
        self.buttonLayout = QHBoxLayout()
        self.buttonLayout.addWidget(self.deleteSelectedFilesBtn)
        self.buttonLayout.addWidget(self.clearBtn)
        self.formLayout.addRow(self.buttonLayout)
        
        self.formLayout.addRow('Library File:', self.libraryPathBtn)
        self.formLayout.addRow('', self.libraryPathLabel)
        
        self.formLayout.addRow('Output Files:', self.figPathBtn)
        self.formLayout.addRow('', self.figPathLabel)

        self.formLayout.addRow('Fragment Tolerance, ppm:', self.ppmToleranceEdit)
        self.formLayout.addRow('Min Matched Peaks:', self.minMatchedPeaksEdit)
        self.formLayout.addRow('Precursor Ion Mass Tolerance,da:', self.precursorIonMassToleranceEdit)
        self.formLayout.addRow('Intensity, 3e3:', self.intensityEdit)
        self.formLayout.addRow('Cosine_threshold,0.1~0.99:', self.cosineEdit)
        
        #self.formLayout.addRow('Number of Scans:', self.scansEdit)
        # Scans range input
        self.scansRangeLabel = QLabel("Scans Range:")
        self.scansLowerSpin = QSpinBox()
        self.scansUpperSpin = QSpinBox()

        # Configure the spin boxes for an example range of 1 to 10000
        self.scansLowerSpin.setMinimum(1)
        self.scansLowerSpin.setMaximum(10000)
        self.scansUpperSpin.setMinimum(1)
        self.scansUpperSpin.setMaximum(10000)
        self.scansLowerSpin.setValue(1)
        self.scansUpperSpin.setValue(1)
       
        # Create a layout for the range inputs
        self.scansRangeLayout = QHBoxLayout()
        self.scansRangeLayout.addWidget(self.scansLowerSpin)
        self.scansRangeLayout.addWidget(QLabel("to"))
        self.scansRangeLayout.addWidget(self.scansUpperSpin)
        
        # Add the range input layout to the form layout
        self.formLayout.addRow(self.scansRangeLabel, self.scansRangeLayout)
        
        # Create a checkbox for plot generation
        self.generatePlotsCheckbox = QCheckBox('Generate Identified Plots')
        self.generatePlotsCheckbox.setChecked(False)  # Default is unchecked (no plots generated)
        self.formLayout.addRow(self.generatePlotsCheckbox)
        
        layout.addLayout(self.formLayout)
        
        # Add the Clear button to the layout
#         layout.addWidget(self.clearBtn)
        
        # Button to start the analysis
        self.startBtn = QPushButton('Start Analysis',self)
        self.startBtn.clicked.connect(self.startAnalysis)
        layout.addWidget(self.startBtn)
        
        #  logging displaying area...
        self.logTextEdit = QTextEdit()  # Log display area
        layout.addWidget(self.logTextEdit)
        
        # Setup logging
        logTextBox = QTextEditLogger(self.logTextEdit)
        logTextBox.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logging.getLogger().addHandler(logTextBox)
        logging.getLogger().setLevel(logging.INFO)
        
        # Redirect stdout and stderr
        sys.stdout = StreamToLogger(logging.getLogger('STDOUT'), logging.INFO)
        sys.stderr = StreamToLogger(logging.getLogger('STDERR'), logging.ERROR)
        
        
        self.resize(700, 600)
        
        self.setLayout(layout)

    def selectLibraryFile(self):
        fileName, _ = QFileDialog.getOpenFileName(self, "Select Library File", "", "All Files (*);;Text Files (*.txt)")
        if fileName:
            self.libraryFilePath = fileName
            self.libraryPathLabel.setText(fileName)  # Update the label to show the selected path
            
    def selectFigFolder(self):
        folder = QFileDialog.getExistingDirectory(self, "Select Folder to Save Figures")
        if folder:
            self.figPath = folder
            self.figPathLabel.setText(folder)  # Update the label to show the selected folder
            
    def selectInputFiles(self):
        files, _ = QFileDialog.getOpenFileNames(self, "Select mzML Files", "", "mzML Files (*.mzML);;mzXML Files (*.mzXML);;All Files (*)")
        if files:
            self.InputFilePaths = files
            self.mzmlFilesList.clear()  # Clear existing items before adding new ones
            for file in files:
                self.mzmlFilesList.addItem(file)  # Add each selected file path to the list widget        
            # Calculate the number of scans from the last selected file
 
    def deleteSelectedFiles(self):
        selectedItems = self.mzmlFilesList.selectedItems()
        if not selectedItems:  # No item selected
            return
        for item in selectedItems:
            itemIndex = self.mzmlFilesList.row(item)
            self.mzmlFilesList.takeItem(itemIndex)
            self.InputFilePaths.pop(itemIndex)  # Ensure this list remains synchronized
            
    def clearAllFiles(self):
        """Function to clear all selected files and reset the interface."""
        self.mzmlFilesList.clear()
        self.libraryPathLabel.setText('No file selected')
        self.figPathLabel.setText('No folder selected')
        self.InputFilePaths = []
        self.libraryFilePath = ''
        self.figPath = ''
    
    def startAnalysis(self):
        
        logging.info("Analysis started.")
        
        try: 
            
            library_loader = LibraryLoadingStrategy(self.libraryFilePath)
            library = library_loader.load_library()
            ppm_tolerance = float(self.ppmToleranceEdit.text())
            minmatchedpeaks = int(self.minMatchedPeaksEdit.text())
            PrecursorIonMassTolerance = float(self.precursorIonMassToleranceEdit.text())
            intensity_threshold = float(self.intensityEdit.text()) 
            cosine_threshold = float(self.cosineEdit.text()) 
            lowerscan = self.scansLowerSpin.value()
            higherscan = self.scansUpperSpin.value() 
            # Retrieve the state of the generate plots checkbox
        
            generate_plots = self.generatePlotsCheckbox.isChecked()
         # set it to Default scans in the input file if nothing is set

            for InputFilePath in self.InputFilePaths:


                analyzer = QueryTargetedSpectrum(InputFilePath,intensity_threshold)

                if higherscan == 1:
                      # Adjust based on your implementation
                    total_scans = analyzer.get_scans()     
                    higherscan = total_scans 

                 # Get the number of scans from the input field

                result_dicts = main_processing_function(lowerscan, higherscan, analyzer, library, PrecursorIonMassTolerance,
                                                        cosine_threshold,ppm_tolerance, minmatchedpeaks, self.figPath, 
                                                        generate_plots = generate_plots
                                                       )
            
                # higherscan = self.scansUpperSpin.value() 
            
                save_results(result_dicts, self.figPath, InputFilePath)
            
            logging.info("Analysis completed successfully.")
            
        except Exception as e:
            
            logging.error(f"An error occurred: {str(e)}")
            

# def main():
#     app = QApplication(sys.argv)
#     ex = Identification_GUI()
#     ex.show()
#     sys.exit(app.exec_())

# if __name__ == '__main__':
#     main()


