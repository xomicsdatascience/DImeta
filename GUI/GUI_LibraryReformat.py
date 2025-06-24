
import sys
import os
import uuid
from PyQt5.QtWidgets import (QApplication, QWidget, QVBoxLayout, QHBoxLayout, QPushButton,QDialog,QSpinBox, QPushButton,
                             QLabel, QTextEdit, QListWidget, QLineEdit, QFileDialog, QMessageBox,
                             QInputDialog)
from PyQt5.QtCore import Qt
from LibraryHandling import LibraryLoadingStrategy, LibraryReformat, LibrarySaveStrategy


class LibraryReformatterGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()
        
    def initUI(self):
        self.setWindowTitle('Library Reformatter')
        self.setGeometry(200, 200, 640, 480)
        
        layout = QVBoxLayout()
        
        # Input Files Label and List
        self.input_label = QLabel("Input Files (.msp/.mgf):")
        layout.addWidget(self.input_label)
        self.listWidget = QListWidget()
        layout.addWidget(self.listWidget)
        

        # Output Directory
        self.output_label = QLabel("Output Directory:")
        layout.addWidget(self.output_label)
        self.output_path_display = QLineEdit()
        layout.addWidget(self.output_path_display)
        

        
        # TopNum Selection
        self.topnum_label = QLabel("Enter the number of top intensities for reformatting:")
        layout.addWidget(self.topnum_label)
        self.topnum_spinBox = QSpinBox(self)
        self.topnum_spinBox.setMinimum(1)
        self.topnum_spinBox.setMaximum(15)
        self.topnum_spinBox.setValue(1)  # Default value
        layout.addWidget(self.topnum_spinBox)
         
        # Log Area
        self.log_label = QLabel("Log:")
        layout.addWidget(self.log_label)
        self.log_area = QTextEdit()
        self.log_area.setReadOnly(True)
        layout.addWidget(self.log_area)
        
        # Buttons
        button_layout = QHBoxLayout()
        self.open_files_button = QPushButton('Select Libraries')
        self.open_files_button.clicked.connect(self.open_files)
        button_layout.addWidget(self.open_files_button)
         # Buttons
        self.save_directory_button = QPushButton('Select Output Directory')
        self.save_directory_button.clicked.connect(self.select_output_directory)
        button_layout.addWidget(self.save_directory_button)
        self.process_files_button = QPushButton('Process and Reformat Files')
        self.process_files_button.clicked.connect(self.process_files)
        button_layout.addWidget(self.process_files_button)
        
        layout.addLayout(button_layout)
        
        self.setLayout(layout)
        
        self.filepaths = []
        self.output_directory = ''

        
    def open_files(self):
        options = QFileDialog.Options()
        files, _ = QFileDialog.getOpenFileNames(self, "Open Mass Spectra Files", "",
                                                "Supported Files (*.mgf *.msp);;All Files (*)", options=options)
        if files:
            self.filepaths = files
            self.listWidget.clear()
            self.listWidget.addItems(files)
        else:
            QMessageBox.warning(self, "File Selection", "No files were selected.")
            
    def select_output_directory(self):
        directory = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if directory:
            self.output_directory = directory
            self.output_path_display.setText(directory)
        else:
            QMessageBox.warning(self, "Directory Selection", "No directory was selected.")
    
    def log(self, message):
        self.log_area.append(message)
        
    def process_files(self):
        if not self.filepaths:
            QMessageBox.warning(self, "Error", "Please select one or more input files.")
            return
        if not self.output_directory:
            QMessageBox.warning(self, "Error", "Please select an output directory.")
            return

        topnum = self.topnum_spinBox.value()
        self.log(f"Starting file processing with top {topnum} peaks in library mass spectra...")
        # 这里添加文件处理逻辑
        
        try:
            self.log("Combining libraries...")
            combined_library = LibraryLoadingStrategy.combine_libraries(self.filepaths)
            self.log("Libraries combined successfully.")

            self.log(f"Reformatting library for top {topnum} peaks in library mass spectra...")
            reformat = LibraryReformat(topnum)
            reformatted_library = reformat.reformat_library(combined_library)
            self.log("Library reformatted successfully.")

            random_filename = f"reformatted_library_{uuid.uuid4()}.msp"
            output_file_path = os.path.join(self.output_directory, random_filename)
            self.log("Saving reformatted library...")
            LibrarySaveStrategy.save_library_to_msp_class(reformatted_library, output_file_path)
            self.log(f"Library saved successfully to {output_file_path}.")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred: {str(e)}")
            self.log(f"An error occurred: {str(e)}")

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = LibraryReformatterGUI()
    ex.show()
    sys.exit(app.exec_())


