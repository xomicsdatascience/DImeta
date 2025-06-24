import sys
from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QVBoxLayout, QFileDialog, QLabel, QStatusBar, QFormLayout, QTextEdit
import os
import pandas as pd
from meta_quan_merge import Meta_df_Merge




class Meta_Quantification_GUI(QWidget):
    
    def __init__(self):
        super().__init__()
        self.initUI()
    
    def initUI(self):
#         self.setWindowTitle('Metabolites quan')
        layout = QVBoxLayout()
        self.formLayout = QFormLayout()

        self.btn_select_input = QPushButton('Select Input Folder', self)
        self.btn_select_input.clicked.connect(self.select_input_folder)
        layout.addWidget(self.btn_select_input)

        self.input_folder_label = QLabel('Input Folder: Not selected', self)
        layout.addWidget(self.input_folder_label)

        self.btn_select_output = QPushButton('Select Output Folder', self)
        self.btn_select_output.clicked.connect(self.select_output_folder)
        layout.addWidget(self.btn_select_output)

        self.output_folder_label = QLabel('Output Folder: Not selected', self)
        layout.addWidget(self.output_folder_label)

        self.btn_process_files = QPushButton('Process Files in Folder', self)
        self.btn_process_files.clicked.connect(self.process_files)
        layout.addWidget(self.btn_process_files)

        # Text edit for log messages
        self.log_window = QTextEdit(self)
        self.log_window.setReadOnly(True)
        layout.addWidget(self.log_window)

        self.resize(700, 600)
        
        self.setLayout(layout)
        
    def log(self, message):
        
        self.log_window.append(message)
    
    def select_input_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Select Folder")
        if folder:
            self.input_folder = folder
            self.input_folder_label.setText(f'Input Folder: {folder}')
            self.log(f"Selected input folder: {folder}")

    def select_output_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Select Folder")
        if folder:
            self.output_folder = folder
            self.output_folder_label.setText(f'Output Folder: {folder}')
            self.log(f"Selected output folder: {folder}")

    def process_files(self):
        if hasattr(self, 'input_folder') and hasattr(self, 'output_folder'):
            try:
                self.log("Starting file processing...")
                processor = Meta_df_Merge(self.input_folder)
                final_df = processor.merge_dfs()
                message = processor.save_final_df(final_df, self.output_folder)
                self.log("Process finished")
            except Exception as e:
                self.log(f"Error: {str(e)}")
        else:
            self.status_bar.showMessage("Please select both input and output folders.")
            

# if __name__ == '__main__':
    
#     app = QApplication(sys.argv)
#     ex = Meta_Quantification_GUI()
#     ex.show()
#     sys.exit(app.exec_())
# def main():
#     app = QApplication(sys.argv)
#     mainWindow = Main_GUI()
#     mainWindow.show()
#     sys.exit(app.exec_()) 
