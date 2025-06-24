

from pyteomics import mgf
import os
from pyteomics import mzml
import pyteomics
import matplotlib.pyplot as plt
import pandas as pd
import glob
import heapq


class LibraryLoadingStrategy:
    
    '''class for loading Libraries'''
    def __init__(self, file_path):
        self.file_path = file_path

    def load_library(self):
        """Automatically load library based on file extension."""
        _, file_extension = os.path.splitext(self.file_path)
        if file_extension.lower() == '.mgf':
            return self._load_mgf()
        elif file_extension.lower() == '.msp':
            return self._load_msp()
        else:
            raise ValueError(f"Unsupported file format: {file_extension}")

    def _load_mgf(self)-> list:
        """Load a library from a .mgf file."""
        # Assuming mgf.IndexedMGF is a class from an external library
        spectra_library = mgf.IndexedMGF(self.file_path, index_by_scans=True)
        return spectra_library

    def _load_msp(self)-> list:
        """Load a library from a .msp file."""
        
        with open(self.file_path, 'r', encoding='utf-8') as file:
            spectra = []
            spectrum = {} # Start with an empty dictionary
            for line in file:
                line = line.strip().lower() # remove leading and trailing whitespace characters, format uppercase to lowercase
                if not line:  # Skip empty lines
                    continue
                if line.startswith('name:'):   
                    if 'mz' in spectrum:         # Check if the current spectrum already has 'mz' key
                        spectra.append(spectrum)
                        spectrum = {}            # Reset for new spectrum
                if ':' in line:                 # Meta data line
                    key, value = line.split(':', 1)
                    spectrum[key.strip()] = value.strip()
                else:                          # Spectrum data line
                    mz, intensity = line.split()
                    spectrum.setdefault('mz', []).append(float(mz))
                    spectrum.setdefault('intensity', []).append(float(intensity))
            if 'mz' in spectrum:       # append the last spectrum
                spectra.append(spectrum)
            return spectra
    
    @classmethod
    def combine_libraries(cls, file_paths) -> list:
        """Combine all spectrums in all input libraries from the provided file paths into a new library."""
        combined_library = []
        for file_path in file_paths:
            loader = cls(file_path)
            combined_library.extend(loader.load_library())
        return combined_library


class LibraryReformat:
    
    def __init__(self, topnum):
        self.topnum = topnum

    def reformat_spectrum(self, spectrum):
        if len(spectrum['mz']) <= self.topnum:
            return spectrum.copy()

        # Use a min heap to keep track of top elements
        paired = zip(spectrum['mz'], spectrum['intensity'])
        top_pairs = heapq.nlargest(self.topnum, paired, key=lambda pair: pair[1])

        top_mz, top_intensity = zip(*top_pairs) if top_pairs else ([], [])

        new_spectrum = {**spectrum, 'mz': list(top_mz), 'intensity': list(top_intensity)}
        new_spectrum['num peaks'] = str(len(top_mz))

        return new_spectrum

    def reformat_library(self, library):

        reformatted_library = []
        for spectrum in library:
            if 'precursor_type' in spectrum:
                spectrum['precursortype'] = spectrum.pop('precursor_type')
            if 'precursormz' in spectrum:
                reformatted_library.append(self.reformat_spectrum(spectrum))
        return reformatted_library
    
class LibrarySaveStrategy:

    @classmethod
    def save_library_to_msp_class(cls, library, output_file_path):
        """Save the library to a new .msp file using a class method."""
        with open(output_file_path, 'w', encoding='utf-8') as file:
            for spectrum in library:
                for key, value in spectrum.items():
                    if key not in ['mz', 'intensity']:
                        file.write(f"{key.capitalize()}: {value}\n")

                mz_intensities = zip(spectrum.get('mz', []), spectrum.get('intensity', []))
                for mz, intensity in mz_intensities:
                    file.write(f"{mz} {intensity}\n")

                file.write("\n")    

def read_path(dir)-> list:
    return([os.path.join(dir, file) for file in os.listdir(dir)])



