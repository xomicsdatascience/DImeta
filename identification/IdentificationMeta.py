
import pyteomics.mzml
from pyteomics import mzxml, mzml
import os
import pandas as pd
import glob
import numpy as np
import heapq
from collections import defaultdict


class QueryTargetedSpectrum:
    """in this class, read query spectrum from either mzxml or mzml files,
    produce a real time library depending on the precusor ion range from a specific scan number,
    read how many scans in a specific input file,
    label all spectra in the real-time library and name it as target spectrum"""

    def __init__(self, filepath,intensity_threshold=3000):
        
        self.filepath = filepath
        _, file_extension = os.path.splitext(filepath)
        self.intensity_threshold = intensity_threshold
        
        if file_extension.lower() == '.mzml':
            self.file_type = 'mzml'
            self.tmp = pyteomics.mzml.read(filepath, use_index=True)
        elif file_extension.lower() == '.mzxml':
            self.file_type = 'mzxml'
            self.tmp = pyteomics.mzxml.read(filepath, use_index=True)
        else:
            raise ValueError(f"Unsupported file format: {file_extension}")

            
    def get_query_spectrum(self, scan) -> list:
        if self.file_type == 'mzml':
            return self._get_query_spectrum_mzml(scan)
        elif self.file_type == 'mzxml':
            return self._get_query_spectrum_mzxml(scan)
    
    def _get_query_spectrum_mzml(self, scan) -> list:
        # Implementation for mzML files
        ms_level = self.tmp.get_by_index(scan)['ms level']
        mzs, intens, label = [], [], []

        if ms_level == 2:
            mz = self.tmp.get_by_index(scan)['m/z array']
            inten = self.tmp.get_by_index(scan)['intensity array']
            filter_mz_inten = [(x, y) for x, y in zip(mz, inten) if y > self.intensity_threshold]   # filter the input spectrum intensity 
            if filter_mz_inten:
                mzs, intens = zip(*filter_mz_inten)
            else:
                mzs, intens = [], []
            label = [str(scan)] * len(mzs)
        else:
            label = []

        return list(zip(mzs, intens, label))

    def _get_query_spectrum_mzxml(self, scan) -> list:
        ms_level = self.tmp[scan]['msLevel']
        mzs, intens, label = [], [], []

        if ms_level == 2:
            mz = self.tmp[scan]['m/z array']
            inten = self.tmp[scan]['intensity array']
            filter_mz_inten = [(x, y) for x, y in zip(mz, inten) if y > self.intensity_threshold]  # filter the input spectrum intensity 
            if filter_mz_inten:
                mzs, intens = zip(*filter_mz_inten)
            else:
                mzs, intens = [], []
            label = [str(scan)] * len(mzs)
        else:
            label = []

        return list(zip(mzs, intens, label))

    
    def get_realtime_lib(self, scan, library, PIMT) -> list:
        if self.file_type == 'mzml':
            return self._get_realtime_lib_mzml(scan, library, PIMT)
        elif self.file_type == 'mzxml':
            return self._get_realtime_lib_mzxml(scan, library, PIMT)
            

    def _get_realtime_lib_mzml(self, scan, library, PIMT) -> list:
        
        ms_level = self.tmp.get_by_index(scan)['ms level']
        
        realtime_lib = []
        
        if ms_level == 2:
            precursor = float(self.tmp.get_by_index(scan)['precursorList']['precursor'][0]['isolationWindow']['isolation window target m/z'])
            realtime_lib = [item for item in library if (precursor - PIMT < float(item['precursormz']) < precursor + PIMT)]
        
        return realtime_lib

    def _get_realtime_lib_mzxml(self, scan, library, PIMT) -> list:
        
        '''Library should be reformatted by the DIMA LibraryReformatted class and define the maximum peaks in each library mass spectrum,
        PIMT: Precursor Ion Mass Tolerance, defines mass tolerance for MS1 '''
        
        ms_level = self.tmp[scan]['msLevel']
        
        realtime_lib = []
        if ms_level == 2:
            precursor = float(self.tmp[scan]['precursorMz'][0]['precursorMz'])
            realtime_lib = [item for item in library if (precursor - PIMT < float(item['precursormz']) < precursor + PIMT)]
        
        return realtime_lib
    
    
    def get_scans(self):
        if self.file_type == 'mzml':
            return self._countscans_mzml()
        elif self.file_type == 'mzxml':
            return self._countscans_mzxml()
    
    def _countscans_mzml(self)->int:
        # Initialize a counter
        scan_count = 0
        # Read the .mzML file and iterate over its contents
        for spectrum in self.tmp:
        # Increment the counter for each spectrum found
            scan_count += 1
        return scan_count

    def _countscans_mzxml(self)->int:
        # Initialize a counter
        scan_count = 0
        # Read the .mzXML file and iterate over its contents
        for spectrum in self.tmp:
        # Increment the counter for each spectrum found
            scan_count += 1
        return scan_count
    
    
    def get_precusorMZ(self,scan):
        if self.file_type == 'mzml':
            return self._getprecusor_mzml(scan)
        elif self.file_type == 'mzxml':
            return self._getprecusor_mzxml(scan)
    
    def _getprecusor_mzml(self,scan)->float:
        # Initialize a counter
        precusormz= float(self.tmp.get_by_index(scan)['precursorList']['precursor'][0]['isolationWindow']['isolation window target m/z']
                if 'precursorList' in self.tmp.get_by_index(scan) else 'Nan')
        return precusormz

    def _getprecusor_mzxml(self,scan)->float:
        
        precusormz= float(self.tmp.get_by_index(scan)['precursorMz'][0]['precursorMz']
                if 'precursorMz' in self.tmp.get_by_index(scan) else 'Nan')
        return precusormz
    
    def get_compensation_voltage(self,scan):
        
        if self.file_type == 'mzml':
            return self._get_compensation_voltage_mzml(scan)
        
        elif self.file_type == 'mzxml':
            return self._get_compensation_voltage_mzxml(scan)
    
    def _get_compensation_voltage_mzml(self,scan)->str:
        
        comp_vol = float(self.tmp.get_by_index(scan)['FAIMS compensation voltage'] if 'FAIMS compensation voltage' in self.tmp.get_by_index(scan) else 'Nan')
        return comp_vol

    def _get_compensation_voltage_mzxml(self,scan)->str:
        
        comp_vol= float(self.tmp.get_by_index(scan)['compensationVoltage'] if 'compensationVoltage' in self.tmp.get_by_index(scan) else 'Nan')
        return comp_vol
    
      
    def get_target_spectrum(self, scan, library, PIMT) -> list:
        # This method can remain largely unchanged, as it calls get_realtime_lib
        real = self.get_realtime_lib(scan, library, PIMT)
        target_spectrum = []
        for i in range(len(real)):
            target_spectrum.extend(list(zip(real[i]['mz'], real[i]['intensity'], [str(i)] * len(real[i]['mz']))))
        return sorted(target_spectrum, key=lambda x: x[0])


    
    

def within_tolerance_ppm(mz1, mz2, ppm):
    """ Check if mz2 is within the "ppm" tolerance of mz1. """
    tolerance = mz1 * ppm / 1e6
    return abs(mz1 - mz2) <= tolerance

def within_tolerance_da(mz1, mz2, da):
    """ Check if mz2 is within the "Da" tolerance of mz1. """
    tolerance = da
    return abs(mz1 - mz2) <= tolerance

def match_spectrum(query_spectrum, target_spectrum, ppm_tolerance):
    """
    Match a spectrum to a query spectrum within a given PPM tolerance and return the matched peaks along with their three-dimensional values.
    
    :param query_spectrum: List of (m/z, intensity, value) tuples for the query spectrum.
    :param target_spectrum: List of (m/z, intensity, value) tuples for the target spectrum.
    :param ppm_tolerance: PPM tolerance for matching.
    :return: Matched peaks as a list of tuples (query_mz, query_intensity, query_value, target_mz, target_intensity, target_value).
    """
    matches = []
    for query_mz, query_intensity, query_label in query_spectrum:
        for target_mz, target_intensity, target_label in target_spectrum:
            if within_tolerance_ppm(query_mz, target_mz, ppm_tolerance):
                matches.append((query_mz, query_intensity, query_label,target_mz, target_intensity, target_label))
                #break  # Assuming only one match per query peak
            #if within_tolerance_da(query_mz, target_mz, da):
                #matches.append((query_mz, query_intensity, query_label, target_mz, target_intensity, target_label)) 
    return matches

def cosine_similarity(vector1, vector2):
    """
    Calculate the cosine similarity between two vectors.
    
    Args:
    vector1 (np.array): First vector.
    vector2 (np.array): Second vector.
    
    Returns:
    float: Cosine similarity between vector1 and vector2.
    """
    dot_product = np.dot(vector1, vector2)
    norm_vector1 = np.linalg.norm(vector1)
    norm_vector2 = np.linalg.norm(vector2)
    similarity = dot_product / (norm_vector1 * norm_vector2)
    return similarity

def macc_score(numMatchedPeaks, cosineScore)->float:
    """
    Taking both the cosine score and the number of matching fragment peaks into consideration
    """
    return numMatchedPeaks ** (1 / 5) * cosineScore

def custom_sort(item):
    """sort the list according to the macc_score"""
    return macc_score(len(item[1]),item[0])
 
def group_tuples_by_same_value(tuples_list,same_value_position):
    """
    Groups tuples in the provided list by their last value.

    :param tuples_list: List of tuples.
    :return: List of lists, where each sublist contains tuples with the same last value(label).
    """
    grouped_tuples = {}
    for t in tuples_list:
        key = t[same_value_position]  # Last element of the tuple
        grouped_tuples.setdefault(key, []).append(t)

    return list(grouped_tuples.values())

def normalize_to_100(numbers):
    """normalize the spectrum to 0 to 100"""
    if not numbers or max(numbers) == 0:
        return [0] * len(numbers)

    max_value = max(numbers)
    return [(x / max_value) * 100 for x in numbers]

def filter_tuples(tuples)->list:
    """
    Filters a list of tuples to ensure no duplicate first values,
    and for each unique fourth value, keeps the tuple with the highest second value.

    Parameters:
    - tuples: A list of tuples, where each tuple has at least four elements.

    Returns:
    - A list of filtered tuples.
    """   
    # Step 1: Filter based on the fourth value, choosing the tuple with the highest second value
    max_second_by_fourth = defaultdict(lambda: (None, float('-inf'), None, None))
    for t in tuples:
        if t[1] > max_second_by_fourth[t[3]][1]:
            max_second_by_fourth[t[3]] = t

    filtered_by_second = list(max_second_by_fourth.values())

    # Step 2: Ensure no duplicate first values
    unique_by_first = {}
    for t in filtered_by_second:
        # If first value already exists and its second value is lower, replace it
        if t[0] not in unique_by_first or unique_by_first[t[0]][1] < t[1]:
            unique_by_first[t[0]] = t

    # Final list of tuples after filtering
    final_tuples = list(unique_by_first.values())
    
    return final_tuples




