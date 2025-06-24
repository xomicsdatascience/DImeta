#!/usr/bin/env python
# coding: utf-8

import sys
import os
import uuid
import matplotlib.pyplot as plt
import os
import glob
import numpy as np
import heapq
import re

from LibraryHandling import LibraryLoadingStrategy, LibraryReformat, LibrarySaveStrategy, read_path
from IdentificationMeta import QueryTargetedSpectrum, match_spectrum, group_tuples_by_same_value, filter_tuples, cosine_similarity, custom_sort, macc_score, normalize_to_100
import pandas as pd 

#result_dict = defaultdict(list)

    
def get_spectra(analyzer, scan_index, library, PrecursorIonMassTolerance):
    query_spectrum = analyzer.get_query_spectrum(scan_index)
    realtime_library = analyzer.get_realtime_lib(scan_index, library, PrecursorIonMassTolerance)
    target_spectrum = analyzer.get_target_spectrum(scan_index, library, PrecursorIonMassTolerance)
    return query_spectrum, realtime_library, target_spectrum


def match_and_calculate_cosine_similarity(query_spectrum, target_spectrum, ppm_tolerance, minmatchedpeaks):
    matched_peaks = match_spectrum(query_spectrum, target_spectrum, ppm_tolerance)
    cosine_scores = []
    for matched_spectrum in group_tuples_by_same_value(matched_peaks, -1):
        filtered_matches = filter_tuples(matched_spectrum)
        if len(filtered_matches) >= minmatchedpeaks:
            cosine_score = cosine_similarity([t[1] for t in filtered_matches], [t[4] for t in filtered_matches])
            cosine_scores.append((cosine_score, filtered_matches))
    return cosine_scores

def sanitize_filename(filename):
    # Replace any invalid characters with an underscore
    return re.sub(r'[<>:"/\\|?*]', '_', filename)
    
def generate_plot(compound_info, cos, scan_index, fig_path):
# Normalize intensities to 100
    real_query_spectrum = cos[1]
    library_mz = compound_info['mz']
    library_intensity = normalize_to_100(compound_info['intensity'])
    query_mz = [x[0] for x in real_query_spectrum]
    query_intensity = normalize_to_100([x[1] for x in real_query_spectrum])
    plt.rcParams['font.family'] = 'Arial'
    # Plot library spectrum
    plt.vlines(library_mz, [0], library_intensity, colors="red", label='Library Spectrum')
    # Plot query spectrum (inverted)
    plt.vlines(query_mz, [0], [-item for item in query_intensity], colors="blue", label='Query Spectrum')  
    # Plot horizontal line at y=0
    plt.axhline(y=0, color='gray', alpha = 0.8, linewidth=1)
    plt.legend()
    plt.title(f"{compound_info.get('name', 'Unknown Compound')} | Cosine Score: {cos[0]:.2f} | Scan: {scan_index}") 
    # Sanitize the compound name for the filename
    compound_name = compound_info.get('name', 'Unknown Compound')
    sanitized_name = sanitize_filename(compound_name)
    # Create the filename using the sanitized compound name
    filename = f"{sanitized_name}_{scan_index}.svg"
    # Save the figure
    plt.savefig(os.path.join(fig_path, filename), dpi=300, bbox_inches='tight')
    # plt.show()
    plt.close()


def main_processing_function(lowerscan,higherscan, analyzer, library, PrecursorIonMassTolerance, 
                             cosine_threshold, ppm_tolerance, minmatchedpeaks, fig_path, generate_plots=False):
    result_dict = { 'PrecursorMZ': [],'Compensation Voltage': [], 'Cosine_score': [], 'Ion_count': [], 'Scan': [], 
                   'Compound': [], 'CompoundMZ': [], 'Adduct': [], 'Formula': [], 'Macc_score': [], 'Matched_peaks': [] }
    #result_dict = {}
    for scan_index in range(lowerscan,higherscan):
        query_spectrum, realtime_library, target_spectrum = get_spectra(analyzer, scan_index, library, PrecursorIonMassTolerance)
        
        cosine_scores = match_and_calculate_cosine_similarity(query_spectrum, target_spectrum, ppm_tolerance, minmatchedpeaks)
        
        if cosine_scores:           
            # Filter cosine scores based on the threshold
            filtered_scores = [cos for cos in cosine_scores if cos[0] > cosine_threshold]           
            if not filtered_scores:
                # If no scores meet the threshold, select the strongest top one
                filtered_scores = [max(cosine_scores, key=lambda x: x[0])]   
            # Process each selected score
            for cos in filtered_scores:
                number = int(cos[1][0][-1])
                ioncount = round(sum([it[1] for it in cos[1]]), 3)
                
                # Retrieve compound information
                try:
                    compound_info = realtime_library[number]
                except IndexError:
                    # Handle case where index is out of bounds
                    compound_info = {}               
#         if cosine_scores:
#             cos = sorted(cosine_scores, key=lambda x: custom_sort(x), reverse=True)[0]  # here only select the top 1 matched spectrum
#             number = int(cos[1][0][-1])
#             ioncount = round(sum([it[1] for it in cos[1]]), 3)
#             try:
#                 compound_info = realtime_library[number]
            
#             except IndexError:
#                 # If `number` is out of bounds for `realtime_library`
#                 compound_info = {}
            
            #result_dict, compound_info, number = update_results(cos, scan_index, analyzer, realtime_library, result_dict)
            #compound_info = realtime_library.get(number, {})
                result_dict['PrecursorMZ'].append(str(analyzer.get_precusorMZ(scan_index)))
                result_dict['Cosine_score'].append(cos[0])
                result_dict['Ion_count'].append(ioncount)
                result_dict['Scan'].append(scan_index)
                result_dict['Compound'].append(compound_info.get('name', ''))
                result_dict['CompoundMZ'].append(compound_info.get('precursormz', ''))
                result_dict['Adduct'].append(compound_info.get('precursortype', '').upper())
                result_dict['Formula'].append(compound_info.get('formula', '').upper())
                result_dict['Macc_score'].append(macc_score(len(cos[1]), cos[0]))
                result_dict['Matched_peaks'].append(len(cos[1])) 
                result_dict['Compensation Voltage'].append(str(analyzer.get_compensation_voltage(scan_index))) 
                # **Generate plot for each match if enabled**
                if generate_plots:
                    generate_plot(compound_info, cos, scan_index, fig_path)
    return result_dict

    
def save_results(result_dict, fig_path, mzml_file_path):
    
    df = pd.DataFrame(result_dict)
    df.to_excel(os.path.join(fig_path, f"{mzml_file_path.split('/')[-1].split('.')[0]}.xlsx"), index=False)







