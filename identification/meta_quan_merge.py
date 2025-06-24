
import numpy as np 
import pandas as pd
import os




class Meta_df_Merge:
    
    def __init__(self, folder_path):
        self.folder_path = folder_path

    def reformat_df(self, df):
        
        '''Reformat the DataFrame with CV labeled and precursor rounded'''
        df['precursor'] = df['PrecursorMZ'].round(2)
        df['label'] = df['precursor'].astype(str) + '_' + df['Compensation Voltage'].astype(str)
        selected_columns = df[['label', 'Ion_count']]
        selected_columns = selected_columns.sort_values(by='Ion_count', ascending=False).drop_duplicates(subset='label', keep='first')
        selected_columns = selected_columns.set_index('label')
        
        return selected_columns

    def merge_dfs(self):
        '''Merge DataFrames from all Excel files in the specified directory'''
        dfs = []
        
        for filename in os.listdir(self.folder_path):
            if filename.endswith('.xlsx'):
                path = os.path.join(self.folder_path, filename)
                df = pd.read_excel(path)
                df_reformatted = self.reformat_df(df)
                file_label = filename.replace('.xlsx', '')
                df_reformatted.columns = [f"{file_label}" for col in df_reformatted.columns]
                dfs.append(df_reformatted)

        final_df = pd.concat(dfs, axis=1)
        
        return final_df

    def save_final_df(self, final_df, output_folder=None, output_filename='merged_output.csv'):
        
        '''Save the final merged DataFrame to a specified output folder and filename'''
        
        if output_folder is None:
            output_folder = self.folder_path  # Use the initial folder path if no output folder is specified
        output_path = os.path.normpath(os.path.join(output_folder, output_filename))
        
        final_df.to_csv(output_path)
        
        print(f"Saved merged DataFrame to {output_path}")





