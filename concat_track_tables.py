#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@program: concat_track_tables.py
@author: Samuel Piquer-Esteban
@date: 29 May 2023 
@version: v1

"""

#Set program name
__program__ = 'concat_track_tables.py'

#Import moduls
import pandas as pd
from argparse import ArgumentParser
import os 

#Program functions
def check_existence_directory_parameter(dir_path, dir_type, parameter):
    """
    This function will check if the provided directory path exist.

    Parameters
    ----------
    dir_path : str
        Path for the provided directory.
    dir_type : str
        Type of directory (normally, "input" or "output").
    parameter : str
        The name of the parameter.

    Raises
    ------
    Exception
        If the provided directory does not exit raises an exception.

    Returns
    -------
    None.

    """
    if dir_path is None:
        pass
    else:
        if os.path.isdir(dir_path)==False:
            frase = ''.join(['Error! The provided ', dir_type, ' Directory path is not a directory!\n Check the provided ', parameter, ' parameter!'])
            raise Exception(frase)


def treat_output_file_prefix_name(prefix, file_name):
    """
    This function creates the name for the file treating 
    the --output_name_prefix parameter.

    Parameters
    ----------
    prefix : str
        The provided output_name_prefix parameter.
    file_name : str
        The name of the file.

    Returns
    -------
    final_name : str
        The final name of the file.

    """
    if prefix is None:
        final_name = file_name
    else:
        final_name = prefix + '_' + file_name
    return final_name


def treat_duplicated_outfiles(directory, file_extension, file_name):
    """
    This function generates a unique outputfile name for the provided
    file_name and output directory.

    Parameters
    ----------
    directory : str
        The output directory to check.
    file_extension : str
        The extension of the provided output file_name.
    file_name : str
        The name of the provided output file.

    Returns
    -------
    str
        If file_name is not in directory, returns the same provided file_name.
        Else it will return a new_file_name.

    """
    #Get files in directory
    files = os.listdir(directory)
    #Get copy for file name
    new_file_name = file_name
    #IF/ELSE
    if new_file_name in files:
        #Count
        counter = 1
        #While in list try new number combination
        while new_file_name in files:
            new_file_name = ''.join([file_name.split(file_extension)[0],'(',str(counter),')',file_extension])
            counter += 1
        #Return final file name
        return new_file_name
    else:
        #If file not in directory return the given file_name
        return file_name


def treat_output_directory_parameter_outfiles(file_name, file_extension, outputdir_path):
     """
     This function will treat the output_directory parameter and return
     the full path for the output file_name provided.

     Parameters
     ----------
     file_name : str
         Name for the output file.
     file_extension : str
         The extension of the provided output file_name.
     outputdir_path : str
         Output directory path.

     Returns
     -------
     outputfile_path : str
         Path for the output file.

     """
     ##If not path provided create file in current/ else create file in provided path
     if outputdir_path is None:
         #Check if file_name exits in current directory
         #If exits get new name
         out_name = treat_duplicated_outfiles(os.getcwd(), file_extension, file_name)
         #Get final full path
         outputfile_path = os.path.join(os.getcwd(), out_name)
     else:
         #Check if file_name exits in provided output directory
         #If exits get new name
         out_name = treat_duplicated_outfiles(outputdir_path, file_extension, file_name)
         #Get final full path
         outputfile_path = os.path.join(outputdir_path, out_name)
     return outputfile_path


def get_list_files_in_directory(directory, pattern):
    """
    This function gets a list with the files 
    of the provided directory for the provided pattern.

    Parameters
    ----------
    directory : str
        Path to the directory.
    pattern : str
        Extension pattern to recognize files of interest.
        
    Raises
    ------
    Exception
        If there are no files detected with the provided pattern raises an exception.

    Returns
    -------
    files : list
        List of the files with pattern present in the directory.
    
    """
    #Get files in directory
    files = os.listdir(directory)
    
    #Filter to retain only files with extension
    files = [x for x in files if x.endswith(pattern)]
    
    #Check if files were detected
    if len(files) == 0:
        frase0 = 'Error! There are no files with pattern extension "'
        frase1 = '" in the provided Input Directory!'
        raise Exception(''.join([frase0, pattern, frase1]))
    
    return files


def check_track_tables_headers(track_tables_list):
    """
    This function checks if all track tables have the same colnames.

    Parameters
    ----------
    track_tables_list : list
        List with the different pandas dataframe track tables

    Returns
    -------
    None.

    """
    #Set results
    cols_len = []
    set_cols = set()
    
    for track in track_tables_list:
        #Len cols in track
        cols_len.append(len(track.columns))
        #Aupdate set
        set_cols.update(list(track.columns))
        
    #Check results
    if not all(i == len(set_cols) for i in cols_len):
        print('Error! Some Track Tables present different column names! Check for pseudo-NA columns!')
        

#Main Program
def main():
    #Setting Arguments
    parser = ArgumentParser()
    ##Parameter raw_file_path
    parser.add_argument(
            '-i','--input_directory',
            action = 'store',
            required = True,
            help ='Input Directory. Indicate the path to Input Directory. It must contains all track tables using ".tsv" extension.'
    )
    ##Parameter output_directory
    parser.add_argument(
            '-o','--output_directory', 
            action = 'store',
            required = False,
            help = 'Output Directory (Optional). Indicate the path to the Output Directory. Output files will be created in the current directory if not indicated.'
    )
    ##Parameter output_file_prefix
    parser.add_argument(
            '-op','--output_name_prefix', 
            action = 'store',
            required = False,
            help = 'Output Name Prefix (Optional). Indicate prefix name for the output files. '
    )
    
    #Print program name
    print('')
    print(__program__)
    
    #Process common arguments
    args = parser.parse_args()
    input_path = args.input_directory
    outputdir_path = args.output_directory
    output_name_prefix = args.output_name_prefix
    
    #Check if the files exists and load files/ Try-except-else
    try:
        #0) Initial checks and steps
        
        ##Check that provided input directory exist
        check_existence_directory_parameter(input_path, 'Input', '--input_directory')
    
        ##Check that provided output directory exist
        check_existence_directory_parameter(outputdir_path, 'Output', '--output_directory')
        
        ##List tsv files in input dir
        input_files = get_list_files_in_directory(input_path, '.tsv')
        input_files_full = [os.path.join(input_path, i) for i in input_files]
        
        #1) Try to load files
        
        #Set list with track_tables
        track_tables = []
        
        #Message
        print('\n1) Loading files...')
        
        #Try to read each track table
        for file in input_files_full:
            print(file)
            track_tables.append(pd.read_csv(file, sep = '\t'))
        
        #Check that all col_names are the same
        check_track_tables_headers(track_tables)
        
        #2) Creating Track Tables
        ##Message
        print('\n2) Concatenating Track Tables...')
        
        #CONCAT 
        result = pd.concat(track_tables)
        print('(Number of Rows, Number of Cols):',result.shape)
        
        ##Save final file
        ###Get final file name
        out_table = treat_output_file_prefix_name(output_name_prefix, 'concat_track_table.tsv')
        ###Treat output_directory parameter / Get full output file path
        outputfile_table = treat_output_directory_parameter_outfiles(out_table, '.tsv', outputdir_path)
        ###Message
        print('\nSaving Concat Track Table in file:')
        print(outputfile_table)
        result.to_csv(outputfile_table, header = True, index = False, sep = '\t')
    
    except Exception as ex:
        print('\nThe system returned the following exception: ', ex)
    
    finally:
        #Print empty line for aesthetic purposes
        print('')    

if __name__ == '__main__':
    main()

