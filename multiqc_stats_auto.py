#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@program: multiqc_stats_auto.py
@author: Samuel Piquer-Esteban
@date: 21 Apr 2023 
@version: v1

"""

#Set program name
__program__ = 'multiqc_stats_auto.py'

#Import moduls
import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter
from functools import reduce
import os 

#Program Functions
def process_path(path):
    """
    This function gets the full path to the multiqc files of interest.

    Parameters
    ----------
    path : str
        Path for the multiqc_data directory.

    Returns
    -------
    path : str
        Full path to the multiqc file of interest.

    """
    if path.endswith('/'):
        path=path[:-1]
    path=path+'/multiqc_fastqc.txt'
    return path


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


def check_processed_paths_lenght(processed_paths, processed_names):
    """
    This function checks if the --processed_files_names and --processed_files_paths
    parameters have the same lenght.

    Parameters
    ----------
    processed_paths : list
        This to multiqc_data directories for the different qc steps.
    processed_names : list
        Names of the corresponding qc steps.

    Raises
    ------
    Exception
        If this two parameters do not have the same lenght raises exception.

    Returns
    -------
    None.

    """
    if len(processed_paths) != len(processed_names):
        raise Exception('Error! --processed_files_names and --processed_files_paths must have the same length!')


def check_fastq_PAIRED_patterns(fastq_pattern, r1_files_pattern, r2_files_pattern):
    """
    This function checks if the provided Fastq PAIRED patterns are compatible
    with the provided Fastq File Pattern.

    Parameters
    ----------
    fastq_pattern : str
        The provided Fastq File Pattern.
    r1_files_pattern : str
        The provided R1 File Pattern.
    r2_files_pattern : str
        The provided R2 File Pattern.

    Raises
    ------
    Exception
        If the patterns are not compatible raises an exception.

    Returns
    -------
    None.

    """
    if r1_files_pattern.endswith(fastq_pattern) == False or r2_files_pattern.endswith(fastq_pattern) == False:
        raise Exception('Error! The provided Fastq Pattern does not match with some of the PAIRED files patterns!\nCheck the provided patterns parameters!')
    elif r1_files_pattern == r2_files_pattern:
        raise Exception('Error! The provided PAIRED files patterns are identical!\nCheck the provided patterns parameters!')
    elif r1_files_pattern == fastq_pattern or r2_files_pattern == fastq_pattern:
        raise Exception('Error! The provided Fastq Pattern is identical to some of the PAIRED files patterns!\nCheck the provided patterns parameters!')


def process_multiqc_counts(multiqc_table, fastq_pattern, r2_files_pattern, sep, n_sep, col_name):
    """
    This function processes the multiqc file of interest to get the stat counts.

    Parameters
    ----------
    multiqc_table : pandas data.frame
        The multiqc stats file as pandas data.frame.
    fastq_pattern : str
        The provided Fastq File Pattern.
    r2_files_pattern : str
        The provided R2 File Pattern.
    sep : str
        Sample Name separator.
    n_sep : int
        Sample Name separator appereance.
    col_name : str
        Name of the qc step that will be adopted as column name.
    
    Raises
    ------
    Exception
        If duplicated sample names after filtering R2 and removing 
        Fastq Patterns are found raises an Exception.
    
    Returns
    -------
    stats_counts : pandas data.frame
        Table with counts per sample for the qc step.

    """
    #Filter rows by File sufix (R2 pattern) -> Keep only SE and R1 PE
    r2_sample_end = r2_files_pattern.replace(fastq_pattern, '')
    multiqc_table_no_R2 = multiqc_table[~multiqc_table.Sample.str.endswith(r2_sample_end)]
    #Select columns of interest make copy as new pandas df
    stats_counts = multiqc_table_no_R2[['Sample','Total Sequences']].copy()
    #Message number of samples evaluated
    print('\nQC Step:', col_name)
    print('Number of samples evaluated:',len(stats_counts))    
    #Process sample names to keep only basename
    stats_counts['Sample'] = stats_counts.Sample.str.split(sep).apply(lambda x:x[:n_sep]).str.join(sep)
    #Rename 'Total Sequences' column to col_name
    stats_counts = stats_counts.rename(columns={'Total Sequences':col_name})
    #Change float to int
    stats_counts[col_name] = stats_counts[col_name].astype('int')
    #Check if there are duplicates in Samples column 
    if any(stats_counts.duplicated(subset=['Sample'])):
        print('\nSample names:')
        print(list(stats_counts['Sample']))
        raise Exception('Some samples present duplicated names! Check your files!')
    #Return result pandas df
    return stats_counts


def calculate_len_ori_tables(multiqc_table):
    """
    This function gets the min and max lenght of the sequences.

    Parameters
    ----------
    multiqc_table : pandas data.frame
        The multiqc stats file as pandas data.frame.

    Returns
    -------
    result : list
        List with min and max read lenght values in the provided multiqc_table.

    """
    #Get Max and Min lenghts per sample
    ##Set min and max list
    min_list = []
    max_list = []
    ##Iter rows for Sequence length column
    for i,r in multiqc_table.iterrows():
        ##Get temp info
        temp_seq_len = r['Sequence length']
        ##If/else
        if '-' in str(temp_seq_len):
            temp_min = temp_seq_len.split('-')[0]
            min_list.append(temp_min)
            temp_max = temp_seq_len.split('-')[1]
            max_list.append(temp_max)
        else:
            temp_min = temp_seq_len
            min_list.append(temp_min)
            temp_max = temp_seq_len
            max_list.append(temp_max)
    ##Convert list to integers
    min_list = [int(float(i)) for i in min_list]
    max_list = [int(float(i)) for i in max_list]
    #Calculate stats
    min_length_total = min(min_list)
    max_length_total = max(max_list)
    #Put in a list
    result = [min_length_total, max_length_total]
    #Return list
    return result


def calculate_reads_stats(track_table, col_name):
    """
    This function generates QC Stats of interest based on the Track Table
    for the provided column name QC step.

    Parameters
    ----------
    track_table : pandas data.frame
        The provided track table.
    col_name : str
        The name of the QC Step of interest (column name).

    Returns
    -------
    result : list
        List with the calculated stats.

    """
    #Calculate total reads stats
    total_count = track_table[col_name].sum()
    mean_total_count = track_table[col_name].mean()
    min_total_count = track_table[col_name].min()
    max_total_count = track_table[col_name].max()
    #Calculate % reads left
    perc_reads_left = (total_count/track_table['input_raw_reads'].sum())*100
    #Calculate total reads left versus input raw samples(percentages)
    reads_left_per_sample = (track_table[col_name]/track_table['input_raw_reads'])*100
    #Calcula left reads stats
    mean_left_per_sample = reads_left_per_sample.mean()
    sd_left_per_sample = reads_left_per_sample.std()
    min_left_per_sample = reads_left_per_sample.min()
    max_left_per_sample = reads_left_per_sample.max()
    #Put in a list
    result = [total_count, mean_total_count, min_total_count, max_total_count, perc_reads_left, mean_left_per_sample, sd_left_per_sample, min_left_per_sample, max_left_per_sample]
    #Return list
    return result


def read_filter_file(filter_file_path):
    """
    This function tries to read the provided filter file.

    Parameters
    ----------
    filter_file_path : str
        Path to the provided filter file.

    Raises
    ------
    Exception
        If the file is empty raises an exception.

    Returns
    -------
    samples_2_filter : list
        List with the sample names to filter.

    """
    #Open file
    with open(filter_file_path) as filter_file:
        
        #Read file as list
        samples_2_filter = [line.strip() for line in filter_file]
        
        #Remove empty lines
        samples_2_filter = list(filter(None, samples_2_filter))
        
        #Check if the file is empty
        if len(samples_2_filter) == 0:
            raise Exception('Error! All lines are empty in the provided Filter File!')
    #Return list
    return samples_2_filter


def check_samples_name_in_track_table(counts_track_table_df, list_values):
    """
    This function checks if all values in the filter list are in the Sample
    column from the generated track table.

    Parameters
    ----------
    counts_track_table_df : pandas dataframe 
        The provided counts track table to check.
    list_values : list
        The provided list of sample values to filter.

    Returns
    -------
    None.
        
    """
    result = all(elem in list(counts_track_table_df['Sample']) for elem in list_values)
    if result == False:
        raise Exception('Error! Some of the provided values in the Filter File are not in the Counts Track Table! Check your Filter File!')


def get_percentage_track_table(counts_track_table_df):
    """
    This function will generate the percentages track table from 
    the counts track table. 

    Parameters
    ----------
    counts_track_table_df : pandas dataframe
        The provided counts track table.

    Returns
    -------
    percentage_track_table : pandas dataframe
        The resulting percentages track table.

    """
    
    #Set empty pandas dataframe (with Sample and input_raw_reads from counts_table)
    percentage_track_table = counts_track_table_df[['Sample','input_raw_reads']].copy()
    
    #For each column in counts_track_table_df calculate percentage with respect to raw reads
    columns_of_interest = [i for i in counts_track_table_df.columns.values if i not in ['Sample','input_raw_reads']]
    for col_name in columns_of_interest:
        #Calculate total reads left versus input raw samples(percentages)
        temp_reads_left_perc_sample = (counts_track_table_df[col_name]/counts_track_table_df['input_raw_reads'])*100
        #Save temp counts column
        percentage_track_table[col_name] = counts_track_table_df[col_name].copy()

        #Save temp percentage column
        percentage_track_table[col_name + '(%)[reference: raw reads]'] = temp_reads_left_perc_sample
    
    #Return percentages_track_table
    return percentage_track_table
    
    
#Main Program
def main():
    #Setting Arguments
    parser = ArgumentParser(formatter_class=RawTextHelpFormatter)
    ##Parameter raw_file_path
    parser.add_argument(
            '-i1','--raw_file_path',
            action = 'store',
            required = True,
            help ='Raw Samples Data Path. Indicate rute to multiqc_data directory.'
    )
    ##Parameter processed_files_paths
    parser.add_argument(
            '-i2','--processed_files_paths',
            nargs = '+',
            required = True,
            help = 'Processed Samples Data Paths. Indicate rute to multiqc_data directories. Use like: -i2 path1/multiqc_data/ path2/multiqc_data/ path3/multiqc_data/ [...]'
    )
    ##Parameter processed_files_names
    parser.add_argument(
            '-i2_names','--processed_files_names',
            nargs = '+',
            required = True,
            help = 'Processed Samples Data Names. Indicate the corresponding names for the processing steps. Use like: -i2_names name1 name2 name3 [...]. Must be in the same order than -i2 parameter.'
    )
    ##Parameter sample_name_sep
    parser.add_argument(
            '-sep','--sample_name_sep', 
            action = 'store',
            default = '_',
            required = False,
            help = 'Sample Name separator (Optional). Samples names will by separated by the provided character [Default="_"].'
    )
    ##Parameter sample_name_sep_appereance
    parser.add_argument(
            '-n_sep','--sample_name_sep_appereance', 
            action = 'store',
            required = False,
            type = int,
            default = 1,
            help = 'Sample Name separator Appereance (Optional). Indicate by which appereance of the separator the file name can be divided in sample_name + rest [Default=1 appereance].'
    )
    ##Parameter fastq_pattern
    parser.add_argument(
            '-p','--fastq_pattern', 
            action = 'store',
            default = '.fastq.gz',
            required = False,
            help = 'Fastq File Extension (Optional) [Default:".fastq.gz"]. Indicate the extension to identify Fastq files.'
    )
    ##Parameter r1_patterns
    parser.add_argument(
            '-r1','--r1_pattern', 
            action = 'store',
            default = '_1.fastq.gz',
            required = False,
            help = 'R1 File Pattern (Optional) [Default:"_1.fastq.gz"]. Indicate the pattern to identify R1 PAIRED Fastq files.'
    )
    ##Parameter r2_patterns
    parser.add_argument(
            '-r2','--r2_pattern', 
            action = 'store',
            default = '_2.fastq.gz',
            required = False,
            help = 'R2 File Pattern (Optional) [Default:"_2.fastq.gz"]. Indicate the pattern to identify R2 PAIRED Fastq files.'
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
    ##Parameter filter_mode
    parser.add_argument(
            '-f', '--filter_mode',
            action = 'store_true',
            required = False,
            help = 'Filter Mode (Optional). If indicated, it will enable Filter Samples Mode.'
    )
    ##Parameter metadata_table
    parser.add_argument(
            '-t','--filter_file', 
            action = 'store',
            required = False,
            help = 'Filter File [Required for Filter Mode]. Indicate the path to the Filter File [expected format: TXT file]. Each line in the Filter File will correspond to a sample name that you wish to filter out.'
    )
    
    #Print program name
    print('')
    print(__program__)
    
    #Process common arguments
    args = parser.parse_args()
    raw_path = args.raw_file_path
    processed_paths = args.processed_files_paths
    processed_names = args.processed_files_names
    sample_name_separator = args.sample_name_sep
    sep_appereance = args.sample_name_sep_appereance
    fastq_pattern = args.fastq_pattern
    r1_files_pattern = args.r1_pattern
    r2_files_pattern = args.r2_pattern
    outputdir_path = args.output_directory
    output_name_prefix = args.output_name_prefix
    filter_mode = args.filter_mode
    
    #Process specific arguments for filter mode
    if filter_mode == True:
        filter_file_path = args.filter_file
    
    #Process paths(add file of interest)
    raw_path_file=process_path(raw_path)
    processed_path_files=[process_path(i) for i in processed_paths]
    
    #Check if the files exists and load files/ Try-except-else
    try:
        #0) Initial checks
    
        ##Check that provided output directory exist
        check_existence_directory_parameter(outputdir_path, 'Output', '--output_directory')
        
        ##check that --processed_files_names and --processed_files_paths parameters have the same lenght 
        check_processed_paths_lenght(processed_paths, processed_names)
    
        ##Check PAIRED Fastq patterns
        check_fastq_PAIRED_patterns(fastq_pattern, r1_files_pattern, r2_files_pattern)
        
        #1) Try to load files
        
        #Message
        print('\n1) Loading files...')
        
        #Load raw samples stats file
        print('\nLoading Raw Multiqc Stats File:')
        print(raw_path_file)
        stats_raw = pd.read_csv(raw_path_file, sep = '\t')
        
        #Load processed samples stats files for different steps on a list
        print('\nLoading Process Steps Multiqc Stats Files:')
        stats_processed_steps = []
        for file_path in processed_path_files:
            print(file_path)
            stats_processed_steps.append(pd.read_csv(file_path, sep = '\t'))
        
        #Try to load filter file if filter mode is activated
        if filter_mode == True:
            print('\nLoading Filter File:')
            print(filter_file_path)
            samples_2_filter = read_filter_file(filter_file_path)
        
        #2) Creating Track Tables
        ##Message
        print('\n2) Creating Track Tables...')
        
        #2.1) Create counts track table
        ##Get raw samples counts table
        stats_raw_counts = process_multiqc_counts(stats_raw, fastq_pattern, r2_files_pattern, sample_name_separator, sep_appereance, 'input_raw_reads')
        
        ##Get samples counts for different qc steps
        stats_process_counts_list = []
        for i in range(len(stats_processed_steps)):
            temp = process_multiqc_counts(stats_processed_steps[i], fastq_pattern, r2_files_pattern, sample_name_separator, sep_appereance, processed_names[i])
            stats_process_counts_list.append(temp)
        
        ##Merge count tables
        ###Combine lists of pandas data.frames
        dfs_counts_list = [stats_raw_counts] + stats_process_counts_list
        ###Merge data.frames using reduce()
        counts_track_table = reduce(lambda left,right: pd.merge(left,right,on=['Sample'], how='outer'), dfs_counts_list)
        
        #Try to filter the track table if filter mode is activated and additional check (check that the provided samples names are in the generated track table)
        if filter_mode == True:
            ##Message
            print('\nApplying Filter Step:')
            print('\nTotal Samples to filter:', len(samples_2_filter))
            
            
            ##Check that provided values are in the Sample column
            check_samples_name_in_track_table(counts_track_table, samples_2_filter)
            
            ##Filter counts_track_table
            print('\nSamples Before Filtering:', len(counts_track_table))
            counts_track_table = counts_track_table[~counts_track_table['Sample'].isin(samples_2_filter)]
            print('Samples After Filtering:', len(counts_track_table))
            
            #Previous steps for filter stats_raw and stats_processed_steps
            ##NOTE: We re-add the sep, to avoid unspecific matches. For example, sample1 could match with sample1, sample10, sample100, etc.
            ##NOTE2: We add the fastq_pattern to try que build the fastq name (in some SE cases, there will not be a separator and this is the way that I came up with to take into account this cases)  
            
            #Get samples_2_filter + sample_name_separator
            samples_2_filter_sep = [i + sample_name_separator for i in samples_2_filter]
            #Get samples_2_filter + fastq_pattern
            samples_2_filter_file = [i + fastq_pattern for i in samples_2_filter]
            
            ##Filter stats_raw
            stats_raw = stats_raw[~stats_raw['Filename'].str.contains('|'.join(samples_2_filter_sep + samples_2_filter_file))]
            
            ##Filter stats_processed_steps
            for i in range(len(stats_processed_steps)):
                temp_stats = stats_processed_steps[i]
                temp_stats_filtered = temp_stats[~temp_stats['Sample'].str.contains('|'.join(samples_2_filter_sep + samples_2_filter_file))]
                stats_processed_steps[i] = temp_stats_filtered    
        
        ##Save final file
        ###Get final file name
        out_counts_table = treat_output_file_prefix_name(output_name_prefix, 'counts_track_table.tsv')
        ###Treat output_directory parameter / Get full output file path
        outputfile_counts_table = treat_output_directory_parameter_outfiles(out_counts_table, '.tsv', outputdir_path)
        ###Message
        print('\nSaving Counts Track Table in file:')
        print(outputfile_counts_table)
        counts_track_table.to_csv(outputfile_counts_table, header = True, index = False, sep = '\t')
        
        #2.2) Create percentages track table
    
        #Get percentages track table
        percentage_track_table = get_percentage_track_table(counts_track_table)
        
        ##Save final file
        ###Get final file name
        out_percentages_table = treat_output_file_prefix_name(output_name_prefix, 'percentages_track_table.tsv')
        ###Treat output_directory parameter / Get full output file path
        outputfile_percentages_table = treat_output_directory_parameter_outfiles(out_percentages_table, '.tsv', outputdir_path)
        ###Message
        print('\nSaving Percentages Track Table in file:')
        print(outputfile_percentages_table)
        percentage_track_table.to_csv(outputfile_percentages_table, header = True, index = False, sep = '\t')
        
        #3) Create the Steps QC Stats File
        
        ##Message
        print('\n3) Creating QC Stats Table...')
        
        ##Set future columns names(header)
        header_see_MultiQC = ['Adapters', 'N_Content', 'Reads_Quality_Hist', 'Reads_Avg_Quality_Range']
        header_len_ori_tables = ['Reads_Min_Length', 'Reads_Max_Length']
        header_reads = ['Total_Reads', 'Avg_Reads_Per_Sample', 'Min_Reads', 'Max_Reads']
        header_left = ['Total_Reads_Left(%)', 'Avg_Reads_Left_Per_Sample(%)', 'SD_Reads_Left_Per_Sample(%)', 'Min_Reads_Left(%)', 'Max_Reads_Left(%)']
        header = ['Step'] + header_see_MultiQC + header_len_ori_tables + header_reads + header_left
          
        ##Set results list and repetitive empty list
        final_lines = []
        see_MultiQC = ['']*len(header_see_MultiQC)    
        
        ##Calculate QC Stats for Raw Samples
        init_stats = ['Raw Samples'] + see_MultiQC + calculate_len_ori_tables(stats_raw) + calculate_reads_stats(counts_track_table, 'input_raw_reads')
        final_lines.append(init_stats)
        
        ##Calculate QC Stats for the rest of the QC Stesp
        for i in range(len(stats_processed_steps)):
            temp_step = processed_names[i]
            temp_stats = [temp_step] + see_MultiQC + calculate_len_ori_tables(stats_processed_steps[i]) + calculate_reads_stats(counts_track_table, processed_names[i])
            final_lines.append(temp_stats)
        
        ##Save final file
        ###Create pandas data.frame
        qc_stats_results = pd.DataFrame(final_lines, columns = header)
        ###Get final file name
        out_stats_table = treat_output_file_prefix_name(output_name_prefix, 'qc_stats_table.xlsx')
        ###Treat output_directory parameter / Get full output file path
        outputfile_stats_table = treat_output_directory_parameter_outfiles(out_stats_table, '.xlsx', outputdir_path)
        ###Message
        print('\nSaving QC Stats Table in file:')
        print(outputfile_stats_table)
        qc_stats_results.to_excel(outputfile_stats_table, header = True, index = False, float_format = '%.3f')
        
    except Exception as ex:
        print('\nThe system returned the following exception: ', ex)
    
    finally:
        #Print empty line for aesthetic purposes
        print('')    

if __name__ == '__main__':
    main()