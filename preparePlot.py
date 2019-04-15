#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 09:44:34 2019

@author: Kim
"""

# Libraries 
import pandas as pd 
import numpy as np
import csv
from math import fmod

# To access csv with specific line number
from itertools import islice

# To compute running time
import time

# To access a line in a file
import linecache


def len_big_file(file_name):
    
    with open(file_name, newline='') as csvfile:
        file = csv.reader(csvfile, delimiter='\t', quotechar='|')        
        
        a = 0
        b = 100_000
        
        ante_prev_b = a
        prev_b = a
        
        len_is_close = 0
        
        while (a != b):
            name_a = linecache.getline(file_name, a*33 + 2)
            name_b = linecache.getline(file_name, b*33 + 2)
            
    
            # if string is empty == if b is out of range
            if not name_b:
                # Save current lower bound
                    
                ante_prev_b = prev_b
                prev_b = b
                
                      
                if (len_is_close == 1):
                    #b -=1
                    b = (a+b)//2
                else:
                    # From [a, b] to [a, (a+b)//2]
                    b = (a+b)//2
                    if (b == ante_prev_b):
                        print('len_is_close')
                        len_is_close = 1
                        a = b
                        b = prev_b

            # if string has any, if b is in range
            else:
                ante_prev_b = prev_b
                prev_b = b
                
                if (len_is_close == 1):
                    a = b
                    b *=2
                else:
                    # Double the width of the interval [a,b] into [a, 2b]
                    b *= 2
                    if (b == ante_prev_b):
                        len_is_close = 1
                        a = prev_b
            
        total_prot = b+1
        row_count = total_prot*33
    return row_count


def dataframe_prot(file_name, num_prot_start, num_prot_end):
    # index is the list that contains the names of the proteins 
    # as indexes for the dataframe
    index = np.array([])
    
    # primary is the list that contains the names of each amino acid (aa) 
    # that compose the protein
    primary = np.array([])
    
    # tertiary_X is the list that contains the x, y, z coordinates of each aa 
    #  is the number of proteins indexed
    #  is the number of aa in the protein indexed i
    tertiary_X = []
    tertiary_Y = []
    tertiary_Z = []
    tertiary_tmp_ = []
    
    # islice function has index shifted by 1
    start_row_slice = 33*num_prot_start
    end_row_slice = 33*(num_prot_end + 1)
    row_num_ = start_row_slice + 1
    
    # Column number to group by 3 the coordinates x,y,z
    col_num_ = 1
    
    ####################################
    start_time = time.time()
    ####################################

    with open(file_name, newline='') as csvfile:
        file = csv.reader(csvfile, delimiter='\t', quotechar='|')
        
        ####################################
        show_elapsed_time(start_time, 'File is opened: ')
        
        ####################################
        
        
        for row in islice(file, start_row_slice, end_row_slice):
            
            # index
            if (fmod(row_num_, 33) == 2):
                # Remove the symbols [, ], and "
                index = np.append(index, str(row)[2:-2] )
            
            # Primary  
            elif (fmod(row_num_, 33) == 4):
                # Remove the symbols [, ], and "
                primary = np.append(primary, str(row)[2:-2] )
                
            # TODO
            #elif (fmod(row_num_, 33) == 6):
                #mydict['EVOLUTIONARY'].append(row)
            
            # tertiary_X
            elif (fmod(row_num_, 33) == 28): 
                col_ = []
                # i is all the number on this row
                for i in row:
                    col_.append(float(i))
                    # Group by 3
                    if( fmod(col_num_, 3) == 0 ):
                        tertiary_tmp_.append( col_ )
                        col_ = []
                    col_num_ += 1
                
                # Group by protein
                tertiary_X.append( tertiary_tmp_ )
                tertiary_tmp_ = []
    
    
            # tertiary_Y
            elif (fmod(row_num_, 33) == 29): 
                col_ = []
                # i is all the number on this row
                for i in row:
                    col_.append(float(i))
                    # Group by 3
                    if( fmod(col_num_, 3) == 0 ):
                        tertiary_tmp_.append( col_ )
                        col_ = []
                    col_num_ += 1
                # Group by protein
                tertiary_Y.append( tertiary_tmp_ )
                tertiary_tmp_ = []
    
    
            # tertiary_Z
            elif (fmod(row_num_, 33) == 30): 
                col_ = []
                # i is all the number on this row
                for i in row:
                    col_.append(float(i))
                    # Group by 3
                    if( fmod(col_num_, 3) == 0 ):
                        tertiary_tmp_.append( col_ )
                        col_ = []
                    col_num_ += 1
                
                # Group by protein
                tertiary_Z.append( tertiary_tmp_ )
                tertiary_tmp_ = []
    
            # TODO
            #elif (fmod(row_num_, 33) == 32):
                #mydict['MASK'].append(row)
                #row_num_ = row_num_
            
            row_num_ += 1 
    # End for row in file: 
    
    ####################################
    show_elapsed_time(start_time, 'Loop is finished: ')
    
    
    del tertiary_tmp_, row_num_

    #####################################
    
    columns = np.array(['ID', 'Primary'])
    data = np.array([index, primary]).T
    
    df = pd.DataFrame(data, index=index, columns=columns)
    
    df['tertiary_X'] = tertiary_X
    df['tertiary_Y'] = tertiary_Z
    df['tertiary_Z'] = tertiary_Z
    
    return df


# Function to select a color for each amino acid
def from_acid_to_color(myarray):
    for i in range(0, len(myarray), 1):
        output = myarray
        if myarray[i] == 'A':
            output[i] = 'black'
            
        elif myarray[i] == 'C':
            output[i] = 'grey'
            
        elif myarray[i] == 'D' :
            output[i] = 'firebrick'
            
        elif myarray[i] == 'E' :
            output[i] = 'red'
            
        elif myarray[i] == 'F' :
            output[i] = 'sienna'
            
        elif myarray[i] == 'G' :
            output[i] = 'tan'
            
        elif myarray[i] == 'H' :
            output[i] = 'gold'
            
        elif myarray[i] == 'I' :
            output[i] = 'darkkhaki'
            
        elif myarray[i] == 'K' :
            output[i] = 'olivedrab'
            
        elif myarray[i] == 'L': 
            output[i] = 'chartreuse'
            
        elif myarray[i] == 'M' :
            output[i] = 'darkgreen'
            
        elif myarray[i] == 'N' :
            output[i] = 'mediumspringgreen'
            
        elif myarray[i] == 'P' :
            output[i] = 'darkcyan'
            
        elif myarray[i] == 'Q' :
            output[i] = 'deepskyblue'
            
        elif myarray[i] == 'R' :
            output[i] = 'slategrey'
            
        elif myarray[i] == 'S' :
            output[i] = 'royalblue'
            
        elif myarray[i] == 'T' :
            output[i] = 'navy'
            
        elif myarray[i] == 'V' :
            output[i] = 'blue'
            
        elif myarray[i] == 'W' :
            output[i] = 'mediumpurple'
            
        elif myarray[i] == 'Y' :
            output[i] = 'darkorchid'
        
        elif myarray[i] == '' :
            output[i] = 'white'

    return output

# Not used yet
def distance_3D(p1, p2):
    """
    p1 and p2 are lists with 3 elements [x, y, z]
    """
    result = (p1[0] - p2[0])**2
    result = (p1[1] - p2[1])**2 + result
    result = (p1[2] - p2[2])**2 + result
    result = result**(1/2)
    
    return result

def show_elapsed_time(start_time, description):
    """
    start_time: start time of the clock
    description: string foor description purpose
    """
    elapsed_time = time.time() - start_time
    print(description,
          round( elapsed_time//3600), 'hr', 
          round( fmod(elapsed_time//60, 60) ), 'min', 
          round( fmod(elapsed_time, 60) ), 's')
    
