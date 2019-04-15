#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 14:46:18 2019

@author: Kim NOEL
"""

################################
# Import libraries
################################
# To compute running time
import time

# To remove files of previous run
import csv

# To compute modulo for cycles rules
from math import fmod

# To plot
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Import support python file
import preparePlot as fun


################################
# Parameters that the user can change
################################
length_interval = 1

list_files =[
             'testing.csv',
             'training_30.csv', 'training_50.csv', 'training_70.csv',
             'training_90.csv', 'training_95.csv', 'training_100.csv',
             'validation.csv'
             ]

# ========= USER: SELECT THE FILE TO EXCTRACT THE DATA FROM =========
for file_name in ['testing.csv','validation.csv']: #list_files:
    
    print('\n')
    print('File name:', file_name)
    start_time = time.time()
    
    with open(file_name, newline='') as csvfile:
        file = csv.reader(csvfile, delimiter='\t', quotechar='|')
        row_count = fun.len_big_file(file_name)
        
        fun.show_elapsed_time(start_time, 'row_count is extracted: ')
        
        
    # Total number of proteins in the file
    total_prot = row_count//33

    print('row_count =', row_count)
    print('total_prot =', total_prot)
    
    
    ################################
    # Parameters that the user can change
    ################################
    
    num_prot_start = 0
    num_prot_end = length_interval + 0
    
    # Retreive data in df which has 'length_interval' number of rows
    # This avoid to work with a df with 100_000 rows
    while(num_prot_end < total_prot-1):
        
        #num_prot_start = 0
        #num_prot_end = sliding_end  # last included    
        
        last_df_index = num_prot_end - num_prot_start + 1
        
        df = fun.dataframe_prot(file_name, num_prot_start, num_prot_end)
        
        for df_index in range(0,  last_df_index):
    
            # Name of the protein without the # for export purposes
            name_prot = df.iloc[df_index]['ID']
            name_prot = name_prot.replace('#', '-')
            
            # References for the amino acids (aa)
            total_aa = len(df.iloc[df_index]['Primary'])
            list_aa = list( df.iloc[df_index]['Primary'] )
            
            # Print the possibles values of num_prot !== df_index
            print('\n')
            print('Select the ID number of the protein between 0 and ', total_prot-1)
            
            # Print what is planned
            print('The protein with ID number', num_prot_start + df_index, 'and labeled', name_prot, 'is being plotted.')
            print('The protein has', total_aa, 'amino acids.')
            
            
            ################################
            # Remove image from previous run
            ################################
            #for f in glob.glob("protein*.png"):
            #    os.remove(f)
            
            
            ################################
            # Define parameters for plotting
            ################################
            
            # For plotting the protein
            mpl.rcParams['legend.fontsize'] = 10
            dpi = 150
            
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            plt.figure(num=None, figsize=(80, 60), dpi=dpi, facecolor='w', edgecolor='k')
            ax = Axes3D(fig)
            ax.set_title("Protein from " + file_name + " indexed " + str(num_prot_start + df_index) + " and named " + name_prot.replace('-', '#') )
            ax.set_xlabel('x axis [pm]')
            ax.set_ylabel('y axis [pm]')
            ax.set_zlabel('z axis [pm]')
            
            
            bary_center_X = []
            bary_center_Y = []
            bary_center_Z = []
            bary_center = []
            
            list_color = ['dodgerblue', 'red', 'mediumpurple', 'darkorchid', 'orange']
            index_color = 0
            
            ################################
            # Plot amino acids in protein
            ################################
            for i in range(0, total_aa -1 ):
                # Extract 3x3 matrix of coordinates
                x = df.iloc[df_index]['tertiary_X'][i]
                y = df.iloc[df_index]['tertiary_Y'][i]
                z = df.iloc[df_index]['tertiary_Z'][i]
                
                # Compute gravity point of the 9 elementary points that form
                # the 3 segments == representaion of the aa
                xG = (x[0] + x[1] + x[2])/3
                yG = (y[0] + y[1] + y[2])/3
                zG = (z[0] + z[1] + z[2])/3
                
                if(xG != 0):
                    # Add new coordinates of the aa point if it is not [0, 0, 0]
                    bary_center_X.append(xG)
                    bary_center_Y.append(yG)
                    bary_center_Z.append(zG)
                    bary_center.append([xG, yG, zG])
                    
                    # Plot amino acids
                    ax.scatter(xG, yG, zG, c=fun.from_acid_to_color(list_aa)[i], marker='o', s=10)
                    
                    # Plot lines between aa
                    ax.plot(bary_center_X, bary_center_Y, bary_center_Z, '-o', c=list_color[index_color])
                    
                    # Plot the figure of the protein being built
                    #fig.savefig("protein%d.png" % i, dpi=dpi)
                    
                else:
                    # If point with coordinate 0 is found, then start a new polyline
                    # to join the aa
                    bary_center_X = []
                    bary_center_Y = []
                    bary_center_Z = []
                    index_color = int(fmod(index_color+1,5))    
            # End for i        
            #plt.show()
            
            
            
            ################################
            # Export the plots 
            ################################
            # Create file of 
            fig.set_size_inches((8, 6), forward=False)
            fig.savefig("protein_" + file_name[:-4] + "__" + str(num_prot_start + df_index) + "__" + name_prot + ".png", dpi=dpi)
            
            # Plot image and rotate around azimuth
            #for ii in range(0,360,36):
            #    ax.view_init(elev=15., azim=ii)
            #    fig.set_size_inches((8, 6), forward=False)
                
                #fig.savefig("proteindegree%d.png" % ii, dpi=dpi, c='grey')
            
            
            ################################
            # Compute running time in s 
            ################################
            fun.show_elapsed_time(start_time, 'Running time end program: ')
            
            
        #End for df_index   
        del df
        
        num_prot_start = num_prot_end + 1
        num_prot_end = num_prot_end + length_interval
    #end while        
    
            
            
    #End for df_index 
# End for file_name in list_files:
