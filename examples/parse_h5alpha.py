# -*- coding: utf-8 -*-
"""
This program is use to save gene sequences to individual files from the genbank files

Created on Wed Mar  4 21:07:19 2020

@author: taomihog
"""
import pandas as pd

path_genebankfile = '200304_NEB_H5alpha_genbank.txt'
path_fasta =        '200304_NEB_H5alpha.txt'
path_location =     'summary.csv'
output_folder = 'parsed/'






# so tedious!!
def genelocation(path,path_save):
    m=0
    with open(path, 'r') as file:
    
        
        lines=file.readlines()
        direction = 1 # for complement direction = -1
        for i in range(1000000):
            try:
                line = lines[i]
                if '     gene            ' in line:
                    if not 'CDS' in lines[i+2]: continue
                    location = line.replace('     gene            ','')
                    if location.count('complement') ==0:
                        direction = 1
                    else:
                        location = location.replace('complement(','').replace(')','')
                        dirction = -1
                    
                    for k in range(len(location)):
                        k += 1
                        if not location[0:k].isdecimal():
                            a = int(location[0:k-1])
                            b = int(location[k+1:len(location)-1])
    #                        print (k,a,b)
                            break
                    
                    for j in range(20):
                        j += 7 #save some time
                        if '/product=' in lines[j+i]:
                            break
                    gene_name1  = lines[j+i].replace('/product=','').replace('\n','').replace('  ','')
                    if gene_name1.count('\"') == 2:
                        gene_name1=gene_name1.replace('\"','')
                        pass
                    elif gene_name1.count('\"') == 1:
                        gene_name1  = gene_name1.replace('\"','')
                        for w in [1,2,3]:
                            if lines[j+i+w].count('\"') == 0:
                                gene_name1 += lines[j+i+w]
                            else:
                                gene_name1 += lines[j+i+w]
                                break
                        gene_name1 = gene_name1.replace('\"','').replace('\n','').replace('  ','')
                    else:
                        print ('error')#should be situation like this
                    print (i,gene_name1,direction,a,b,b-a+1)
                    data.loc[m] = [i,gene_name1,direction,a,b,b-a+1]
                    m+=1
            except:
                break
    #    print (data)
        data.to_csv(path_save,index = False)





#delete the illegal chars for file names
import string
def clean_file_name(filename):
    valid_chars = "-_.() %s%s" % (string.ascii_letters, string.digits)
    cleaned_filename = ''.join(char if char in valid_chars else '_' for char in filename)
    return cleaned_filename





if __name__ == "__main__":
    
    
    import os
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    data = pd.DataFrame(columns=['line','gene_name','direction','start','end','length'])
    
    #get all genes locations
    genelocation(path_genebankfile, path_location)

    # save individual files
    with open(path_fasta, 'r') as file:
        temp = file.readlines()
        genome = ''
        for seq in temp:
            # print (seq)
            if not ('>' in seq):
                genome += seq.replace('\n','')
                
    data = pd.read_csv(path_location)
    
    
    for i in range(data.count(axis = 0)[0]): 
        # I can only use this slow method, not familiar with dataFrame or string
        filename = output_folder + clean_file_name(data['gene_name'].iloc[i][1:]) + '.txt'
        print(filename)
        with open(filename, 'w') as file:
            file.write(genome[data['start'].iloc[i]:data['end'].iloc[i]])