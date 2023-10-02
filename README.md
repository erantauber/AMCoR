# AMCoR
AMCoR stands for analysis of DNA methylation in DNA coding sequences (CDS). This script would take alignments of CDS and identify pairs of codons (dyads) that form a CpG site. More detailed information is included in the file. Please cite:  Creasey L and Tauber E (2023), No codon is an island: Epigenetic selection for CpG codon dyad.
###################################
### Analyse all data in 'Test data' ###
###################################

# This script will analyse all sequences in the 'test data' folder.
# Currently it is only set up for this specific data set, meaning a FASTA file named in the way the example data is
# final results for individual genes are outputted into the test_results folder, a few general result documents are also generated, these are outputted to the source folder.
# The code could be adapted to look for .txt files or other files with some effort. 

### Instructions!

# This is very simple, put all test data into the test_data folder, then execute this script.
# With our specific data set there is a gene_correspondents file to match ENSG code to gene names. This may lead to errors with other data sets. 
# If you add new sequences that do not have a gene ID then you will likely get an error. I can edit the code so that it gives an 'NA' rather then an error if needed. 


# This code is split into 4 sections 


# First is the function used for identifying dyads
# I recommend against changing this unless you want the general information document to present something different
# If you do add or change data, you will need to extract it in the second section

#Needed for generating the csv
import csv
#os is needed to iterate over the folder so you can do all data points there in 
import os
#A counting package needed further on. 
from collections import Counter


#Function for counting Dyads, more functionality can easily be added if needed. 
def dyad_counter(data, file_ID, species_name, dyad_dic, dyad_freq, dyad_species_name, dyad_ratio_dic, average_length):
        #Create/reset objects
        i = 0
        CpG = 0
        dyad_locations = []
        dyad_seq=[]
        line_num=''
        gene_ID=''
        missing=0
        #A while loop that goes up in codons and checks between the codons for CpGs, it then counts these.
        while i< len(data):
            #counts up the index each loop
            i=i+3
            #checks if the index is above the size of the sequence after the addition. 
            if i>=len(data):
                break
        #checks if both nessecary bases are C and G, if so it adds 1 to the count. 
            if data[i]=='G' and data[i-1]=='C':
                CpG+=1
                codons = data[i-3:i+3]
                #This section updates the 3 dictionaries with the data we need for analysing the codons. 
                #If the dyad position is already in the main dicitonary then it needs to update the relevant dictionaries
                if i in dyad_dic:
                    #These lines add the species name to the species name dictionary, this data is mainly visual for the end
                    temp_name_list = dyad_species_name[i]
                    temp_name_list.append(species_name)
                    temp_name_dic = {i:temp_name_list}
                    dyad_species_name.update(temp_name_dic)
                    #The next dictionary check is to see if the location's sequence is already present, if codons are conserved then is will happen a lot, it adds to the frequency then moves to the next codon
                    if codons in dyad_dic[i]:
                        dyad_freq[i]+=1
                        second_temp=dyad_ratio_dic[i]                                
                        second_temp.append(codons)
                        #print(second_temp)
                        temp_dic = {i:second_temp}
                        dyad_ratio_dic.update(temp_dic)
                        continue
                    #If the codon is not present, then a new codon is added to entry at the given position
                    else:
                        temp_list=[]
                        second_temp=[]
                        for num in dyad_dic.items():
                            if num[0]==i:
                                dyad_freq[i]+=1
                                temp_list=dyad_dic[i]
                                temp_list.append(codons)
                                temp_dic = {i:temp_list}
                                dyad_dic.update(temp_dic)
                                second_temp=dyad_ratio_dic[i]                                
                                second_temp.append(codons)
                                #print(second_temp)
                                temp_dic = {i:second_temp}
                                dyad_ratio_dic.update(temp_dic)
                                
                #If the entry is new, then the first set of data for the three dictionaries is entered.                
                else:
                    dyad_dic[i]=[codons]
                    dyad_freq[i]=1
                    dyad_species_name[i] = [species_name]
                    dyad_ratio_dic[i]=[codons]
                    
        #This code turns the file_ID into a number, it then uses that number to get teh corresponding line in gene_correspondents, which then gives it the gene code which is recorded. 
        if file_ID !='':
            line_num = file_ID.strip('.fasta')
            line_num = int(line_num.strip('gene_'))

            target_line = content[line_num-1]
            target_line = target_line.split()
            gene_ID = target_line[0]
        gap_removed = data.replace('-', '')
        if len(gap_removed)>25:
            average_length = average_length + len(gap_removed)
        else:
            missing=1
        dyad_frequency = CpG/len(gap_removed)*100
        #adds the data to the CSV        
        file_num = 'Gene_'+str(line_num)
        row_data = [file_num, gene_ID, species_name, len(gap_removed), CpG, dyad_frequency]    
        writer.writerow(row_data)
        #returns the data, most the data is re-entered for the next sequence, but must be returned like this to save it.  
        return(dyad_dic, dyad_freq, dyad_species_name, gene_ID, dyad_ratio_dic, average_length, CpG, missing)

#This is the second section (although technically it is where the code starts)
#Generally speaking it is where the user input the data, and then the general document information is gathered and inputted into a csv
#At the end of this section, you have the general document, and then a series of dictionaries that are then used to generate the individual locations document (section 3)

    
#This gets the list of gene names to Gene IDs
with open('gene_correspondence.txt') as name_file:
    content = name_file.readlines()
directory = 'test_data'
#empty object for later use
gen_rows = {}
gene_length=[]
final_gene_ID=[]
conserved_dic={}
early_dic={}
conserved_loc = {}
length_dic ={}
big_dic={}
small_dic={}
averages_dic={}
animal_dic={}
missing_dic={}
# iterate over files in the test_data directory. Could be easily edited, or changed to be a user input
for filename in os.listdir(directory):       
#Creates the CSV with the desired headers  
    with open('test_results/' + str(filename + '_general_results')+'.csv', 'w') as result:
        headers = ['File Number', 'Gene ID', 'Species', 'Gene Length', 'No. of CpG', 'Dyad Freq. Factor']
        writer = csv.writer(result)
        writer.writerow(headers)
        #empty objects needed for later on. 
        dyad_dic={}
        dyad_ratio_dic={}
        seq_list=[]
        dic_info={}
        temp_location_list={}
        num_of_sequences = 0
        dyad_freq={}
        dyad_species_name={}       
        file = filename
        line_num=0
        average_length =0 
        #This is if the user wants to conduct the further test. 
        test_check = 'y'
        #opens the file
        with open('test_data/'+file) as f:
            count=0
            #For loop that goes over every line in the provided document. 
            for all_lines in f:
                #If the line starts with a '>' then it is added as a new entry (as > indicates that it is a new sequence)
                if all_lines[0] =='>':
                    #This section trims the data of \n, then runs the function which looks for dyads and then imports the data onto the first CSV. 
                    if count!=0:
                        data = data.replace('\n','').upper()
                        fun_res = dyad_counter(data, file_ID, species_name, dyad_dic, dyad_freq, dyad_species_name, dyad_ratio_dic, average_length)
                        #These lines extract the data generated in the function above, the data will be needed for the section below. 
                        dyad_dic = fun_res[0]
                        dyad_freq = fun_res[1]
                        dyad_species_name = fun_res[2]
                        temp_gene_ID = fun_res[3]
                        dyad_ratio_dic = fun_res[4]
                        average_length = fun_res[5]
                        CpG=fun_res[6]
                        missing=fun_res[7]
                        num_of_sequences +=1
                        if species_name in animal_dic.keys():
                            animal_dic[species_name]=animal_dic[species_name]+CpG
                            missing_dic[species_name]=missing_dic[species_name]+missing
                        else:
                            animal_dic[species_name]=CpG
                            
                            missing_dic[species_name]=missing
                        
                    #This is for the first line of the document, as there is no data to be run yet. 
                    elif count==0:
                        count = 1
                    #The next few lines get information about the gene that data is going to be collected for (specifically the Gene ID and species name)

                    file_ID = str(file).replace('>','')
                    species_name = str(all_lines.replace('>',''))
                    species_name = str(species_name.replace('\n',''))
                    
                    #Resets data to nothing for the next set of data and then continues onto the first line of data. 
                    data = ''
                    continue
                #Data is just data with the new line
                data = data+all_lines
            #As the end of the loop does not end on a >, the code for calculating dyads must be repeated once again 
            data = data.replace('\n','').upper()
            fun_res = dyad_counter(data, file_ID, species_name, dyad_dic, dyad_freq, dyad_species_name, dyad_ratio_dic, average_length)
            dyad_dic = fun_res[0]
            dyad_freq = fun_res[1]
            dyad_species_name = fun_res[2]                
            num_of_sequences +=1
            average_length = fun_res[5]
            final_average = average_length/num_of_sequences
            row_data = ['Average length of sequence', '', final_average, '', '']
            writer.writerow(row_data)
            averages_dic[temp_gene_ID] = final_average
        
        
        
        total_dyads=0
        for each in dyad_freq.items():
            total_dyads = each[1]+total_dyads
        temp_final = {filename:total_dyads}
        gen_rows.update(temp_final)
        gene_length.append(len(data))
        final_gene_ID.append(temp_gene_ID)

#section 3        
#If the inputted file has more then one sequence of DNA it will do the section below
#This section looks at each location where there was at least 1 dyad and analyse it
#It looks at the number dyads, the codon sequences, the amino acid sequence, the frequency of each codon, and then attempts to give an array of score that grade how likely the dyad pattern found is a result of chance. 


    #This section compared the dyads and looks for dyads that are highly conserved, and gives the species. It also looks at how likely that dyad sequence is conserved.
    if num_of_sequences>1:
            #Dictionary that can assign codons to amino acids
        codon_to_amino = {'Phe':'TTC',
                          'Leu':'CTC',
                          'Ile':'ATC',
                          'Val':['GTT', 'GTC', 'GTA', 'GTG'],
                          'Ser':'TCC',
                          'Pro':'CCC',
                          'Thr':'ACC',
                          'Ala':['GCT', 'GCC', 'GCA', 'GCG'],
                          'Tyr':'TAC',
                          'His':'CAC',
                          'Asn':'AAC',
                          'Asp':['GAT', 'GAC'],
                          'Glu':['GAA', 'GAG'],
                          'Cys':'TGC',
                          'Arg':'CGC',
                          'Ser':'AGC',
                          'Gly':['GGT', 'GGC', 'GGA', 'GGG']

        }
            #dictionary with likelihood that each amino acid is conserved. 
        front_codon_dic = {'Phe': 24.71948300108424,
                             'Leu': 6.231128355613783,
                             'Ile': 27.856867589620553,
                             'Met': 0.0,
                             'Val': 7.107237039222716,
                             'Ser': 21.207978547877214,
                             'Pro': 17.20568849034703,
                             'Thr': 14.975463866491085,
                             'Ala': 20.088245011013637,
                             'Tyr': 32.991588512980336,
                             'His': 39.75633601626467,
                             'Gln': 0.0,
                             'Asn': 29.174469941162602,
                             'Lys': 0.0,
                             'Asp': 28.986902793614437,
                             'Glu': 0.0,
                             'Cys': 27.218103207070428,
                             'Trp': 0.0,
                             'Arg': 9.957816426216132,
                             'Gly': 18.35457996927749}
        #Check to see if the further is wanted. 
        if test_check == 'Y' or test_check == 'y':
        #This writes a new CSV with the shared dyad data. 
            with open('test_results/' + str(file + '_dyad_results')+'.csv', 'w') as dyad_file:
                #Sorts the core dictionary so it is index order
                dyad_dic = dict(sorted(list(dyad_dic.items())))
                dyad_ratio_dic = dict(sorted(dyad_ratio_dic.items()))
                header = ['Gene ID','Dyad Index', 'Number of appearances', 'Dyad Sequences', 'Amino Acids', 'how frequent each amino acid appears','Ratio on how conserved the gene is.','Difference between expected and observed', 'Random Score', 'Mammal corrected score','Total Number of Sequences',  'Species present']
                writer = csv.writer(dyad_file)
                writer.writerow(header)
                finaltemp=[]
                templist=[]
                count2=-1
                conserved_dic[file]=0
                early_dic[file]=0
                for each in dyad_ratio_dic.items():
                    temp = list(each[1])
                    tempor = Counter(temp)
                    fulltemp=[]
                    for all in tempor.items():
                        fulltemp.append(all[1])
                    finaltemp.append(fulltemp)
                        

                #For loop that goes over each core dictionary entry and extracts the relevant data. It also gives the shared index so that data cane be extracted from the other dictionaries.
                for each_dyad in dyad_dic.items():
                    #redundancy_rate is set to 1 to start with, this is used later when observing how likely the codon pattern should be conserved. 
                    redundancy_rate = 1
                    count2+=1
                    amino_list = []
                    list_of_redundancies=[]
                    #Extracting the information from dyad_dic to get the locations and the associated codon sequences
                    dyad_loc = each_dyad[0]
                    dyad_seqs = each_dyad[1]
                    #This looks into the dyad_freq dictionary and gets the final number of dyads at the associated location. 
                    final_freq = dyad_freq[each_dyad[0]]
                    random_score = 0
                    #Gets the species name in the same way the final frequency is collected. 
                    final_name = dyad_species_name[each_dyad[0]]
                    gene_code = final_gene_ID[-1]
                    if each_dyad[0] in conserved_loc.keys():
                        
                        conserved_loc[each_dyad[0]] = conserved_loc[each_dyad[0]] + final_freq
                    else:
                        conserved_loc[each_dyad[0]] = final_freq
                    if final_freq>200:
                        conserved_dic[file]+=1
                        if dyad_loc<201:
                            early_dic[file]+=1
                        ####Here is new code 
                        #if each_dyad[0] in conserved_loc.keys():
                         #  conserved_loc[each_dyad[0]] = conserved_loc[each_dyad[0]] + 1
                            
                        #else:
                         #   conserved_loc[each_dyad[0]] = 1
                        fix_data = data.strip('-')
                        temp_per = int((each_dyad[0]/len(fix_data))*100)
                        if temp_per in length_dic:
                            length_dic[temp_per]+=1
                        else:
                            length_dic[temp_per]=1
                        if final_average>2000:
                            if temp_per in big_dic:
                                big_dic[temp_per]+=1
                            else: 
                                big_dic[temp_per]=1
                        else:
                            if temp_per in small_dic:
                                small_dic[temp_per]+=1
                            else:
                                small_dic[temp_per]=1
                    #The next for loop looks at all the codon pairs in dyad_seqs made above.
                    for codon_pair in dyad_seqs:
                        #First the codon pair is split up to make the front and back codons
                        front_codon = codon_pair[0:3]
                        back_codon = codon_pair[3:6]
                        #The next two for loops search for the codon sequences in the codon to amino dicitonary above and converts it to the relevant amino acid. 
                        for each_amino in codon_to_amino.items():
                            #Some amino acids have multiple corisponding codons, in these cases a new for loop goes over each item and checks if front codon is in that list. 
                            if len(each_amino[1])>1:
                                for list_amino in each_amino[1]:
                                    if front_codon == list_amino:               
                                        #If the front codon is found in an amino's codon list then front amino is assigned to the corrisonding amino acid. 
                                        front_amino = each_amino[0]
                            #If there is only one associated codon to an amino, then it just needs to do a basic check to see if it is present, then assign accordingly. 
                            if front_codon in each_amino:
                                front_amino = each_amino[0]
                        #This is a repeat of the above but looking at the back codon instead. 
                        for each_amino in codon_to_amino.items():
                            if len(each_amino[1])>1:
                                for list_amino in each_amino[1]:
                                    if back_codon == list_amino:                
                                        back_amino = each_amino[0]
                            if back_codon in each_amino:
                                back_amino = each_amino[0]
                        if [front_amino, back_amino] not in amino_list:
                            amino_list.append([front_amino, back_amino])
                        #A conserved ratio that indicates how conserved the gene is. 
                        conserved_ratio = final_freq/len(amino_list)
                        for every in range(len(finaltemp[count2])):                            
                            list_of_redundancies.append(front_codon_dic[front_amino])
                    #The final set of if statements calculate the likelihood that all the CpGs are conserved by chance. It also calculates a final score which is an arbitrar score but does indicate if a location has a hgih chance ot be interesting. 
                    if final_freq== 1:
                        redundancy_rate = 'NA'
                        random_score = 0
                        mammal_score=0 
                    else:
                        redundancy_rate = float((sum(list_of_redundancies)/100)/len(list_of_redundancies))
                        random_score = final_freq-(num_of_sequences*redundancy_rate)
                        mammal_score = final_freq-(num_of_sequences*(float(redundancy_rate/10)))
                    codon_freq = list(dyad_ratio_dic.items())        
                    rows = [gene_code, dyad_loc, final_freq, dyad_seqs,amino_list, finaltemp[count2], conserved_ratio, redundancy_rate, random_score, mammal_score, num_of_sequences, final_name]
                    writer.writerow(rows)
                

#Section 4 
# The section generates the general documents 
# The first is a general document that has all genes, their respective number of CpGs, length of the gene, average number of CpGs per species, and a ratio value (number of CpGs/length of the gene)
# The second document shows each dyad location, and how many conserved dyads are present (as in number of genes with a conserved dyad at that location)
# The third document uses normalised data. Dyad location is transfomred into percentage along the gene length, frequency at each percentage is recorded, 
# and frequency is also recorded for small genes (<2000bp) and large genes (>2000) at each percentage point. 

#First is the general summary file                    
counter=0
#Initially, a file is generated
with open('General_summary_update.csv', 'w') as gen:
    headers = ['File name', 'Gene ID', 'Number of CpGs', 'Average number of CpGs', 'Length of Gene', 'Relative dyad concentration', 'How many highly conserved Ratio', 'Conservation against general presence','Number of conserved dyads','Number of early (<200) conserved dyads']
    writer = csv.writer(gen)
    writer.writerow(headers)
    
    #For loop for each gene, with associated number of dyads. 
    for each_gene in gen_rows.items():
        # ratio = relative dyad concentration
        ratio = each_gene[1]/gene_length[counter]
        num_conserved_ratio = (conserved_dic[each_gene[0]]/gene_length[counter])*100
        average_gene_length = averages_dic[final_gene_ID[counter]]
        if num_conserved_ratio == 0:
            ratio_ratio=0
        else:
            #Ratio_ratio = relative dyad concentration divided by the relative number of conserved dyads
            # essnetially, there may be many dyads or not many conserved, or there are few dyads but they are conserved. 
            ratio_ratio = ratio/num_conserved_ratio
        #This is a list of data that forms one line, which is then written into the document made. 
        fin = [each_gene[0], final_gene_ID[counter], each_gene[1], (each_gene[1]/261), average_gene_length, ratio, num_conserved_ratio, ratio_ratio,conserved_dic[each_gene[0]], early_dic[each_gene[0]]]
        writer.writerow(fin)
        counter+=1

#The document detailed locations of conserved dyads is generated
with open('conserved_locations.csv', 'w') as conserve:
    headers = ['location', 'no.of conserved dyads']
    writer = csv.writer(conserve)
    writer.writerow(headers)
    #The conserved locations dictionaries first needs to be sorted into numerical order
    conserved_loc = dict(sorted(conserved_loc.items()))
    #for loop going over each conserved location and writing a line into the document show the index (ind) and the total number of conserved dyads (tot)
    for each_conserved_loc in conserved_loc.items():
        ind = each_conserved_loc[0]
        tot = each_conserved_loc[1]
        fin = [ind, tot]
        writer.writerow(fin)

#Final document has data on the positions of conserved dyads along a normalised percentage gradient (percentage long the gene length)
#It also looks if it is a large or small gene and assigns it accordingly (allowing comparisons between the two sizes) 
with open('percentage_locations.csv', 'w') as per_con:
    headers = ['Percent', 'Freq', 'big', 'small']
    writer = csv.writer(per_con)
    writer.writerow(headers)
    #Dictionaries are sorted into percentage order, from 0-99
    length_dic = dict(sorted(length_dic.items()))
    big_dic = dict(sorted(big_dic.items()))
    small_dic = dict(sorted(small_dic.items()))
    #for loop that looks at each percentage point and records the number of dyads at each point
    for each_length_percentage in length_dic.items():
        ind = each_length_percentage[0]
        tot = each_length_percentage[1]
        #Errors can be thrown if there are 0 dyads for a certain size, so an If statement is needed to catch said errors
        if each_length_percentage[0] not in big_dic.keys():
            big_val=0
        else:
            big_val = big_dic[each_length_percentage[0]]
        if each_length_percentage[0] not in small_dic.keys():
            small_val=0
        else:
            small_val = small_dic[each_length_percentage[0]]
        fin = [ind, tot, big_val, small_val]
        writer.writerow(fin)
        
with open('animal__frequency.csv', 'w') as ani_freq:
    headers = ['species', 'Number of CpGs', 'Number of missing genes']
    writer= csv.writer(ani_freq)
    writer.writerow(headers)
    for each in animal_dic.items():
        animal = each[0]
        freq = each[1]
        miss = missing_dic[each[0]]
        fin = [animal, freq, miss]
        writer.writerow(fin)
