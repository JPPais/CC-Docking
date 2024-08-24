import os
import glob

# Define the base path
base_path = '/home/afortuna/DprE1_Docking/VINA/*/'

# Get a list of directories matching the pattern
analysis_folders = glob.glob(os.path.join(base_path, '02_analysis*'))

##Defining path
#input_path = '/home/afortuna/DprE1_Docking/VINA/03_docking_vina/02_analysis_4FDN_22'

# Loop through each directory
for input_path in analysis_folders:
    # Defining files:
    input_file1 = input_path + '/output.txt'
    best_file1 = input_path + '/best_distNO.txt'
    best_file6 = input_path + '/best_distcenter.txt'
    best_file2 = input_path + '/best_score_distNO.txt'
    best_file7 = input_path + '/best_score_distcenter.txt'
    best_file3 = input_path + '/best_score.txt'
    best_file4 = input_path + '/short_S.txt'
    best_file5 = input_path + '/short_FAD.txt'

    def best_dist_calc(input_file, best_file):
        # Initialize dictionaries to store data
        shortest_combined_distancesNO = {}
        tableNO = []

        # Read data from the input file and populate dictionaries
        with open(input_file, 'r') as file:
            next(file)  # Skip the header line
            for line in file:
                columns = line.split()
                if columns[3].isdigit() == False:
                    identifier_A, identifier_B, distance_SNO, distance_FADNO = columns[0], columns[1], float(columns[6]), float(columns[7])

                    # Calculate combined distance
                    combined_distanceNO = distance_SNO + distance_FADNO

                    # Check and update shortest distances
                    if identifier_B not in shortest_combined_distancesNO or combined_distanceNO < shortest_combined_distancesNO[identifier_B]:
                        shortest_combined_distancesNO[identifier_B] = combined_distanceNO
                    # Append the line data to the table
                    tableNO.append(columns)
                
        # Write the results to an output file
        with open(best_file, 'w') as best_dist_file:
            best_dist_file.write("File_Name\tIdentifier\tCompound_Name\tNNumber\tMIC\tLog(MIC)\tDistNO_to_S\tDistNO_to_FAD\tDistC_to_S\tDistC_to_FAD\tScore\n")

            # Read data from the input file and find matching entries
            with open(input_file, 'r') as out:
                next(out)  # Skip the header line
                for line in out:
                    columns = line.split()
                    File_Name, Identifier = columns[0], columns[1]

                    # Find the matching entries and write them to the output file
                    if Identifier in shortest_combined_distancesNO and shortest_combined_distancesNO[Identifier] == float(columns[6]) + float(columns[7]):
                        best_dist_file.write('\t'.join(columns) + '\n')

    def best_dist_center(input_file, best_file):
        # Initialize dictionaries to store data
        shortest_combined_distancesC = {}
        tableC = []

        # Read data from the input file and populate dictionaries
        with open(input_file, 'r') as file:
            next(file)  # Skip the header line
            for line in file:
                columns = line.split()
                identifier_A, identifier_B, distance_SC, distance_FADC = columns[0], columns[1], float(columns[8]), float(columns[9])

                # Calculate combined distance
                combined_distanceC = distance_SC + distance_FADC

                # Check and update shortest distances     
                if identifier_B not in shortest_combined_distancesC or combined_distanceC < shortest_combined_distancesC[identifier_B]:
                    shortest_combined_distancesC[identifier_B] = combined_distanceC
                tableC.append(columns)

        # Write the results to an output file
        with open(best_file, 'w') as best_dist_file:
            best_dist_file.write("File_Name\tIdentifier\tCompound_Name\tNNumber\tMIC\tLog(MIC)\tDistNO_to_S\tDistNO_to_FAD\tDistC_to_S\tDistC_to_FAD\tScore\n")

            # Read data from the input file and find matching entries
            with open(input_file, 'r') as out:
                next(out)  # Skip the header line
                for line in out:
                    columns = line.split()
                    File_Name, Identifier = columns[0], columns[1]

                    # Find the matching entries and write them to the output file
                    if Identifier in shortest_combined_distancesC and shortest_combined_distancesC[Identifier] == float(columns[8]) + float(columns[9]):
                        best_dist_file.write('\t'.join(columns) + '\n')



    #First usage, output treatment
    results1 = best_dist_calc(input_file1, best_file1)
    resultscenter1 = best_dist_center(input_file1, best_file6)

    # Create a dictionary to store the lowest scores for each identifier
    lowest_scores = {}

    # Open the input file for reading
    with open(input_file1, "r") as file1:
        # Read and skip the header line
        header = file1.readline()
        
        # Loop through the remaining lines
        for line in file1:
            # Split the line into columns using tab as the delimiter
            columns = line.strip().split('\t')
            
            # Extract the identifier and score values
            identifier = columns[1]
            score = float(columns[10])
            
            # Check if this identifier is already in the dictionary
            if identifier in lowest_scores:
                # Compare the score with the lowest score for this identifier
                if score < lowest_scores[identifier][0]:
                    # Update the lowest score and store the line
                    lowest_scores[identifier] = (score, [line])
                elif score == lowest_scores[identifier][0]:
                    # If the score is the same, add the line to the list
                    lowest_scores[identifier][1].append(line)
            else:
                # If this is the first entry for the identifier, initialize the lowest score
                lowest_scores[identifier] = (score, [line])

        # Open a new file for writing the results
        with open(best_file3, "w") as output_file:
            # Write the header to the output file
            output_file.write(header)
            
            # Write the lines with the lowest scores for each identifier
            for identifier, (_, lines) in lowest_scores.items():
                output_file.write("".join(lines))

    results2 = best_dist_calc(best_file3, best_file2)
    resultscenter2 = best_dist_center(best_file3, best_file7)

    shortest_dist_S = {}
    shortest_dist_FAD = {}

    # Open the input file for reading
    with open(input_file1, "r") as file1:
        # Read and skip the header line
        header = file1.readline()
        
        # Loop through the remaining lines
        for line in file1:
            # Split the line into columns using tab as the delimiter
            columns = line.strip().split('\t')
            
            # Extract the identifier and score values
            identifier = columns[1]
            dist_S = float(columns[5])
            dist_FAD = float(columns[6])
            
            # Check if this identifier is already in the dictionary
            if identifier in shortest_dist_S:
                
                if dist_S < shortest_dist_S[identifier][0]:
                    
                    shortest_dist_S[identifier] = (dist_S, [line])
                elif dist_S == shortest_dist_S[identifier][0]:
                    
                    shortest_dist_S[identifier][1].append(line)
            else:
                # If this is the first entry for the identifier, initialize the lowest dist_S
                shortest_dist_S[identifier] = (dist_S, [line])

            # Check if this identifier is already in the dictionary
            if identifier in shortest_dist_FAD:
                
                if dist_FAD < shortest_dist_FAD[identifier][0]:
                    
                    shortest_dist_FAD[identifier] = (dist_FAD, [line])
                elif dist_FAD == shortest_dist_FAD[identifier][0]:
                    
                    shortest_dist_FAD[identifier][1].append(line)
            else:
                # If this is the first entry for the identifier, initialize the lowest dist_S
                shortest_dist_FAD[identifier] = (dist_FAD, [line])            

    # Open a new file for writing the results
    with open(best_file4, "w") as output_file:
        # Write the header to the output file
        output_file.write(header)
        
        # Write the lines with the lowest scores for each identifier
        for identifier, (_, lines) in shortest_dist_S.items():
            output_file.write("".join(lines))

    # Open a new file for writing the results
    with open(best_file5, "w") as output_file:
        # Write the header to the output file
        output_file.write(header)
        
        # Write the lines with the lowest scores for each identifier
        for identifier, (_, lines) in shortest_dist_FAD.items():
            output_file.write("".join(lines))

    print("Shortest combined distances and best score files writen")
