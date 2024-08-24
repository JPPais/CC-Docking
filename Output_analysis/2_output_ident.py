import math
import os
import glob

# Define the base path
base_path = '/path/*/'

# Get a list of directories matching the pattern
analysis_folders = glob.glob(os.path.join(base_path, 'folder_*'))

# Loop through each directory
for input_path in analysis_folders:
    # Input files
    smiles_file = 'compounds.txt'

    dist_file = input_path + '/dist_results.txt'
    center_file = input_path + '/withoutNO2.txt'

    # Output files
    out_file = input_path + '/output.txt'

    # Initialize dictionaries to store data from File A and File C
    distance_data = {}
    aromatic_data = {}
    # Read data from File A and store it in the distance_data dictionary
    with open(dist_file, 'r') as file_a:
        next(file_a)
        for line in file_a:
            columns = line.strip().split()
            identifier_a = columns[0]
            distance_a = columns[2]
            distance_b = columns[3]
            score = columns[4]
            n_number = columns[1]
            distance_c = columns[5]
            distance_d = columns[6]
            name = identifier_a + "_" + n_number + "."
            distance_data[name] = (identifier_a, distance_a, distance_b, score, distance_c, distance_d)

    # Read data from File C and update the distance_data dictionary
    with open(center_file, 'r') as file_c:
        next(file_c)
        for line in file_c:
            columns = line.strip().split()
            identifier_c = columns[0]
            aromatic_score = columns[1]
            aromatic_dist_FAD = columns[11]
            aromatic_dist_S = columns[12]
    #        name_center = identifier_c.split('_')[0]
            aromatic_data[identifier_c] = (identifier_c, aromatic_score, aromatic_dist_FAD, aromatic_dist_S)
            
    # Open the output file for writing
    with open(out_file, 'w') as output_file:
        output_file.write("File_Name\tIdentifier\tCompound_Name\tNNumber\tMIC\tLog(MIC)\tDistance_to_S\tDistance_to_FAD\tDCenter_to_S\tDCenter_to_FAD\tScore\n")

        # Read data from File B and merge it with data from File A and File C
        with open(smiles_file, 'r') as file_b:
            next(file_b)
            for line in file_b:
                columns = line.strip().split()
                identifier = columns[0]
                chemical_name = columns[1]
                mic = columns[3]
                logmic = math.log10(float(mic))

                # Search for matching names in distance_data
                matching_names = [name for name in distance_data if f"{identifier}_" in name]

                for name in matching_names:
                    identifier_a, distance_a, distance_b, score, distance_c, distance_d = distance_data[name]

                    # Extract the N"number" from the name
                    n_number = name.split('_')[-1].split('.')[0]

                    output_file.write(f"{identifier_a}\t{identifier}\t{chemical_name}\t{n_number}\t{mic}\t{logmic:.5f}\t{distance_a}\t{distance_b}\t{distance_d}\t{distance_c}\t{score}\n")

                # Search for matching identifiers in aromatic_data
                matching_aromatic = [aromatic_key for aromatic_key in aromatic_data if f"{identifier}_" in aromatic_key]

                for aromatic_key in matching_aromatic:
                    identifier_c, aromatic_score, aromatic_dist_FAD, aromatic_dist_S = aromatic_data[aromatic_key]
                    output_file.write(f"{identifier_c}\t{identifier}\t{chemical_name}\t0\t{mic}\t{logmic:.5f}\t0\t0\t{aromatic_dist_S}\t{aromatic_dist_FAD}\t{aromatic_score}\n")

    print("Merging and writing to output.txt complete.")
