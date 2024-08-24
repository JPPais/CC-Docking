import os 
import math
import numpy as np
import glob

# Define the base path
base_path = '/path/*/'

# Get a list of directories matching the pattern
analysis_folders = glob.glob(os.path.join(base_path, 'folder_*'))

# Loop through each directory
for folder_path in analysis_folders:

    # File path for the "receptor" file
    receptor_file_path = '/path/receptor.pdbqt'

    coord_file = folder_path + '/results_coordinates.txt'
    dista_file = folder_path + '/dist_results.txt'

    out_file = folder_path + '/aromatic_coordinates.txt'
    new_output_file = folder_path + '/withoutNO2.txt' 

    # Function to parse a .pdbqt file and extract atom coordinates
    def parse_pdbqt(file_path):
        atom_coordinates = {}
        n_count = 1  # Counter for unique names for matching N atoms

        with open(file_path, 'r') as pdbqt_file:
            for line in pdbqt_file:

                if line.startswith("ATOM"):  # Change the marker to "ATOM" for .pdbqt format
                    atom_name = line[11:16].strip()  # Adjust the column positions for atom_name
                    atom_charge = line[69:76].strip()  # Adjust the column positions for atom_type
                    x, y, z = float(line[32:39]), float(line[39:47]), float(line[47:54])  # Adjust the column positions for coordinates
                    if atom_name == 'N' and (atom_charge == '0.062' or atom_charge == '0.067'):
                        Lig_N = f'N{n_count}'  # Generate a unique name
                        atom_coordinates[Lig_N] = (x, y, z)
                        n_count += 1

        return atom_coordinates

    # Function to parse a .pdbqt file and extract score
    def get_scores(file_path):
        score = None
        with open(file_path, 'r') as pdbqt_file:
            for line in pdbqt_file:
                if line.startswith("REMARK VINA RESULT:"):
                    fields = line.split()
                    if len(fields) >= 4:
                        score = fields[3]
                        break
        return score

    # Function to retrieve "S" and "FAD" atom coordinates from the receptor file
    def retrieve_s_and_fad_coordinates(receptor_file):
        s_coordinates = {}
        fad_coordinates = {}
        with open(receptor_file, 'r') as receptor:
            for line in receptor:
                line = line.strip()
                if line.startswith("ATOM"):
                    atom_name = line[11:16].strip()  # Adjust the column positions
                    res_num = line[23:26].strip()  # Adjust the column positions
                    res_type = line[17:21].strip()  # Adjust the column positions
                    if atom_name == 'SG' and res_num == '382':
                        x, y, z = float(line[32:39]), float(line[39:47]), float(line[47:54])
                        s_coordinates['S'] = (x, y, z)
                    elif atom_name == 'N5' and res_type == 'FAD':
                        x, y, z = float(line[32:39]), float(line[39:47]), float(line[47:54])
                        fad_coordinates['FAD'] = (x, y, z)
        
        return s_coordinates, fad_coordinates

    # Function to calculate the Euclidean distance between two 3D points (atoms)
    def calculate_distance(atom1, atom2):
        x1, y1, z1 = atom1
        x2, y2, z2 = atom2
        return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)


    def extract_aromatic_atoms(line):
        atom_type = line[77:79].strip()
        if atom_type == 'A':
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            return x, y, z
        return None

    def process_qt_file(pdbqt_file, output_file):
        aromatic_atoms = []

        with open(pdbqt_file, 'r') as file:
            for line in file:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    atom_coords = extract_aromatic_atoms(line)
                    if atom_coords and len(aromatic_atoms) < 6:
                        aromatic_atoms.append(atom_coords)

        if len(aromatic_atoms) == 6:
            center = tuple(np.mean(aromatic_atoms, axis=0))

            with open(output_file, 'a') as result_file:
                result_file.write(f"{os.path.splitext(os.path.basename(pdbqt_file))[0]}\t")
                score = get_scores(pdbqt_file)
                result_file.write(f"{score}\t")
                result_file.write('\t'.join(map(str, np.concatenate(aromatic_atoms).tolist() + list(center))) + '\n')

    # List all .pdbqt files in the folder
    pdbqt_files = [f for f in os.listdir(folder_path) if f.endswith('.pdbqt')]

    # Parse each .pdbqt file and extract the coordinates of matching N atoms
    all_atom_coordinates = {}

    for file_name in pdbqt_files:
        file_path = os.path.join(folder_path, file_name)
        atom_coordinates = parse_pdbqt(file_path)
        all_atom_coordinates[file_name] = atom_coordinates

    all_scores = {}

    for file_name in pdbqt_files:
        file_path = os.path.join(folder_path, file_name)
        scores = get_scores(file_path)
        all_scores[file_name] = scores

    # Retrieve "S" and "FAD" coordinates from the receptor file
    s_coordinates, fad_coordinates = retrieve_s_and_fad_coordinates(receptor_file_path)

    # Write the results to a file named "results.txt"
    with open(coord_file, 'w') as results_file:
        # Write the coordinates of N atoms from all files
        for file_name, atom_coords in all_atom_coordinates.items():
            results_file.write(f"File: {file_name}\n")
            for atom_name, coordinates in atom_coords.items():
                results_file.write(f"  Atom {atom_name} Coordinates: {coordinates}\n")

        # Write the coordinates for "S" and "FAD" atoms if found
        if s_coordinates:
            results_file.write(f"Coordinates for 'S' atom: {s_coordinates}\n")
        else:
            results_file.write("No 'S' atom found with the specified criteria.\n")

        if fad_coordinates:
            results_file.write(f"Coordinates for 'FAD' atom: {fad_coordinates}\n")
        else:
            results_file.write("No 'FAD' atom found with the specified criteria.\n")

    print("Results have been written to 'results.txt'.")

    # Calculate distances and write results to "dist_results.txt"
    with open(dista_file, 'w') as dist_results_file:
        dist_results_file.write("Filename\tN_Atom_Name\tDistance_to_S\tDistance_to_FAD\tScore\tCenter_to_FAD\tCenter_to_S\n")

        # Calculate distances between "S" and "FAD" and each N atom in the folder
        for file_name in pdbqt_files:
            file_path = os.path.join(folder_path, file_name)
            atom_coordinates = parse_pdbqt(file_path)
            for atom_name, coordinates in atom_coordinates.items():
                if atom_name.startswith('N'):
                    distance_to_S = calculate_distance(s_coordinates['S'], coordinates)
                    distance_to_FAD = calculate_distance(fad_coordinates['FAD'], coordinates)

                    filename = os.path.splitext(file_name)[0]
                    score = get_scores(file_path)  # Extract the score value
                    result_line = f"{filename}\t{atom_name}\t{distance_to_S:.5f}\t{distance_to_FAD:.5f}\t{score}\t"

                    # Calculate distances between the center and FAD/S
                    center = np.mean(list(atom_coordinates.values()), axis=0)
                    center_to_FAD = calculate_distance(center, fad_coordinates['FAD'])
                    center_to_S = calculate_distance(center, s_coordinates['S'])

                    result_line += f"{center_to_FAD:.5f}\t{center_to_S:.5f}\n"
                    dist_results_file.write(result_line)

    print("Results written to 'dist_results.txt'.")

    with open(out_file, 'w') as file:
        file.write("File\tScore\t" + "\t".join([f"Atom{i}X\tAtom{i}Y\tAtom{i}Z" for i in range(1, 7)]) + "\t")
        file.write("CenterX\tCenterY\tCenterZ\n")

    pdbqt_files = [f for f in os.listdir(folder_path) if f.endswith(".pdbqt")]

    for pdbqt_file in pdbqt_files:
        process_qt_file(os.path.join(folder_path, pdbqt_file), out_file)

    # Load compound names from "aromatic_coordinates" and "dist_results"
    with open(out_file, 'r') as aromatic_file:
        aromatic_compounds = set(line.split('\t')[0] for line in aromatic_file.readlines()[1:])

    with open(dista_file, 'r') as dist_results_file:
        dist_results_compounds = set(line.split('\t')[0] for line in dist_results_file.readlines()[1:])

    # Find compounds that exist in "aromatic_coordinates" but not in "dist_results"
    missing_compounds = aromatic_compounds - dist_results_compounds

    # Write the missing compounds and their distances to the new output file
    with open(new_output_file, 'w') as new_output:
        new_output.write("Compound\tScore\tCenterX\tCenterY\tCenterZ\tFAD_X\tFAD_Y\tFAD_Z\tS_X\tS_Y\tS_Z\tCenter_to_FAD\tCenter_to_S\n")

        for compound in missing_compounds:
            # Retrieve relevant information for the missing compound from "aromatic_coordinates"
            with open(out_file, 'r') as aromatic_file:
                for line in aromatic_file:
                    if line.startswith(compound):
                        data = line.split('\t')
                        score = data[1]
                        # Extract coordinates from the last three columns, considering potential variations
                        center_coords = float(data[20]), float(data[21]), float(data[22])

                        # Coordinates for FAD and S  - Note: I need to modify the sign of X coordinate, not sure why
                        fad_coords = [-float(fad_coordinates['FAD'][0]), float(fad_coordinates['FAD'][1]), float(fad_coordinates['FAD'][2])]
                        s_coords = [-float(s_coordinates['S'][0]), float(s_coordinates['S'][1]), float(s_coordinates['S'][2])]

                        # Calculate distances between the center and FAD/S
                        center_to_FAD = calculate_distance(center_coords, fad_coords)
                        center_to_S = calculate_distance(center_coords, s_coords)

                        # Write the results to the new output file
                        new_output.write(f"{compound}\t{score}\t{center_coords[0]:.5f}\t{center_coords[1]:.5f}\t{center_coords[2]:.5f}\t"
                                        f"{fad_coords[0]:.5f}\t{fad_coords[1]:.5f}\t{fad_coords[2]:.5f}\t"
                                        f"{s_coords[0]:.5f}\t{s_coords[1]:.5f}\t{s_coords[2]:.5f}\t"
                                        f"{center_to_FAD:.5f}\t{center_to_S:.5f}\n")

    print("New output file written to", new_output_file)
