import os
import shutil

source_dir = '/path/'
dest_dir = '/path/results/'
target_files = ['best_distNO.txt', 'best_distcenter.txt', 'best_score.txt', 'best_score_distNO.txt', 'best_score_distcenter.txt', 'output.txt']  # Add the specific file names you want to copy and rename

def result_rename(source_directory, destination_directory, target_files):
    if not os.path.exists(destination_directory):
        os.makedirs(destination_directory)  # Create the destination directory if it doesn't exist

    for root, dirs, files in os.walk(source_directory):
        for file in files:
            if file in target_files:
                original_file_path = os.path.join(root, file)
                relative_path = os.path.relpath(original_file_path, source_directory)
                renamed_path = os.path.join(destination_directory, relative_path.replace(os.path.sep, '_'))
                shutil.copy2(original_file_path, renamed_path)

    print("Files copied and renamed successfully.")

result_rename(source_dir, dest_dir, target_files)

def process_file(input_file, output_folder):
    # Create the output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Create the output file path
    output_file = os.path.join(output_folder, os.path.basename(input_file))

    # Open the input and output files
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Iterate through each line in the input file
        for line in infile:
            # Split the line into columns
            columns = line.strip().split()
#            outfile.write(line)
            # Check if the second column contains the desired name (adjust the index if columns are 0-indexed)
            if len(columns) > 1 and columns[1] in ['Identifier', 'A1', 'A2' ]: #Replace A1, A2 by a list of your molecular identifiers to select specific compounds or compound lists to further analysis
                # Write the entire line to the output file
                outfile.write(line)

def process_folder(input_folder, output_folder, search_text):
    # Get a list of all files in the input folder
    all_files = os.listdir(input_folder)

    # Filter files based on the search text
    filtered_files = [f for f in all_files if search_text in f]

    # Process each file in the filtered list
    for file_name in filtered_files:
        # Check if the file is a .txt file
        if file_name.endswith('.txt'):
            # Get the full path of the input file
            input_file = os.path.join(input_folder, file_name)
            # Process the file and copy relevant lines to the output folder
            process_file(input_file, output_folder)

            input_file = os.path.join(output_folder, file_name)
            sumcalc(input_file, output_folder_calc)


def sumcalc(input_file, output_folder_calc):
    
        # Create the output folder if it doesn't exist
    if not os.path.exists(output_folder_calc):
        os.makedirs(output_folder_calc)

    # Create the output file path
    output_file = os.path.join(output_folder_calc, os.path.basename(input_file))

    # Open the input file for reading 
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # Define the header for the new columns
    header = lines[0].strip() + "sumNO\tsumCenter\t" 

    # Initialize an empty list to store modified lines
    modified_lines = []

    for line in lines[1:]:  # Skip the header line
        # Split the line into columns
        columns = line.strip().split('\t')
        
        # Calculate the sum of columns 7 and 8
        sum_no = float(columns[6]) + float(columns[7])
        
        # Calculate the sum of columns 9 and 10
        sum_center = float(columns[8]) + float(columns[9])
        
        # Add the sums as new columns to the line
        modified_line = '\t'.join(columns + [f'{sum_no:.5f}', f'{sum_center:.5f}'])
        
        # Append the modified line to the list
        modified_lines.append(modified_line)

    # Write the modified data to the output file
    with open(output_file, 'w') as f:
        # Write the header to the output file
        f.write(header + "\n")
        
        # Write the modified lines to the output file
        for line in modified_lines:
            f.write(line + "\n")

input_folder = '/path/results/'
output_folder = '/path/results/select'
output_folder_calc = '/path/results/select'
search_text1 = "distNO"  # Change this any text you want to search for selecting the files for data extraction.
#search_text2 = "distcenter"

# Call the function to process the entire folder
process_folder(input_folder, output_folder, search_text1)
#process_folder(input_folder, output_folder, search_text2)
