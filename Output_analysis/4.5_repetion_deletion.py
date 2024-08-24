#It may be necessary to delete repeated values, depending on certain aspects of the data you are treating.
#As such, this script can be skipped if the prepared data does not contain repeated values

import os

# Define the input folder containing .txt files and the output folder
input_folder = "/path/results/select"
output_folder = "/path/results/select" #It will overwrite the files in the folder unless you provide an alternative path

# Create the output folder if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Get a list of all .txt files in the input folder
txt_files = [f for f in os.listdir(input_folder) if f.endswith('.txt')]

# Iterate through each .txt file
for txt_file in txt_files:
    # Define the input and output file paths for the current file
    input_file = os.path.join(input_folder, txt_file)
    output_file = os.path.join(output_folder, txt_file)

    # Open the input file for reading
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # Initialize an empty list to store modified lines
    modified_lines = []

    # Iterate through the lines of the file
    for i in range(len(lines) - 1):
        current_line = lines[i].strip().split('\t')
        next_line = lines[i + 1].strip().split('\t')
        
        # Check if the text in the column of the current line is equal to the text in the column of the next line
        if current_line[1] != next_line[1]:
            modified_lines.append(lines[i])

    # Add the last line to the modified lines list since it won't be compared with another line
    modified_lines.append(lines[-1])

    # Write the modified data to the output file
    with open(output_file, 'w') as f:
        for line in modified_lines:
            f.write(line)
    columns_to_check = [11, 12]
    count_values_above_cutoff(input_file, output_rep_filename, columns_to_check, cutoff_value)

