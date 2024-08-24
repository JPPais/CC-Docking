import pandas as pd
import os
import numpy as np

input_folder = "/home/afortuna/DprE1_Docking/VINA/Results_All/select"  # Replace with the path to your input folder
output_folder = "/home/afortuna/DprE1_Docking/VINA/Results_All/Results_analysis_select"  # Replace with the path to your output folder

# Cutoff values
x_cutoff = -0.35 #define your LogMIC cutoff values for binary classification
y_cutoff_range = np.arange(2, 15.2, 0.2)  # Range from 2 to 15 with a step of 0.2
y_cutoff_range_sum = np.arange(5, 20.2, 0.2)  # Range from 2 to 15 with a step of 0.2
score_cutoff_range = np.arange(-17, 0, 0.2)

# Ensure the output folder exists, create it if necessary
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

def assign_quadrant(row, x_cutoff, y_cutoff):
    try:
        x = float(row['x'])
        y = float(row['y'])
    except ValueError:
        return 'Invalid'

    if x >= x_cutoff and y >= y_cutoff:
        return '+x+y'
    elif x >= x_cutoff and y < y_cutoff:
        return '+x-y'
    elif x < x_cutoff and y >= y_cutoff:
        return '-x+y'
    elif x < x_cutoff and y < y_cutoff:
        return '-x-y'
    else:
        return 'Invalid'

# Function to process the input file
def process_input_file(input_file, y_col, y_cutoff_range):
    df = pd.read_csv(input_file, header=None, usecols=[5, y_col], delimiter='\t')  # Assuming values are tab-separated
    
    df.columns = ['x', 'y']

    output_file = os.path.join(output_folder, os.path.basename(input_file).replace('.txt', '_analysis.csv'))

    # Iterate over the range of Y cutoff values
    for y_cutoff in y_cutoff_range:
        # Apply the function to create a new column 'quadrant'
        df['quadrant'] = df.apply(lambda row: assign_quadrant(row, x_cutoff, y_cutoff), axis=1)

        # Count the number of inhibitors in each quadrant
        quadrant_counts = df['quadrant'].value_counts()

        # Example usage:
        TP = quadrant_counts.get('-x-y', 0)
        TN = quadrant_counts.get('+x+y', 0)
        FP = quadrant_counts.get('+x-y', 0)
        FN = quadrant_counts.get('-x+y', 0)

        if (TP + TN + FP + FN) == 0:
            accuracy = np.nan
        else:
            accuracy = (TP + TN) / (TP + TN + FP + FN)

        if (TP + FN) == 0:
            sensitivity = np.nan
        else:
            sensitivity = TP / (TP + FN) #sometimes called or recall

        if (TP + FP) == 0:
            precision = np.nan
        else:
            precision = TP / (TP + FP)

        if (TP + FP) == 0:
            specificity = np.nan
        else:
            specificity = TN / (TN + FP)
        
        if np.isnan(precision) or np.isnan(sensitivity) or (precision + sensitivity) == 0:
            f1_score = np.nan
        else:
            f1_score = 2 * (precision * sensitivity) / (precision + sensitivity)

        if ((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)) == 0:
            mcc = np.nan
        else:
            mcc = (TP * TN - FP * FN) / ((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)) ** 0.5 #Matthews Correlation Coefficient (MCC)

        # Write data to a CSV file
        with open(output_file, 'a') as file:
            # Check if the header has been written
            if os.stat(output_file).st_size == 0:
                # Write header
                file.write(f'File Name,Cutoff X,Cutoff Y,Total Inhibitors,True Negative,False Negative,False Positive,True Positive,Specificity,Sensitivity,Accuracy,Precision,F1_Score,MCC\n')
        
            # Write data
            file.write(f'{os.path.basename(input_file)},{x_cutoff},{round(y_cutoff, 1)},{len(df)},{TN},{FN},{FP},{TP},{specificity},{sensitivity},{accuracy},{precision},{f1_score},{mcc}\n')

# List of input files
input_files = [os.path.join(input_folder, file) for file in os.listdir(input_folder) if file.endswith(".txt")]

# Process each input file
for input_file in input_files:
    process_input_file(input_file, 6, y_cutoff_range)
    process_input_file(input_file, 7, y_cutoff_range)
    process_input_file(input_file, 8, y_cutoff_range)
    process_input_file(input_file, 9, y_cutoff_range)
    process_input_file(input_file, 10, score_cutoff_range)
    process_input_file(input_file, 11, y_cutoff_range_sum)
    process_input_file(input_file, 12, y_cutoff_range_sum)

# Function to find the lowest and highest values and write to a new file
def find_and_write_extremes(file_list):
    # Create a new file for the results
    result_file = "Results_analysis_select.csv"

    with open(result_file, 'w') as result_file:
        result_file.write('File Name,Cutoff X,Cutoff Y,Total Inhibitors,True Negative,False Negative,False Positive,True Positive,Specificity,Sensitivity,Accuracy,Precision,F1_Score,MCC,Distance Type\n')

        # Process each output file
        for output_file in file_list:
            # Read the output file into a DataFrame
            df = pd.read_csv(output_file)

            # Find all rows with the minimum value in the 'Total_False' column
            max_lines = df[df['MCC'] == df['MCC'].max()]

            # Iterate over each row with the minimum value
            for _, max_line in max_lines.iterrows():
                # Determine the distance type based on the rules provided
                if 8 <= max_line.name <= 37:
                    distance_type = "NDist_S"
                elif 76 <= max_line.name <= 108:
                    distance_type = "NDist_FAD"
                elif 138 <= max_line.name <= 175:
                    distance_type = "CDist_S"
                elif 206 <= max_line.name <= 231:
                    distance_type = "CDist_FAD"
                elif 267 <= max_line.name <= 339:
                    distance_type = "Score"
                elif 350 <= max_line.name <= 420:
                    distance_type = "Sum_dist_NO"
                elif 430 <= max_line.name <= 500:
                    distance_type = "Sum_dist_Center"

                # Write the data to the result file
                result_file.write(','.join(map(str, max_line)))
                result_file.write(f',{distance_type}\n')

find_and_write_extremes([os.path.join(output_folder, file.replace('.txt', '_analysis.csv')) for file in os.listdir(input_folder) if file.endswith(".txt")])
