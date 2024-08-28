In this folder we have two scripts to run our Machine Learning model
We are using unpublished results as input and therefore the input data is not present in this folder, however you can addapt to your work.

first you have to run: 
**data_prep_vis.py**

- Data Preparation: this script reads molecular data from a CSV file and standardizes SMILES strings. Molecular descriptors are calculated using the Mordred library and stored in a new DataFrame.

- Feature Selection: it calculates and filters molecular descriptors based on variance thresholds, and removes highly correlated features.

- Categorization: Categorizes compounds based on a defined LogMIC cutoff, creating "Strong" and "Weak" categories. In this example we are interest in distinguish compounds with weak and strong antitubercular activity but you can divide it in three or more categories.

- Visualization: Generates various plots (e.g., box plots, histograms, scatter plot matrix) to visualize the distribution and relationships of molecular descriptors. It is extremelly important to analyse or data before running the ML part, specially to check how each feature relates with others.

- Final Outputs: Saves processed data and visualizations as CSV and PNG files


