
# Load libraries
import pandas as pd 
from pandas import read_csv
from pandas.plotting import scatter_matrix
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem import PandasTools
from rdkit.Chem import MolStandardize
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdRGroupDecomposition
from rdkit.Chem import rdDepictor
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Draw import MolsToGridImage
from rdkit.Chem.Draw import IPythonConsole
from mordred import Calculator, descriptors
from sklearn.feature_selection import VarianceThreshold
import seaborn as sns

#prepare descriptors_data
inputdata = pd.read_csv('ML_List.csv',sep=",")

PandasTools.AddMoleculeColumnToFrame(inputdata, smilesCol='SMILES')
inputdata.rename(columns = {'ROMol':'mol'}, inplace = True)

pd.set_option("display.max_rows", 12)
inputdata.to_csv('output.csv', index=False)

#Perform molecule standardization
SMILES_nstd = inputdata.SMILES

def SMILES_standardization(SMILES_nstd):
    uncharger = rdMolStandardize.Uncharger()
    try:
        # Combining neutralization with the parent fragment identification
        m = Chem.MolFromSMILES(SMILES_nstd)
        std = uncharger.uncharge(rdMolStandardize.FragmentParent(m))
        SMILES = Chem.MolToSMILES(std)
    except:
        SMILES = np.nan
    return SMILES
 
calc = Calculator(descriptors, ignore_3D=True)
mol = inputdata.mol
descriptors = calc.pandas(mol)

MW = descriptors['MW'] # Exact molecular weight
nAtom = descriptors['nAtom'] # Number of all atoms
nHeavyAtom = descriptors['nHeavyAtom'] # Number of heavy atoms
nAcid = descriptors['nAcid'] # Acidic group count
nBase = descriptors['nBase'] # Basic group count
nRot = descriptors['nRot'] # Rotatable bonds count
nAromBond = descriptors['nAromBond'] # Aromatic bonds count
nRing = descriptors['nRing'] # Ring count
nAroRing = descriptors['naRing'] # Aromatic ring count
HBA = descriptors['nHBAcc'] # Number of hydrogen bond acceptor
HBD = descriptors['nHBDon'] # Number of hydrogen bond donnor
LogP = descriptors['SLogP'] # Wildman-Crippen LogP
TopoPSA = descriptors['TopoPSA'] # Topological polar surface area
Polarizability = descriptors['bpol'] # Bond polarizability
FractionCSP3 = descriptors['FCSP3'] # The fraction of C atoms that are SP3 hybridized
Lipinski = descriptors['Lipinski'] # Lipinski rule of five 

descriptors_list = pd.DataFrame(inputdata.Identifier)

descriptors_list['MW'] = MW
descriptors_list['nAtom'] = nAtom
descriptors_list['nHeavyAtom'] = nHeavyAtom
descriptors_list['nAcid'] = nAcid
descriptors_list['nBase'] = nBase
descriptors_list['nRot'] = nRot
descriptors_list['nAromBond'] = nAromBond
descriptors_list['nRing'] = nRing
descriptors_list['nAroRing'] = nAroRing
descriptors_list['HBA'] = HBA
descriptors_list['HBD'] = HBD
descriptors_list['LogP'] = LogP
descriptors_list['TopoPSA'] = TopoPSA
descriptors_list['Polarizability'] = Polarizability
descriptors_list['FractionCSP3'] = FractionCSP3
descriptors_list['Lipinski'] = Lipinski

descriptors_list.to_csv('descriptors_list.csv', index=False)  # Saves the DataFrame to 'descriptors_list.csv'

# Load the descriptors_datas
descriptors_data1 = pd.read_csv('ML_List.csv')
descriptors_data2 = pd.read_csv('descriptors_list.csv')

# Extract the last four columns of descriptors_data1
last_four_columns = descriptors_data1.iloc[:, -4:]

# Concatenate the last four columns of descriptors_data1 with descriptors_data2
merged_descriptors_data = pd.concat([descriptors_data2, last_four_columns], axis=1)

# Save the updated descriptors_data2 to a new CSV file
merged_descriptors_data.to_csv('descriptors_list.csv', index=False) 

# Read the CSV file and select only the columns of interest
new_data = pd.read_csv('descriptors_list.csv', usecols=['Identifier', 'LogMIC', 'vina_alphafold_22', 'Autodock_4P8L_22', 'Vinardo_4P8L_22'])

# Merge the new data with the existing descriptors_list DataFrame
descriptors_list = descriptors_list.merge(new_data, on='Identifier', how='left')

# Optionally, handle missing values
descriptors_list.fillna(0, inplace=True)  # Replace missing values with 0

# Calculate the variance along each column
var_descriptors_list = descriptors_list.drop(columns=['Identifier'])
column_variances = var_descriptors_list.var()

# Create a DataFrame to display variances of each feature
variance_df = pd.DataFrame({'Feature': column_variances.index, 'Variance': column_variances.values})
pd.set_option('display.max_rows', None)
print(variance_df)

# Define a variance threshold to feature selection
sel = VarianceThreshold(threshold=0.2)
descriptors_filtered = sel.fit_transform(descriptors_list)

# Get the indices of the selected features
selected_indices = sel.get_support(indices=True)

# Get the selected column names
selected_column_names = descriptors_list.columns[selected_indices]

# Create a DataFrame with only the selected features
descriptors_list = pd.DataFrame(descriptors_filtered, columns=selected_column_names)

descriptors_list.to_csv('descriptors_list_selected.csv', index=False)  # Saves the DataFrame


# List to store correlation coefficients
correlation_coefficients = []
descriptors_data = descriptors_list.drop(columns=['Identifier'])

# Iterate over pairs of columns to calculate the correlation coefficient
for i in range(len(descriptors_data.columns)):
    for j in range(i + 1, len(descriptors_data.columns)):
        column1 = descriptors_data.columns[i]
        column2 = descriptors_data.columns[j]
        # Calculate correlation coefficient between the columns
        correlation_coefficient = descriptors_data[column1].corr(descriptors_data[column2])
        color = '\x1b[32m' if correlation_coefficient < 0.8 else '\x1b[31m'  # Color the output based on correlation strength
        # Print the correlation coefficient with appropriate color
        print(f"Pearson's correlation coefficient between '{column1}' and '{column2}': {color}{correlation_coefficient:.2f}\x1b[0m")

descriptors_data = descriptors_data.drop(columns=['nAtom', 'nHeavyAtom', 'nRing', 'nAromBond', 'HBA'])

# shape
print(descriptors_data.shape)

# descriptions
print(descriptors_data.describe())

# Define the cutoff - 3 categories
#cutoff_strong = -0.7
#cutoff_weak = 0.7
#descriptors_data['strong'] = descriptors_data['LogMIC'] < cutoff_strong
#descriptors_data['moderate'] = (descriptors_data['LogMIC'] <= cutoff_weak) & (descriptors_data['LogMIC'] >= cutoff_strong)
#descriptors_data['weak'] = descriptors_data['LogMIC'] < cutoff_weak

#descriptors_data['Category'] = descriptors_data.apply(
#    lambda row: 'Strong' if row['strong'] else 'Moderate' if row['moderate'] else 'Weak',
#    axis=1
#)

# Define the cutoff - 2 categories
cutoff = -0.35
descriptors_data['strong'] = descriptors_data['LogMIC'] < cutoff
descriptors_data['weak'] = descriptors_data['LogMIC'] > cutoff

descriptors_data['Category'] = descriptors_data.apply(
    lambda row: 'Strong' if row['strong'] else 'Weak',
    axis=1
)

# Plot the counts for each category
plt.figure(figsize=(8, 6))  # Optional: Set the figure size
sns.countplot(x='Category', data=descriptors_data)
plt.title('Counts of Categories')  # Optional: Add a title to your plot

# Save the plot as a PNG file
plt.savefig('category_counts.png', dpi=300, bbox_inches='tight')

#descriptors_data = descriptors_data.drop(columns=['strong','moderate','weak'])
descriptors_data = descriptors_data.drop(columns=['strong','weak'])
descriptors_data.to_csv('descriptors_data.csv', index=False)
print(descriptors_data.groupby('Category').size())

plot_data = descriptors_data.drop(columns=['Lipinski', 'Category'])

# Box plot
plot_data.plot(kind='box', subplots=True, layout=(3,4), sharex=False, sharey=False)
plt.savefig('box_plot.png')  # Save the box plot

# Histograms
plot_data.hist()
plt.savefig('histograms.png')  # Save the histograms

# Scatter plot matrix
scatter_matrix(plot_data, alpha=0.8, figsize=(11, 11), diagonal='kde')

# Save the scatter plot matrix
plt.savefig('scatter_matrix.png')

#category vs logMIC

file_path = 'descriptors_data.csv'  # Replace with your file path
data = pd.read_csv(file_path)

# Define color mapping for categories
color_map = {
    'Strong': 'red',
    'Moderate': 'gray',
    'Weak': 'green'
}

# Create a scatter plot
plt.figure(figsize=(10, 6))
for category, color in color_map.items():
    subset = data[data['Category'] == category]
    plt.scatter(subset['LogMIC'], subset['Category'], color=color, label=category)

# Add labels and title
plt.xlabel('LogMIC')
plt.ylabel('Category')
plt.title('LogMIC vs Category')
plt.legend(title='Category')

# Save the plot to a file
output_file = 'logmic_vs_category_plot.png'  # Specify the file name and format
plt.savefig(output_file)
