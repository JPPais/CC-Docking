# CC-Docking_Analysis
This repository is dedicated to the analysis of docking results to extract and treat ligand-protein distances and score, while processing them for predictive modeling  

Before performing the molecular docking run, ligand and receptor preparation is necessary. Ligand preparation is included, requiring only a list componds in SMILES format, with additional data necessary for the output analysis. For receptor preparation, AutoDockTools was employed, following "traditional" methodologies.

The docking runs were performed using Autodock Vina, hence all the output data analysis scripts are developed accordingly. Furthermore, docking execution scripts are not provided, since they are higly dependent of the computing setup you have access to. Also, the system employed delivered an *.out file containing the results for all the compounds docked and an initial preparation was performed. The *.out file was parsed and was splited in individual "compound_name.pdbqt" files for each compound.

After this, the sequence of scripts for analysis is provided, and leads to the assessment of a predictive model based in the docking results, with a binary classification system (active / inactive). 
