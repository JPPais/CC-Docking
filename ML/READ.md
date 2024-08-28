In this folder, we have two scripts to run our Machine Learning model.
We are using unpublished results as input so the input data is not present in this folder, but you can adapt it to your work.

first, you have to run: 
**data_prep_vis.py**

- Data Preparation: this script reads molecular data from a CSV file and standardizes SMILES strings. Molecular descriptors are calculated using the Mordred library and stored in a new DataFrame.

- Feature Selection: it calculates and filters molecular descriptors based on variance thresholds, and removes highly correlated features.

- Categorization: categorizes compounds based on a defined LogMIC cutoff, creating "Strong" and "Weak" categories. In this example, we are interested in distinguishing compounds with weak and strong antitubercular activity but you can divide it into three or more categories.

- Visualization: generates various plots (e.g., box plots, histograms, scatter plot matrices) to visualize the distribution and relationships of molecular descriptors. It is extremely important to analyze our data before running the ML part, especially to check how each feature relates to others.

- Final Outputs: saves processed data and visualizations as CSV and PNG files

After having the processed data outputs, you can run: **ML.py**

- Feature Types: this script converts features to floating-point numbers for compatibility with machine learning models.
  
- Model Selection: The script evaluates several models, including Logistic Regression, Decision Trees, SVM, Random Forest, Gradient Boosting, and others.

- Metric Selection: It uses accuracy, Matthews correlation coefficient (MCC), F1 score, and precision to assess model performance.

- Model Evaluation: Each model is trained and validated using stratified cross-validation over multiple runs to ensure robust results.

- Misclassification Analysis: The script tracks misclassified samples to identify model weaknesses.

- Output: The performance metrics, classification reports, and misclassification details for each model and run are recorded in a text file. After all runs, the script calculates and records the average performance metrics for each model.

- Result Summary: Summarizes key metrics and misclassification statistics for easy comparison of models.




