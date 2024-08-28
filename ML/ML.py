import numpy as np
from pandas import read_csv
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.metrics import accuracy_score, precision_score, f1_score, matthews_corrcoef, classification_report
from sklearn.linear_model import LogisticRegression
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, AdaBoostClassifier, ExtraTreesClassifier
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
import matplotlib.pyplot as plt

# Load the dataset
inputfile = 'descriptors_data.csv'
names = ['MW','nRot','nAroRing','HBD','LogP','TopoPSA','Polarizability','Lipinski','LogMIC','vina_alphafold_22','Autodock_4P8L_22','Vinardo_4P8L_22', 'Category']
dataset = read_csv(inputfile, names=names, header=0)
array = dataset.values
X = dataset.iloc[:, 0:12].astype(float).values  # Convert features to float
y = dataset.iloc[:, 12].values  # Target (last column, Category)

# Spot Check Algorithms
models = []
models.append(('LR', LogisticRegression(solver='liblinear', multi_class='ovr', class_weight='balanced')))
models.append(('LDA', LinearDiscriminantAnalysis()))
models.append(('KNN', KNeighborsClassifier()))
models.append(('CART', DecisionTreeClassifier(class_weight='balanced')))
models.append(('NB', GaussianNB()))
models.append(('SVM', SVC(gamma='auto', class_weight='balanced')))
models.append(('RF', RandomForestClassifier(class_weight='balanced', random_state=1)))
models.append(('GBM', GradientBoostingClassifier(random_state=1)))
models.append(('QDA', QuadraticDiscriminantAnalysis()))
models.append(('ADA', AdaBoostClassifier(random_state=1)))
models.append(('ETC', ExtraTreesClassifier(class_weight='balanced', random_state=1)))

# List of metrics
metrics = ['accuracy', 'mcc', 'f1', 'precision']
metric_funcs = {
    'accuracy': accuracy_score,
    'mcc': matthews_corrcoef,
    'f1': lambda y_true, y_pred: f1_score(y_true, y_pred, average='weighted'),
    'precision': lambda y_true, y_pred: precision_score(y_true, y_pred, average='weighted')
}

# Number of runs
n_runs = 5

# File to write results
output_file = "model_evaluation_results.txt"

# Open the file to write
with open(output_file, "w") as f:
    # Initialize a dictionary to store cumulative metrics
    cumulative_metrics = {name: {metric: [] for metric in metrics} for name, _ in models}

    # Perform the evaluation over multiple runs
    for run in range(n_runs):
        f.write(f"Run {run+1}/{n_runs}\n")
        X_train, X_validation, Y_train, Y_validation = train_test_split(X, y, test_size=0.20, random_state=run, shuffle=True)

        for name, model in models:
            f.write(f"\nEvaluating {name} on the validation set...\n")
            model.fit(X_train, Y_train)
            predictions = model.predict(X_validation)

            # Calculate and print all metrics
            for metric_name in metrics:
                metric_value = metric_funcs[metric_name](Y_validation, predictions)
                f.write(f"{metric_name.capitalize()}: {metric_value:.4f}\n")
                cumulative_metrics[name][metric_name].append(metric_value)

            # Write classification report
            f.write(classification_report(Y_validation, predictions) + "\n")

            # Print the number of predictions assigned to each category
            unique, counts = np.unique(predictions, return_counts=True)
            f.write("Number of predictions assigned to each category:\n")
            for category, count in zip(unique, counts):
                f.write(f"{category}: {count}\n")

            # Analyzing Misclassifications
            misclassified = Y_validation != predictions
            f.write(f"Number of misclassified samples: {misclassified.sum()}\n")

            # Inspect misclassified samples if needed
            misclassified_indices = np.where(misclassified)[0]
            f.write(f"Misclassified sample indices: {misclassified_indices.tolist()}\n")

        f.write("\n" + "="*50 + "\n")

    # Calculate and write the averages
    f.write("\nAverages over all runs:\n")
    for name in cumulative_metrics:
        f.write(f"\n{name} Averages:\n")
        for metric_name in metrics:
            avg_value = np.mean(cumulative_metrics[name][metric_name])
            f.write(f"Average {metric_name.capitalize()}: {avg_value:.4f}\n")

print(f"Results written to {output_file}")
