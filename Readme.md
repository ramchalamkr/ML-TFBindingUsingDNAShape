# TF binding sites prediction using DNAShapeR and Machine Learning Approaches.
The code has been written in Python
#### First execute the comp561_DNAshape_allChromosomes.py to generate the positive DnaShape values. make sure the data provided for the project is in the same folder or specify the relative path in the file. Store the results as a pickle file.
#### Next execute the comp561_DNAshape_negative.py to generate the negative Dnashape values. Store the results as a pickle file.
#### Execute each of the files separately generatePositiveFeaturesMax.py, generatePositiveFeaturesUAK21.py, generatePositiveFeaturesZNF263.py,
generatePositiveFeaturesE2F4.py. They use the output of 1 to generate final positive Feature set for each TF. These final feature sets are stored as another pickle file.
#### Execute the file generateNegativeFeatures.py. This uses output of step 2 to generate the final feature set for negative values. This result is stored as another pickle file
#### Final step is to execute the ML_PredictionForTF_Final.py. This file uses the output of step 3 and 4 to run the classifiers. Based on the TF length small modification is required in the code to run it for different TFs.

