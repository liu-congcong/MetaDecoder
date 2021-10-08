import numpy
from sklearn.ensemble import IsolationForest


def isolation_forest(sequences, x, weight, threads, outlier, random_number):
    '''
    Parameters:
        sequences: an array of sequences.
        x: input data.
        weight: an array of sequence lengths.
        threads: number of threads.
        outlier: outlier ratio 0 - 1.
        random_number: the random number.
    Return:
        (inliers, outliers)
    '''
    isolationforest = IsolationForest(
        max_samples = sequences.shape[0],
        contamination = outlier,
        n_jobs = threads,
        random_state = random_number
    )
    isolationforest.fit(x, sample_weight = weight)
    prediction = isolationforest.predict(x)
    return (sequences[prediction == 1], sequences[prediction == -1])