import numpy as np


def all_vs_all_similarity_matrix(spectrums, similarity_function,
                                 filename=None, safety_points=None):
    """Calculate similarity matrix of all spectrums vs all spectrums.

    Args:
    ----
    pectrums, similarity_function,
                                 filename=None, safety_points=None
    """

    if safety_points is not None:
        # Save matrix along process
        total_num_calculations = int((len(spectrums)**2)/2 + 0.5 * len(spectrums))
        safety_interval = int(total_num_calculations/safety_points)

    similarities = np.zeros((len(spectrums), len(spectrums)))
    num_matches = np.zeros((len(spectrums), len(spectrums)))

    count = 0
    for i in range(len(spectrums)):
        for j in range(i, len(spectrums)):
            score, matches = similarity_function(spectrums[i], spectrums[j])
            similarities[i, j] = score
            num_matches[i, j] = matches
            count += 1
            # Show progress
            if (count+1) % 10000 == 0:
                print("\r", "About {:.3f}% of similarity scores calculated.".format(100 * count/total_num_calculations), end="")

            # Create safety points
            if filename is not None and safety_points is not None:
                if (count+1) % safety_interval == 0:
                    safety_filename = filename.split(".")[0] + "safety"
                    np.save(safety_filename + ".npy", similarities)
                    np.save(safety_filename + "_matches.npy", num_matches)

    # Symmetric matrix --> fill
    for i in range(1, len(spectrums)):
        for j in range(i):
            similarities[i, j] = similarities[j, i]
            num_matches[i, j] = num_matches[j, i]

     # Save final results
    if filename is not None:
        np.save(filename, similarities)
        np.save(filename.split(".")[0] + "_matches.npy", num_matches)

    return similarities, num_matches
