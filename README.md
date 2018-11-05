#tn

## Objective
To train different clustering methods on a subset of the TN data, and then validate (evaluate) the trained methods using the remaining TN data to determine which method is the most effective.

* A clustering method defines a partition on a database population.
* We refer to each element of a partition as a subset. Usually a clustering method will only label some subsets as "clusters", so we want to have a more general term.
* We want to evaluate a partition on the following criteria:
* The number or proportion of new cases contained within the smallest number of subsets.
* The size of the subsets containing new cases.
* The predictive value of subsets containing new cases - e.g., do the cases currently contained in the subset (not new) make that subset distinctive from other subsets a consistent way, such that we are able to predict those subsets as having the next cases?
* Each clustering method can produce an infinite number of partitions based on its clustering thresholds. For example, TN93 has a distance cutoff that can vary from 0 to infinity.
* We want to split the database into training and validation portions. We want to fit a clustering method to the training portion by right-censoring the cases with respect to their sample collection dates, and then validate on the remaining portion.

## End points
* Are clustering methods effective at all for predicting the next cases?
* If so, which clustering method is the most effective?
