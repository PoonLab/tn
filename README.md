# genSpaceMAUP
*Optimizing pairwise distance based component clustering through the Modifiable Areal Unit Problem*

The script genSpaceMAUP.R takes in tn93 distance data from a set of population sequences and selects the most optimal cutoff threshold for the purposes of predicting component growth in a network created from this data. This can be used to help parameterize molecular clustering techniques for outbreak detection and growth-forcasting in public health, and was designed to illustrate the effect that parameters such as the tn93 distance cutoff, can have on the overall information within a set of clusters. 

Our code does this by using data from earlier time to establish a log-linked model of cluster growth (discrete network components are defined here as "clusters"). Using the latest time point, we can then validate that model of growth, and across different cutoff thresholds, this model will perform differently overall. The measure of performance used here is Generalized AIC (GAIC), which represents the difference in performance between a null and proposed model. In this case, our proposed model is based on recency (ie. two cases close together in time are more likely to connect to eachother in a network) and is corrected by the overall edge density of cases from a given year. Because modulating the tn93 cutoff distance used to create networks generates a trade-off between case coverage and random error, a selectable optimum GAIC should exist between extremes and this will be represented as a minimum, negative spike in GAIC.


### NOTE

* Fasta file headers must be of the form ID_Year for example: "K113H63_2005", before tn93 analysis is run on them. Year may be collection or diagnostic, models will work better with a diagnostic year
  
* This method performs best on large data sets (>1000 sequences). It can also be ineffective if the newest time-point has a lack of data (<100 sequences). Currently, it will automatically disregard cases from the newest year if there are fewer than 63 in that set, treating the next newest year as the new year.

* Other analysis scripts ("_An" or "_cmpr" suffix) will have varying levels of accessability. Users will need to consult the inline comments for actual usage


### USAGE:

`Rscript genSpaceMAUP.R -f ../ExampleData/Example_tn93.txt -o Example`

The "Example" directory contains example output and data-set for user testing. 


### REQUIREMENTS:

* R.utils [https://cran.r-project.org/web/packages/R.utils/index.html)

* dplyr [https://cran.r-project.org/web/packages/dplyr/index.html]
 
* tn93 [https://github.com/veg/tn93]


### OPTIONS:

* **f**: The filepath to a tn93 output file. By default, this takes the standard input.

* **o**: The output filepath (do not include extensions). This defaults to the input filename without extensions

* **g**: The file path to saved graphical info. If you have already run this and saved a graph, this will save the trouble of making one.


### OUTPUT:
* A PDF image with the suffix "_VS.pdf" summarizing the proposed model's fit, the null model's fit and the difference between them.

* A ".rds" file with the suffix "_GD.rds" representing graph and model information at 50 different cutoff distances. 

* A ".rds" file with the suffix "_Optimum.rds" representing graph and model information at the optimum cutoff distance. 

  * **v**: A data frame of vertices. Each vertex represents a dated sequence

  * **e**: A data frame of edges between vertices. This will cover a completely connected graph, made from **$v**

  * **f**: Used in the creation of our predictive model. Records the minimum retrospective edge from each case and it's time-lag (summarizes as positives for a binomial model).

  * **g**: Growth of clusters. Each new case linked to a given cluster is considered an instance of growth. Note that new cases are only linked to their closest retrospective neighbour and only if that distance is below the cutoff distance.

  * **c**: Size of clusters (before growth)

  * **gaic**: The GAIC measurement. Measured from $nullFit and **$ageFit**

  * **nullFit**: The fit of a null model of growth. In this case, every case is considered equally likely to connect to new cases

  * **ageFit**: The fit of our proposed model. In this case, cases closer to the newest cases in time are considered more likely to connect to them, thus weigthing their clusters higher

  * **ageMod**: The model generated from our data set which informs ageFit. Currently this is a binomial model based on **$f**
