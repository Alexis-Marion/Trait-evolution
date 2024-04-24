# Discrete trait evolution - binary trait \#1


***
#### Contents
* [Introduction to discrete character evolution](https://github.com/dsilvestro/PyRate/blob/master/tutorials/pyrate_tutorial_1.md#generate-pyrate-input-file-option-1)
* [Data preparation](https://github.com/dsilvestro/PyRate/blob/master/tutorials/pyrate_tutorial_1.md#generate-pyrate-input-file-option-1)  
* [The equal rate model](https://github.com/dsilvestro/PyRate/blob/master/tutorials/pyrate_tutorial_1.md#generate-pyrate-input-file-option-2)  
* [Beyond the equal rate model](https://github.com/dsilvestro/PyRate/blob/master/tutorials/pyrate_tutorial_1.md#estimation-of-speciation-and-extinction-rates-through-time)
* [Plotting ancestral state esimates](https://github.com/dsilvestro/PyRate/blob/master/tutorials/pyrate_tutorial_1.md#estimation-of-speciation-and-extinction-rates-through-time)

* [Return to Index](https://github.com/dsilvestro/PyRate/tree/master/tutorials#pyrate-tutorials---index) 
***



# Introduction to discrete character evolution
Understanding the tempo and mode of character evolution across the Tree of Life is a long-standing topic in macroevolution.
However, with

As phylogenetic data became increasingly available over th 
new methods are be implemented

 a plethora of discrete models of trait evolution have flourished. Such models are highly versatile and allow for a large range of hypothesis testing:

- Do my traits evolve at the same rates?
- What is the ancestral character state for my clade?
- Are there one or convergent appearances for my trait?
- Does Dollo's law apply to my character set?
- Are the evolution of two (or more) traits correlated?

Throughout the following tutorials, we will be able to address these questions.
In this first tutorial, we will be introducing the main properties of discrete models of trait evolution, their implementation, and how we can use them to answer macroevolutionary questions.


# Data preparation



| Taxon\_name   | Discrete_Trait |
| ------------- |:-------------:|
sp_1 | 2 
sp_2 | 2 
sp_3 | 2 
sp_4 | 2 
sp_5 | 2 
sp_6 | 1 
sp_7 | 1 
sp_8 | 1 
sp_9 | 1 
sp_10 | 1 


2. **Launch R** as explained above.  
3. **Load the *pyrate_utilities.r* file** as explained above.  
4. **Parse the raw data and generate PyRate input file**. Type in R:   
  `extract.ages(file="…/PyRate/example_files/Ursidae.txt", replicates=10)` This function includes `replicates` and `cutoff` options as the `extract.ages.pbdb()` function described above (option 1). 


# Estimation of speciation and extinction rates through time

## Defining the preservation model
**Non-homogeneous Poisson process of preservation (NHPP)**. By default, PyRate assumes an NHPP model in which preservation rates change during the lifespan of each lineage following a bell-shaped distribution (Silvestro et al. 2014 Syst Biol). 

**Homogeneous Poisson process (HPP)**. Alternatively, an HPP model, in which the preservation rate is constant through time, is also available, using the command `-mHPP`:

`python PyRate.py .../Canis_pbdb_data_PyRate.py -mHPP`

**Time-variable Poisson process (TPP)**. An alternative model of preservation assumes that preservation rates are constant within a predefined time frame, but can vary across time frames (e.g. geological epochs). This model is particularly useful if we expect rate heterogeneity to occur mostly through time, rather than among lineages.

We can set up a model in which preservation rates are estimated independently within each geological epoch by using the command `-qShift` and providing a file containing the times that delimit the different epochs (an example is provided in the example file _PyRate-master/example\_files/epochs\_q.txt_):  

`python PyRate.py .../Canis_pbdb_data_PyRate.py -qShift .../epochs_q.txt`  

The default prior on the vector of preservation rates is a single gamma distribution with shape = 1.5 and rate = 1.5. These parameter values can be changed using the command `-pP`, for example `-pP 2 0.1` changes the shape to 2 and the rate to 0.1.
Alternatively, the rate parameter of the prior can be estimated from the data to reduce the subjectivity of the choice and allow the prior to better adapt to the data. This option is available by setting the rate parameter to 0, e.g.: 

`python PyRate.py .../Canis_pbdb_data_PyRate.py -qShift .../epochs_q.txt -pP 1.5 0`

When the rate is set to zero, PyRate assigns a vague exponential hyper-prior to the rate and samples the rate along with all other model parameters. This approach is recommended whenever multiple preservation rates are estimated


**Gamma model of rate heterogeneity**. NHPP, HPP, and TPP models can all be coupled with a Gamma model of rate heterogeneity, which enables us to account for heterogeneity in the preservation rate across lineages. To set the Gamma model we add the flag `-mG`:

`python PyRate.py .../Canis_pbdb_data_PyRate.py -mG` [NHPP model]  
`python PyRate.py .../Canis_pbdb_data_PyRate.py -mHPP -mG` [HPP model]
`python PyRate.py .../Canis_pbdb_data_PyRate.py -qShift .../epochs_q.txt -mG`  


**Saving per-lineage relative preservation rates**. When combining a TPP model with a Gamma model you can log to a file the estimated relative preservation rate for each lineage. This is done by adding the flag `-log_sp_q_rates` to the command and will save one additional log file containing one column for each lineage in the data set with the estimated relative rate of preservation (note that the actual rate will depend on the mean preservation rate and it variation through time). 

`python PyRate.py .../Canis_pbdb_data_PyRate.py -qShift .../epochs_q.txt -mG -log_sp_q_rates` 

### Model testing across preservation models
A maximum likelihood test across (described [here](https://www.cambridge.org/core/journals/paleobiology/article/improved-estimation-of-macroevolutionary-rates-from-fossil-data-using-a-bayesian-framework/334F08A74A6C92F1FEAD91A71FE59A1C)) is available to assess which of NHPP, HPP, or TPP is best supported by the data. To run the test you need to provide the input file and the file providing the times of rate shift for the TPP model and add the flag `-PPmodeltest`

``python PyRate.py .../Canis_pbdb_data_PyRate.py -qShift .../epochs_q.txt -PPmodeltest``

Note that the Gamma model is not tested here, but we recommend to add it to the best model among NHPP, HPP, or TPP as selected by the model testing. 

## Analysis setup
Here we describe the main settings for a standard analysis of fossil occurrence data using PyRate. The analysis will estimate:

1. **origination and extinction times** of each lineage  
2. **preservation rate** and its level of **heterogeneity**  
3. **speciation and extinction rates** through time. 

Temporal rate variation is introduced by rate shifts. The number and temporal placement of shifts are estimated from the data using the RJMCMC algorithm (default) or the BDMCMC algorithm. **A tutorial specific to the RJMCMC is available [here](https://github.com/dsilvestro/PyRate/blob/master/tutorials/pyrate_tutorial_3.md#pyrate-tutorial-3).** Both BDMCMC and RJMCMC estimate birth-death models with rate shifts, but **the RJMCMC algorithm is recommended as it generally provides more accurate results (see open access paper [here](https://www.cambridge.org/core/journals/paleobiology/article/improved-estimation-of-macroevolutionary-rates-from-fossil-data-using-a-bayesian-framework/334F08A74A6C92F1FEAD91A71FE59A1C))**.


The analysis requires a **PyRate input file** generated by the R function described above. The first argument we need to provide is the input file:

`python PyRate.py .../Canis_pbdb_data_PyRate.py -A 2`

In most operating systems (including Mac OS, Windows, and Ubuntu) you can drag and drop the input file onto the terminal window to paste the full path and file name. The additional command `-A 2` is used to run the BDMCMC algorithm, whereas using `-A 4` (or omitting the flag) will set the algorithm to RJMCMC. 

PyRate includes default settings for all parameters except for the input file. While several parameters should be changed only when experiencing convergence issues, there are a few that are very important as they change the basic model assumptions or the length of the analyses - see below.

Since the input file generated in the previous steps included 10 randomized replicates of the fossil ages, we can specify which replicate we want to analyze. **Ideally, we should analyze multiple randomized replicates and combine the results to incorporate dating uncertainties** in our rate estimates. To specify the which replicate we want to analyze, we use the flag `-j` followed by the replicate number. For instance using:

`python PyRate.py .../Canis_pbdb_data_PyRate.py -A 2 -mHPP -mG -j 1`

we set the analysis to consider the first replicate.

We can (and in some cases should) change the number of MCMC iterations and the sampling frequency. By default PyRate will run 10,000,000 iterations and sample the parameters every 1,000 iterations. Depending on the size of the data set you may have to **increase the number iterations to reach convergence** (in which case it might be a good idea to sample the chain less frequently to reduce the size of the output files). This is done using the commands `-n` and `-s`:

`python PyRate.py .../Canis_pbdb_data_PyRate.py -mG -n 20000000 -s 5000`

Under these settings PyRate will run for 20 million iterations sampling every 5,000. Thus the resulting log files (see below) will include 4,000 posterior samples.  


## Output files
**Note that the output from RJMCMC analyses is by default different, as described [here](https://github.com/dsilvestro/PyRate/blob/master/tutorials/pyrate_tutorial_3.md#rjmcmc-output).**

The PyRate analysis described above produces three output files, saved in a folder named *pyrate\_mcmc\_logs* in the same directory as the input file:

#### 1) sum.txt  
   Text file providing the complete list of settings used in the analysis.  
#### 2) mcmc.log   
Tab-separated table with the MCMC samples of the posterior, prior, likelihoods of the preservation process and of the birth-death (indicated by *PP_lik* and *BD_lik*, respectively), the preservation rate (*q_rate*), the shape parameter of its gamma-distributed heterogeneity (*alpha*), the number of sampled rate shifts (*k_birth*, *k_death*), the time of origin of the oldest lineage (*root_age*), the total branch length (*tot_length*), and the times of speciation and extinction of all taxa in the data set (*\*_TS* and *\*_TE*, respectively). When using the TPP model of preservation, the preservation rates between shifts are indicated as *q_0, q_1, ... q_n* (from older to younger).
#### 3) marginal_rates.log 
Tab-separated table with the posterior samples of the marginal rates of speciation, extinction, and net diversification, calculated within 1 time unit (typically Myr). 


## Summarize the results
The log files can be opened in the program [**Tracer**](https://github.com/beast-dev/tracer/releases/tag/v1.7.1) to check if the MCMC has converged, for example looking the the Effective Sample Sizes (ESS), and determine the proportion of burnin. 
When running multiple replicates the ESS values should be computed from the log files of each replicate (i.e. not from a [combined log file](https://github.com/dsilvestro/PyRate/blob/master/tutorials/pyrate_tutorial_1.md#combine-log-files-across-replicates)), because each replicate is based on a different randomized dataset and therefore might converge to slightly different solutions.
These differences however can be used to incorporate the uncertainties in the dating of the fossil occurrences.

The **mcmc.log** file can be used to calculate the sampling frequencies of birth-death models with different number of rate shifts. This is done by using the PyRate command `-mProb` followed by the log file:

`python PyRate.py -mProb .../Canis_pbdb_data_mcmc.log -b 200`

where the flag `-b 200` indicates that the first 200 samples will be removed (i.e. the first 200,000 iterations, if the sampling frequency was set to 1,000). This command will provide a table with the relative probabilities of birth-death models with different number of rate shifts. 


The **marginal_rates.log** file can be used to generate rates-through-time plots using the function `-plot`:

`python PyRate.py -plot .../Canis_pbdb_data_marginal_rates.log -b 200`

This will generate an R script and a PDF file with the RTT plots showing speciation, extinction, and net diversification through time. A slightly different flavor of the RTT plot can be obtained using the flag `-plot2` instead of `-plot`. 



#### Combine log files across replicates
To combine log files from different replicates into one you can use the command:

`PyRate.py -combLog path_to_your_log_files -tag mcmc -b 100`

where `path_to_your_log_files` specifies the directory where the log files are; `-tag mcmc` sets PyRate to combine all files that contain _mcmc.log_ in the file name (you can use different tags if you need); and `-b 100` specifies that the first 100 samples from each log file should be excluded as burnin – the appropriate number of burnin samples to be excluded should be determined after inspecting the log files, e.g. using Tracer.

To generate a rates-through-time plot that combines all replicates, you can use the command:

`PyRate.py -plot path_to_your_log_files -tag Canis_pbdb -b 100`

This will combine all the _marginal\_rates.log_ files which include `Canis_pbdb` in the file name and combine the results in a single plot.​ Different tags can be used to determine which files are to be combined.    


#### Plot preservation rates through time
Preservation rates estimated using the TPP model can be summarized and plotted using the command `-plotQ` followed by the `mcmc.log` file (including path to file):

`PyRate.py -plotQ .../Canis_pbdb_data_mcmc.log -qShift epochs.txt -b 100`

Note that the file with the times of rate shift (e.g. `PyRate-master/example_files/epochs_q.txt`) should be provided as well.


***
## Speciation and extinction rates within fixed time bins
#### Analysis setup
PyRate can also fit birth-death models in which the **number and temporal placement of rate shifts is fixed a priori**, e.g. based on geological epochs. In this case a file with the predefined times of rate shifts must be provided using the command `-fixShift`. The format of this file is very simple and an example is available here: `.../PyRate-master/example_files/epochs.txt`. This model assumes half-Cauchy prior distributions for speciation and extinction rates between shifts, with a hyper-prior on the respective scale parameter to reduce the risk of over parameterization. 
To enforce fixed times of rate shifts we use the following command:

`python PyRate.py .../Canis_pbdb_data_PyRate.py -fixShift .../epochs.txt`

The other options described above to set preservation model, length of the MCMC, and sampling frequency are also available in this case.

#### Summarize the results
Running PyRate with fixed times of rate shifts produces the same 3 output files described in the previous analysis. The main difference is in the *mcmc.log* file where we will no longer have the estimate number of rate shifts (*k_birth*, *k_death*) as these are fixed. However, the log file now includes speciation/ extinction rates between shifts (named *lambda_0, lambda_1,* ... and  *mu_0, mu_1, ...*, respectively), and the estimated scale parameters of the half-Cauchy prior distributions assigned to speciation and extinction rates, indicated by *hypL* and *hypM*, respectively.

RTT plots can be generated as in the previous analysis using the command `-plot` (or `-plot2`) followed by the path to the *marginal_rates.log* file and setting the number of samples be discarded as burnin (e.g. `-b 100`).



## Setting fixed shifts at the boundaries, while searching for rate shifts between them

Sometimes fossil data sets are truncated by max and min boundaries (for instance because occurrences are only available within a time window). This can cause apparent rate shifts at the edges of the time range, which reflect the sampling bias. In this case, you can fix _a priori_ times of rate shift based on the known temporal boundaries of the data set and estimate the rates within the time window, ignoring what happens beyond the boundaries.  This feature can be combined with the [RJMCMC algorithm](https://github.com/dsilvestro/PyRate/blob/master/tutorials/pyrate_tutorial_3.md#pyrate-tutorial-3) (`-A 4`) to infer rate shifts within the allowed time window:

`python PyRate.py <data_set> -A 4 -edgeShift 18 2`

where -A 4 sets the RJMCMC algorithm (that looks for rate shifts), and `-edgeShift 18 2` sets fixed times of shifts at times 18 and 2. With these settings the algorithm will only search for shifts within this time range.
You can also set only a max age boundary shift using:

`python PyRate.py <data_set> -A 4 -edgeShift 18 0`

or a min age boundary shift using:

`python PyRate.py <data_set> -A 4 -edgeShift inf 0`

When summarizing the results, only rates estimated within the time window should be considered, e.g. using the [`plotRJ` command](https://github.com/dsilvestro/PyRate/blob/master/tutorials/pyrate_tutorial_3.md#plot-rates-through-time-and-rate-shifts) with the flags `-root_plot` and `-min_age_plot`.
