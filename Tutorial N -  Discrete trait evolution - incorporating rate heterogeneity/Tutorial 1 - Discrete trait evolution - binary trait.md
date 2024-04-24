# Discrete trait evolution - binary trait \#1


***
#### Contents
* [Introduction to discrete character evolution](Introduction-to-discrete-character-evolution)
* [Data preparation](Data-preparation)  
* [The equal rate model](The-equal-rate-model)  
* [Beyond the equal rate model](Beyond-the-equal-rate-model)
* [Plotting ancestral state estimates](Plotting-ancestral-state-esimates)

* [Return to Index](https://github.com/Alexis-Marion/Trait-evolution/) 
***



# Introduction to discrete character evolution
Understanding the tempo and mode of character evolution across the Tree of Life is a long-standing topic in macroevolution.
The first attempts to outline major character transition and radiation were mostly performed on small trees and relied on relatively simple frameworks (_e.g_ maximum parsimony). Yet, such analyses, although useful, most likely lack resolution and precision and probably were prone to methodological biases. Such issues, likely stem from the relatively small size of the examined phylogeny and were also constrained by computational limits.
However, as both phylogenetic data and powerful computers became increasingly available over the last decades, the possibility of evaluating more complex scenarios of trait evolution became more tangible.
Thus through rapid methodological improvements, a plethora of discrete models of trait evolution have flourished.
Such models are highly diverse and versatile and allow for a wide range of hypothesis testing:

- Do my traits evolve at the same rates?
- What is the ancestral character state for my clade?
- Are there one or convergent appearances for my trait?
- Does Dollo's law apply to my character set?
- Are the evolution of two (or more) traits correlated?

Throughout the following tutorials, we will be able to address these questions.
Specifically, in this first tutorial, we will be introducing the main properties of discrete models of trait evolution, their implementation, and how we can use them to answer macroevolutionary questions.

# Data preparation

The first step before performing any phylogenetic comparative analyses is to prepare the data for the analyses. Usually, comparative biologists use dataframe to perform their analyses. Conventionally, the first column is the taxon name, and the second is the examined trait. In the case of a binary character state, the dataframe may look like this :

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

Throughout all the following tutorials, we will follow the same scheme (although with some variations).

Although we assume that we possess a well-formatted dataframe, we are not impervious to potential issues. First of all, some species present in the datasets might be absent from the phylogeny. Likewise, some species sampled in the phylogeny might be absent in the dataset. It is also possible that the user misspecified a character state for one, or more, taxa.

Thus to prepare the data for any following analyses, we will run the script `data_prepartion_for_comparative_analyses.r`. This script allows the user to identify incompatibilities between their phylogenies and their dataset, and to highlight possible misspecifications by plotting the character states on the phylogeny. This script also includes the option to automatically fill the dataset with NA taxa included in the phylogeny, but not in the dataframe.

# The equal rate model

The most simple model for discrete character evolution is called the Markov model (MK). First implemented for discrete trait data by Pagel (1994), this model assumes that transitions between character states follow a Markov process (thus the name). The MK model is a variation of the famous Jukes & Cantor model for nucleotide transition (also known as JC69).
The assumptions of the MK model in its purest form are rather simple and almost identical to JC69. The MK model assumes for k unorder character states, stationary character distribution π where π equals [1/k 1/k ... 1/k], equal transition rates (q) between any k character states (qab = qba). Consequently, the MK model in this form estimates only one shared transition parameter. For this reason, this model is also known as the equal rate model (ER).

# Beyond the equal rate model

Although the ER is already very useful, assuming equal transition rates among all character states might not be biologically accurate. Thus, similarly to DNA model of sequence evolution, new models relaxed the assumption, by allowing rates to vary across character states. The two most well-known unequal rate models are the Symmetrical rate (SYM; Paradis et al., 2004) and the All Rates Differ (ARD; Paradis et al., 2004). The SYM model assumes that for a K-states trait, with K > 2, qab = qba but qab \neq qca. In the case of a binary state, qab and qba are the only possible transtion and thus the SYM model collapses into the equal rate model. Conversely, the ARD model assumes that all rates may be different, thus qab \neq qba. Thus
Generally, assuming a trait with K unordered states, the number of possible transition rates will be q = ((K-1)*K) for an ARD model and ((K-1)*K)/2 for a SYM model. Thus, as the number of possible states increases, so do the potential transition rates. Consequently, a major limitation arises from these two models: overparametrization. Overparametrization is the tendency for one model (usually maximum-likelihood) to have more parameters than can be estimated from the data. Usually, such limits are obtainable when the number of estimated parameters comes close to the number of data, which is in our case the number of tips.
Namely, if one works on a relatively small-sized phylogeny, let's ~50 species, and want to estimate transition rates for a 7-state character on an ARD framework, given the equation found above, q = 42, which is very close to what. Consequently, the model will most likely be overparametrized, and thus model selection is de facto flawed.


# Plotting ancestral state estimate



