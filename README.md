Time Sensitive Influence Maximisation
=====================================

Project for the CS229 (machine learning) course at KAUST.

## Abstract

In a diffusion network, the estimated likelihood of influence of one node on another depends
on the transmission rate of the influencer node. Existing approaches learn the transmission
rates and influence likelihoods from data, but consider the transmission rate to be time-
invariant. This project proposes a modification by modeling the transmission rate as a distribution
over time; specifically, over the 24 hours of a day. New influence **estimation** and **maximisation**
algorithms are proposed in this model, and implemented as modifications of **ConTinEst** and **InfluMax**,
respectively.

## Implementation

This project modifies the influence propogation model to consider transmission rates
of nodes that not static, but a distribution over the 24 hours of a day. This approach
is then composed of two parts:

   1. Influence estimation, built on [ConTinEst](http://www.cc.gatech.edu/~ndu8/DuSonZhaMan-NIPS-2013.html).
   2. Influence maximisation, built on [InfluMax](http://people.tuebingen.mpg.de/manuelgr/influmax/).

## Data

Both ConTinEst and InfluMax provide toy data and synthetic data generators.
I will evaluate my performance on these toy/synthetic data sets, and probably
on a real dataset from Twitter if time permits.
