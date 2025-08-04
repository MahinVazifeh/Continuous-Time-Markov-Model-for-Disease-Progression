## Multi-State Markov Modeling of Mondino Foundation Data
Overview

This R script performs preprocessing, analysis, and modeling of clinical data using Markov models and Hidden Markov Models (HMMs). It primarily focuses on modeling transitions between disease severity states in patients, based on neurological and clinical assessment data, using packages like msm, momentuHMM, and markovchain.
# Key Objectives

    Discretize EDSS clinical scores into severity states.

    Preprocess and filter patient data for modeling.

    Visualize distributions and correlations between clinical features.

    Fit Markov and Hidden Markov models to examine disease progression.

    Test the influence of covariates (e.g., Age, Visual, Deambulation) on transition dynamics.

    Visualize transition matrices and model goodness-of-fit.

# Main Components
1. Libraries Used

    Statistical modeling: msm, momentuHMM, MASS

    Markov modeling: markovchain, msmtools

    Visualization: ggplot2, ggpubr, corrplot, diagram

    Data manipulation: Hmisc, foreign, cluster, gclus

2. Discretization

    Function create_discret() transforms continuous EDSS scores into three discrete states:

        1 = Mild (EDSS 0–4)

        2 = Moderate (EDSS 4.5–7.5)

        3 = Severe (EDSS 8–10)

3. Data Preprocessing

    Loads clinical CSV data (first_data.csv).

    Rounds feature values and assigns column names.

    Filters out patients with only a single observation (required for Markov modeling).

    Saves cleaned data to updated_data.csv.

4. Data Visualization

    Uses ggplot2 to plot histograms and bar charts for:

        Disease features (e.g., Pyramidal, Visual, Mental)

        Patient characteristics (e.g., Age, Sex)

    Correlation heatmaps created using corrplot.

5. Statistical Modeling

    Ordinal regression (polr) to examine predictors of severity state.

    Markov Chain Estimation:

        Transition matrices estimated via markovchainFit().

        Plots transition diagrams using transitionPlot().

6. Hidden Markov Model (HMM)

    create_hmm_model() function wraps the HMM fitting process using msm():

        Outputs include Q-matrix, transition probability matrix, sojourn times, and hazard rates.

    Model fitted with and without covariates like:

        Age, Visual, Deambulation, Sensitive, etc.

    Some covariates (e.g., Pyramidal, Mental) caused numerical errors due to overfitting or instability.

7. Graphical Outputs

    Transition probabilities and Q-matrix plotted using plotmat().

    Visual comparison of estimated vs. real state prevalence via prevplot().

# Files Required

    first_data.csv: Raw patient-level input data.

    final_data.csv: Cleaned and corrected version for HMM modeling.
