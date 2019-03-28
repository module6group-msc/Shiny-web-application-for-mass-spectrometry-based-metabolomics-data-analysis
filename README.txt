Web application for time series analysis based on mass spectrometry data using R.

This is a web-based Shiny application, built during Module 6 of MSc Bioinformatics at University of Birmingham.
The main purpose of that application is to perform time series analysis for mass spectrometry-based metabolomics data.
The application gives the opportunity to perform:
-Normalization of dataset
-Parametric statistical test (repeated measures ANOVA)
-Non-parametric statistical test (Friedman test)
-Post hoc test for ANOVA models)
-(Post hoc test for Freidman tests is under construction...)
-Plots and graphs for significantly changed metabolites collected from ANOVA
-Plots and graphs for significantly changed metabolites collected from Friedman test

Demo dataset
We bulit the application making the stastistical analysis for a given metabolomics dataset.
The demo dataset was provided by Professor Warwick Dunn at University of Birmingham and was an untargeted metabolomics dataset. 
The biological study investigated 38 adults in three different age ranges (<30, 30-50 and >50) who undertook an exercise-related 
intervention with muscle biopsies collected at six different time points, 3 time points before a 20-week exercise intervention and
3 time points after the 20-week exercise intervention. The 3 time points were fasted state, immediately after exercise and 1.5 hours 
after feeding. Muscle biopsies were analysed applying ultra-performance liquid chromatography-mass spectrometry (LC/MS) and raw data 
was processed applying the open source software XCMS to generate the data matrix applied for the demonstration.
For our purposes, the factors that taken into consideration were TIMEPOINTS and different SUBJECTS.
The demo dataset was at .csv format with the following details:
•	First row: the names of the subjects
•	Second row: the names of the time points
•	Rest of the rows: each metabolite in each row
•	First column: “Subject” (first row) and “Timepoint” (second row) titles and the names of the metabolites measured (rest of the rows)
•	Rest of the columns: Intensity scores after LC/MS and using XCMS software.
PLEASE NOTE that the application works properly with datasets in the same format.

The application was built using "shiny" R package. (R version 3.5.3). To use the application it is not necessery to run RStudio or use R in any kind of
software enviroment.
The code of the statistical analysis and the application are available in GitHub:
https://github.com/module6group-msc/Shiny-web-application-for-mass-spectrometry-based-metabolomics-data-analysis

No programming skills required to use the application.
