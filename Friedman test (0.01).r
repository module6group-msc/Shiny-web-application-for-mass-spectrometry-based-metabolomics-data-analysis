#Required packages
if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)
if (!require("clusterSim")) install.packages("clusterSim")
library(clusterSim)
if (!require("psych")) install.packages("psych")
library(psych)
if (!require("ggplot2")) install.packages("ggplot2")
library(ggplot2)
if (!require("gplots")) install.packages("gplots")
library(gplots)

#Set your working directory (where our datafile is)
#Load data
#without_missing.csv is the original Dummy file that Rick gave us, WITHOUT QC and SUBJECTS WITH MISSING VALUES
dt <- read.csv("without_missing.csv", header = T, sep = ",", quote = "", dec = ".")
#Make first column the rownames of "data"
data <- data.frame(dt[,-1], row.names=dt[,1])

#This will create a new directory in your working directory with the name "Friedman"
#At the end, everything will be inside the Friedman directory
Friedman_0.01 <- paste("Friedman_0.01") #set the name
dir.create(Friedman_0.01) #create directory
setwd(Friedman_0.01) #set the new directory as working directory

###Normalization###
length_data <- nrow(data) #number of rows in data
for_nomr=data[3:length_data,] #collect only the region of data that contains scores

for(i in 1:length(for_nomr)) { #for every column in for_nomr dataframe (actually is for every subject for every timepoint)
  for_nomr[,i] <-as.numeric(as.character(for_nomr[,i]))} #convert every score to numeric

#Normalization function
#formula: ((score-mean of column)/sqrt(sum((score-mean of column)^2))) 
normalized=data.Normalization (for_nomr,type="n12",normalization="column")
#other option: a2=data.Normalization (for_nomr,type="n12",normalization="row")
#we can not do it by row, beacause we get warnings

for(i in 1:length(normalized)) { #for every column in normalized dataframe (actually is for every subject for every timepoint)
  normalized[,i] <-as.factor(normalized[,i])} #convert every score to factor

first_two_rows <- data[1:2,] #get the Subject and Timepoint rows from data
new_data <- rbind(first_two_rows, normalized) #make them one dataframe

###Friedman's Test###
significant_meta=data.frame(matrix(nrow = nrow(new_data) , ncol = 2)) #empty dataframe to save the significant metabolites with their pvalues
for (m in 3:nrow(new_data)) { #for every metabolite for the whole datafile of normalized metabolites
  #start counting from 3rd row, because in 1st row are the Subjects and in 2nd are the Timepoints
  idx=m #idx is the actual number of row in new_data dataframe
  Score <- idx
  metabolite = data.frame(y=new_data[Score,]) #create a new dataframe for the scores of each metabolite (of each row)
  
  metabolite[nrow(metabolite)+1,] = new_data[1,] # +Subject row
  metabolite[nrow(metabolite)+1,] = new_data[2,] # +Timepoint row
  
  metabolite_rev=t(metabolite) #reverse the dataframe, in order to have them in columns
  
  # t() command creates table
  metabolite_rev=as.data.frame(metabolite_rev) #make it dataframe
  
  names(metabolite_rev)[1] <- "Score" #the colname of the score column was the name of the metabolite. So, in order to be easier to work with, we rename the colname as "Score"
  
  metabolite_rev$Score <- as.numeric(as.character(metabolite_rev$Score)) #convert the Score column to numeric
  
  #FRIEDMAN model
  model = friedman.test(Score~Timepoint | Subject, data = metabolite_rev)
  
  #The summary(model) is a list of elements. In order to be able to pick one particular element, we need to unlist them
  sum_model = unlist(model) 
  
  if (sum_model["p.value"]<0.01) { # check if this particular metabolite significant p-value
    significant_meta[m,1] <- idx # keep the idx number of it (the row number in the new_data dataframe)
    significant_meta[m,2] <- sum_model["p.value"] # keep the corresponding p-value
  }
}

#significant_meta has all idxs of the significant metabolites and for no significant has NAs
new_sign <- na.omit(significant_meta, cols="X2") #remove the rows with NAs (the idxs of no significant metabolites)
#new_sign has the idxs only of the significant metabolites

all_significant=data.frame(matrix(nrow = nrow(new_sign) , ncol = 2)) #empty dataframe to save the actual names of significant metabolites with their pvalues
for (s in 1:nrow(new_sign)) { #for every row in new_sign dataframe (for every significant metabolite)
  idx2=new_sign[s,1] #idx2 is the number of row in new_sign dataframe
  sign_metabolite = data.frame(y=new_data[idx2,]) #create a new dataframe for each significant metabolite (for every row in new_data dataframe and each row is choosen by the idx2 number)
  all_significant[s,1] <- row.names(sign_metabolite) #keep the name of the significant metabolite
  all_significant[s,2] <- new_sign[s,2]} #keep its corresponding pvalue

colnames(all_significant)[1] <- "Sign_metabolite" #change the name of the 1st column
colnames(all_significant)[2] <- "p-value" #change the name of the 2nd column
write.csv(all_significant, file = "Friedman.Significant.csv", row.names = FALSE) #make a .csv file with the significant metabolites and their pvalues
  
###Graphs###
for (o in 1:nrow(new_sign)) { #for every significant metabolite (for every row in new_sign dataframe ---- remember the new_sign dataframe contains the idx for each metabolite (the number of row in data dataframe))
  #we need to get the original values of each metabolite, before the normalization
  idx3=new_sign[o,1] #idx3 is the number of row in new_sign dataframe
  metabolite_original = data.frame(y=data[idx3,]) #create a new dataframe for each significant metabolite (for every row in data dataframe and each row is choosen by the idx3 number)
  
  metabolite_original[nrow(metabolite_original)+1,] = data[1,] # +Subject row
  metabolite_original[nrow(metabolite_original)+1,] = data[2,] # +Timepoint row
  
  #reverse the dataframe, in order to have them in columns
  metabolite_original_rev=t(metabolite_original)
  
  # t() command creates table
  metabolite_original_rev=as.data.frame(metabolite_original_rev) #make it dataframe
  
  names(metabolite_original_rev)[1] <- "Score" #the colname of the score column was the name of the metabolite. So, in order to be easier to work with, we rename the colname as "Score"
  
  metabolite_original_rev$Score <- as.numeric(as.character(metabolite_original_rev$Score)) #convert the Score column to numeric
  
  #We need to transform the scores because we have to big values. Otherwise, the graphs do not seem to be informative.
  #Log transformation of score (use log10 to transform the values)
  transformations <- as.data.frame(log10(metabolite_original_rev$Score))
  #Bind Subject and Timepoint rows
  first_two_rows_again <- as.data.frame(t(metabolite_original[2:3,])) #pick the Subject and Timepoint rows
  trans_data <- cbind(first_two_rows_again, transformations) #make them one dataframe
  names(trans_data)[3] <- "LogScore" #rename the column with the log10 scores
  
  subject_num <- length(levels(trans_data$Subject)) #number of different levels of Subject factor, here is 38
  timepoints_num <- length(levels(trans_data$Timepoint)) #number of different levels of Timepoint factor, here is 6
  
  newdir <- paste0(row.names(metabolite_original)[1]) #pick the name of the metabolite
  dir.create(newdir) #and create a new directory with that name
  cwd <- getwd() #keep "in mind" your MAIN working directory
  setwd(newdir) # "go inside" the directory that you just created and make it your temporary working directory
  mypath <- file.path(getwd() ,paste("Bp_",row.names(metabolite_original)[1], ".png", sep = "")) #path to save the graph
  png(file=mypath, width=862, height=392) #format of saving graph
  mytitle = paste(row.names(metabolite_original)[1],": LogScore distribution across all timepoints and subjects") #main title of the graph
  ytitle = paste(row.names(metabolite_original)[1]) #label of y-axis
  #make the graph
  boxplot(trans_data$LogScore, 
          horizontal = TRUE, 
          main= mytitle,
          las=2,
          xlab= "Log10_transformed_intensity_score", ylab= ytitle,
          col = "blue")
  dev.off() #necessary command to be able to open the .png file that we just created
  
  mypath2 <- file.path(getwd() ,paste("Bp_timepoints_", row.names(metabolite_original)[1], ".png", sep = "")) #path to save the graph
  png(file=mypath2, width=862, height=392) #format of saving graph
  mytitle2 = paste(row.names(metabolite_original)[1],": Boxplot comparing LogScores of six Timepoints") #main title of the graph
  #make the graph
  boxplot(trans_data$LogScore~trans_data$Timepoint, 
          main=mytitle2, 
          col= rainbow(timepoints_num), 
          las=2,
          par(mar = c(4, 8, 1, 1)+ 0.4),
          xlab= "Log10_transformed_intensity_score",
          horizontal = TRUE)
  dev.off() #necessary command to be able to open the .png file that we just created
  
  mypath3 <- file.path(getwd() ,paste("Bp_subjects_timepoints_", row.names(metabolite_original)[1], ".png", sep = "")) #path to save the graph
  png(file=mypath3, width=862, height=1000) #format of saving graph
  mytitle3 = paste( row.names(metabolite_original)[1],": Boxplot comparing LogScores of all subjects for all timepoints") #main title of the graph
  #make the graph
  boxplot(trans_data$LogScore~trans_data$Subject, 
          main=mytitle3, 
          xlab= "Log10_transformed_intensity_score",
          col= rainbow(subject_num), 
          horizontal = TRUE,
          las=2)
  dev.off() #necessary command to be able to open the .png file that we just created
  
  mypath4 <- file.path(getwd() ,paste("Mean_",  row.names(metabolite_original)[1], ".png", sep = "")) #path to save the graph
  png(file=mypath4, width=862, height=392) #format of saving graph
  mytitle4 = paste( row.names(metabolite_original)[1],": Mean LogScore across Timepoints") #main title of the graph
  #make the graph
  plotmeans(trans_data$LogScore~trans_data$Timepoint,
            digits=2, 
            ccol="red", 
            mean.labels=T, 
            main=mytitle4,
            xlab= "Timepoints",
            ylab= "Mean_Log10_transformed_intensity_score")
  dev.off() #necessary command to be able to open the .png file that we just created
  
  means<- round(tapply(trans_data$LogScore, trans_data$Timepoint, mean), digits=2) #mean for each timepoint
  mypath5 <- file.path(getwd() ,paste("Bp_timepoints_meandotted_",  row.names(metabolite_original)[1], ".png", sep = "")) #path to save the graph
  png(file=mypath5, width=862, height=392) #format of saving graph
  mytitle5 = paste(row.names(metabolite_original)[1],": Boxplot comparing LogScores of six Timepoints (blackdotted mean)") #main title of the graph
  #make the graph
  boxplot(trans_data$LogScore ~ trans_data$Timepoint, 
          main=mytitle5,
          xlab="Timepoints", ylab="Log10_transformed_intensity_score", col=rainbow(timepoints_num))
  points(means, col="black", pch=16, cex=1.5) #add dots with means
  dev.off() #necessary command to be able to open the .png file that we just created
  
  mypath6 <- file.path(getwd() ,paste("Vp_", row.names(metabolite_original)[1], ".png", sep = "")) #path to save the graph
  png(file=mypath6, width=862, height=900) #format of saving graph
  title_violin1 = paste(row.names(metabolite_original)[1],": Violin plots comparing LogScores of six Timepoints") #main title of the graph
  #make the graph
  print(ggplot(trans_data, aes(x=trans_data$Timepoint, y=trans_data$LogScore)) + 
          geom_violin() +
          labs(title=title_violin1,x="Timepoints", y = "Log10_transformed_intensity_score")+
          coord_flip()+
          theme_classic())
  dev.off() #necessary command to be able to open the .png file that we just created
  
  mypath7 <- file.path(getwd() ,paste("Vp_withbp_", row.names(metabolite_original)[1], ".png", sep = "")) #path to save the graph
  png(file=mypath7, width=862, height=900) #format of saving graph
  title_violin2 = paste(row.names(metabolite_original)[1],": Violin plots (with boxplots) comparing LogScores of six Timepoints") #main title of the graph
  #make the graph
  print(ggplot(trans_data, aes(x=trans_data$Timepoint, y=trans_data$LogScore)) +
          geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
          labs(title=title_violin2,x="Timepoints", y = "Log10_transformed_intensity_score")+
          coord_flip()+
          geom_boxplot(width=0.1) + theme_classic())
  dev.off() #necessary command to be able to open the .png file that we just created
  
  mypath8 <- file.path(getwd() ,paste("Vp_withbp_coloured_", row.names(metabolite_original)[1], ".png", sep = "")) #path to save the graph
  png(file=mypath8, width=862, height=900) #format of saving graph
  #make the graph
  print(ggplot(trans_data, aes(x=trans_data$Timepoint, y=trans_data$LogScore, fill=trans_data$Timepoint)) + 
          geom_violin(trim=FALSE, show.legend=FALSE)+
          geom_boxplot(width=0.1, fill="white")+
          labs(title=title_violin2,x="Timepoints", y = "Log10_transformed_intensity_score")+ 
          coord_flip()+
          scale_fill_brewer(palette="Dark2") + theme_minimal())
  dev.off() #necessary command to be able to open the .png file that we just created
  
  setwd(cwd) #return to the directory that you kept "in mind" and make it the working directory again
  #Actually, it returns to the MAIN working directory which is the "Friedman_0.01" in order to make a new directory for the next metabolite
}

