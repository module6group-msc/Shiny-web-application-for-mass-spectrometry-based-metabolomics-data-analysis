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

#This will create a new directory in your working directory with the name "Anova"
#At the end, everithing is going to be inside Anova directory
anova <- paste("Anova_with_post_hoc") #set the name
dir.create(anova) #create directory
setwd(anova) #set the new directory as working directory
main <- getwd()
###Normalization###
length_data <- nrow(data) #number of rows in data
for_nomr=data[3:length_data,] #colect only the region of data that contains scores

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
write.csv(new_data, file = "Normalized.csv", row.names = FALSE) #make a .csv file with the normalized values of raw data

###Repeated measures ANOVA###
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
  
  #ANOVA model
  model = aov(Score~Timepoint+Error(Subject/Timepoint),metabolite_rev)
  
  #The summary(model) is a list of elements. In order to be able to pick one particular element, we need to unlist them
  sum_model = unlist(summary(model)) 
  
  if (sum_model["Error: Subject:Timepoint.Pr(>F)1"]<0.01) { #check if the metabolite that is under investgation has significant p-value
    significant_meta[m,1] <- idx #keep the idx number of it (the row number in the new_data dataframe)
    significant_meta[m,2] <- sum_model["Error: Subject:Timepoint.Pr(>F)1"] #keep the corresponding pvalue
  }
}

#significant_meta has all idxs of the significant metabolites and for no significant has NAs
new_sign <- na.omit(significant_meta) #remove the rows with NAs (the idxs of no significant metabolites)
#new_sign has the idxs only of the significant metabolites

all_significant=data.frame(matrix(nrow = nrow(new_sign) , ncol = 2)) #empty dataframe to save the actual names of significant metabolites with their pvalues
for (s in 1:nrow(new_sign)) { #for every row in new_sign dataframe (for every significant metabolite)
  idx2=new_sign[s,1] #idx2 is the number of row in new_sign dataframe
  sign_metabolite = data.frame(y=new_data[idx2,]) #create a new dataframe for each significant metabolite (for every row in new_data dataframe and each row is choosen by the idx2 number)
  all_significant[s,1] <- row.names(sign_metabolite) #keep the name of the significant metabolite
  all_significant[s,2] <- new_sign[s,2]} #keep its corresponding pvalue

colnames(all_significant)[1] <- "Sign_metabolite" #change the name of the 1st column
colnames(all_significant)[2] <- "p-value" #change the name of the 2nd column
write.csv(all_significant, file = "Significant.csv", row.names = FALSE) #make a .csv file with the significant metabolites and their pvalues

#run a post hoc to get the combinations
model_fph <- aov(Score~Timepoint+Subject,metabolite_rev)
post_hoc <-TukeyHSD(model_fph,"Timepoint",conf.level=0.99, data=metabolite_rev)
#post_hoc_anova <- unlist(post_hoc)
#list with all possible combinations
list <- post_hoc$Timepoint
#dataframe with all possible combinations
list_df <- as.data.frame(list)
#how many combinations we have
nrow(list_df)
#all possible combinations
all_comb <- row.names(list_df)
all_comb <- as.matrix(all_comb)
#create a dataframe with all possible combinations (with valid names)
all_comb_valid_names=data.frame(matrix(nrow = nrow(all_comb) , ncol = 1))
for (n in 1:nrow(all_comb)) {
  validn <- gsub("-", "_", all_comb[n,1])
  all_comb_valid_names[n,1] <- validn
  names(all_comb_valid_names)[1] <- "Combinations"}

#This will create a new directory in your working directory with the name "Post_hoc_test"
#At the end, everithing is going to be inside Anova directory
post_hoc_test <- paste("Post_hoc_test") #set the name
dir.create(post_hoc_test) #create directory
setwd(post_hoc_test) #set the new directory as working directory

#create .csv file for each combination
for (i in 1:nrow(all_comb_valid_names)) {
  write.csv(all_comb_valid_names$Combinations[i], file=paste0(all_comb_valid_names$Combinations[i],".csv"), quote = FALSE)
}


#Post Hoc test
for (fph in 1:nrow(new_sign)) { #for every significant metabolite (for every row in new_sign dataframe ---- remember the new_sign dataframe contains the idx for each metabolite (the number of row in data dataframe))
  #we need to get the original values of each metabolite, before the normalization
  idx3=new_sign[fph,1] #idx3 is the number of row in new_sign dataframe
  metabolite_fph = data.frame(y=new_data[idx3,]) #create a new dataframe for each significant metabolite (for every row in data dataframe and each row is choosen by the idx3 number)
  
  metabolite_fph[nrow(metabolite_fph)+1,] = new_data[1,] # +Subject row
  metabolite_fph[nrow(metabolite_fph)+1,] = new_data[2,] # +Timepoint row
  
  #reverse the dataframe, in order to have them in columns
  metabolite_fph_rev=t(metabolite_fph)
  
  # t() command creates table
  metabolite_fph_rev=as.data.frame(metabolite_fph_rev) #make it dataframe
  
  names(metabolite_fph_rev)[1] <- "Score" #the colname of the score column was the name of the metabolite. So, in order to be easier to work with, we rename the colname as "Score"
  
  metabolite_fph_rev$Score <- as.numeric(as.character(metabolite_fph_rev$Score)) #convert the Score column to numeric
#necessery model for post hoc test
model_fph <- aov(Score~Timepoint+Subject,metabolite_fph_rev)
#post hoc test
post_hoc <-TukeyHSD(model_fph,"Timepoint",data=metabolite_fph_rev)
#get the column with the pvalues
df_for_pvalue <- as.data.frame(post_hoc$Timepoint)
padj <- df_for_pvalue[4]
#get the significant combinations
sign_timepoints <- data.frame(matrix(nrow = nrow(padj) , ncol = 2))
for (ph in 1:nrow(padj)) {
  if (padj[ph,1] < 0.05) {
    sign_timepoints[ph,1] <- row.names(padj)[ph]
    sign_timepoints[ph,2] <- padj[ph,1]
  }
#exclude the NAs
significant_timepoints <- na.omit(sign_timepoints)
#only the significant combinations (without the p-value)
comb <- significant_timepoints[,1]
#significant combinations
comb <- as.matrix(comb)
#create a dataframe with significant combinations (with valid names)
comb_valid_names=data.frame(matrix(nrow = nrow(comb) , ncol = 1))
}
for (n in 1:nrow(comb)) {
  validn_comp <- gsub("-", "_", comb[n,1])
  comb_valid_names[n,1] <- validn_comp
  names(comb_valid_names)[1] <- "significant_Combinations"}

#get the list with the .csv files that we created
temp = list.files(pattern="*.csv")
#keep the names without the ".csv"
names_of_the_csv <- data.frame(matrix(nrow = length(temp) , ncol = 1))
for (name in 1:length(temp)) {
  names_of_the_csv[name,1] <- tools::file_path_sans_ext(basename(temp[name]))
}
#get only the first column of sign_timepoints
first_of_sign_timepoints <- comb_valid_names[,1]
#write into corresponding .csv the name of the metabolite
for (csv in 1:(nrow(names_of_the_csv)-1)) {
    if (names_of_the_csv[csv,1]%in%comb_valid_names[,1] == TRUE) {
      write.table( row.names(metabolite_fph)[1],  
                   file=temp[csv], 
                   append = T, 
                   sep=',', 
                   row.names=F, 
                   col.names=F )
    }}

}

setwd(main) #to return to "Anova" folder

#This will create a new directory in your working directory with the name "Plots"
#At the end, everithing is going to be inside Anova directory
plots <- paste("Plots") #set the name
dir.create(plots) #create directory
setwd(plots) #set the new directory as working directory

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
  mytitle = paste(row.names(metabolite_original)[1],": Log10 transformed intensity score distribution across all timepoints and subjects") #main title of the graph
  ytitle = paste(row.names(metabolite_original)[1]) #label of y-axis
  #make the graph
  boxplot(trans_data$LogScore, 
          horizontal = TRUE, 
          main= mytitle,
          las=2,
          xlab= "Log10 transformed intensity score", ylab= ytitle,
          col = "blue")
  dev.off() #necessary command to be able to open the .png file that we just created
  
  mypath2 <- file.path(getwd() ,paste("Bp_timepoints_", row.names(metabolite_original)[1], ".png", sep = "")) #path to save the graph
  png(file=mypath2, width=862, height=392) #format of saving graph
  mytitle2 = paste(row.names(metabolite_original)[1],": Boxplot comparing Log10 transformed intensity scores of six Timepoints") #main title of the graph
  #make the graph
  boxplot(trans_data$LogScore~trans_data$Timepoint, 
          main=mytitle2, 
          col= rainbow(timepoints_num), 
          las=2,
          par(mar = c(4, 8, 1, 1)+ 0.4),
          xlab= "Log10 transformed intensity scores",
          horizontal = TRUE)
  dev.off() #necessary command to be able to open the .png file that we just created
  
  mypath3 <- file.path(getwd() ,paste("Bp_subjects_timepoints_", row.names(metabolite_original)[1], ".png", sep = "")) #path to save the graph
  png(file=mypath3, width=862, height=1000) #format of saving graph
  mytitle3 = paste( row.names(metabolite_original)[1],": Boxplot comparing Log10 transformed intensity scores of all subjects for all timepoints") #main title of the graph
  #make the graph
  boxplot(trans_data$LogScore~trans_data$Subject, 
          main=mytitle3, 
          xlab= "Log10 transformed intensity scores",
          col= rainbow(subject_num), 
          horizontal = TRUE,
          las=2)
  dev.off() #necessary command to be able to open the .png file that we just created
  
  mypath4 <- file.path(getwd() ,paste("Mean_",  row.names(metabolite_original)[1], ".png", sep = "")) #path to save the graph
  png(file=mypath4, width=862, height=392) #format of saving graph
  mytitle4 = paste( row.names(metabolite_original)[1],": Mean Log10 transformed intensity score across Timepoints") #main title of the graph
  #make the graph
  plotmeans(trans_data$LogScore~trans_data$Timepoint,
            digits=2, 
            ccol="red", 
            mean.labels=T, 
            main=mytitle4,
            xlab= "Timepoints",
            ylab= "Mean Log10 transformed intensity score")
  dev.off() #necessary command to be able to open the .png file that we just created
  
  means<- round(tapply(trans_data$LogScore, trans_data$Timepoint, mean), digits=2) #mean for each timepoint
  mypath5 <- file.path(getwd() ,paste("Bp_timepoints_meandotted_",  row.names(metabolite_original)[1], ".png", sep = "")) #path to save the graph
  png(file=mypath5, width=862, height=392) #format of saving graph
  mytitle5 = paste(row.names(metabolite_original)[1],": Boxplot comparing Log10 transformed intensity scores of six Timepoints (blackdotted mean)") #main title of the graph
  #make the graph
  boxplot(trans_data$LogScore ~ trans_data$Timepoint, 
          main=mytitle5,
          xlab="Timepoints", ylab="Log10 transformed intensity scores", col=rainbow(timepoints_num))
  points(means, col="black", pch=16, cex=1.5) #add dots with means
  dev.off() #necessary command to be able to open the .png file that we just created

  mypath8 <- file.path(getwd() ,paste("Vp_withbp_coloured_", row.names(metabolite_original)[1], ".png", sep = "")) #path to save the graph
  png(file=mypath8, width=862, height=900) #format of saving graph
  title_violin2 = paste(row.names(metabolite_original)[1],": Violin plots (with boxplots) comparing Log10 transformed intensity scores of six Timepoints") #main title of the graph
  #make the graph
  print(ggplot(trans_data, aes(x=trans_data$Timepoint, y=trans_data$LogScore, fill=trans_data$Timepoint)) + 
          geom_violin(trim=FALSE, show.legend=FALSE)+
          geom_boxplot(width=0.1, fill="white")+
          labs(title=title_violin2,x="Timepoints", y = "Log10 transformed intensity scores")+ 
          coord_flip()+
          scale_fill_brewer(palette="Dark2") + theme_minimal())
  dev.off() #necessary command to be able to open the .png file that we just created
  
  setwd(cwd) #return to the directory that you kept "in mind" and make it the working directory again
  #Actually, it returns to the MAIN working directory which is the "Anova" in order to make a new directory for the next metabolite
}
