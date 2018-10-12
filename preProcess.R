library(igraph)

#Takes in a dataframe to split into 2 partitions
sampleSplit <- function(df, divide) {
  #@param df: A dataframe to be randomly divided in two
  #@param divide: The ratio of that divide (ie. 50/50, 80/20, ect)
  #@return: A list of indexes that will point to random rows in a determined portion of df 
  
  partition <- floor(divide*nrow(df))
  index <- sample(nrow(df), size=partition)
  return (index)
}

#Expecting the output from a tn93 run formatted to a csv file. with 2 ID columns and a pairwise distance
#The ID columns should be in an ID_Date format.
args = commandArgs(trailingOnly = T)
input <- read.csv(args[1], stringsAsFactors = F)

#Iterates through tn93 output to create columns specifying collection year and id as their own columns
for(row in 1:nrow(input)) {
  splitID1 <- strsplit(input[row,"ID1"], "_")[[1]]
  splitID2 <- strsplit(input[row,"ID2"], "_")[[1]]
  input[row, "ID1"] <- splitID1[1] 
  input[row, "Date1"] <- as.integer(splitID1[2])
  input[row, "ID2"] <- splitID2[1]
  input[row, "Date2"] <- as.integer(splitID2[2])
  
  #INQUIRY: Sot sure if there is a faster way to censor the edgelist for future years
  input[row, "maxDate"] <- as.integer(max(input[row, "Date2"],  input[row, "Date1"]))
}

#Initializing 2 data frames to point to 2 partitions of the data. 
train <- data.frame(ID1 = character(), ID2 = character(), Date1 = integer(), Date2=integer(), minDate=integer())
test <- data.frame(ID1 = character(), ID2 = character(), Date1 = integer(), Date2=integer(), minDate=integer()) 

#Populates the dataframes with even splits of each year, iterating through subsets of the input defined by the lowest collection date
year = min(input$maxDate)

while (year<=max(input$maxDate)){
  
  #The inputted data frame excluding all cases beyond a current year and all cases already inputted copied into the test or train partitions
  inputAtYear <- input[input$maxDate==year, ] 
  index <- sampleSplit(inputAtYear, 0.5)
  
  #adding the current year's data into the partitions
  train <- rbind(inputAtYear[index, ],train)
  test <- rbind(inputAtYear[-index, ], test)
  
  #graph <- graph_from_data_frame(train, T, NULL)
  #graph <- graph_from_data_frame(test, T, NULL)
  
  year=year+1
}

#Optional data and output testing
if (is.na(args[2])==F){
  str(train)
  str(test)
  str(input)
  pdf("dataSummary.pdf")
  plot(hist(input$maxDate))
  plot(graph_from_data_frame(input, T, NULL))
  dev.off() 
}