#USAGE: Rscript preProcess.R tn93output.csv year

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
year = as.integer(args[2])

#Alters tn93 output to create columns specifying collection year and id as their own columns
temp <- sapply(input$ID1, function(x) strsplit(x, '_')[[1]])
input$ID1 <- temp[1,]
input$Date1 <- temp[2,]

temp <- sapply(input$ID2, function(x) strsplit(x, '_')[[1]])
input$ID2 <- temp[1,]
input$Date2 <- temp[2,]

input$maxDate <- apply(cbind(input$Date1, input$Date2), 1, max)

#The inputted data frame excluding all cases beyond a current year and all cases already inputted copied into the test or train partitions
inputToYear <- input[input$maxDate<=year, ] 
index <- sampleSplit(inputToYear, 0.5)

#adding the current year's data into the partitions
train <-inputToYear[index, ]
test <- inputToYear[-index, ]