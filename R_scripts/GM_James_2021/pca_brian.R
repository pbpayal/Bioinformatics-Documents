# Second analysis:  2018-07-13
# Analysis params
# Remove B11-68 (duplicate) and B11-69, B15-149, B15-151
# WT were all grouped as one (not segregated by age)
# Ancova corrected by both GENDER and BIOPSY_SITE
# Noise cut-off = 3


library(stats)
library(gplots)
library(rgl)
library(multcomp)
library(effects)
library(car)
library(compute.es)
library(ggplot2)
library(pastecs)
library(stringr)

###############
#### Setup ####
###############


rm(list=ls())
options(object.size=Inf)

# Set current working directory
# setwd("/Users/uapinyoyingpb/Notebook/Archive/2018/lama2_project/2018-05-14_tuxedo2_analysis/s6_assembled_transcripts/gtfs/extracted_tpms")
setwd("/Users/uapinyoyingpb/Notebook/Archive/2018/lama2_project/2018-05-14_tuxedo2_analysis/s6_assembled_transcripts/expr_for_rsem")


##########################
#### Custom functions ####
##########################

image.scale <- function (z, col, x, y = NULL, size = NULL, digits = 3, labels = c("breaks", 
                                                                                  "ranges")) 
{ 
  # sort out the location 
  n <- length(col) 
  usr <- par("usr") 
  mx <- mean(usr[1:2]); my <- mean(usr[3:4]) 
  dx <- diff(usr[1:2]); dy <- diff(usr[3:4]) 
  if (missing(x)) 
    x <- mx + 1.05*dx/2 # default x to right of image 
  else if (is.list(x)) { 
    if (length(x$x) == 2) 
      size <- c(diff(x$x), -diff(x$y)/n) 
    y <- x$y[1] 
    x <- x$x[1] 
  } else x <- x[1] 
  if (is.null(size)) 
    if (is.null(y)) { 
      size <- 0.618*dy/n # default size, golden ratio 
      y <- my + 0.618*dy/2 # default y to give centred scale 
    } else size <- (y-my)*2/n 
  if (length(size)==1) 
    size <- rep(size, 2) # default square boxes 
  if (is.null(y)) 
    y <- my + n*size[2]/2 
  # draw the image scale 
  i <- seq(along = col) 
  rect(x, y - i * size[2], x + size[1], y - (i - 1) * size[2], 
       col = rev(col), xpd = TRUE) 
  # sort out the labels 
  rng <- range(z, na.rm = TRUE) 
  bks <- seq(from = rng[2], to = rng[1], length = n + 1) 
  bks <- formatC(bks, format="f", digits=digits) 
  labels <- match.arg(labels) 
  if (labels == "breaks") 
    ypts <- y - c(0, i) * size[2] 
  else { 
    bks <- paste(bks[-1], bks[-(n+1)], sep = " - ") 
    ypts <- y - (i - 0.5) * size[2] 
  } 
  text(x = x + 1.2 * size[1], y = ypts, labels = bks, adj = 
         ifelse(size[1]>0, 0, 1), xpd = TRUE) 
} 

#####################
#### import data ####
#####################

# Data is from stringtie's special output for expression tables meant for RSEM
# dont use the GTF files output for ballgown

LAMA2_15_148_expr <- read.table("LAMA2-15-148_expression.txt", skip=1)
# LAMA2_15_149_expr <- read.table("LAMA2-15-149_expression.txt", skip=1) # biopsy looks mostly normal?
# LAMA2_15_151_expr <- read.table("LAMA2-15-151_expression.txt", skip=1) # Partial merosin
LAMA2_3935_expr <- read.table("LAMA2-3935_expression.txt", skip=1)
LAMA2_3946_expr <- read.table("LAMA2-3946_expression.txt", skip=1)
LAMA2_4002_expr <- read.table("LAMA2-4002_expression.txt", skip=1)
LAMA2_4007_expr <- read.table("LAMA2-4007_expression.txt", skip=1)
LAMA2_4010_expr <- read.table("LAMA2-4010_expression.txt", skip=1)
LAMA2_4027_expr <- read.table("LAMA2-4027_expression.txt", skip=1)
LAMA2_4043_expr <- read.table("LAMA2-4043_expression.txt", skip=1)
LAMA2_4120_expr <- read.table("LAMA2-4120_expression.txt", skip=1)
LAMA2_4244_expr <- read.table("LAMA2-4244_expression.txt", skip=1)
LAMA2_4363_expr <- read.table("LAMA2-4363_expression.txt", skip=1)
LAMA2_4400_expr <- read.table("LAMA2-4400_expression.txt", skip=1)
LAMA2_4414_expr <- read.table("LAMA2-4414_expression.txt", skip=1)
LAMA2_4422_expr <- read.table("LAMA2-4422_expression.txt", skip=1)
# LAMA2_B11_68_expr <- read.table("LAMA2-B11-68_expression.txt", skip=1) # repeated LAMA2_15_148
# LAMA2_B11_69_expr <- read.table("LAMA2-B11-69_expression.txt", skip=1) # outlier - unknown reason
LAMA2_B14_14_expr <- read.table("LAMA2-B14-14_expression.txt", skip=1)
LAMA2_B16_21_expr <- read.table("LAMA2-B16-21_expression.txt", skip=1)
M13_14_expr <- read.table("M13-14_expression.txt", skip=1)
M13_24_expr <- read.table("M13-24_expression.txt", skip=1)
M13_2_expr <- read.table("M13-2_expression.txt", skip=1)
M13_32_expr <- read.table("M13-32_expression.txt", skip=1)
M13_37_expr <- read.table("M13-37_expression.txt", skip=1)
M13_38_expr <- read.table("M13-38_expression.txt", skip=1)
M13_42_expr <- read.table("M13-42_expression.txt", skip=1)
M13_49_expr <- read.table("M13-49_expression.txt", skip=1)
M13_8_expr <- read.table("M13-8_expression.txt", skip=1)

#############################
#### Assign column names ####
#############################

nnn <- c("Gene","X","Chr","Strand","Start","Stop","Coverage","FPKM","TPM")

dimnames(LAMA2_15_148_expr)[[2]] <- nnn
# dimnames(LAMA2_15_149_expr)[[2]] <- nnn
# dimnames(LAMA2_15_151_expr)[[2]] <- nnn
dimnames(LAMA2_3935_expr)[[2]] <- nnn
dimnames(LAMA2_3946_expr)[[2]] <- nnn
dimnames(LAMA2_4002_expr)[[2]] <- nnn
dimnames(LAMA2_4007_expr)[[2]] <- nnn
dimnames(LAMA2_4010_expr)[[2]] <- nnn
dimnames(LAMA2_4027_expr)[[2]] <- nnn
dimnames(LAMA2_4043_expr)[[2]] <- nnn
dimnames(LAMA2_4120_expr)[[2]] <- nnn
dimnames(LAMA2_4244_expr)[[2]] <- nnn
dimnames(LAMA2_4363_expr)[[2]] <- nnn
dimnames(LAMA2_4400_expr)[[2]] <- nnn
dimnames(LAMA2_4414_expr)[[2]] <- nnn
dimnames(LAMA2_4422_expr)[[2]] <- nnn
# dimnames(LAMA2_B11_68_expr)[[2]] <- nnn
# dimnames(LAMA2_B11_69_expr)[[2]] <- nnn
dimnames(LAMA2_B14_14_expr)[[2]] <- nnn
dimnames(LAMA2_B16_21_expr)[[2]] <- nnn
dimnames(M13_14_expr)[[2]] <- nnn
dimnames(M13_24_expr)[[2]] <- nnn
dimnames(M13_2_expr)[[2]] <- nnn
dimnames(M13_32_expr)[[2]] <- nnn
dimnames(M13_37_expr)[[2]] <- nnn
dimnames(M13_38_expr)[[2]] <- nnn
dimnames(M13_42_expr)[[2]] <- nnn
dimnames(M13_49_expr)[[2]] <- nnn
dimnames(M13_8_expr)[[2]] <- nnn


#####################################################################################
### Assigning unique identifiers as rownames by constructing one out of bed info ####
#####################################################################################

temp1 <- LAMA2_15_148_expr # copy the current data.table
temp1 <- temp1[,c(1,3,5,6,4)] # rearrage columns
temp1[,5] <- ifelse(temp1[,5]=="-","Minus","Plus") # change - and + strand to words
temp1 <- apply(temp1,1,paste,collapse="_") # Turn the table into a single column of items and delimit by '_' (e.g. col1_col2_col3)
temp1 <- str_replace_all(temp1, " ", "") # remove empty spaces
dimnames(LAMA2_15_148_expr)[[1]] <- temp1 # now put this newly constructed list of row names 

# temp1 <- LAMA2_15_149_expr
# temp1 <- temp1[,c(1,3,5,6,4)]
# temp1[,5] <- ifelse(temp1[,5]=="-","Minus","Plus")
# temp1 <- apply(temp1,1,paste,collapse="_")
# temp1 <- str_replace_all(temp1, " ", "")
# dimnames(LAMA2_15_149_expr)[[1]] <- temp1

# temp1 <- LAMA2_15_151_expr
# temp1 <- temp1[,c(1,3,5,6,4)]
# temp1[,5] <- ifelse(temp1[,5]=="-","Minus","Plus")
# temp1 <- apply(temp1,1,paste,collapse="_")
# temp1 <- str_replace_all(temp1, " ", "")
# dimnames(LAMA2_15_151_expr)[[1]] <- temp1

temp1 <- LAMA2_3935_expr
temp1 <- temp1[,c(1,3,5,6,4)]
temp1[,5] <- ifelse(temp1[,5]=="-","Minus","Plus")
temp1 <- apply(temp1,1,paste,collapse="_")
temp1 <- str_replace_all(temp1, " ", "")
dimnames(LAMA2_3935_expr)[[1]] <- temp1

temp1 <- LAMA2_3946_expr
temp1 <- temp1[,c(1,3,5,6,4)]
temp1[,5] <- ifelse(temp1[,5]=="-","Minus","Plus")
temp1 <- apply(temp1,1,paste,collapse="_")
temp1 <- str_replace_all(temp1, " ", "")
dimnames(LAMA2_3946_expr)[[1]] <- temp1

temp1 <- LAMA2_4002_expr
temp1 <- temp1[,c(1,3,5,6,4)]
temp1[,5] <- ifelse(temp1[,5]=="-","Minus","Plus")
temp1 <- apply(temp1,1,paste,collapse="_")
temp1 <- str_replace_all(temp1, " ", "")
dimnames(LAMA2_4002_expr)[[1]] <- temp1

temp1 <- LAMA2_4007_expr
temp1 <- temp1[,c(1,3,5,6,4)]
temp1[,5] <- ifelse(temp1[,5]=="-","Minus","Plus")
temp1 <- apply(temp1,1,paste,collapse="_")
temp1 <- str_replace_all(temp1, " ", "")
dimnames(LAMA2_4007_expr)[[1]] <- temp1

temp1 <- LAMA2_4010_expr
temp1 <- temp1[,c(1,3,5,6,4)]
temp1[,5] <- ifelse(temp1[,5]=="-","Minus","Plus")
temp1 <- apply(temp1,1,paste,collapse="_")
temp1 <- str_replace_all(temp1, " ", "")
dimnames(LAMA2_4010_expr)[[1]] <- temp1

temp1 <- LAMA2_4027_expr
temp1 <- temp1[,c(1,3,5,6,4)]
temp1[,5] <- ifelse(temp1[,5]=="-","Minus","Plus")
temp1 <- apply(temp1,1,paste,collapse="_")
temp1 <- str_replace_all(temp1, " ", "")
dimnames(LAMA2_4027_expr)[[1]] <- temp1

temp1 <- LAMA2_4043_expr
temp1 <- temp1[,c(1,3,5,6,4)]
temp1[,5] <- ifelse(temp1[,5]=="-","Minus","Plus")
temp1 <- apply(temp1,1,paste,collapse="_")
temp1 <- str_replace_all(temp1, " ", "")
dimnames(LAMA2_4043_expr)[[1]] <- temp1

temp1 <- LAMA2_4120_expr
temp1 <- temp1[,c(1,3,5,6,4)]
temp1[,5] <- ifelse(temp1[,5]=="-","Minus","Plus")
temp1 <- apply(temp1,1,paste,collapse="_")
temp1 <- str_replace_all(temp1, " ", "")
dimnames(LAMA2_4120_expr)[[1]] <- temp1

temp1 <- LAMA2_4244_expr
temp1 <- temp1[,c(1,3,5,6,4)]
temp1[,5] <- ifelse(temp1[,5]=="-","Minus","Plus")
temp1 <- apply(temp1,1,paste,collapse="_")
temp1 <- str_replace_all(temp1, " ", "")
dimnames(LAMA2_4244_expr)[[1]] <- temp1

temp1 <- LAMA2_4363_expr
temp1 <- temp1[,c(1,3,5,6,4)]
temp1[,5] <- ifelse(temp1[,5]=="-","Minus","Plus")
temp1 <- apply(temp1,1,paste,collapse="_")
temp1 <- str_replace_all(temp1, " ", "")
dimnames(LAMA2_4363_expr)[[1]] <- temp1

temp1 <- LAMA2_4400_expr
temp1 <- temp1[,c(1,3,5,6,4)]
temp1[,5] <- ifelse(temp1[,5]=="-","Minus","Plus")
temp1 <- apply(temp1,1,paste,collapse="_")
temp1 <- str_replace_all(temp1, " ", "")
dimnames(LAMA2_4400_expr)[[1]] <- temp1

temp1 <- LAMA2_4414_expr
temp1 <- temp1[,c(1,3,5,6,4)]
temp1[,5] <- ifelse(temp1[,5]=="-","Minus","Plus")
temp1 <- apply(temp1,1,paste,collapse="_")
temp1 <- str_replace_all(temp1, " ", "")
dimnames(LAMA2_4414_expr)[[1]] <- temp1

temp1 <- LAMA2_4422_expr
temp1 <- temp1[,c(1,3,5,6,4)]
temp1[,5] <- ifelse(temp1[,5]=="-","Minus","Plus")
temp1 <- apply(temp1,1,paste,collapse="_")
temp1 <- str_replace_all(temp1, " ", "")
dimnames(LAMA2_4422_expr)[[1]] <- temp1

# temp1 <- LAMA2_B11_68_expr
# temp1 <- temp1[,c(1,3,5,6,4)]
# temp1[,5] <- ifelse(temp1[,5]=="-","Minus","Plus")
# temp1 <- apply(temp1,1,paste,collapse="_")
# temp1 <- str_replace_all(temp1, " ", "")
# dimnames(LAMA2_B11_68_expr)[[1]] <- temp1

# temp1 <- LAMA2_B11_69_expr
# temp1 <- temp1[,c(1,3,5,6,4)]
# temp1[,5] <- ifelse(temp1[,5]=="-","Minus","Plus")
# temp1 <- apply(temp1,1,paste,collapse="_")
# temp1 <- str_replace_all(temp1, " ", "")
# dimnames(LAMA2_B11_69_expr)[[1]] <- temp1

temp1 <- LAMA2_B14_14_expr
temp1 <- temp1[,c(1,3,5,6,4)]
temp1[,5] <- ifelse(temp1[,5]=="-","Minus","Plus")
temp1 <- apply(temp1,1,paste,collapse="_")
temp1 <- str_replace_all(temp1, " ", "")
dimnames(LAMA2_B14_14_expr)[[1]] <- temp1

temp1 <- LAMA2_B16_21_expr
temp1 <- temp1[,c(1,3,5,6,4)]
temp1[,5] <- ifelse(temp1[,5]=="-","Minus","Plus")
temp1 <- apply(temp1,1,paste,collapse="_")
temp1 <- str_replace_all(temp1, " ", "")
dimnames(LAMA2_B16_21_expr)[[1]] <- temp1

temp1 <- M13_14_expr
temp1 <- temp1[,c(1,3,5,6,4)]
temp1[,5] <- ifelse(temp1[,5]=="-","Minus","Plus")
temp1 <- apply(temp1,1,paste,collapse="_")
temp1 <- str_replace_all(temp1, " ", "")
dimnames(M13_14_expr)[[1]] <- temp1

temp1 <- M13_24_expr
temp1 <- temp1[,c(1,3,5,6,4)]
temp1[,5] <- ifelse(temp1[,5]=="-","Minus","Plus")
temp1 <- apply(temp1,1,paste,collapse="_")
temp1 <- str_replace_all(temp1, " ", "")
dimnames(M13_24_expr)[[1]] <- temp1

temp1 <- M13_2_expr
temp1 <- temp1[,c(1,3,5,6,4)]
temp1[,5] <- ifelse(temp1[,5]=="-","Minus","Plus")
temp1 <- apply(temp1,1,paste,collapse="_")
temp1 <- str_replace_all(temp1, " ", "")
dimnames(M13_2_expr)[[1]] <- temp1

temp1 <- M13_32_expr
temp1 <- temp1[,c(1,3,5,6,4)]
temp1[,5] <- ifelse(temp1[,5]=="-","Minus","Plus")
temp1 <- apply(temp1,1,paste,collapse="_")
temp1 <- str_replace_all(temp1, " ", "")
dimnames(M13_32_expr)[[1]] <- temp1

temp1 <- M13_37_expr
temp1 <- temp1[,c(1,3,5,6,4)]
temp1[,5] <- ifelse(temp1[,5]=="-","Minus","Plus")
temp1 <- apply(temp1,1,paste,collapse="_")
temp1 <- str_replace_all(temp1, " ", "")
dimnames(M13_37_expr)[[1]] <- temp1

temp1 <- M13_38_expr
temp1 <- temp1[,c(1,3,5,6,4)]
temp1[,5] <- ifelse(temp1[,5]=="-","Minus","Plus")
temp1 <- apply(temp1,1,paste,collapse="_")
temp1 <- str_replace_all(temp1, " ", "")
dimnames(M13_38_expr)[[1]] <- temp1

temp1 <- M13_42_expr
temp1 <- temp1[,c(1,3,5,6,4)]
temp1[,5] <- ifelse(temp1[,5]=="-","Minus","Plus")
temp1 <- apply(temp1,1,paste,collapse="_")
temp1 <- str_replace_all(temp1, " ", "")
dimnames(M13_42_expr)[[1]] <- temp1

temp1 <- M13_49_expr
temp1 <- temp1[,c(1,3,5,6,4)]
temp1[,5] <- ifelse(temp1[,5]=="-","Minus","Plus")
temp1 <- apply(temp1,1,paste,collapse="_")
temp1 <- str_replace_all(temp1, " ", "")
dimnames(M13_49_expr)[[1]] <- temp1

temp1 <- M13_8_expr
temp1 <- temp1[,c(1,3,5,6,4)]
temp1[,5] <- ifelse(temp1[,5]=="-","Minus","Plus")
temp1 <- apply(temp1,1,paste,collapse="_")
temp1 <- str_replace_all(temp1, " ", "")
dimnames(M13_8_expr)[[1]] <- temp1


#------------------------------------------------------------------------------
# Define Unique Union Row Names
#------------------------------------------------------------------------------

temp1 <- dimnames(LAMA2_15_148_expr)[[1]]

# temp1 <- c(temp1,dimnames(LAMA2_15_149_expr)[[1]])
# temp1 <- unique(sort(temp1))
# temp1 <- c(temp1,dimnames(LAMA2_15_151_expr)[[1]])
# temp1 <- unique(sort(temp1))

temp1 <- c(temp1,dimnames(LAMA2_3935_expr)[[1]])
temp1 <- unique(sort(temp1))
temp1 <- c(temp1,dimnames(LAMA2_3946_expr)[[1]])
temp1 <- unique(sort(temp1))
temp1 <- c(temp1,dimnames(LAMA2_4002_expr)[[1]])
temp1 <- unique(sort(temp1))
temp1 <- c(temp1,dimnames(LAMA2_4007_expr)[[1]])
temp1 <- unique(sort(temp1))
temp1 <- c(temp1,dimnames(LAMA2_4010_expr)[[1]])
temp1 <- unique(sort(temp1))
temp1 <- c(temp1,dimnames(LAMA2_4027_expr)[[1]])
temp1 <- unique(sort(temp1))
temp1 <- c(temp1,dimnames(LAMA2_4043_expr)[[1]])
temp1 <- unique(sort(temp1))
temp1 <- c(temp1,dimnames(LAMA2_4120_expr)[[1]])
temp1 <- unique(sort(temp1))
temp1 <- c(temp1,dimnames(LAMA2_4244_expr)[[1]])
temp1 <- unique(sort(temp1))
temp1 <- c(temp1,dimnames(LAMA2_4363_expr)[[1]])
temp1 <- unique(sort(temp1))
temp1 <- c(temp1,dimnames(LAMA2_4400_expr)[[1]])
temp1 <- unique(sort(temp1))
temp1 <- c(temp1,dimnames(LAMA2_4414_expr)[[1]])
temp1 <- unique(sort(temp1))
temp1 <- c(temp1,dimnames(LAMA2_4422_expr)[[1]])
temp1 <- unique(sort(temp1))

## Omit sample B11_68
# temp1 <- c(temp1,dimnames(LAMA2_B11_68_expr)[[1]])
# temp1 <- unique(sort(temp1))
# temp1 <- c(temp1,dimnames(LAMA2_B11_69_expr)[[1]])
# temp1 <- unique(sort(temp1))

temp1 <- c(temp1,dimnames(LAMA2_B14_14_expr)[[1]])
temp1 <- unique(sort(temp1))
temp1 <- c(temp1,dimnames(LAMA2_B16_21_expr)[[1]])
temp1 <- unique(sort(temp1))
temp1 <- c(temp1,dimnames(M13_14_expr)[[1]])
temp1 <- unique(sort(temp1))
temp1 <- c(temp1,dimnames(M13_24_expr)[[1]])
temp1 <- unique(sort(temp1))
temp1 <- c(temp1,dimnames(M13_2_expr)[[1]])
temp1 <- unique(sort(temp1))
temp1 <- c(temp1,dimnames(M13_32_expr)[[1]])
temp1 <- unique(sort(temp1))
temp1 <- c(temp1,dimnames(M13_37_expr)[[1]])
temp1 <- unique(sort(temp1))
temp1 <- c(temp1,dimnames(M13_38_expr)[[1]])
temp1 <- unique(sort(temp1))
temp1 <- c(temp1,dimnames(M13_42_expr)[[1]])
temp1 <- unique(sort(temp1))
temp1 <- c(temp1,dimnames(M13_49_expr)[[1]])
temp1 <- unique(sort(temp1))
temp1 <- c(temp1,dimnames(M13_8_expr)[[1]])
temp1 <- unique(sort(temp1))

###############################
#### Meld Final Data Table ####
###############################

running.final <- rep(99999,length(temp1)) # add values to a list to be replaced

temp2 <- LAMA2_15_148_expr # copy whole table
temp3 <- dimnames(temp2)[[1]]  # grab just the rownames / unique identifiers
temp4 <- intersect(temp3,temp1) # Determine the common rownames between the total list of transcripts and current sample
temp5 <- temp2[temp4,9] # Populate the TPM values from the sample from current samples according to common row names
names(temp5) <- temp4 # now name the new table of values from the common rownames of current sample
if(length(temp5)<length(temp1)) {  # if current sample rows must be shorter than master list
  temp6 <- setdiff(temp1,temp5) # Get only unique rows only between the two samples, no dups
  temp7 <- rep(0,length(temp6)) # determine length of unique rows and populate with zeros (fill in for transcripts that are not in some samples)
  names(temp7) <- temp6 # take the names fromt he unique rows and put it in the newly populated 0'ed rows
  temp5 <- c(temp5,temp7) # Now merge rows from the unique to current sample transcript list
}
temp5 <- temp5[temp1] # Include all rownames from the master list into the curr sample transcript list
running.final <- cbind(running.final,temp5) # now add it to the growing master table

# temp2 <- LAMA2_15_149_expr
# temp3 <- dimnames(temp2)[[1]]
# temp4 <- intersect(temp3,temp1)
# temp5 <- temp2[temp4,9]
# names(temp5) <- temp4
# if(length(temp5)<length(temp1)) {
#   temp6 <- setdiff(temp1,temp5)
#   temp7 <- rep(0,length(temp6))
#   names(temp7) <- temp6
#   temp5 <- c(temp5,temp7)
# }
# temp5 <- temp5[temp1]
# running.final <- cbind(running.final,temp5)

# temp2 <- LAMA2_15_151_expr
# temp3 <- dimnames(temp2)[[1]]
# temp4 <- intersect(temp3,temp1)
# temp5 <- temp2[temp4,9]
# names(temp5) <- temp4
# if(length(temp5)<length(temp1)) {
#   temp6 <- setdiff(temp1,temp5)
#   temp7 <- rep(0,length(temp6))
#   names(temp7) <- temp6
#   temp5 <- c(temp5,temp7)
# }
# temp5 <- temp5[temp1]
# running.final <- cbind(running.final,temp5)

temp2 <- LAMA2_3935_expr
temp3 <- dimnames(temp2)[[1]]
temp4 <- intersect(temp3,temp1)
temp5 <- temp2[temp4,9]
names(temp5) <- temp4
if(length(temp5)<length(temp1)) {
  temp6 <- setdiff(temp1,temp5)
  temp7 <- rep(0,length(temp6))
  names(temp7) <- temp6
  temp5 <- c(temp5,temp7)
}
temp5 <- temp5[temp1]
running.final <- cbind(running.final,temp5)

temp2 <- LAMA2_3946_expr
temp3 <- dimnames(temp2)[[1]]
temp4 <- intersect(temp3,temp1)
temp5 <- temp2[temp4,9]
names(temp5) <- temp4
if(length(temp5)<length(temp1)) {
  temp6 <- setdiff(temp1,temp5)
  temp7 <- rep(0,length(temp6))
  names(temp7) <- temp6
  temp5 <- c(temp5,temp7)
}
temp5 <- temp5[temp1]
running.final <- cbind(running.final,temp5)

temp2 <- LAMA2_4002_expr
temp3 <- dimnames(temp2)[[1]]
temp4 <- intersect(temp3,temp1)
temp5 <- temp2[temp4,9]
names(temp5) <- temp4
if(length(temp5)<length(temp1)) {
  temp6 <- setdiff(temp1,temp5)
  temp7 <- rep(0,length(temp6))
  names(temp7) <- temp6
  temp5 <- c(temp5,temp7)
}
temp5 <- temp5[temp1]
running.final <- cbind(running.final,temp5)

temp2 <- LAMA2_4007_expr
temp3 <- dimnames(temp2)[[1]]
temp4 <- intersect(temp3,temp1)
temp5 <- temp2[temp4,9]
names(temp5) <- temp4
if(length(temp5)<length(temp1)) {
  temp6 <- setdiff(temp1,temp5)
  temp7 <- rep(0,length(temp6))
  names(temp7) <- temp6
  temp5 <- c(temp5,temp7)
}
temp5 <- temp5[temp1]
running.final <- cbind(running.final,temp5)

temp2 <- LAMA2_4010_expr
temp3 <- dimnames(temp2)[[1]]
temp4 <- intersect(temp3,temp1)
temp5 <- temp2[temp4,9]
names(temp5) <- temp4
if(length(temp5)<length(temp1)) {
  temp6 <- setdiff(temp1,temp5)
  temp7 <- rep(0,length(temp6))
  names(temp7) <- temp6
  temp5 <- c(temp5,temp7)
}
temp5 <- temp5[temp1]
running.final <- cbind(running.final,temp5)

temp2 <- LAMA2_4027_expr
temp3 <- dimnames(temp2)[[1]]
temp4 <- intersect(temp3,temp1)
temp5 <- temp2[temp4,9]
names(temp5) <- temp4
if(length(temp5)<length(temp1)) {
  temp6 <- setdiff(temp1,temp5)
  temp7 <- rep(0,length(temp6))
  names(temp7) <- temp6
  temp5 <- c(temp5,temp7)
}
temp5 <- temp5[temp1]
running.final <- cbind(running.final,temp5)

temp2 <- LAMA2_4043_expr
temp3 <- dimnames(temp2)[[1]]
temp4 <- intersect(temp3,temp1)
temp5 <- temp2[temp4,9]
names(temp5) <- temp4
if(length(temp5)<length(temp1)) {
  temp6 <- setdiff(temp1,temp5)
  temp7 <- rep(0,length(temp6))
  names(temp7) <- temp6
  temp5 <- c(temp5,temp7)
}
temp5 <- temp5[temp1]
running.final <- cbind(running.final,temp5)

temp2 <- LAMA2_4120_expr
temp3 <- dimnames(temp2)[[1]]
temp4 <- intersect(temp3,temp1)
temp5 <- temp2[temp4,9]
names(temp5) <- temp4
if(length(temp5)<length(temp1)) {
  temp6 <- setdiff(temp1,temp5)
  temp7 <- rep(0,length(temp6))
  names(temp7) <- temp6
  temp5 <- c(temp5,temp7)
}
temp5 <- temp5[temp1]
running.final <- cbind(running.final,temp5)

temp2 <- LAMA2_4244_expr
temp3 <- dimnames(temp2)[[1]]
temp4 <- intersect(temp3,temp1)
temp5 <- temp2[temp4,9]
names(temp5) <- temp4
if(length(temp5)<length(temp1)) {
  temp6 <- setdiff(temp1,temp5)
  temp7 <- rep(0,length(temp6))
  names(temp7) <- temp6
  temp5 <- c(temp5,temp7)
}
temp5 <- temp5[temp1]
running.final <- cbind(running.final,temp5)

temp2 <- LAMA2_4363_expr
temp3 <- dimnames(temp2)[[1]]
temp4 <- intersect(temp3,temp1)
temp5 <- temp2[temp4,9]
names(temp5) <- temp4
if(length(temp5)<length(temp1)) {
  temp6 <- setdiff(temp1,temp5)
  temp7 <- rep(0,length(temp6))
  names(temp7) <- temp6
  temp5 <- c(temp5,temp7)
}
temp5 <- temp5[temp1]
running.final <- cbind(running.final,temp5)

temp2 <- LAMA2_4400_expr
temp3 <- dimnames(temp2)[[1]]
temp4 <- intersect(temp3,temp1)
temp5 <- temp2[temp4,9]
names(temp5) <- temp4
if(length(temp5)<length(temp1)) {
  temp6 <- setdiff(temp1,temp5)
  temp7 <- rep(0,length(temp6))
  names(temp7) <- temp6
  temp5 <- c(temp5,temp7)
}
temp5 <- temp5[temp1]
running.final <- cbind(running.final,temp5)

temp2 <- LAMA2_4414_expr
temp3 <- dimnames(temp2)[[1]]
temp4 <- intersect(temp3,temp1)
temp5 <- temp2[temp4,9]
names(temp5) <- temp4
if(length(temp5)<length(temp1)) {
  temp6 <- setdiff(temp1,temp5)
  temp7 <- rep(0,length(temp6))
  names(temp7) <- temp6
  temp5 <- c(temp5,temp7)
}
temp5 <- temp5[temp1]
running.final <- cbind(running.final,temp5)

temp2 <- LAMA2_4422_expr
temp3 <- dimnames(temp2)[[1]]
temp4 <- intersect(temp3,temp1)
temp5 <- temp2[temp4,9]
names(temp5) <- temp4
if(length(temp5)<length(temp1)) {
  temp6 <- setdiff(temp1,temp5)
  temp7 <- rep(0,length(temp6))
  names(temp7) <- temp6
  temp5 <- c(temp5,temp7)
}
temp5 <- temp5[temp1]
running.final <- cbind(running.final,temp5)

##  Not doing this repeat sample
# temp2 <- LAMA2_B11_68_expr
# temp3 <- dimnames(temp2)[[1]]
# temp4 <- intersect(temp3,temp1)
# temp5 <- temp2[temp4,9]
# names(temp5) <- temp4
# if(length(temp5)<length(temp1)) {
#   temp6 <- setdiff(temp1,temp5)
#   temp7 <- rep(0,length(temp6))
#   names(temp7) <- temp6
#   temp5 <- c(temp5,temp7)
# }
# temp5 <- temp5[temp1]
# running.final <- cbind(running.final,temp5)

# temp2 <- LAMA2_B11_69_expr
# temp3 <- dimnames(temp2)[[1]]
# temp4 <- intersect(temp3,temp1)
# temp5 <- temp2[temp4,9]
# names(temp5) <- temp4
# if(length(temp5)<length(temp1)) {
#   temp6 <- setdiff(temp1,temp5)
#   temp7 <- rep(0,length(temp6))
#   names(temp7) <- temp6
#   temp5 <- c(temp5,temp7)
# }
# temp5 <- temp5[temp1]
# running.final <- cbind(running.final,temp5)

temp2 <- LAMA2_B14_14_expr
temp3 <- dimnames(temp2)[[1]]
temp4 <- intersect(temp3,temp1)
temp5 <- temp2[temp4,9]
names(temp5) <- temp4
if(length(temp5)<length(temp1)) {
  temp6 <- setdiff(temp1,temp5)
  temp7 <- rep(0,length(temp6))
  names(temp7) <- temp6
  temp5 <- c(temp5,temp7)
}
temp5 <- temp5[temp1]
running.final <- cbind(running.final,temp5)

temp2 <- LAMA2_B16_21_expr
temp3 <- dimnames(temp2)[[1]]
temp4 <- intersect(temp3,temp1)
temp5 <- temp2[temp4,9]
names(temp5) <- temp4
if(length(temp5)<length(temp1)) {
  temp6 <- setdiff(temp1,temp5)
  temp7 <- rep(0,length(temp6))
  names(temp7) <- temp6
  temp5 <- c(temp5,temp7)
}
temp5 <- temp5[temp1]
running.final <- cbind(running.final,temp5)

temp2 <- M13_14_expr
temp3 <- dimnames(temp2)[[1]]
temp4 <- intersect(temp3,temp1)
temp5 <- temp2[temp4,9]
names(temp5) <- temp4
if(length(temp5)<length(temp1)) {
  temp6 <- setdiff(temp1,temp5)
  temp7 <- rep(0,length(temp6))
  names(temp7) <- temp6
  temp5 <- c(temp5,temp7)
}
temp5 <- temp5[temp1]
running.final <- cbind(running.final,temp5)

temp2 <- M13_24_expr
temp3 <- dimnames(temp2)[[1]]
temp4 <- intersect(temp3,temp1)
temp5 <- temp2[temp4,9]
names(temp5) <- temp4
if(length(temp5)<length(temp1)) {
  temp6 <- setdiff(temp1,temp5)
  temp7 <- rep(0,length(temp6))
  names(temp7) <- temp6
  temp5 <- c(temp5,temp7)
}
temp5 <- temp5[temp1]
running.final <- cbind(running.final,temp5)

temp2 <- M13_2_expr
temp3 <- dimnames(temp2)[[1]]
temp4 <- intersect(temp3,temp1)
temp5 <- temp2[temp4,9]
names(temp5) <- temp4
if(length(temp5)<length(temp1)) {
  temp6 <- setdiff(temp1,temp5)
  temp7 <- rep(0,length(temp6))
  names(temp7) <- temp6
  temp5 <- c(temp5,temp7)
}
temp5 <- temp5[temp1]
running.final <- cbind(running.final,temp5)

temp2 <- M13_32_expr
temp3 <- dimnames(temp2)[[1]]
temp4 <- intersect(temp3,temp1)
temp5 <- temp2[temp4,9]
names(temp5) <- temp4
if(length(temp5)<length(temp1)) {
  temp6 <- setdiff(temp1,temp5)
  temp7 <- rep(0,length(temp6))
  names(temp7) <- temp6
  temp5 <- c(temp5,temp7)
}
temp5 <- temp5[temp1]
running.final <- cbind(running.final,temp5)

temp2 <- M13_37_expr
temp3 <- dimnames(temp2)[[1]]
temp4 <- intersect(temp3,temp1)
temp5 <- temp2[temp4,9]
names(temp5) <- temp4
if(length(temp5)<length(temp1)) {
  temp6 <- setdiff(temp1,temp5)
  temp7 <- rep(0,length(temp6))
  names(temp7) <- temp6
  temp5 <- c(temp5,temp7)
}
temp5 <- temp5[temp1]
running.final <- cbind(running.final,temp5)

temp2 <- M13_38_expr
temp3 <- dimnames(temp2)[[1]]
temp4 <- intersect(temp3,temp1)
temp5 <- temp2[temp4,9]
names(temp5) <- temp4
if(length(temp5)<length(temp1)) {
  temp6 <- setdiff(temp1,temp5)
  temp7 <- rep(0,length(temp6))
  names(temp7) <- temp6
  temp5 <- c(temp5,temp7)
}
temp5 <- temp5[temp1]
running.final <- cbind(running.final,temp5)

temp2 <- M13_42_expr
temp3 <- dimnames(temp2)[[1]]
temp4 <- intersect(temp3,temp1)
temp5 <- temp2[temp4,9]
names(temp5) <- temp4
if(length(temp5)<length(temp1)) {
  temp6 <- setdiff(temp1,temp5)
  temp7 <- rep(0,length(temp6))
  names(temp7) <- temp6
  temp5 <- c(temp5,temp7)
}
temp5 <- temp5[temp1]
running.final <- cbind(running.final,temp5)

temp2 <- M13_49_expr
temp3 <- dimnames(temp2)[[1]]
temp4 <- intersect(temp3,temp1)
temp5 <- temp2[temp4,9]
names(temp5) <- temp4
if(length(temp5)<length(temp1)) {
  temp6 <- setdiff(temp1,temp5)
  temp7 <- rep(0,length(temp6))
  names(temp7) <- temp6
  temp5 <- c(temp5,temp7)
}
temp5 <- temp5[temp1]
running.final <- cbind(running.final,temp5)

temp2 <- M13_8_expr
temp3 <- dimnames(temp2)[[1]]
temp4 <- intersect(temp3,temp1)
temp5 <- temp2[temp4,9]
names(temp5) <- temp4
if(length(temp5)<length(temp1)) {
  temp6 <- setdiff(temp1,temp5)
  temp7 <- rep(0,length(temp6))
  names(temp7) <- temp6
  temp5 <- c(temp5,temp7)
}
temp5 <- temp5[temp1]
running.final <- cbind(running.final,temp5)

# Remove the first column with 9's that was unused
running.final <- running.final[,-c(1)]

# Rename the colnames of the running.final table of TPM and transcripts

dimnames(running.final)[[2]] <- c("LAMA2_15_148", 
                                  # "LAMA2_15_149", 
                                  # "LAMA2_15_151", 
                                  "LAMA2_3935", 
                                  "LAMA2_3946", 
                                  "LAMA2_4002", 
                                  "LAMA2_4007", 
                                  "LAMA2_4010", 
                                  "LAMA2_4027", 
                                  "LAMA2_4043", 
                                  "LAMA2_4120", 
                                  "LAMA2_4244", 
                                  "LAMA2_4363", 
                                  "LAMA2_4400", 
                                  "LAMA2_4414", 
                                  "LAMA2_4422", 
                                  # "LAMA2_B11_68", 
                                  # "LAMA2_B11_69", 
                                  "LAMA2_B14_14", 
                                  "LAMA2_B16_21", 
                                  "M13_14", 
                                  "M13_24", 
                                  "M13_2", 
                                  "M13_32", 
                                  "M13_37", 
                                  "M13_38", 
                                  "M13_42", 
                                  "M13_49", 
                                  "M13_8")

#### OUTPUT ####
write.table(running.final,"2018-07-13_LAMA2_Cross_Sample_Gene_Level_Expression_rm_outliers.txt",sep="\t")
################



###########################################################
#### Setup sample information for statistical analysis ####
###########################################################


sample.info <- dimnames(running.final)[[2]] # Get all the column / sample names from master table
sample.info <- rbind(sample.info,c(rep("LAMA2",16), rep("WT",9))) # Create a second header to label the sample groups as LAMA2 or WT

# specify their age groups
sample.info <- rbind(sample.info,c('young', 
                                   # 'old', # LAMA2-15-149
                                   # 'old', # LAMA2-15-151
                                   'between', 
                                   'between', 
                                   'between',
                                   'young',
                                   'young',
                                   'young',
                                   'old', 
                                   'old',
                                   'between', 
                                   'between', 
                                   'between',
                                   'old',
                                   'young',
                                   # 'young', # B11-68
                                   # 'young',  # B11-69
                                   'young',
                                   'young',
                                   'young',
                                   'between',
                                   'young',
                                   'old',
                                   'young',
                                   'young',
                                   'between', 
                                   'old',
                                   'young'))

# Specify their genders
sample.info <- rbind(sample.info,c("F", 
                                   # "NA", # LAMA2-15-149
                                   # "F",  # LAMA2-15-151 
                                   "M", 
                                   "M", 
                                   "F", 
                                   "M", 
                                   "F", 
                                   "F", 
                                   "F", 
                                   "M", 
                                   "M", 
                                   "M", 
                                   "F", 
                                   "M", 
                                   "M",
                                   # "F", # LAMA2-B11-68
                                   # "M", # LAMA2-B11-69
                                   "M", 
                                   "M", 
                                   "M", 
                                   "F", 
                                   "M", 
                                   "M", 
                                   "M", 
                                   "M", 
                                   "M", 
                                   "F", 
                                   "M"))

# Specify their biopsy sites
sample.info <- rbind(sample.info, c("quad", 
                                    # "quad", # LAMA2-15-149
                                    # "quad", # LAMA2-15-151
                                    "bicep", 
                                    "bicep", 
                                    "bicep", 
                                    "bicep", 
                                    "bicep", 
                                    "bicep", 
                                    "bicep", 
                                    "bicep",
                                    "bicep", 
                                    "bicep", 
                                    "bicep", 
                                    "bicep", 
                                    "bicep", 
                                    # "quad", # LAMA2-B11-68
                                    # "quad", # LAMA2-B11-69
                                    "quad", 
                                    "quad", 
                                    "quad", 
                                    "quad", 
                                    "quad", 
                                    "quad", 
                                    "quad", 
                                    "quad", 
                                    "quad", 
                                    "quad", 
                                    "quad"))

# Age in weeks
sample.info <- rbind(sample.info, c(9, 
                                    # NA, # LAMA2-15-149
                                    # 36, # LAMA2-15-151
                                    17, 
                                    NA, 
                                    19, 
                                    6, 
                                    6, 
                                    8, 
                                    61, 
                                    38, 
                                    14, 
                                    17, 
                                    24, 
                                    37, 
                                    11, 
                                    # 9, # LAMA2-B11-68
                                    # 2, # LAMA2-B11-69
                                    5, 
                                    2, 
                                    10, 
                                    12, 
                                    5, 
                                    21, 
                                    4, 
                                    4, 
                                    15, 
                                    28, 
                                    7))

# Target1 - lumps all WT samples into a single group
sample.info <- rbind(sample.info, c("lama2_young",
                                    # "lama2_old", # LAMA2-15-149
                                    # "lama2_old", # LAMA2-15-151
                                    "lama2_between",
                                    "lama2_between",
                                    "lama2_between",
                                    "lama2_young",
                                    "lama2_young",
                                    "lama2_young",
                                    "lama2_old",
                                    "lama2_old",
                                    "lama2_between",
                                    "lama2_between",
                                    "lama2_between",
                                    "lama2_old",
                                    "lama2_young",
                                    # "lama2_young", # LAMA2-B11-68
                                    # "lama2_young", # LAMA2-B11-69
                                    "lama2_young",
                                    "lama2_young",
                                    "wt",
                                    "wt",
                                    "wt",
                                    "wt",
                                    "wt",
                                    "wt",
                                    "wt",
                                    "wt",
                                    "wt"))

# Target2 - Separates WT samples by age group
sample.info <- rbind(sample.info, c("lama2_young",
                                    # "lama2_old", # LAMA2-15-149
                                    # "lama2_old", # LAMA2-15-151
                                    "lama2_between",
                                    "lama2_between",
                                    "lama2_between",
                                    "lama2_young",
                                    "lama2_young",
                                    "lama2_young",
                                    "lama2_old",
                                    "lama2_old",
                                    "lama2_between",
                                    "lama2_between",
                                    "lama2_between",
                                    "lama2_old",
                                    "lama2_young",
                                    # "lama2_young", # LAMA2-B11-68
                                    # "lama2_young", # LAMA2-B11-69
                                    "lama2_young",
                                    "lama2_young",
                                    "wt_young",
                                    "wt_between",
                                    "wt_young",
                                    "wt_old",
                                    "wt_young",
                                    "wt_young",
                                    "wt_between",
                                    "wt_old",
                                    "wt_young"))

# Assign a color to each of the groups LAMA2-Young, LAMA2-Between, LAMA2-old, Controls (mixed)
sample.info <- rbind(sample.info,c("yellow", 
                                   # "red", # LAMA2-15-149
                                   # "red", # LAMA2-15-151
                                   "orange", 
                                   "orange", 
                                   "orange", 
                                   "yellow", 
                                   "yellow", 
                                   "yellow",
                                   "red", 
                                   "red",
                                   "orange", 
                                   "orange", 
                                   "orange", 
                                   "red", 
                                   "yellow", 
                                   # "green",  # LAMA2-B11-68
                                   # "yellow", # LAMA2-B11-69
                                   "yellow", 
                                   "yellow", 
                                   "blue", 
                                   "blue", 
                                   "blue",
                                   "blue", 
                                   "blue", 
                                   "blue", 
                                   "blue", 
                                   "blue", 
                                   "blue"))

# Shapes to signify gender, Females 19, NA 17, Males 15
sample.info <- rbind(sample.info,c("19", 
                                   # "17", # LAMA2-15-149
                                   # "19", # LAMA2-15-151
                                   "15", 
                                   "15", 
                                   "19", 
                                   "15", 
                                   "19", 
                                   "19", 
                                   "19", 
                                   "15", 
                                   "15", 
                                   "15", 
                                   "19", 
                                   "15", 
                                   "15",
                                   # "19", # LAMA2-B11-68
                                   # "15", # LAMA2-B11-69
                                   "15", 
                                   "15", 
                                   "15", 
                                   "19", 
                                   "15", 
                                   "15", 
                                   "15", 
                                   "15", 
                                   "15", 
                                   "19", 
                                   "15"))



sample.info <- t(sample.info) # Transform the sample.info table from long to tall format
dimnames(sample.info)[[1]] <- sample.info[,1] # Turn the first column/sample names into rownames
sample.info <- sample.info[,-1] # now remove that first column

dimnames(sample.info)[[2]] <- c("Group", "Class","Gender","Biopsy_Site", "Age_Months", "Target1", "Target2", "Color","Shape")
sample.info <- as.data.frame(sample.info)

# Set order of colors
sample.info$Color <- factor(sample.info$Color, levels = c("blue", "yellow", "orange", "red"))

# Sort samples by group, then by age
sample.info <- sample.info[order(sample.info$Color),]

# Copy the master data table
sample.data <- running.final

# Sort the columns of the data table based on the order of the sample.info
sample.data <- sample.data[,row.names(sample.info)]


##########################################################
#### Generate cross sample box plot pre-normalization ####
##########################################################

# Before normalization, we want to see the log2 TPM + 2 expression per sample
# This illustrates how different each sample's expression values are,
# Differences may be due to sequencing depth, preprocessing and alignment
# We need to normalize the data before we can actually make comparisons between samples


#dev.off() # clear old plot devices
dim(sample.data) # original has 58,866 rows, 28 cols
sample.data <- sample.data+2  # Add 2 to all values to increase floor from 0
sample.data <- log2(sample.data) # take the log2 of each value
removenodetect <- apply(sample.data,1,max)>1  # if any items in a row have higher value than 1 set to TRUE, else FALSE (list of T/F for each unique col)
sample.data <- sample.data[removenodetect,] # now basically keep only the rows that are TRUE in the removenodetect list
dim(sample.data) # results in 42,548 rows and 28 cols; removed 16,318 transcripts

col.2.use <- as.character(unlist(sample.info$Color)) # grab the color column and turn into a char vector
WorkingList <- vector("list", dim(sample.data)[2]) # Turn each column into a list of expression values that are easily accessible by index
for ( i in 1:dim(sample.data)[2]) { # For every column (sample) in sample.data
  WorkingList[[i]] <- as.numeric(sample.data[,i]) # turn each item (row) in column into a numeric value (was character)
}
boxplot(WorkingList,col=col.2.use,outline=FALSE, main="Cross Sample Box Plot Pre-normalization", 
        xlab="Sample Index (n=28)",ylab="Expression (Log2(TPM+2))") # Edit colors here as needed or use col=8


########################################
#### Perform Quantile Normalization ####
########################################

# Normalization method that brings all samples to the same total expression/read count
# so that we can compare them. This is one of many normalization methods

# dev.off()
plot(sort(sample.data[,1]),type="l")
running.sort <- NULL
for (i in 1:dim(sample.data)[2]) {
  temp <- sort(sample.data[,i])
  running.sort <- cbind(running.sort,temp)
  lines(sort(sample.data[,i]),col=col.2.use[i]) # Edit colors here as needed
}
running.median <- apply(running.sort,1,median)
running.normalized <- NULL
for (i in 1:dim(sample.data)[2]) {
  temp <- sample.data[,i]
  temp <- sort(temp,index.return = TRUE)
  temp <- temp$ix
  names(temp) <- as.numeric(unlist(running.median))
  temp <- sort(temp)
  temp <- as.numeric(unlist(names(temp)))
  running.normalized <- cbind(running.normalized,temp)
  rm(temp)
}
#dev.list()

dimnames(running.normalized) <- dimnames(sample.data)

####### VARIABLE !!!!
write.table(running.normalized,"2018-07-13_rm_outliers_running.normalized.txt",sep="\t")
#######

#------------------------------------------------------------------------------
# Generate Cross Sample Box Plot - Post-Normalization
#------------------------------------------------------------------------------

# Visualize the normalized samples

dev.off()
WorkingList <- vector("list", dim(running.normalized)[2])
for ( i in 1:dim(running.normalized)[2]) {
  WorkingList[[i]] <- as.numeric(running.normalized[,i])
}
boxplot(WorkingList,col=col.2.use,outline=FALSE,main="Cross Sample Box Plot Post-normalization",
        xlab="Sample Index (n=32) ",ylab="Expression (Quantile(Log2(TPM+2)))") # Edit colors here as needed or use col=8


#------------------------------------------------------------------------------
# EDA - cov-based PCA (2D) - Covariance based 
#------------------------------------------------------------------------------

# We do the covariance based PCA to see how the avg magnitude of gene expression differ between samples

dev.off()
par(mfrow=c(1,2))
shapes.2.use <- as.numeric(unlist(sample.info$Shape))
shapes.2.use <- ifelse(shapes.2.use==1,15,19)
shapes.2.out <- ifelse(shapes.2.use==15,22,21)
pca.results <- princomp(running.normalized,cor=F)
pca.loadings <- loadings(pca.results)
totalvar <- (pca.results$sdev^2)
variancePer <- round(totalvar/sum(totalvar)*100,1)
variancePer <- variancePer[1:3]
palette("default")
# x.range <- range(pca.loadings[,1])
# 0.1954101 0.2030271
x.range <- range(0.1945, 0.204)

# y.range <- range(pca.loadings[,2])
# -0.2157176  0.3258338
y.range = range(-0.2157176, 0.3258338)
xlab <- paste(c("Principal Component 1 (",variancePer[1],"%)"),collapse="")	
ylab <- paste(c("Principal Component 2 (",variancePer[2],"%)"),collapse="")
plot(pca.loadings[,1],pca.loadings[,2],type='n',xlab=xlab,ylab=ylab,bg="gray",
     xlim=x.range, ylim=y.range)
points(pca.loadings[,1],pca.loadings[,2],lwd=0.5,col=col.2.use,cex=2,pch=shapes.2.use)
points(pca.loadings[,1],pca.loadings[,2],lwd=0.5,col=1,cex=2,pch=shapes.2.out)
text(pca.loadings[,1],pca.loadings[,2],rownames(sample.info), cex=0.6,pos=1)


# x2.range <- range(pca.loadings[,2])
# -0.2157176  0.3258338
x2.range <- range(-0.23, 0.34)
# y2.range <- range(pca.loadings[,3])
y2.range <- range(-0.2679105, 0.5559832)
xlab <- paste(c("Principal Component 2 (",variancePer[2],"%)"),collapse="")	
ylab <- paste(c("Principal Component 3 (",variancePer[3],"%)"),collapse="")
plot(pca.loadings[,2],pca.loadings[,3],type='n',xlab=xlab,ylab=ylab,bg="gray",
     xlim=x2.range, ylim=y2.range)
points(pca.loadings[,2],pca.loadings[,3],lwd=0.5,col=col.2.use,cex=2,pch=shapes.2.use)
points(pca.loadings[,2],pca.loadings[,3],lwd=0.5,col=1,cex=2,pch=shapes.2.out)
text(pca.loadings[,2],pca.loadings[,3],rownames(sample.info), cex=0.6,pos=1)

#------------------------------------------------------------------------------
# EDA cor-based Heat Map (unclustered)
#------------------------------------------------------------------------------

# Here we use the correlation to determine the avg directionality of genes (Up or down-reg) in each sample
# In otherwords we are comparing the patterns of expression between samples

dev.off()
cdat <- running.normalized
cdat <- cor(cdat)
par("mar" = c(8, 8, 5, 10) + 0.1) # wide righthand margin 
greys <- rgb(r=(3:0)/3, g=(3:0)/3, b=(3:0)/3, names=paste("greys",3:0,sep="."))
x <- (1:nrow(cdat))
y <- (1:ncol(cdat))
image(x, y, cdat, col=greys, labels=FALSE, tick=FALSE, xlab="", ylab="")
image.scale(as.numeric(cdat),col=greys,digits = 3)
axis(1, at=seq(1, nrow(cdat), by=1), labels=FALSE, tick=FALSE)
text(x=seq(1, ncol(cdat), by=1), par('usr')[3]- 0.80, labels=colnames(cdat), srt=90, pos=2, xpd = TRUE, adj=5,0, cex=0.85)

axis(2, at=seq(1, ncol(cdat), by=1), labels=FALSE, tick=FALSE)
text(y=seq(1, nrow(cdat), by=1), par('usr')[1]- 0.80, labels=rownames(cdat), srt=0, pos=2, xpd = TRUE, adj=5,0, cex=0.85)
box()

for (i in 1:dim(sample.info)[1]) {
  points((i:dim(sample.info)[1]), (i:dim(sample.info)[1]),col=col.2.use[i],pch=15,cex=1.5)
}

# At this stage, you may have found specific samples you want to remove from the analysis if they are big outliers
# Then you would have to go upstream and remove some of the samples and redo everything up to this point.

#------------------------------------------------------------------------------
# Perform noise modeling
#------------------------------------------------------------------------------

# For each group of samples calculate the mean, SD and CV 
# we want to plot the CV vs mean quantile(log2(TPM+2))
# This will allow us to see where there is too much variation (noise/ high CV) compared to the mean expression

group1 <- running.normalized[,col.2.use=="blue"] # Take all the normals and create sub-table
g1.mean <- apply(group1,1,mean) # get the average for each transcript across all samples within the group
g1.sd <- apply(group1,1,sd) # get standard deviation
g1.cv <- (g1.sd/g1.mean)*100 # get coefficient of variance
lowess.g1 <- lowess(g1.cv~g1.mean,f=1/64) # lowess is a gplot function to create smooth lines for plotting

group2 <- running.normalized[,col.2.use=="yellow"] # all young lama2-cmds
g2.mean <- apply(group2,1,mean)
g2.sd <- apply(group2,1,sd)
g2.cv <- (g2.sd/g2.mean)*100
lowess.g2 <- lowess(g2.cv~g2.mean,f=1/64)

group3 <- running.normalized[,col.2.use=="orange"] # all between lama2-cmds
g3.mean <- apply(group3,1,mean)
g3.sd <- apply(group3,1,sd)
g3.cv <- (g3.sd/g3.mean)*100
lowess.g3 <- lowess(g3.cv~g3.mean,f=1/64)

group4 <- running.normalized[,col.2.use=="red"] # all old lama2-cmds
g4.mean <- apply(group4,1,mean)
g4.sd <- apply(group4,1,sd)
g4.cv <- (g4.sd/g4.mean)*100
lowess.g4 <- lowess(g4.cv~g4.mean,f=1/64)

#-------------------------------------------------------------------------------
# Generate Noise Model Plot
#-------------------------------------------------------------------------------

# The graph shows that there is an exponential amount of noise (high CV) at < mean of 3
# But then it begins to level out and linearize at around a mean of 3 and above
# This is the sweet spot that we want to set our noise cut-off

dev.off()

xrange <- range(c(1:15))
yrange <- range(c(0:150))
plot(xrange,yrange,type='n',xlab="Mean (Quantile(log2(TPM+2))",ylab="Coeficient of Variation")
points(g1.mean,g1.cv,pch=21,col=8)
points(g2.mean,g2.cv,pch=21,col=8)
points(g3.mean,g3.cv,pch=21,col=8)
points(g4.mean,g4.cv,pch=21,col=8)
lines(lowess.g1,col="blue",lwd=2)
lines(lowess.g2,col="yellow",lwd=2)
lines(lowess.g3,col="orange",lwd=2)
lines(lowess.g4,col="red",lwd=2)
abline(h=0)
box()

## Zoomed in graph
# xrange <- range(c(1:6))
# yrange <- range(c(0:50))
# plot(xrange,yrange,type='n',xlab="Mean (Quantile(log2(TPM+2))",ylab="Coeficient of Variation")
# points(g1.mean,g1.cv,pch=21,col=8)
# points(g2.mean,g2.cv,pch=21,col=8)
# points(g3.mean,g3.cv,pch=21,col=8)
# points(g4.mean,g4.cv,pch=21,col=8)
# lines(lowess.g1,col="blue",lwd=2)
# lines(lowess.g2,col="yellow",lwd=2)
# lines(lowess.g3,col="orange",lwd=2)
# lines(lowess.g4,col="red",lwd=2)
# abline(h=0)
# box()

noise.pick <- 3
abline(v=noise.pick,lty=2,lwd=2)

#-------------------------------------------------------------------------------
# Remove Noise & Floor
#-------------------------------------------------------------------------------

dim(running.normalized)
rrr <- running.normalized
rrr <- ifelse(rrr<noise.pick,noise.pick,rrr)
temp <- apply(rrr,1,mean)>noise.pick
Final.working.data <- rrr[temp,]
dim(Final.working.data)

#### Variable
write.table(Final.working.data,"2018-07-13_Final.working.data_rm_outliers.txt",sep="\t")
####

#------------------------------------------------------------------------------
# Repeat EDA -> cov-based PCA (2D) post-noise
#------------------------------------------------------------------------------

# We want to repeat the covariance based PCA to make sure that the filtration process didn't
# drastically alter the way the samples group or the avg magnitude of gene expression of each sample


dev.off()
par(mfrow=c(1,2))
shapes.2.use <- as.numeric(unlist(sample.info$Shape))
shapes.2.use <- ifelse(shapes.2.use==1,15,19)
shapes.2.out <- ifelse(shapes.2.use==15,22,21)
pca.results <- princomp(Final.working.data,cor=F)
pca.loadings <- loadings(pca.results)
totalvar <- (pca.results$sdev^2)
variancePer <- round(totalvar/sum(totalvar)*100,1)
variancePer <- variancePer[1:3]
palette("default")
# x.range <- range(pca.loadings[,1])
# x.range <- range(0.1835,0.194)
# y.range = range(-0.25,0.32)
x.range <- range(0.1945, 0.2055)
y.range = range(-0.25, 0.3258338)
xlab <- paste(c("Principal Component 1 (",variancePer[1],"%)"),collapse="")	
ylab <- paste(c("Principal Component 2 (",variancePer[2],"%)"),collapse="")
plot(pca.loadings[,1],pca.loadings[,2],type='n',xlab=xlab,ylab=ylab,bg="gray",
     xlim=x.range, ylim=y.range)
points(pca.loadings[,1],pca.loadings[,2],lwd=0.5,col=col.2.use,cex=2,pch=shapes.2.use)
points(pca.loadings[,1],pca.loadings[,2],lwd=0.5,col=1,cex=2,pch=shapes.2.out)
text(pca.loadings[,1],pca.loadings[,2],rownames(sample.info), cex=0.6,pos=1)
# x2.range <- range(-0.27, 0.32)
# y2.range <- range(-0.35, 0.5209382)
x2.range <- range(-0.24, 0.34)
y2.range <- range(-0.31, 0.5559832)

xlab <- paste(c("Principal Component 2 (",variancePer[2],"%)"),collapse="")	
ylab <- paste(c("Principal Component 3 (",variancePer[3],"%)"),collapse="")
plot(pca.loadings[,2],pca.loadings[,3],type='n',xlab=xlab,ylab=ylab,bg="gray",
     xlim=x2.range, ylim=y2.range)
points(pca.loadings[,2],pca.loadings[,3],lwd=0.5,col=col.2.use,cex=2,pch=shapes.2.use)
points(pca.loadings[,2],pca.loadings[,3],lwd=0.5,col=1,cex=2,pch=shapes.2.out)
text(pca.loadings[,2],pca.loadings[,3],rownames(sample.info), cex=0.6,pos=1)


#------------------------------------------------------------------------------
# EDA cor-based Heat Map (unclustered)
#------------------------------------------------------------------------------

# Here you do expect the pattern of expression to change a bit.  The "noise" was a
# common element between samples that was removed.  Therefore, samples within groups
# should now become more similar and samples between groups should become more different
# The result is that the range of correlation values should increase. Some of the
# differences may have come from removing any outliers you decided to take out of the analysis


dev.off()
cdat <- Final.working.data
cdat <- cor(cdat)
par("mar" = c(8, 8, 5, 10) + 0.1) # wide righthand margin 
greys <- rgb(r=(3:0)/3, g=(3:0)/3, b=(3:0)/3, names=paste("greys",3:0,sep="."))
x <- (1:nrow(cdat))
y <- (1:ncol(cdat))
image(x, y, cdat, col=greys, labels=FALSE, tick=FALSE, xlab="", ylab="")
image.scale(as.numeric(cdat),col=greys,digits = 3)
axis(1, at=seq(1, nrow(cdat), by=1), labels=FALSE, tick=FALSE)
text(x=seq(1, ncol(cdat), by=1), par('usr')[3]- 0.80, labels=colnames(cdat), srt=90, pos=2, xpd = TRUE, adj=5,0, cex=0.85)

axis(2, at=seq(1, ncol(cdat), by=1), labels=FALSE, tick=FALSE)
text(y=seq(1, nrow(cdat), by=1), par('usr')[1]- 0.80, labels=rownames(cdat), srt=0, pos=2, xpd = TRUE, adj=5,0, cex=0.85)
box()

for (i in 1:dim(sample.info)[1]) {
  points((i:dim(sample.info)[1]), (i:dim(sample.info)[1]),col=col.2.use[i],pch=15,cex=1.5)
}

#------------------------------------------------------------------------------
# Post noise filter cor-based Heat Map (clustered)
#------------------------------------------------------------------------------

# Now redo the heatmap using hierarchical clustering. This may take a long time to draw the graph and annotations
# because there are so many genes and samples (litterally like 5 mins or more to complete).
# Again, the heatmap is here to show you the pattern of expression or directionality of each gene (up or down-reg)
# With this information, we can see how samples are grouping up based on their expression and decide if some of them
# are outliers (e.g. are some of the disease samples more similar to control samples? etc). Removing outliers may
# improve the power of our analysis, but also reduce our sample size.

cc <- greenred(64)
dd <- Final.working.data
ccc <- col.2.use
heatmap.2(as.matrix(dd),trace="none",density.info="none",col=cc,ColSideColors=ccc, margins=c(12,8),srtCol=45)



#------------------------------------------------------------------------------
# ANCOVA testing - Statistical testing between groups + corrections 
#------------------------------------------------------------------------------

# Here we start to do some of the statistical testing to see how samples differ
# between groups.  However, there are other possible confounding factors that may
# need to be corrected such as gender and biopsy site.

working.gender <- sample.info$Gender
working.site <- sample.info$Biopsy_Site
working.target <- sample.info$Target1
working05subset <- Final.working.data

running.posthoc <- NULL
running.pred.data <- NULL

#------------------------------------------------------------------------
# Test out without the loop first
#------------------------------------------------------------------------

i <- 1
print(i)
print(dim(working05subset))

# Create a small data frame with relevant groups for analysis
perm.temp <- data.frame(DATA=as.numeric(t(working05subset[i,])),
                        GENDER=factor(working.gender),
                        SITE=factor(working.site),
                        TARGET=factor(working.target))
dimnames(perm.temp)[[2]][1] <- "DATA"
c1.init <- aov(DATA~TARGET+GENDER+SITE,data=perm.temp)
c1.step <- step(c1.init,trace=0) # Step challenge to find meaningful factor combinations
c1.keep <- summary(c1.step)
ccc <- c1.keep[[1]]
ccc <- dimnames(ccc)[[1]]
ccc <- sub("\\s+", "", ccc)
ccc.check <- ccc=="TARGET"
if(sum(ccc.check)==0) {
  ccc <- c("TARGET",ccc)
}

ccc.report <- ccc

# Set up the contrasts / comparisons between groups
if (deviance(c1.step) > sqrt(.Machine$double.eps)) {
  ccc.report <- ccc.report[-c(length(ccc.report))]
  ccc.report <- paste(ccc.report,collapse="+")
  #levels(perm.temp$TARGET) ## check levels for order
  contrasts(perm.temp$TARGET) <- cbind(c(0,0,1,-1), # WT vs lama2_young
                                       c(1,0,0,-1), # WT vs lama2_between
                                       c(0,1,0,-1)) # WT vs lama2_old
                                       # c(1,0,-1,0), # lama2_young vs lama2_between
                                       # c(0,1,-1,0), # lama2_young vs lama2_old
                                       # c(-1,1,0,0), # lama2_between vs lama2_old
                                       # c(1,1,1,-3)) # WT vs all lama2
}

ccc <- paste(c("c1.final <- aov(DATA~",ccc.report,",data=perm.temp)"),collapse="")
eval(parse(text=ccc)) # produces c1.final
temp <- Anova(c1.final,type="III",singular.ok=TRUE) # Set to type III Sum of squares
temp1 <- temp[[1]][2]
temp2 <- temp[[2]][2]
temp3 <- temp[[3]][2]
temp4 <- temp[[4]][2]
out2 <- c(ccc.report,temp1,temp2,temp3,temp4)
names(out2) <- c("ANCOVA_Model","ANCOVA_TARGET_SumSq","ANCOVA_TARGET_DF","ANCOVA_TARGET_Fstat","ANCOVA_TARGET_UncorrectedP")
outp <- predict(c1.final,data=perm.temp)
running.pred.data <- rbind(running.pred.data,outp)
m.summary <- summary.lm(c1.final)
m.adj <- effect("TARGET",c1.final) 

unadj.control.m <- mean(perm.temp[,1][perm.temp[,4]=="wt"])
unadj.target.m <- mean(perm.temp[,1][perm.temp[,4]=="lama2_young"])
out1 <- c(unadj.control.m,unadj.target.m)

names(out1) <- c("UnadjMean_WT","UnadjMean_lama2_young")
out3 <- c(m.adj[[5]][2],m.adj[[5]][1])
names(out3) <- c("ANCOVA_AdjMean_WT","ANCOVA_AdjMean_lama2_young")
postHocs <- glht(c1.final,linfct=mcp(TARGET="Tukey"))
postHocs.s <- summary(postHocs)
postHocs.c <- confint(postHocs)
ph.c <- postHocs.s[[10]]$coefficients[1]
ph.p <- postHocs.s[[10]]$pvalues[1]
ph.t <- postHocs.s[[10]]$tstat[1]
ph.s <- postHocs.s[[10]]$sigma[1]
ph.ci <- postHocs.c[[10]]
compterms <- dimnames(ph.ci)[[1]]
for (j in 1:length(compterms)) {
  compterms[j] <- sub(" - ", "_vs_",compterms[j])
  compterms[j] <- paste(c("PostHoc_",compterms[j]),collapse="")
}
out4 <- NULL

for (j in 1:length(compterms)) {
  out.temp <- c(ph.ci[j,],ph.s[j],ph.t[j],ph.p[j])
  lin.fc <- out.temp[1]
  lin.fc <- 2^lin.fc
  lin.fc <- ifelse(lin.fc<1,-1/lin.fc,lin.fc)
  out.temp <- c(out.temp,lin.fc)
  names(out.temp) <- c(paste(c(compterms[j],"_AdjMeanDiff"),collapse=""),
                       paste(c(compterms[j],"_AdjMeanDiff_lwr"),collapse=""),
                       paste(c(compterms[j],"_AdjMeanDiff_upr"),collapse=""),
                       paste(c(compterms[j],"_stder"),collapse=""),
                       paste(c(compterms[j],"_Fstat"),collapse=""),
                       paste(c(compterms[j],"_Pvalue"),collapse=""),
                       paste(c(compterms[j],"_LinearFoldChange"),collapse=""))
  out4 <- c(out4,out.temp)
}

running.posthoc <- rbind(running.posthoc,c(out1,out2,out3,out4))





####### FULL RUN #################

options(warn = -1)
running.posthoc <- NULL
running.pred.data <- NULL
for (i in 1:dim(working05subset)[[1]]) {
# for (i in 1:20) {
  print(i)
  print(dim(working05subset))
  
  perm.temp <- data.frame(DATA=as.numeric(t(working05subset[i,])),
                          GENDER=factor(working.gender),
                          SITE=factor(working.site),
                          TARGET=factor(working.target))
  
  dimnames(perm.temp)[[2]][1] <- "DATA"
  c1.init <- aov(DATA~TARGET+GENDER+SITE,data=perm.temp)
  c1.step <- step(c1.init,trace=0) # Step challenge to find meaningful factor combinations
  c1.keep <- summary(c1.step)
  ccc <- c1.keep[[1]]
  ccc <- dimnames(ccc)[[1]]
  ccc <- sub("\\s+", "", ccc)
  ccc.check <- ccc=="TARGET"
  if(sum(ccc.check)==0) {
    ccc <- c("TARGET",ccc)
  }
  ccc.report <- ccc
  if (deviance(c1.step) > sqrt(.Machine$double.eps)) {
    ccc.report <- ccc.report[-c(length(ccc.report))]
    ccc.report <- paste(ccc.report,collapse="+")
    #levels(perm.temp$TARGET)
    contrasts(perm.temp$TARGET) <- cbind(c(0,0,1,-1), # WT vs lama2_young
                                         c(1,0,0,-1), # WT vs lama2_between
                                         c(0,1,0,-1)) # WT vs lama2_old
                                         # c(1,0,-1,0), # lama2_young vs lama2_between
                                         # c(0,1,-1,0), # lama2_young vs lama2_old
                                         # c(-1,1,0,0), # lama2_between vs lama2_old
                                         # c(1,1,1,-3)) # WT vs all lama2
    ccc <- paste(c("c1.final <- aov(DATA~",ccc.report,",data=perm.temp)"),collapse="")
    eval(parse(text=ccc)) # produces c1.final
    temp <- Anova(c1.final,type="III",singular.ok=TRUE)
    temp1 <- temp[[1]][2]
    temp2 <- temp[[2]][2]
    temp3 <- temp[[3]][2]
    temp4 <- temp[[4]][2]
    out2 <- c(ccc.report,temp1,temp2,temp3,temp4)
    names(out2) <- c("ANCOVA_Model","ANCOVA_TARGET_SumSq","ANCOVA_TARGET_DF","ANCOVA_TARGET_Fstat","ANCOVA_TARGET_UncorrectedP")
    outp <- predict(c1.final,data=perm.temp)
    running.pred.data <- rbind(running.pred.data,outp)
    m.summary <- summary.lm(c1.final)
    m.adj <- effect("TARGET",c1.final) 
    unadj.wt.m <- mean(perm.temp[,1][perm.temp[,4]=="wt"])
    unadj.lama2y.m <- mean(perm.temp[,1][perm.temp[,4]=="lama2_young"])
    unadj.lama2b.m <- mean(perm.temp[,1][perm.temp[,4]=="lama2_between"])
    unadj.lama2o.m <- mean(perm.temp[,1][perm.temp[,4]=="lama2_old"])
    out1 <- c(unadj.wt.m,unadj.lama2y.m,unadj.lama2b.m,unadj.lama2o.m)
    names(out1) <- c("UnadjMean_WT","UnadjMean_Lama2young","UnadjMean_Lama2between","UnadjMean_Lama2old")
    out3 <- c(m.adj[[5]][4],m.adj[[5]][3],m.adj[[5]][1],m.adj[[5]][2])
    names(out3) <- c("ANCOVA_AdjMean_WT","ANCOVA_AdjMean_Lama2young","ANCOVA_AdjMean_Lama2between","ANCOVA_AdjMean_Lama2old")
    postHocs <- glht(c1.final,linfct=mcp(TARGET="Tukey"))
    postHocs.s <- summary(postHocs)
    postHocs.c <- confint(postHocs)
    ph.c <- postHocs.s[[10]]$coefficients
    ph.p <- postHocs.s[[10]]$pvalues
    ph.t <- postHocs.s[[10]]$tstat
    ph.s <- postHocs.s[[10]]$sigma
    ph.ci <- postHocs.c[[10]]
    compterms <- dimnames(ph.ci)[[1]]
    for (j in 1:length(compterms)) {
      compterms[j] <- sub(" - ", "_vs_",compterms[j])
      compterms[j] <- paste(c("PostHoc_",compterms[j]),collapse="")
    }
    out4 <- NULL
    for (j in 1:length(compterms)) {
      out.temp <- c(ph.ci[j,],ph.s[j],ph.t[j],ph.p[j])
      lin.fc <- out.temp[1]
      lin.fc <- 2^lin.fc
      lin.fc <- ifelse(lin.fc<1,-1/lin.fc,lin.fc)
      out.temp <- c(out.temp,lin.fc)
      names(out.temp) <- c(paste(c(compterms[j],"_AdjMeanDiff"),collapse=""),
                           paste(c(compterms[j],"_AdjMeanDiff_lwr"),collapse=""),
                           paste(c(compterms[j],"_AdjMeanDiff_upr"),collapse=""),
                           paste(c(compterms[j],"_stder"),collapse=""),
                           paste(c(compterms[j],"_Fstat"),collapse=""),
                           paste(c(compterms[j],"_Pvalue"),collapse=""),
                           paste(c(compterms[j],"_LinearFoldChange"),collapse=""))
      out4 <- c(out4,out.temp)
    }
    running.posthoc <- rbind(running.posthoc,c(out1,out2,out3,out4))
  } else {
    # determine length first then set it in rbind command
    # length(c(out1,out2,out3,out4)) # 55
    running.posthoc <- rbind(running.posthoc,c(rep(NA,55)))
    running.pred.data <- rbind(running.pred.data,rep(NA,dim(working05subset)[2]))
    
  }
}
dimnames(running.posthoc)[[1]] <- dimnames(working05subset)[[1]]
dimnames(running.pred.data)[[1]] <- dimnames(working05subset)[[1]]

write.table(running.posthoc,"running.posthoc_rm_extra_contrasts.txt",sep="\t")
write.table(running.pred.data,"running.pred.datarm_extra_contrasts.txt",sep="\t")

colnames(running.posthoc) # find the [#] "ANCOVA_TARGET_UncorrectedP" column number
temp1 <- as.numeric(unlist(running.posthoc[,9])) # select the unadjusted ancova p-val column

temp1 <- ifelse(is.na(temp1),1,temp1)
temp1 <- p.adjust(temp1,method="BH") # FDR MCC Benjamini–Hochberg
running.posthoc <- cbind(running.posthoc,temp1)

# change c(#) to the last column of the ancova output table and name it the corrected p-val
colnames(running.posthoc) # pick the column number corresponding to "temp1"
dimnames(running.posthoc)[[2]][c(56)] <- c("ANCOVA_TARGET_CorrectedP")

write.table(running.posthoc,"2018-07-13_running.posthoc_correctedP_rm_extra_contrasts.txt",sep="\t")

#------------------------------------------------------------------------------
# Select significantly different
#------------------------------------------------------------------------------



# class.results.capture <- running.posthoc
# coded.fc <- ifelse(abs(as.numeric(class.results.capture[,16]))>=1.5,1,0)
# coded.pv <- ifelse(as.numeric(class.results.capture[,17])<0.05,1,0)
# coded.pass <- coded.fc+coded.pv
# coded.pass <- ifelse(is.na(coded.pass),0,coded.pass)
# names(coded.pass) <- dimnames(class.results.capture)[[1]]
# coded.keepers <- names(coded.pass[coded.pass==2])
# length(coded.keepers)
# 
# barplot(c(dim(working05subset)[1],dim(working05subset)[1]-length(coded.keepers),length(coded.keepers)),col=1)
# abline(h=0)
# c(dim(working05subset)[1],dim(working05subset)[1]-length(coded.keepers),length(coded.keepers))

# [1] "UnadjMean_WT"                                         
# [2] "UnadjMean_Lama2young"                                 
# [3] "UnadjMean_Lama2between"                               
# [4] "UnadjMean_Lama2old"                                   
# [5] "ANCOVA_Model"                                         
# [6] "ANCOVA_TARGET_SumSq"                                  
# [7] "ANCOVA_TARGET_DF"                                     
# [8] "ANCOVA_TARGET_Fstat"                                  
# [9] "ANCOVA_TARGET_UncorrectedP"                           
# [10] "ANCOVA_AdjMean_WT"                                    
# [11] "ANCOVA_AdjMean_Lama2young"                            
# [12] "ANCOVA_AdjMean_Lama2between"                          
# [13] "ANCOVA_AdjMean_Lama2old"                              

temp.filter <- running.posthoc[running.posthoc[,56]<0.05,]
temp.fcs <- cbind(temp.filter[,c(20,27,34,41,48,55)])
temp.fcs2 <- matrix(as.numeric(unlist(temp.fcs)),dim(temp.fcs)[1],dim(temp.fcs)[2])
dimnames(temp.fcs2) <- dimnames(temp.fcs)

temp.pvs <- cbind(temp.filter[,c(19,26,33,40,47,54)])
temp.pvs2 <- matrix(as.numeric(unlist(temp.pvs)),dim(temp.pvs)[1],dim(temp.pvs)[2])
dimnames(temp.pvs2) <- dimnames(temp.pvs)

temp.fcs3 <- as.matrix(temp.fcs2)
temp.pvs3 <- as.matrix(temp.pvs2)
temp.fcs.coded <- ifelse(abs(temp.fcs3)>=1.5,1,0)
temp.pvs.coded <- ifelse(temp.pvs3<0.05,1,0)

cx.coded <- temp.fcs.coded+temp.pvs.coded
cx.coded2 <- ifelse(cx.coded==2,1,0)
apply(cx.coded2,2,sum)
pdf("Changers_Plot.pdf")
barplot(apply(cx.coded2,2,sum),col=1,ylab="# Differentially Identified Transcripts",xlab="Comparison")
abline(h=0)
box()
dev.off()
cx.coded3 <- as.matrix(cx.coded2)
keepers <- dimnames(cx.coded2)[[1]][as.logical(apply(cx.coded3,1,max)==1)]


#------------------------------------------------------------------------------
# Post Selection cov-based PCA
#------------------------------------------------------------------------------

dev.off()
par(mfrow=c(1,2))
shapes.2.use <- as.numeric(unlist(sample.info$Shape))
shapes.2.use <- ifelse(shapes.2.use==1,15,19)
shapes.2.out <- ifelse(shapes.2.use==15,22,21)
pca.results <- princomp(Final.working.data[keepers,],cor=F)
pca.loadings <- loadings(pca.results)
totalvar <- (pca.results$sdev^2)
variancePer <- round(totalvar/sum(totalvar)*100,1)
variancePer <- variancePer[1:3]
palette("default")
# x.range <- range(pca.loadings[,1])
x.range <- range(0.185, 0.215)
y.range = range(-0.2270517, 0.3258239)
xlab <- paste(c("Principal Component 1 (",variancePer[1],"%)"),collapse="")	
ylab <- paste(c("Principal Component 2 (",variancePer[2],"%)"),collapse="")
plot(pca.loadings[,1],pca.loadings[,2],type='n',xlab=xlab,ylab=ylab,bg="gray",xlim=x.range, ylim=y.range)
#plot(pca.loadings[,1],pca.loadings[,2],type='n',xlab=xlab,ylab=ylab,bg="gray")
points(pca.loadings[,1],pca.loadings[,2],lwd=0.5,col=col.2.use,cex=2,pch=shapes.2.use)
points(pca.loadings[,1],pca.loadings[,2],lwd=0.5,col=1,cex=2,pch=shapes.2.out)
text(pca.loadings[,1],pca.loadings[,2],rownames(sample.info), cex=0.6,pos=1)
x2.range <- range(-0.2370517, 0.3358239)
y2.range <- range(-0.4336696, 0.5102311)
xlab <- paste(c("Principal Component 2 (",variancePer[2],"%)"),collapse="")	
ylab <- paste(c("Principal Component 3 (",variancePer[3],"%)"),collapse="")
plot(pca.loadings[,2],pca.loadings[,3],type='n',xlab=xlab,ylab=ylab,bg="gray",xlim=x2.range, ylim=y2.range)
# plot(pca.loadings[,2],pca.loadings[,3],type='n',xlab=xlab,ylab=ylab,bg="gray")
points(pca.loadings[,2],pca.loadings[,3],lwd=0.5,col=col.2.use,cex=2,pch=shapes.2.use)
points(pca.loadings[,2],pca.loadings[,3],lwd=0.5,col=1,cex=2,pch=shapes.2.out)
text(pca.loadings[,2],pca.loadings[,3],rownames(sample.info), cex=0.6,pos=1)

#------------------------------------------------------------------------------
# Post Selection cor-based Heat Map (Clustered)
#------------------------------------------------------------------------------

# Actual expression
dev.off()
c <- greenred(64)
dd <- Final.working.data[keepers,]
ccc <- col.2.use
heatmap.2(as.matrix(dd),trace="none",density.info="none",col=cc,ColSideColors=ccc, margins=c(12,8),srtCol=45)

# -----------------------------------
# corrected expression
# -----------------------------------


dev.off()
c <- greenred(64)
dd <- running.pred.data[keepers,]
ccc <- col.2.use
heatmap.2(as.matrix(dd),trace="none",density.info="none",col=cc,ColSideColors=ccc, margins=c(12,8),srtCol=45)

# Focus Expression Box Plot

lookie <- "NM_001007176_chr8_49072341_49076083_Plus"
dev.off()
par(mfrow=c(1,3))
temp1 <- Final.working.data[lookie,]
ppp <- sample.info$Gender
ppp <- ifelse(ppp=="M",19,ppp)
ppp <- ifelse(ppp=="1",17,ppp)
ppp <- ifelse(ppp==3,15,ppp)
plot(temp1,col=col.2.use,pch=ppp,cex=2)
WorkingList <- vector("list", 4)
temp1 <- Final.working.data[lookie,]
WorkingList[[1]] <- as.numeric(temp1[col.2.use=="blue"])
WorkingList[[2]] <- as.numeric(temp1[col.2.use=="yellow"])
WorkingList[[3]] <- as.numeric(temp1[col.2.use=="orange"])
WorkingList[[4]] <- as.numeric(temp1[col.2.use=="red"])
boxplot(WorkingList,col=c("blue","yellow","orange","red"),outline=FALSE,xlab="Index",ylab="ANCOVA Corrected Expression",main=lookie)
WorkingList <- vector("list", 4)
temp1 <- running.pred.data[lookie,]
WorkingList[[1]] <- as.numeric(temp1[col.2.use=="blue"])
WorkingList[[2]] <- as.numeric(temp1[col.2.use=="yellow"])
WorkingList[[3]] <- as.numeric(temp1[col.2.use=="orange"])
WorkingList[[4]] <- as.numeric(temp1[col.2.use=="red"])
boxplot(WorkingList,col=c("blue","yellow","orange","red"),outline=FALSE,xlab="Index",ylab="ANCOVA Corrected Expression",main=lookie)

lookie <- "NM_003356_chr11_74000281_74009237_Minus"
dev.off()
par(mfrow=c(1,3))
temp1 <- Final.working.data[lookie,]
ppp <- sample.info$Biopsy_Site
ppp <- ifelse(ppp=="bicep",19,ppp)
ppp <- ifelse(ppp=="2",17,ppp)
plot(temp1,col=col.2.use,pch=ppp,cex=2)
WorkingList <- vector("list", 4)
temp1 <- Final.working.data[lookie,]
WorkingList[[1]] <- as.numeric(temp1[col.2.use=="blue"])
WorkingList[[2]] <- as.numeric(temp1[col.2.use=="yellow"])
WorkingList[[3]] <- as.numeric(temp1[col.2.use=="orange"])
WorkingList[[4]] <- as.numeric(temp1[col.2.use=="red"])
boxplot(WorkingList,col=c("blue","yellow","orange","red"),outline=FALSE,xlab="Index",ylab="ANCOVA Corrected Expression",main=lookie)
WorkingList <- vector("list", 4)
temp1 <- running.pred.data[lookie,]
WorkingList[[1]] <- as.numeric(temp1[col.2.use=="blue"])
WorkingList[[2]] <- as.numeric(temp1[col.2.use=="yellow"])
WorkingList[[3]] <- as.numeric(temp1[col.2.use=="orange"])
WorkingList[[4]] <- as.numeric(temp1[col.2.use=="red"])
boxplot(WorkingList,col=c("blue","yellow","orange","red"),outline=FALSE,xlab="Index",ylab="ANCOVA Corrected Expression",main=lookie)
