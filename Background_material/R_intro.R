# Welcome in Rstudio! ----

# Comments are placed after a hash tag
# The script is used to store code, the console to execute code
# If you are on a specific line, press "ctrl+enter" to execute

# Variables ----

# Multiple types of data in R
1
"a"
TRUE

# Store in variables
a <- 5
# Show what is stored in a variable (also look at environment tab)
a
# Reuse the variable
a + a

# Ex: Store a number of your choice in a new variable
b <- 3
# Ex: Multiply a and b
a * b
# Ex: Store your name in a variable 'name'
name <- "Sofie"

# Vectors ----

# Multiple values -> vector
b <- 1:10
# Concatenation with the c() function
c(b, b)
# Operations are often applied per value
b + b

# Ex: concatenate a and b
c(a, b)
# Ex: concatenate your name and a -> all values in vector same type
c(name, a)
# Ex: add a to b -> shorter vector is repeated
a + b

# Functions ----

# Example of a function
mean(b)
length(b)
head(b)

# Functions can take multiple arguments
head(b, n = 3)

# Help page
?head

# More structured data ----

## Matrices ----

m <- readRDS("Data/exampleMatrix.RDS")
m
nrow(m)
ncol(m)
m[1,]
m[, 1]
m[, "CD45"]
m[, c(1,2)]

# Ex: Show the first two rows
m[c(1,2), ]
# Ex: Show the CD19 value of the third cell
m[3, "CD19"]

## Lists ----

myList <- list(numbers = c(1,2,4),
               letters = c("a","b"))
myList$numbers
myList[["numbers"]]

# Ex: Create a list with your name, your job and your number of pets
person <- list(name = "Sofie",
               job = "Bioinformatician",
               pets = "1")

# Packages ----

# install.packages("Rtsne")
library("Rtsne")
tsne <- Rtsne(m, perplexity = 2)

library("ggplot2")
to_plot <- data.frame(m,
                      tsne_1 = tsne$Y[,1],
                      tsne_2 = tsne$Y[,2])
ggplot(to_plot) +
  geom_point(aes(x = tsne_1, y = tsne_2, col = FSC.A)) +
  scale_color_distiller(palette = "RdYlBu") +
  theme_minimal()
ggsave("tsne_test.pdf", width = 4, height = 4)

writexl::write_xlsx(tsne$Y, file = "tsne_result.xlsx")
imported_tsne <- readxl::read_xlsx("tsne_result.xlsx")

# Ex: Create a plot of FSC vs SSC, colored by CD19 and save it to a png
ggplot(to_plot) +
  geom_point(aes(x = FSC.A, y = SSC.A, col = CD19)) +
  scale_color_distiller(palette = "RdYlBu") +
  theme_minimal()
ggsave("FSCSSCCD19.png", width = 4, height = 4)


# Cytometry data ----

library(flowCore) # Main library for handling cytometry data

file <- "Data/Raw/Tube01_WT_Day1.fcs" # Path relative to working directory
ff <- read.FCS(file) # Read the fcs file into a flowframe
head(exprs(ff)) # Look at the expression matrix
head(keyword(ff)) # Look at all metadata
keyword(ff, "SPILL") # Look at one specific keyword

pData(parameters(ff)) # Look at the channels and marker names
FlowSOM::GetChannels(ff, c("CD3","CD11b"), exact = FALSE)

# Ex: How many events are measured in tube 5?
ff5 <- read.FCS("Data/Raw/Tube05_KO_Day1.fcs")
nrow(exprs(ff5))
nrow(ff5)

# Ex: On which $DATE was tube 5 measured?
keyword(ff, "$DATE")
names(keyword(ff))

# Ex: What is measured in PE-Cy5-A
FlowSOM::GetMarkers(ff, "PE-Cy5-A")
