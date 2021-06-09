library(ggplot2)
library(dplyr)
library(tidyr)

# Create data

# gene <- c("MYH9","MYH9","MYH9","MYH9","MYH9","MYH9")
cell_type <- c('PT', 'LOH', 'DCT', 'CD-PC', 'CD-IC', 'CD-Trans')
expression <- c('5.022',"10.436","12.243","17.356","12.955","31.818")

df <- data.frame(cell_type,expression)
df

p <- ggplot(df, aes(x=cell_type, y=factor(expression), fill=cell_type)) 
    + geom_violin() + geom_boxplot(width=0.1, color="grey", alpha=0.2)  
    + ggtitle("MYH9")
p

ggplot(df, aes(x=cell_type, y=expression, fill=cell_type)) +
  geom_violin(width=1.4) +
  geom_boxplot(width=2, color="grey", alpha=0.4) +
  theme(legend.position="none",plot.title = element_text(size=11)) +
  ggtitle("MYH9")