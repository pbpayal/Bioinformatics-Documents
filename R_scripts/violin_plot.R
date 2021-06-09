library(ggplot2)
library(dplyr)
library(tidyr)

# MYH9
data <- data.frame(
  cell_type=c("PT","LOH","DCT","CD-PC","CD-IC","CD-Trans"),
  expression=as.numeric(c('5.022',"10.436","12.243","17.356","12.955","31.818"))
)


x <- data[order(df$expression),] # sort by expression
data$cell_type <- factor(data$cell_type) # it must be a factor


dotchart(x$expression,labels=x$cell_type,cex=.7,
         main="MYH9",
         xlab="MYH9 Expression ")


p<-ggplot(x, aes(x=expression, y=cell_type)) + 
  geom_dotplot(binaxis='y', stackdir='center')
p

# # Change dotsize and stack ratio
# q <- ggplot(x, aes(x=expression, y=cell_type)) + 
#   geom_dotplot(binaxis='y', stackdir='center',
#                stackratio=1.5, dotsize=1.5)

# Rotate the dot plot
p + coord_flip()

# Change color by groups
dp <-ggplot(x, aes(x=expression, y=cell_type, fill=expression)) + 
  geom_dotplot(binaxis='y', stackdir='center')+
  labs(title="MYH9 Expression in Kidney Cells",x="Expression", y = "Cell Type")
dp + theme_classic() 


# Sequential color scheme
dp <-ggplot(x, aes(x=expression, y=cell_type, color=expression)) + 
  geom_dotplot(binaxis='y', stackdir='center')+
  labs(title="MYH9 Expression in Kidney Cells",x="Expression", y = "Cell Type")
dp+scale_color_gradient(low="blue", high="red")


# MYH10
myh10 <- data.frame(
  m10_cell_type=c("PT","LOH","DCT","CD-PC","CD-IC","CD-Trans"),
  m_10_expression=as.numeric(c('2.568',"2.277","4.471","25.747","3.875","20.909"))
)
myh10

myh10_ordered <- myh10_ordered[order(myh10$m_10_expression),] # sort by expression
myh10_ordered$m10_cell_type <- factor(myh10_ordered$m10_cell_type) # it must be a factor
myh10_ordered

dotchart(myh10_ordered$m_10_expression,labels=m10_cell_type,cex=.7,
         main="MYH10",
         xlab="MYH10 Expression ")


p<-ggplot(myh10_ordered, aes(x=m_10_expression, y=cell_type)) + 
  geom_dotplot(binaxis='y', stackdir='center')
p

# # Change dotsize and stack ratio
# q <- ggplot(x, aes(x=expression, y=cell_type)) + 
#   geom_dotplot(binaxis='y', stackdir='center',
#                stackratio=1.5, dotsize=1.5)

# Rotate the dot plot
p + coord_flip()

# Change color by groups
dp <-ggplot(x, aes(x=expression, y=cell_type, fill=expression)) + 
  geom_dotplot(binaxis='y', stackdir='center')+
  labs(title="MYH9 Expression in Kidney Cells",x="Expression", y = "Cell Type")
dp + theme_classic() 
