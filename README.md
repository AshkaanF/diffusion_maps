# Diffusion map
#### A test inspired by the recent paper by [Barter and Gross, 2019](https://royalsocietypublishing.org/doi/full/10.1098/rspa.2018.0615)
Contains R scripts that recreate diffusion map using a subset of the [Human Microbiome Project 2](https://portal.hmpdacc.org/) metagenomic species-level characterization. Species were quantified using the marker gene-based [MetaPhlAn2](https://bitbucket.org/biobakery/metaphlan2) software.

#### Example

```R
##
## let's make a quick example
##
## point to your directory
setwd('path/to/directory')

## load the functions we'll need
source('./R/accessory_functions.R')

##
## load some data for testing
##
# load the mapping file from the human microbiome project
meta <- read.csv('./data/hmp_map.csv', header = T, row.names = 1, colClasses = 'character')

## load hmp species-level feature array
m <- read.csv('./data/hmp_species.csv', header = T, row.names = 1)

##
## NOTE: m is a sample x feature array
##
## get the normalized laplacian from the feature table
Lij <- m %>%
  norm.mat() %>%                ## normalize matrix
  get.euc() %>%                 ## get euclidean similarities
  threshold(., top_k = 25) %>%  ## threshold the distance matrix
  get.laplac()                  ## calculate normalized laplacian

## calculate eigenvals/vecs
eig <- Lij %>%
  eigen()
  
## get eigenvalues
evl <- eig %$%
  values %>%
  round(., 10)

## get eigenvectors
evc <- eig %$%
  vectors %>%
  round(., 10)

## get eigenvectors for 2 smallest non-zero evs
dim.1 <- evc[, which(order(evl, decreasing = F) == 2)]
dim.2 <- evc[, which(order(evl, decreasing = F) == 3)]

## append to data
meta$dim.1 <- dim.1
meta$dim.2 <- dim.2

##
## plot
##
ggplot(meta, aes(x = dim.1, y = dim.2, fill = env_2)) +
  theme_classic() +
  xlab('Dimension 1') +
  ylab('Dimension 2') +
  geom_hline(yintercept = 0, linetype = 1, size = 0.5, colour = '#959595') +
  geom_vline(xintercept = 0, linetype = 1, size = 0.5, colour = '#959595') +
  geom_point(shape = 21, alpha = 0.5, size = 2) +
  scale_fill_manual(values = brewer.pal(n = 4, name = 'RdYlBu'), name = NULL) +
  theme(axis.title = element_text(size = 11, colour = '#000000'),
        axis.text = element_text(size = 8, colour = '#000000'),
        legend.position = c(1, 0),
        legend.background = element_blank(),
        legend.justification = c(1, 0),
        legend.key.width = unit(0.2, 'cm'),
        legend.key.height = unit(0.5, 'cm'), 
        legend.text = element_text(size = 5),
        legend.text.align = 0,
        legend.title = element_text(size = 9))

```

Plotting the 2 smallest non-zero eigenvalues shows that the HMP metagenomes are roughly ordered from Gut, to Oral, Skin, and Vaginal microbial communities.


![HMP Example](figures/hmp.png)






