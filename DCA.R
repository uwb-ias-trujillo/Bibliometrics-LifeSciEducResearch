##### Install packages ###
install.packages(c("tidyverse","devtools","rlang")) # to install the bibliometrix most recent version from GITHUB we use the devtools
devtools::install_github("massimoaria/bibliometrix") # next we pull the code from github to install the package.
library('tidyverse') # load tidyverse package
library('bibliometrix')

#### UPLOAD DATA ###
#Sample provided by biblionetrix package
scientometrics_text <- readFiles('http://www.bibliometrix.org/datasets/scientometrics.txt')
data(scientometrics_text)
scient_df <- isi2df(scientometrics_text)

#Working with CBE-LSE data. 
isi<-readFiles('rawdata/CBE_LSE_2008-2019/CBE_LSErecords1-733.txt')
data <- convert2df(isi,dbsource = "wos", format = "plaintext")

