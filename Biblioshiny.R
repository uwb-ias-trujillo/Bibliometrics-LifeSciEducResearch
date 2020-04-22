#This script launches the biblioshiny script (previously named DCA.R)
##### Install packages ###
install.packages(c("tidyverse","devtools","rlang")) # to install the bibliometrix most recent version from GITHUB we use the devtools
devtools::install_github("massimoaria/bibliometrix") # next we pull the code from github to install the package.
library('tidyverse') # load tidyverse package
library('bibliometrix')

biblioshiny()
