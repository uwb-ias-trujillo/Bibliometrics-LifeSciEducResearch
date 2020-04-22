#This script launches the biblioshiny script (previously named DCA.R)
##### Install packages ###
install.packages(c("tidyverse","devtools","rlang", "bibliometrix")) # to install the bibliometrix most recent version from GITHUB we use the devtools
library('tidyverse') # load tidyverse package
library('bibliometrix')

biblioshiny()
