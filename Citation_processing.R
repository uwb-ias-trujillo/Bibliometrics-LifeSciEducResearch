
setwd("~/Workspace/WebOfScience-CoCitation/")


# load -----------------------------------------------------------------------------------
# load libraries 
install.packages(c("stringdist", "stringr", "stringi", "dplyr", "magrittr", "broom", "lazyeval", "reshape2"))

library(stringdist)
library(stringr)
library(stringi)
library(plyr)
library(dplyr)
library(magrittr)
library(broom)
library(lazyeval)
library(tidyr)
library(reshape2)

# load functions 
  source("citation_functions.R")
  
# load raw data from WOS as a csv format
    raw_data <- read.csv("sampleWOSdata.csv", stringsAsFactors = FALSE)
    abb <- read.csv("WOS_journalAbbrev.csv")

# Extract citations from WOS list 
    work_data <- extract_citation(raw_data, "CR")

# Standardize the data-------------------------------------------------------------------- 
   std_data <- standardize(work_data)

# show authors who are longer than 50 characters
   author_temp <- std_data %>% 
        select(authors) %>% 
        filter(nchar(authors) > 50) 
    
# clean and standardize author names -----------------------------------------------------    

# Run fuzzy match function on a .2 threshold 

    std_data$authorBlock <- fuzzy_match(std_data$authors, .25, std_data$authors)
    
    # look at author clusters  
    author_merges <- std_data %>% 
      group_by(authorBlock) %>% 
      summarise(count = n()) %>% 
      filter(count > 15) %>% 
      arrange(desc(count))
    
    
    # Run the clustering algorithm again to recluster top merge
    # get the top merged author  
    name <- author_merges[1, ][[1]]
    
    std_data2 <- std_data %>% 
      filter(authorBlock == name) %>% 
      mutate(authorBlock = fuzzy_match(authors, .20, authors))
    
    std_data3 <- std_data %>% 
      filter(authorBlock != name) %>% 
      bind_rows(std_data2)
    
    std_data4 <- std_data3 %>% 
      filter(authorBlock == name) %>% 
      mutate(authorBlock = fuzzy_match(authors, .15, authors)) 
     
    std_data5 <- std_data3 %>% 
      filter(authorBlock != name) %>% 
      bind_rows(std_data4)

#Add specific criteria for problematic author names 
    std_data5$authorBlock[grep("assaraf", std_data5$authorBlock)]


# Run clustering algorithm----------------------------------------------------------------          
    books <- std_data5 %>% 
        group_by(authorBlock) %>% 
        filter(type == "book") %>% 
        mutate(mergedLabel = fuzzy_match(documentName, .15, cleanLabel))
    
    
    articles <- std_data5 %>% 
        group_by(authorBlock, documentYear) %>% 
        filter(type == "article") %>% 
        mutate(mergedLabel = fuzzy_match(documentName, .25, cleanLabel))
    
        clean_data <- rbind(books, articles)


# create the flat data file  
    flat_data <- clean_data %>% 
        mutate( L1 = as.numeric(L1)) %>% 
        group_by(L1) %>% 
        arrange(mergedLabel) %>% 
        do(data.frame(citations = str_c(.$mergedLabel, collapse = "| ")))
    
# make sure to unwrap this again to check if it was wrapped properly 
    flat_data2 <- cbind(raw_data, flat_data[2])
    

# find the similarity  of the merges------------------------------------------------------ 

    clean_data <- clean_data %>% 
        mutate(jwSimilarity = calc_jw(cleanLabel, mergedLabel))


# plot------------------------------------------------------------------------------------ 

    merged <- clean_data %>% 
        filter(jwSimilarity < 1)
    
    
    hist(round(merged$jwSimilarity, 3),
         xlab = "jw similarity distance",
         main = "Freq Merges by similarity distance \n with year and document type boundary")
    



# save the data------------------------------------------------------------------------------------     
    write.csv(clean_data, "WOS_nodeList.csv", row.names = FALSE)
    write.csv(flat_data2, "WOS_clean.csv", row.names = FALSE)
