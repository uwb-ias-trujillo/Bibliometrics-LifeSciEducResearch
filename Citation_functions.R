require(lazyeval)

extract_citation <- function(x, column){

  citations <- str_split(x[, column], "; ")

  names(citations) <- row.names(x)

  citations2<-melt(citations)
  identical(names(citations[[1]]), names(citations[[2]]) )
  citations2<-cbind(citations2[2],citations2[1])
  # create a label
  names(citations2)[2] <- "label"
  names(citations2)[1] <- "L1"

  # create the citing author index
  citations2$index <- row.names(citations2)
  citations2$index <- citations2$index %>%
    str_replace(., "[\\.][0-9]+", "")

  return(citations2)

}




remove_doi <- function(data){

    doi_pattern <- "doi+"

    # find location of DOI and remove it
    location <- str_locate(data, doi_pattern)

    # subract 2 from string to get 2 places before start of DOI
    location2 <- location[, 1] - 2
    lengths <- str_length(data)

    # replace NAs with length of string for those that don't have a DOI
    location3 <- ifelse(is.na(location2), lengths, location2)

    # extract 1 through length of string
    label <- str_sub(data, 1, location3)

    return(label)
}



standardize <- function(data){


    # show authors who are longer than 50 characters

    # eliminate punctuation except for commas
    pun_pattern <- "[\\.\\-'\\[\\]_]+"
    # include a comma anchor at the end to avoid capturing a year within a group of numbers
    year_pattern <- "[1-2]{1}[0-9]{3}[,]"
    # create a label without author, year, volume, page number
    author_pattern <- "^[a-z][^,0-9]+"


   standard <- data %>%
       mutate(label = str_replace_all(label, pun_pattern, "")) %>%
       mutate(label = str_trim(label, side = "both")) %>%
       mutate(label = tolower(label)) %>%
       mutate(label = remove_doi(label)) %>%
       mutate(label = ifelse(nchar(label) < 3, "",  label)) %>%
       mutate(label = ifelse(is.na(label), "", label)) %>%
       mutate(documentYear = str_extract(label, year_pattern)) %>%
       mutate(documentYear = str_replace(documentYear, ",", "")) %>%
       mutate(documentYear = as.numeric(documentYear)) %>%
       mutate(volume = str_extract(label, "v[0-9]+")) %>%
       mutate(page = str_extract(label, "p[0-9]+")) %>%
       mutate(authors = str_extract(label, author_pattern)) %>%
       mutate(documentName = str_replace(label, author_pattern, "")) %>%
       mutate(documentName = str_replace(documentName, year_pattern, "")) %>%
       mutate(documentName = str_replace(documentName, "in press", "")) %>%
       mutate(documentName = str_replace(documentName, "v[0-9]+", "")) %>%
       mutate(documentName = str_replace(documentName, "p[0-9]+", "" )) %>%
       mutate(documentName = str_replace(documentName, ",", "")) %>%
       mutate(documentName = str_replace(documentName, "[0-9]", "")) %>%
       mutate(documentName = str_replace_all(documentName, "[,]", "")) %>%
       mutate(type = ifelse(is.na(volume), "book", "article")) %>%
       mutate(cleanLabel = str_replace(label, ",$", ""))

   return(standard)

}





fuzzy_match <- function(input_data, threshold, label){

    input_data <- as.character(input_data)

    # replace NA with a character "NA"
    index <- which(is.na(label))
    label[index] <- "NA"

    index2 <- which(is.na(input_data))
    input_data[index2] <- "NA"

    if(length(input_data) >= 2){

    # compute the similarity match
    jw_dist <- stringdistmatrix(input_data, input_data, method = "jw")

    ## use single linkage method as it suits the problem
    hc <- hclust(as.dist(jw_dist), method = 'single')

    #plot(hc, labels = input_data)
    #abline(h = threshold, col = "blue")

    ## based on threshold, get clusters and make labels
    cluster <- cutree(hc, h = threshold)

    cluster_labels <- sapply(unique(cluster), function(y) label[min(which(cluster == y))])
    (output_data <- cluster_labels[cluster])

    return(output_data)

    }else{

    return(label)

    }

}


calc_jw <- function(a, b){1 - stringdist(a, b, method = "jw")}
