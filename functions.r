conceptualStructure <- function (M, field = "ID", ngrams = 1, method = "MCA", quali.supp = NULL, 
          quanti.supp = NULL, minDegree = 2, clust = "auto", k.max = 5, 
          stemming = FALSE, labelsize = 10, documents = 2, graph = TRUE, 
          remove.terms = NULL, synonyms = NULL) 
{
  cbPalette <- c(brewer.pal(9, "Set1")[-6], brewer.pal(8, "Set2")[-7], 
                 brewer.pal(12, "Paired")[-11], brewer.pal(12, "Set3")[-c(2, 
                                                                          8, 12)])
  if (!is.null(quali.supp)) {
    QSUPP = data.frame(M[, quali.supp])
    names(QSUPP) = names(M)[quali.supp]
    row.names(QSUPP) = tolower(row.names(M))
  }
  if (!is.null(quanti.supp)) {
    SUPP = data.frame(M[, quanti.supp])
    names(SUPP) = names(M)[quanti.supp]
    row.names(SUPP) = tolower(row.names(M))
  }
  binary = FALSE
  if (method == "MCA") {
    binary = TRUE
  }
  switch(field, ID = {
    CW <- cocMatrix(M, Field = "ID", type = "matrix", sep = ";", 
                    binary = binary, remove.terms = remove.terms, synonyms = synonyms)
    CW = CW[, colSums(CW) >= minDegree]
    CW = CW[, !(colnames(CW) %in% "NA")]
    CW = CW[rowSums(CW) > 0, ]
  }, DE = {
    CW <- cocMatrix(M, Field = "DE", type = "matrix", sep = ";", 
                    binary = binary, remove.terms = remove.terms, synonyms = synonyms)
    CW = CW[, colSums(CW) >= minDegree]
    CW = CW[rowSums(CW) > 0, ]
    CW = CW[, !(colnames(CW) %in% "NA")]
  }, ID_TM = {
    M = termExtraction(M, Field = "ID", remove.numbers = TRUE, 
                       stemming = stemming, language = "english", remove.terms = remove.terms, 
                       synonyms = synonyms, keep.terms = NULL, verbose = FALSE)
    CW <- cocMatrix(M, Field = "ID_TM", type = "matrix", 
                    sep = ";", binary = binary)
    CW = CW[, colSums(CW) >= minDegree]
    CW = CW[, !(colnames(CW) %in% "NA")]
    CW = CW[rowSums(CW) > 0, ]
  }, DE_TM = {
    M = termExtraction(M, Field = "DE", remove.numbers = TRUE, 
                       stemming = stemming, language = "english", remove.terms = remove.terms, 
                       synonyms = synonyms, keep.terms = NULL, verbose = FALSE)
    CW <- cocMatrix(M, Field = "DE_TM", type = "matrix", 
                    sep = ";", binary = binary)
    CW = CW[, colSums(CW) >= minDegree]
    CW = CW[, !(colnames(CW) %in% "NA")]
    CW = CW[rowSums(CW) > 0, ]
  }, TI = {
    M = termExtraction(M, Field = "TI", remove.numbers = TRUE, 
                       stemming = stemming, language = "english", remove.terms = remove.terms, 
                       synonyms = synonyms, keep.terms = NULL, verbose = FALSE, 
                       ngrams = ngrams)
    CW <- cocMatrix(M, Field = "TI_TM", type = "matrix", 
                    sep = ";", binary = binary)
    CW = CW[, colSums(CW) >= minDegree]
    CW = CW[, !(colnames(CW) %in% "NA")]
    CW = CW[rowSums(CW) > 0, ]
  }, AB = {
    M = termExtraction(M, Field = "AB", remove.numbers = TRUE, 
                       stemming = stemming, language = "english", remove.terms = remove.terms, 
                       synonyms = synonyms, keep.terms = NULL, verbose = FALSE, 
                       ngrams = ngrams)
    CW <- cocMatrix(M, Field = "AB_TM", type = "matrix", 
                    sep = ";", binary = binary)
    CW = CW[, colSums(CW) >= minDegree]
    CW = CW[rowSums(CW) > 0, ]
    CW = CW[, !(colnames(CW) %in% "NA")]
  })
  colnames(CW) = tolower(colnames(CW))
  rownames(CW) = tolower(rownames(CW))
  p = dim(CW)[2]
  quali = NULL
  quanti = NULL
  if (!is.null(quali.supp)) {
    ind = which(row.names(QSUPP) %in% row.names(CW))
    QSUPP = as.data.frame(QSUPP[ind, ])
    CW = cbind(CW, QSUPP)
    quali = (p + 1):dim(CW)[2]
    names(CW)[quali] = names(M)[quali.supp]
  }
  if (!is.null(quanti.supp)) {
    ind = which(row.names(SUPP) %in% row.names(CW))
    SUPP = as.data.frame(SUPP[ind, ])
    CW = cbind(CW, SUPP)
    quanti = (p + 1 + length(quali)):dim(CW)[2]
    names(CW)[quanti] = names(M)[quanti.supp]
  }
  results <- factorial(CW, method = method, quanti = quanti, 
                       quali = quali)
  res.mca <- results$res.mca
  df <- results$df
  docCoord <- results$docCoord
  df_quali <- results$df_quali
  df_quanti <- results$df_quanti
  if ("TC" %in% names(M) & method != "MDS") {
    docCoord$TC = as.numeric(M[toupper(rownames(docCoord)), 
                               "TC"])
  }
  km.res = hclust(dist(df), method = "average")
  if (clust == "auto") {
    clust = min((length(km.res$height) - which.max(diff(km.res$height)) + 
                   1), k.max)
  }
  else {
    clust = max(2, min(as.numeric(clust), k.max))
  }
  km.res$data = df
  km.res$cluster = cutree(km.res, k = clust)
  km.res$data.clust = cbind(km.res$data, km.res$cluster)
  names(km.res$data.clust)[3] = "clust"
  centers <- km.res$data.clust %>% group_by(.data$clust) %>% 
    summarise(Dim.1 = mean(.data$Dim.1), Dim.2 = mean(.data$Dim.2)) %>% 
    as.data.frame()
  km.res$centers = centers[, c(2, 3, 1)]
  b = fviz_cluster(km.res, stand = FALSE, data = df, labelsize = labelsize, 
                   repel = TRUE) + theme_minimal() + scale_color_manual(values = cbPalette[1:clust]) + 
    scale_fill_manual(values = cbPalette[1:clust]) + labs(title = paste("Conceptual Structure Map - method: ", 
                                                                        method, collapse = "", sep = "")) + geom_point() + geom_hline(yintercept = 0, 
                                                                                                                                      linetype = "dashed", color = adjustcolor("grey40", alpha.f = 0.7)) + 
    geom_vline(xintercept = 0, linetype = "dashed", color = adjustcolor("grey40", 
                                                                        alpha.f = 0.7)) + theme(text = element_text(size = labelsize), 
                                                                                                axis.title = element_text(size = labelsize, face = "bold"), 
                                                                                                plot.title = element_text(size = labelsize + 1, face = "bold"), 
                                                                                                panel.background = element_rect(fill = "white", colour = "white"), 
                                                                                                axis.line.x = element_line(color = "black", size = 0.5), 
                                                                                                axis.line.y = element_line(color = "black", size = 0.5), 
                                                                                                panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  if (method != "MDS") {
    b = b + xlab(paste("Dim 1 (", round(res.mca$eigCorr$perc[1], 
                                        2), "%)", sep = "")) + ylab(paste("Dim 2 (", round(res.mca$eigCorr$perc[2], 
                                                                                           2), "%)", sep = ""))
  }
  else {
    b = b + xlab("Dim 1") + ylab("Dim 2")
  }
  if (!is.null(quali.supp)) {
    s_df_quali = df_quali[(abs(df_quali[, 1]) >= quantile(abs(df_quali[, 
                                                                       1]), 0.75) | abs(df_quali[, 2]) >= quantile(abs(df_quali[, 
                                                                                                                                2]), 0.75)), ]
    names(s_df_quali) = c("x", "y")
    s_df_quali$label = row.names(s_df_quali)
    x = s_df_quali$x
    y = s_df_quali$y
    label = s_df_quali$label
    b = b + geom_point(aes(x = x, y = y), data = s_df_quali, 
                       colour = "red", size = 1) + geom_label_repel(aes(x = x, 
                                                                        y = y, label = label, size = 1), data = s_df_quali)
  }
  if (!is.null(quanti.supp)) {
    names(df_quanti) = c("x", "y")
    df_quanti$label = row.names(df_quanti)
    x = df_quanti$x
    y = df_quanti$y
    label = df_quanti$label
    b = b + geom_point(aes(x = x, y = y), data = df_quanti, 
                       colour = "blue", size = 1) + geom_label_repel(aes(x = x, 
                                                                         y = y, label = label, size = 1), data = df_quanti) + 
      geom_segment(data = df_quanti, aes(x = 0, y = 0, 
                                         xend = x, yend = y), size = 1.5, arrow = arrow(length = unit(0.3, 
                                                                                                      "cm")))
  }
  b = b + theme(legend.position = "none")
  coord_b <- plotCoord(b)
  if (isTRUE(graph)) {
    plot(b)
  }
  b_dend <- fviz_dend(km.res, rect = TRUE, k = clust, cex = labelsize/20, 
                      main = "Topic Dendrogram", k_colors = cbPalette[1:clust]) + 
    theme(plot.title = element_text(size = labelsize + 1, 
                                    face = "bold"), axis.title = element_text(size = labelsize, 
                                                                              face = "bold"), panel.background = element_rect(fill = "white", 
                                                                                                                              colour = "white"), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  coord <- plotCoord(b_dend, side = "u")
  
  if (isTRUE(graph)) {
    plot(b_dend)
  }
  if (method != "MDS") {
    if (documents > dim(docCoord)[1]) {
      documents = dim(docCoord)[1]
    }
    centers = data.frame(dim1 = km.res$centers[, 1], dim2 = km.res$centers[, 
                                                                           2])
    centers$color = cbPalette[1:dim(centers)[1]]
    row.names(centers) = paste("cluster", as.character(1:dim(centers)[1]), 
                               sep = "")
    A = euclDist(docCoord[, 1:2], centers)
    docCoord$Cluster = A$color
    A$color = cbPalette[A$color]
    A$contrib <- docCoord$contrib
    A <- A %>% mutate(names = row.names(A)) %>% group_by(.data$color) %>% 
      top_n(n = documents, wt = .data$contrib) %>% select(!"contrib") %>% 
      as.data.frame()
    row.names(A) <- A$names
    A <- A[, -4]
    names(centers) = names(A)
    A = rbind(A, centers)
    x = A$dim1
    y = A$dim2
    A[, 4] = row.names(A)
    names(A)[4] = "nomi"
    df_all = rbind(as.matrix(df), as.matrix(A[, 1:2]))
    rangex = c(min(df_all[, 1]), max(df_all[, 1]))
    rangey = c(min(df_all[, 2]), max(df_all[, 2]))
    b_doc <- ggplot(aes(x = .data$dim1, y = .data$dim2, label = .data$nomi), 
                    data = A) + geom_point(size = 2, color = A$color) + 
      labs(title = "Factorial map of the documents with the highest contributes") + 
      geom_label_repel(box.padding = unit(0.5, "lines"), 
                       size = (log(labelsize * 3)), fontface = "bold", 
                       fill = adjustcolor(A$color, alpha.f = 0.6), color = "white", 
                       segment.alpha = 0.5, segment.color = "gray") + 
      scale_x_continuous(limits = rangex, breaks = seq(round(rangex[1]), 
                                                       round(rangex[2]), 1)) + scale_y_continuous(limits = rangey, 
                                                                                                  breaks = seq(round(rangey[1]), round(rangey[2]), 
                                                                                                               1)) + geom_hline(yintercept = 0, linetype = "dashed", 
                                                                                                                                color = adjustcolor("grey40", alpha.f = 0.7)) + geom_vline(xintercept = 0, 
                                                                                                                                                                                           linetype = "dashed", color = adjustcolor("grey40", 
                                                                                                                                                                                                                                    alpha.f = 0.7)) + theme(text = element_text(size = labelsize), 
                                                                                                                                                                                                                                                            axis.title = element_text(size = labelsize, face = "bold"), 
                                                                                                                                                                                                                                                            plot.title = element_text(size = labelsize + 1, face = "bold"), 
                                                                                                                                                                                                                                                            panel.background = element_rect(fill = "white", colour = "white"), 
                                                                                                                                                                                                                                                            axis.line.x = element_line(color = "black", size = 0.5), 
                                                                                                                                                                                                                                                            axis.line.y = element_line(color = "black", size = 0.5), 
                                                                                                                                                                                                                                                            panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    if (method != "MDS") {
      b_doc = b_doc + xlab(paste("Dim 1 (", round(res.mca$eigCorr$perc[1], 
                                                  2), "%)", sep = "")) + ylab(paste("Dim 2 (", 
                                                                                    round(res.mca$eigCorr$perc[2], 2), "%)", sep = ""))
    }
    else {
      b_doc = b_doc + xlab("Dim 1") + ylab("Dim 2")
    }
    xl <- c(rangex[2] - 0.02 - diff(rangex) * 0.125, rangex[2] - 
              0.02)
    yl <- c(rangey[1], rangey[1] + diff(rangey) * 0.125) + 
      0.02
    
    if (isTRUE(graph)) {
      (plot(b_doc))
    }
    docCoord = docCoord[order(-docCoord$TC), ]
    B = euclDist(docCoord[, 1:2], centers)
    B$color = cbPalette[B$color]
    B$TC <- docCoord$TC
    B <- B %>% mutate(names = row.names(B)) %>% group_by(.data$color) %>% 
      top_n(n = documents, wt = .data$TC) %>% select(!"TC") %>% 
      as.data.frame()
    row.names(B) <- B$names
    B <- B[, -4]
    B = rbind(B, centers)
    x = B$dim1
    y = B$dim2
    B[, 4] = row.names(B)
    names(B)[4] = "nomi"
    df_all_TC = rbind(as.matrix(df), as.matrix(B[, 1:2]))
    rangex = c(min(df_all_TC[, 1]), max(df_all_TC[, 1]))
    rangey = c(min(df_all_TC[, 2]), max(df_all_TC[, 2]))
    b_doc_TC = ggplot(aes(x = .data$dim1, y = .data$dim2, 
                          label = .data$nomi), data = B) + geom_point(size = 2, 
                                                                      color = B$color) + labs(title = "Factorial map of the most cited documents") + 
      geom_label_repel(box.padding = unit(0.5, "lines"), 
                       size = (log(labelsize * 3)), fontface = "bold", 
                       fill = adjustcolor(B$color, alpha.f = 0.6), color = "white", 
                       segment.alpha = 0.5, segment.color = "gray") + 
      scale_x_continuous(limits = rangex, breaks = seq(round(rangex[1]), 
                                                       round(rangex[2]), 1)) + scale_y_continuous(limits = rangey, 
                                                                                                  breaks = seq(round(rangey[1]), round(rangey[2]), 
                                                                                                               1)) + xlab(paste("Dim 1 (", round(res.mca$eigCorr$perc[1], 
                                                                                                                                                 2), "%)", sep = "")) + ylab(paste("Dim 2 (", round(res.mca$eigCorr$perc[2], 
                                                                                                                                                                                                    2), "%)", sep = "")) + geom_hline(yintercept = 0, 
                                                                                                                                                                                                                                      linetype = "dashed", color = adjustcolor("grey60", 
                                                                                                                                                                                                                                                                               alpha.f = 0.7)) + geom_vline(xintercept = 0, 
                                                                                                                                                                                                                                                                                                            linetype = "dashed", color = adjustcolor("grey60", 
                                                                                                                                                                                                                                                                                                                                                     alpha.f = 0.7)) + theme(text = element_text(size = labelsize), 
                                                                                                                                                                                                                                                                                                                                                                             axis.title = element_text(size = labelsize, face = "bold"), 
                                                                                                                                                                                                                                                                                                                                                                             plot.title = element_text(size = labelsize + 1, face = "bold"), 
                                                                                                                                                                                                                                                                                                                                                                             panel.background = element_rect(fill = "white", colour = "white"), 
                                                                                                                                                                                                                                                                                                                                                                             axis.line.x = element_line(color = "black", size = 0.5), 
                                                                                                                                                                                                                                                                                                                                                                             axis.line.y = element_line(color = "black", size = 0.5), 
                                                                                                                                                                                                                                                                                                                                                                             panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    xl <- c(rangex[2] - 0.02 - diff(rangex) * 0.125, rangex[2] - 
              0.02)
    yl <- c(rangey[1], rangey[1] + diff(rangey) * 0.125) + 
      0.02
    if (isTRUE(graph)) {
      plot(b_doc_TC)
    }
    semanticResults = list(net = CW, res = res.mca, km.res = km.res, 
                           graph_terms = b, graph_dendogram = b_dend, graph_documents_Contrib = b_doc, 
                           graph_documents_TC = b_doc_TC, docCoord = docCoord)
  }
  else {
    semanticResults = list(net = CW, res = res.mca, km.res = km.res, 
                           graph_terms = b, graph_dendogram = b_dend, graph_documents_Contrib = NULL, 
                           graph_documents_TC = NULL, docCoord = NULL)
  }
  return(semanticResults)
}


authorProdOverTime <- function (M, k = 10, graph = TRUE)
{
  if (!("DI" %in% names(M))) {
    M$DI = "NA"
    }

  M$TC <- as.numeric(M$TC)
  M$PY <- as.numeric(M$PY)
  M <- M[!is.na(M$PY), ]
  Y <- as.numeric(substr(Sys.time(), 1, 4))
  listAU <- (strsplit(M$AU, ";"))
  nAU <- lengths(listAU)

  df <- data.frame(AU = trimws(unlist(listAU)), SR = rep(M$SR,nAU))
  AU <- df %>% group_by(.data$AU) %>% count() %>% arrange(desc(.data$n)) %>% 
    ungroup()
  k <- min(k, nrow(AU))
  AU <- AU %>% slice_head(n = k)

  df <- df %>% right_join(AU, by = "AU") %>% left_join(M, by = "SR") %>% 
    select(.data$AU.x, .data$PY, .data$TI, .data$SO, .data$DI, .data$TC) %>% 
    mutate(TCpY = .data$TC/(Y - .data$PY + 1)) %>% 
    group_by(.data$AU.x) %>% 
    mutate(n = length(.data$AU.x)) %>% 
    ungroup() %>% 
    rename(Author = .data$AU.x, year = .data$PY, DOI = .data$DI) %>% 
    arrange(desc(.data$n), desc(.data$year)) %>% 
    select(-.data$n)

  df2 <- dplyr::group_by(df, .data$Author, .data$year) %>% 
    dplyr::summarise(freq = length(.data$year), TC = sum(.data$TC), TCpY = sum(.data$TCpY)) %>% 
    as.data.frame()

  df2$Author <- factor(df2$Author, levels = AU$AU[1:k])
  x <- c(0.5, 1.5 * k/10)
  y <- c(min(df$year), min(df$year) + diff(range(df2$year)) * 0.125)
  #data("logo", envir = environment())
  #logo <- grid::rasterGrob(logo, interpolate = TRUE)

  g <- ggplot(df2, aes(x = .data$Author, y = .data$year, 
                       text = paste("Author: ", .data$Author,
                                    "\nYear: ", .data$year,
                                    "\nN. of Articles: ", .data$freq, 
                                    "\nTotal Citations per Year: ", round(.data$TCpY, 2)))) + 
    geom_point(aes(alpha = .data$TCpY, size = .data$freq), color = "dodgerblue4") + 
    scale_size(range = c(2, 6)) + 
    scale_alpha(range = c(0.3, 1)) + 
    scale_y_continuous(breaks = seq(min(df2$year), max(df2$year), by = 2)) + 
    guides(size = guide_legend(order = 1, "Number of\nArticles"), 
           alpha = guide_legend(order = 2, "Times Cited\nper Year")) + 
    theme(legend.position = "right", 
          text = element_text(color = "#444444"), 
          panel.background = element_rect(fill = "#FFFFFF"), 
          plot.title = element_text(size = 24), 
          axis.title = element_text(size = 14, color = "#555555"), 
          axis.title.x = element_text(hjust = 0.95), 
          axis.title.y = element_text(vjust = 1, angle = 90), 
          axis.text.x = element_text(face = "bold", angle = 90), 
          axis.text.y = element_text(face = "bold"), 
          axis.line.x = element_line(color = "grey50", size = 0.5), 
          panel.grid.major.x = element_blank(), 
          panel.grid.major.y = element_line(size = 0.2, color = "grey90")) + 
    labs(title = "Top-Authors' Production over Time", x = "Author", y = "Year") + 
    geom_line(data = df2, aes(x = .data$Author, y = .data$year, group = .data$Author), 
              size = 1, color = "firebrick4", alpha = 0.3) + 
    scale_x_discrete(limits = rev(levels(df2$Author))) + 
    coord_flip() #+ annotation_custom(logo, xmin = x[1], xmax = x[2], ymin = y[1], ymax = y[2])
  df$DOI = as.character(df$DOI)
  res <- list(dfAU = df2, dfPapersAU = df, graph = g)
  if (isTRUE(graph)) {
    plot(g)
  }
  return(res)
}
