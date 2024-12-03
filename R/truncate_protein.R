### truncate_protein
#' Truncate large proteins to be able to draw them along with shorter ones 
#' Or truncate to focus on the specific regions of proteins
#' "...trunc..." text will be placed on the regions that are truncated
#'
#' @param data2Parse - data containing the proteins obtained using get_features function from drawProteins package
#' @param prot2Trunc - accession of the protein to be truncated
#' @param indices - specifies the regions to be truncated e.g. list(c(1, 500), c(600, 1000))
#' @param spaceBetween - Dimensions of the filler space to be placed instead of the regions that were truncated
#'
#' @return A data.frame containing chain structure of the truncated protein and the rest of the proteins in the input
#'
#' @examples
#' truncate_protein(data2Parse = data, prot2Trunc = "Q8WZ42", indices = list(c(100, 2000), c(2500, 10000), c(12000, 20000)), spaceBetween = 200)
#'
truncate_protein <- function(data2Parse, prot2Trunc, indices, spaceBetween = 200){
  if(!is.list(indices)){
    stop("Indices must be supplied in a list format, with begin and end points to truncate specified. e.g. list(c(1, 500), c(600, 1000))")
  }
  data2adjust <- data2Parse[Accession == prot2Trunc,]
  dataRetained <- data2Parse[Accession != prot2Trunc,]
  CHAIN <- data2adjust[type == 'CHAIN'][1]
  
  for(i in 1:length(indices)){
    beginInd <- indices[[i]][1]
    endInd <- indices[[i]][2]
    
    if(i == 1){
      data2adjust[type == 'CHAIN']$description <- 'Truncation'
      data2adjust[type == 'CHAIN']$Source <- 'Truncation'
    }
    data2adjust[type == 'CHAIN']$end[i] <- beginInd
    data2adjust[type == 'CHAIN']$length[i] <- data2adjust[type == 'CHAIN']$end[i] - data2adjust[type == 'CHAIN']$begin[i] + 1
    
    data2adjust <- data2adjust[(begin < beginInd) | (end > endInd)]
    data2adjust_N <- data2adjust[begin < beginInd]
    data2adjust_C <- data2adjust[end > endInd]
    
    if(nrow(data2adjust_N) == 0 & nrow(data2adjust_C) == 0){
      stop("Check your indices, no data retained after truncation.")
    }
    
    if(nrow(data2adjust_N) == 0){
      if(nrow(data2adjust_C[begin < endInd]) > 0){
        data2adjust_C[begin < endInd]$description <- NA
        data2adjust_C[begin < endInd]$begin <- endInd
        data2adjust_C$length <- data2adjust_C$end - data2adjust_C$begin + 1
      }
      data2adjust_C$begin <- data2adjust_C$begin - endInd + beginInd + spaceBetween
      data2adjust_C$end <- data2adjust_C$end - endInd + beginInd + spaceBetween
      data2adjust <- data2adjust_C
    }
    if(nrow(data2adjust_C) == 0){
      if(nrow(data2adjust_N[end > beginInd]) > 0){
        data2adjust_N[end > beginInd]$description <- NA
        data2adjust_N[end > beginInd]$end <- beginInd
        data2adjust_N$length <- data2adjust_N$end - data2adjust_N$begin + 1
      }
      data2adjust <- data2adjust_N
    }
    if(nrow(data2adjust_N) > 0 & nrow(data2adjust_C) > 0){
      if(nrow(data2adjust_C[begin < endInd]) > 0){
        data2adjust_C[begin < endInd]$description <- NA
        data2adjust_C[begin < endInd]$begin <- endInd
        data2adjust_C$length <- data2adjust_C$end - data2adjust_C$begin + 1
      }
      if(nrow(data2adjust_N[end > beginInd]) > 0){
        data2adjust_N[end > beginInd]$description <- NA
        data2adjust_N[end > beginInd]$end <- beginInd
        data2adjust_N$length <- data2adjust_N$end - data2adjust_N$begin + 1
      }
      data2adjust_C$begin <- data2adjust_C$begin - endInd + beginInd + spaceBetween
      data2adjust_C$end <- data2adjust_C$end - endInd + beginInd + spaceBetween
      data2adjust <- rbind(data2adjust_N, data2adjust_C)
    }
    
    newCHAIN <- CHAIN
    newCHAIN$description <- 'Truncation'
    newCHAIN$Source <- 'Truncation'
    newCHAIN$begin <- beginInd + spaceBetween
    newCHAIN$end <- newCHAIN$end - endInd + beginInd + spaceBetween
    newCHAIN$length <- newCHAIN$end - newCHAIN$begin + 1
    newCHAIN$ShortID <- paste(newCHAIN$protName, 'Truncation', newCHAIN$begin, newCHAIN$end, sep = '_')
    newCHAIN$FunFamsID <- ''
    
    truncInfo <- CHAIN
    truncInfo$Source <- 'Truncation'
    truncInfo$type <- 'TRUNCATION'
    truncInfo$begin <- beginInd + 1 + (spaceBetween/2)
    truncInfo$end <- beginInd + spaceBetween - 1
    truncInfo$length <- truncInfo$end - truncInfo$begin + 1
    truncInfo$ShortID <- paste(truncInfo$protName, 'Truncation', (beginInd+1), (beginInd+spaceBetween-1), sep = '_')
    truncInfo$FunFamsID <- ''
    truncInfo$description <- paste('...trunc...\n', (endInd - beginInd + 1), 'aa', sep = '')
    
    data2adjust <- rbind(data2adjust, truncInfo, newCHAIN)
    CHAIN <- newCHAIN
  }
  data2Parse <- rbind(data2adjust, dataRetained)
  data2Parse <- data2Parse[order(begin),]
  return(data2Parse)
}