#' Create ggplot2 object with protein chains from feature database
#' Contains personal changes to the original package from brennanpincardiff/drawProteins because they were more convenient in my own specific project.

### draw_canvas
#' Difference from the original package function (brennanpincardiff/drawProteins): 
#' Larger x and y limits in the canvas
#'
draw_canvas <- function(data){
  begin=end=NULL
  p <- ggplot2::ggplot()
  p <- p + ggplot2::ylim(0.25, max(data$order)+0.25)
  p <- p + ggplot2::xlim(-max(data$end, na.rm=TRUE)*0.05,
                         max(data$end, na.rm=TRUE) + max(data$end, na.rm=TRUE)*0.05)
  p <- p + ggplot2::labs(x = "Amino acid number") # label x-axis
  # p <- p + ggplot2::labs(x = "") # label x-axis
  p <- p + ggplot2::labs(y = "") # label y-axis

  return(p)
}

### draw_chains
#' Differences from the original package function (brennanpincardiff/drawProteins): 
#' Larger font size to label chains, and also larger rectangle lines
#' Chains are outlined for better visualisation
#' Modified to be able to draw large proteins along with smaller proteins. "TRUNCATION" marks truncated proteins
#'
draw_chains <- function(p,
                        data = data,
                        outline = "black",
                        fill = NA,
                        label_chains = TRUE,
                        label_trunc = TRUE,
                        labels_chain = data[data$type == "CHAIN",]$entryName,
                        labels_trunc = data[data$type == "TRUNCATION",]$description,
                        size = 2.0,
                        alpha = 1.0,
                        label_size = 8){
  
  begin=end=NULL
  p <- p + ggplot2::geom_rect(data = data[data$type == "CHAIN",],
                              mapping=ggplot2::aes(xmin=begin,
                                                   xmax=end,
                                                   ymin=order-0.2,
                                                   ymax=order+0.2),
                              colour = outline,
                              fill = fill,
                              size = size,
                              alpha = alpha)
  
  if(label_chains == TRUE){
    p <- p +
      ggplot2::annotate("text", x = -25,
                        y = data[data$type == "CHAIN",]$order,
                        label = labels_chain,
                        hjust = 1,
                        size = label_size)
  }
  if(label_trunc == TRUE){
    if(nrow(data[type == 'TRUNCATION'] > 0)){
      p <- p +
        ggplot2::annotate("text", x = data[data$type == "TRUNCATION",]$begin,
                          y = data[data$type == "TRUNCATION",]$order,
                          label = labels_trunc,
                          size = label_size)
    }
  }
  return(p)
}

### draw_domains
#' Differences from the original package function (brennanpincardiff/drawProteins): 
#' Larger font size to label domains
#' colorPicker column in the data.frame to color domains independently
#' Domains are labeled on an outlined white background to improve visibility on both dark and light coloured domains.
#' Domains text is rotated 90 degrees for a different visualisation
#' @import shadowtext
#' @import statebins
#'
draw_domains <- function(p,
                         data = data,
                         label_domains = TRUE,
                         label_size = 8,
                         alpha = 1.0,
                         show.legend = TRUE){
  begin=end=description=NULL
  tmpData <- data[data$type == 'DOMAIN',]
  colorPicker <- tmpData$colorPicker
  p <- p + statebins:::geom_rrect(data= tmpData,
                              mapping=ggplot2::aes(xmin=begin,
                                                   xmax=end,
                                                   ymin=order-0.25,
                                                   ymax=order+0.25
                                                   ),
                              colour = 'black',
                              fill=colorPicker,
                              alpha = alpha,
                              show.legend = show.legend)
  
  if(label_domains == TRUE){
    p <- p + shadowtext::geom_shadowtext(data = tmpData,
                                         ggplot2::aes(x = (begin + (end-begin)/2 - 5),
                                                      y = order,
                                                      label = description),
                                         size = label_size, angle = 90,
                                         bg.colour = "white",  # Color of the outline
                                         bg.r = 0.15,          # Size of the outline
                                         colour = "black"      # Color of the text
    )
  }
  return(p)
}

### draw_phospho
#' Same as the original package function (brennanpincardiff/drawProteins)
#'
draw_phospho <- function(p, data = data,
                        size = 2,
                        fill = "yellow",
                        alpha = 1.0,
                        show.legend = FALSE){
    begin=end=description=NULL
    p <- p + ggplot2::geom_point(data = phospho_site_info(data),
                                ggplot2::aes(x = begin,
                        y = order+0.25),
                        shape = 21,
                        colour = "black",
                        fill = fill,
                        size = size,
                        alpha = alpha,
                        show.legend = show.legend)
    return(p)
}

### draw_regions
#' Differences from the original package function (brennanpincardiff/drawProteins): 
#' colorPicker column in the data.frame to color regions independently
#' Regions are labeled on an outlined white background to improve visibility on both dark and light coloured domains.
#'
draw_regions <- function(p,
                         data = data,
                         label_regions = TRUE,
                         label_size = 4,
                         alpha = 1.0,
                         show.legend = TRUE,
                         type = "REGION"){
  begin=end=description=NULL
  tmpData <- data[data$type == 'REGION',]
  colorPicker <- tmpData$colorPicker
  p <- p + statebins:::geom_rrect(data= tmpData,
                                  mapping=ggplot2::aes(xmin=begin,
                                                       xmax=end,
                                                       ymin=order-0.2,
                                                       ymax=order+0.2
                                  ),
                                  fill=colorPicker,
                                  alpha = alpha,
                                  show.legend = show.legend)
  
  if(label_regions == TRUE){
    p <- p + shadowtext::geom_shadowtext(data = tmpData,
                                         ggplot2::aes(x = (begin + (end-begin)/2 - 5),
                                                      y = order,
                                                      label = description),
                                         size = label_size, angle = 90,
                                         bg.colour = "white",  # Color of the outline
                                         bg.r = 0.15,          # Size of the outline
                                         colour = "black"      # Color of the text
    )
  }
  return(p)
}

### draw_motif
#' Same as the original package function (brennanpincardiff/drawProteins)
#'
draw_motif <- function(p, data = data, alpha = 1.0, show.legend = TRUE){
    begin=end=description=NULL
    ## plot motifs fill by description
    p <- p + ggplot2::geom_rect(data= data[data$type == "MOTIF",],
                                mapping=ggplot2::aes(xmin=begin,
                                xmax=end,
                                ymin=order-0.25,
                                ymax=order+0.25,
                                fill=description),
                                alpha = alpha,
                                show.legend = show.legend)

    return(p)
}

### draw_repeat
#' Differences from the original package function (brennanpincardiff/drawProteins): 
#' Repeats are labeled on an outlined white background to improve visibility on both dark and light coloured domains.
#' Repeat text is rotated 90 degrees for a different visualisation
#'
draw_repeat <- function(p, data = data,
                        label_size = 0.75,
                        outline = "dimgrey",
                        fill = "dimgrey",
                        alpha = 1.0,
                        label_repeats = TRUE,
                        show.legend = TRUE){
  begin=end=description=NULL
  tmpData <- data[data$type == "REPEAT",]
  ## step 6 plot repeats fill by description
  p <- p + ggplot2::geom_rect(data= tmpData,
                              mapping=ggplot2::aes(xmin=begin,
                                                   xmax=end,
                                                   ymin=order-0.25,
                                                   ymax=order+0.25),
                              colour = outline,
                              fill = fill,
                              alpha = alpha,
                              show.legend = show.legend)
  
  if(label_repeats == TRUE){
    # label repeats (for this they are ANK but remove digits)
    p <- p + shadowtext::geom_shadowtext(data = tmpData,
                                         ggplot2::aes(x = (begin + (end-begin)/2),
                                                      y = order,
                                                      label = gsub("\\d", "", description)),
                                         size = label_size, angle = 90,
                                         bg.colour = "white",  # Color of the outline
                                         bg.r = 0.15,          # Size of the outline
                                         colour = "black"      # Color of the text
    )
  }
  return(p)
}

### draw_recept_dom
#' Same as the original package function (brennanpincardiff/drawProteins)
#'
draw_recept_dom <- function(p,
                            data = data,
                            alpha = 1.0,
                            label_domains = FALSE,
                            label_size = 4,
                            show.legend = TRUE){
    begin=end=description=NULL

    p <- p + ggplot2::geom_rect(data= data[data$type == "TOPO_DOM",],
                            mapping=ggplot2::aes(xmin=begin,
                                xmax=end,
                                ymin=order-0.25,
                                ymax=order+0.25,
                                fill=description),
                            alpha = alpha,
                            show.legend = show.legend)

    p <- p + ggplot2::geom_rect(data= data[data$type == "TRANSMEM",],
                            mapping=ggplot2::aes(xmin=begin,
                                xmax=end,
                                ymin=order-0.25,
                                ymax=order+0.25,
                                fill=description),
                            alpha = alpha,
                            show.legend = show.legend)

    if(label_domains == TRUE){
        p <- p + ggplot2::geom_label(data = data[data$type == "TOPO_DOM", ],
                            ggplot2::aes(x = begin + (end-begin)/2,
                                y = order,
                                label = description),
                                size = label_size)

        p <- p + ggplot2::geom_label(data = data[data$type == "TRANSMEM", ],
                            ggplot2::aes(x = begin + (end-begin)/2,
                                y = order,
                                label = "TM"),
                                size = label_size)
    }

    return(p)
}

### draw_folding
#' Differences from the original package function (brennanpincardiff/drawProteins):
#' Fixed colours for strands, helices and turns
#' Slightly adjusted y coordinates for the location of the rectangles pointing folds
#'
draw_folding <- function(p,
                         data = data,
                         show.legend = FALSE,
                         show_strand = TRUE,
                         show_helix = TRUE,
                         show_turn = TRUE,
                         alpha = 1.0){
  begin=end=description=type=NULL
  # STRAND first
  if(show_strand == TRUE){
    p <- p + ggplot2::geom_rect(data= data[data$type == "STRAND",],
                                mapping=ggplot2::aes(xmin=begin,
                                                     xmax=end,
                                                     ymin=order-0.32,
                                                     ymax=order-0.27),
                                fill='#E66100',
                                alpha = alpha,
                                show.legend = show.legend)
  }
  # then HELIX
  if(show_helix == TRUE){
    p <- p + ggplot2::geom_rect(data= data[data$type == "HELIX",],
                                mapping=ggplot2::aes(xmin=begin,
                                                     xmax=end,
                                                     ymin=order-0.32,
                                                     ymax=order-0.27),
                                fill='#5D3A9B',
                                alpha = alpha,
                                show.legend = show.legend)
  }
  # finally TURN
  if(show_turn == TRUE){
    
    p <- p + ggplot2::geom_rect(data= data[data$type == "TURN",],
                                mapping=ggplot2::aes(xmin=begin,
                                                     xmax=end,
                                                     ymin=order-0.32,
                                                     ymax=order-0.27),
                                fill='#1E88E5',
                                alpha = alpha,
                                show.legend = show.legend)
  }
  return(p)
}

### draw_variants
#' New function to label variants on the plot using geom_point
#' Keyword "VARIANT" for selecting variants in the data and keyword "ClinVar_ClinSign" for the class of the variant
#' Colour pathogenic, benign, and unknown variants differently
#'
draw_variants <- function(p, data = data,
                         size = 2,
                         alpha = 1.0,
                         show.legend = FALSE){
  begin=end=description=NULL
  p <- p + ggplot2::geom_point(data = data[type == "VARIANT" & ClinVar_ClinSign == 'Pathogenic',],
                               ggplot2::aes(x = begin,
                                            y = order-0.25),
                               shape = 21,
                               colour = "black",
                               fill = 'firebrick',
                               size = size,
                               alpha = alpha,
                               show.legend = show.legend)
  p <- p + ggplot2::geom_point(data = data[type == "VARIANT" & ClinVar_ClinSign == 'Benign',],
                               ggplot2::aes(x = begin,
                                            y = order+0.25),
                               shape = 21,
                               colour = "black",
                               fill = 'forestgreen',
                               size = size,
                               alpha = alpha,
                               show.legend = show.legend)
  p <- p + ggplot2::geom_point(data = data[type == "VARIANT" & ClinVar_ClinSign == 'Unknown',],
                               ggplot2::aes(x = begin,
                                            y = order+0.25),
                               shape = 21,
                               colour = "black",
                               fill = '#FFD858',
                               size = size,
                               alpha = alpha,
                               show.legend = show.legend)
  return(p)
}

### draw_disulphide
#' New function to label disulphide bonds on the plot using geom_point
#' Keyword "DISULFID" for selecting disulphide bonds in the data
#'
draw_disulphide <- function(p, data = data,
                          size = 3,
                          fill = '#E4D00A',
                          alpha = 1.0,
                          show.legend = FALSE){
  begin=end=description=NULL
  p <- p + ggplot2::geom_point(data = data[data$type == "DISULFID",],
                               ggplot2::aes(x = begin,
                                            y = order+0.25),
                               shape = 21,
                               colour = "black",
                               fill = fill,
                               size = size,
                               alpha = alpha,
                               show.legend = show.legend)
  return(p)
}

### draw_disorder
#' New function to label disordered regions on the plot using geom_rrect from statebins package
#' Keyword "DISORDER" for selecting disordered regions in the data
#' colorPicker column in the data.frame to color disordered regions independently
#'
draw_disorder <- function(p,
                          data = data,
                          label_size = 6,
                          alpha = 1.0,
                          label_disorder = TRUE,
                          show.legend = TRUE){
  begin=end=description=NULL
  tmpData <- data[data$type == 'DISORDER',]
  colorPicker <- tmpData$colorPicker
  
  p <- p + statebins:::geom_rrect(data= tmpData,
                                  mapping=ggplot2::aes(xmin=begin,
                                                       xmax=end,
                                                       ymin=order-0.2,
                                                       ymax=order+0.2
                                  ),
                                  colour = 'black',
                                  fill=colorPicker,
                                  alpha = alpha,
                                  show.legend = show.legend)
  if(label_disorder){
    p <- p + shadowtext::geom_shadowtext(data = tmpData,
                                         ggplot2::aes(x = (begin + (end-begin)/2 - 5),
                                                      y = order,
                                                      label = 'Disorder'),
                                         size = label_size, angle = 90,
                                         bg.colour = "white",  # Color of the outline
                                         bg.r = 0.15,          # Size of the outline
                                         colour = "black"      # Color of the text
    )
  }
  return(p)
}