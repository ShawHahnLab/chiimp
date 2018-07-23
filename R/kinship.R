# Cross-reference genotyping results with kinship information for samples of
# known identity.

#' Match genotyping results with kinship data
#'
#' Given a genotyping results summary data frame, known genotypes, and
#' per-individual kinship data, create a data frame matching called genotypes
#' with known kinship for every row in the input genotyping data frame.
#'
#' @details
#'
#' If the \code{genotypes.known} argument is not given, \code{results_summary}
#' will be used in its place, meaning that the consistency between individuals'
#' genotypes within a single data set will be checked against the kinship data
#' from \code{individual_attrs}.
#'
#' The first two output columns are the Mother and Father of each individual, or
#' NA if not provided.  Each following output column is a factor providing a
#' fixed vocabulary of possible categories, or NA if information is missing for
#' the called genotype and/or known kinship. The two latter pairs of columns are
#' largely redundant but are all provided for convenience.
#'
#' All columns:
#'
#'  * Mother, Father: identifying names (corresponding to the Name columns in
#'    \code{genotypes.known} and \code{results_summary}) for the mother and
#'    father of each individual.
#'  * Allele1Source, Allele2Source: the possible source of each allele:
#'    "Mother," "Father," "Either," or "Neither."
#'  * MotherAllele, FatherAllele: the allele(s) matching each parent:
#'    "Allele1,"  "Allele2," "Either," or "Neither."
#' @md
#'
#' @param results_summary cross-sample summary data frame as produced by
#'   \code{\link{analyze_dataset}}.
#' @param individual_attrs data frame of per-individual attributes including
#'   Mother and Father.
#' @param genotypes.known data frame of known genotypes that should be compared
#'   to the observed genotypes in the results, as loaded by
#'   \code{\link{load_genotypes}}.  If not specified, the same data frame as
#'   \code{results_summary} will be used, effectively just cross-referencing the
#'   determined genotypes with the known kinship data.
#'
#' @return data frame with columns for parentage, the parent(s) matching each
#'   allele, and the allele(s) matching each parent.  See Details section for
#'   more information.
#'
#' @export
match_kin <- function(results_summary,
                      individual_attrs,
                      genotypes.known=results_summary) {
  # The matching Mother and Father for each row in results_summary.
  .idx <- match(results_summary$Name, individual_attrs$Name)
  kinship <- individual_attrs[.idx, c("Mother", "Father")]
  # Data frame of individual and parent allele sequences for each row in
  # results_summary.
  id_m <- paste(kinship[["Mother"]], results_summary$Locus)
  id_f <- paste(kinship[["Father"]], results_summary$Locus)
  id_gtypes <- paste(genotypes.known$Name, genotypes.known$Locus)
  am <- genotypes.known[match(id_m, id_gtypes),
                        c("Allele1Seq", "Allele2Seq")]
  af <- genotypes.known[match(id_f, id_gtypes),
                        c("Allele1Seq", "Allele2Seq")]
  ac <- results_summary[, c("Allele1Seq", "Allele2Seq")]
  a <- cbind(ac, am, af)
  # Create the new data frames for alleles-> kin and kin -> alleles.
  alleles_to_kin <- match_kin_alleles(a)
  a2 <- cbind(alleles_to_kin, match(id_m, id_gtypes), match(id_f, id_gtypes))
  kin_to_alleles <- match_alleles_kin(a2)
  result <- cbind(kinship,
                  alleles_to_kin,
                  kin_to_alleles)
  return(result)
}

# TODO add helper function to convert kinship data into structure ready to be
# passed into kinship2::pedigree().  This should be a data frame with columns
# id, dadid, momid (these apparently can be character vectors), sex (seems OK
# with character M and F for this), and any others we want to include via the
# "affected" matrix of 0/1/NA values.  But how to handle attributes with many
# possible values, like alleles?  These can be extra arguments at plot time,
# like to "col"; maybe these can be shuttled along with the pedigree object as
# attributes and then passed into plot().
#
# Test this whole idea with a set of Gombe samples for one locus, calling
# match_kin against the individual attributes (need to include the sex info),
# then through the helper function, and then along to the plot function.  We
# should be able to re-create Hannah's family tree from the paper.

# Messy First Pass --------------------------------------------------------


# Create a longer-format data frame of relationships for pairs of individuals,
# more suitable as input to graph functions.  Each row is a single
# relationship from one individual to another.  There is redundancy in the
# relationships such as A being Mother to B and B being Child to A; use the
# Relationship column to filter as needed. Input should be a combined
# results_summary/known genotype and kinship data frame.
#
# These rows can represent edges in a graph structure between individuals as
# nodes.  Node (individual) attributes can then be referenced using the "Name"
# and "Other" columns from this data frame.
tabulate_kinship <- function(kinship_data) {
  # This could totally be done with a reshape call.
  d1 <- kinship_data[, c("Name", "Mother", "MotherAllele", "Locus")]
  d1$Relationship <- "Mother"
  colnames(d1) <- c("Name", "Other", "OtherAllele", "Locus", "Relationship")
  d2 <- kinship_data[, c("Name", "Father", "FatherAllele", "Locus")]
  d2$Relationship <- "Father"
  colnames(d2) <- c("Name", "Other", "OtherAllele", "Locus", "Relationship")

  # Add entries from the mate relationships.  This will be a connection for each
  # pair of parents, in both directions, and without duplication across loci.
  d3 <- data.frame(
    Name = c(kinship_data$Mother, kinship_data$Father),
    Other = c(kinship_data$Father, kinship_data$Mother),
    OtherAllele = NA,
    Locus = NA,
    Relationship = "Mate",
    stringsAsFactors = FALSE)
  d3 <- unique(d3)

  # Add entries for the child relationships.  This will just be the inverse of
  # both parent relationships, and without duplication across loci.
  d4 <- data.frame(
    Name = c(d1$Other, d2$Other),
    Other = c(d1$Name, d2$Name),
    OtherAllele = NA,
    Locus = NA,
    Relationship = "Child",
    stringsAsFactors = FALSE)
  d4 <- unique(d4)

  # Bind together and remove any blanks.
  d <- rbind(d1, d2, d3, d4)
  d$Name[d$Name == ""] <- NA
  d$Other[d$Other == ""] <- NA
  d <- subset(d, !is.na(Name) & !is.na(Other))

  d
}

# Subset the structure created by tabulate_kinship, looking at all relationships
# at most "dist" connections away from the given node.
kinship_subset <- function(d, name, dist=1) {
  nms <- name
  for (i in 1:dist) {
    d_sub <- subset(d, Name %in% nms | Other %in% nms)
    nms <- unique(c(d_sub$Name, d_sub$Other))
  }
  d_sub
}

# Remove locus-specific info and dereplicate.
kinship_delocus <- function(d) {
  unique(d[, -match("Locus", colnames(d))])
}

# quick example of a styled graph display with igraph.
# http://igraph.org/r/doc/plot.common.html
# With this same idea we could style edges based on genotyping
# matches/mismatches, number of loci, etc.
example_graph <- function(d) {
  g <- igraph::graph_from_data_frame(d)
  ### Match edge attributes
  # Mapping of relationship to display settings.
  edge_attr_map <- data.frame(
    arrow.mode = c(">", ">", "-", "-"),
    color      = c("red", "blue", "green", "black"),
    width      = c(1, 1, 0, 0),
    arrow.width = c(1, 1, 0, 0),
    stringsAsFactors = FALSE)
  rownames(edge_attr_map) <- c("Mother", "Father", "Child", "Mate")
  igraph::E(g)$arrow.mode <- edge_attr_map[d$Relationship, "arrow.mode"]
  igraph::E(g)$color <- edge_attr_map[d$Relationship, "color"]
  igraph::E(g)$width <- edge_attr_map[d$Relationship, "width"]
  igraph::E(g)$arrow.width <- edge_attr_map[d$Relationship, "arrow.width"]
  igraph::E(g)$arrow.size <- 0.3
  #igraph::E(g)$curved <- TRUE
  # Match nodes to attributes
  node_attrs <- chimp_attrs[match(names(V(g)), chimp_attrs$Name), ]
  node_shape = c(Male="square", Female="circle", Unknown="none")
  igraph::V(g)$shape <- node_shape[ifelse(is.na(node_attrs$Sex),
                                         "Unknown",
                                         node_attrs$Sex)]
  # Plot
  igraph::plot.igraph(g)
  return(g)
}

# graph across loci, color-coding by allele match/mismatch.
example_graph_loci <- function(d) {
  g <- igraph::graph_from_data_frame(d)
  ### Match edge attributes
  # Mapping of relationship to display settings.
  edge_attr_map <- data.frame(
    arrow.mode = c(">", ">", "-", "-"),
    color      = c("red", "blue", "green", "black"),
    width      = c(1, 1, 0, 0),
    arrow.width = c(1, 1, 0, 0),
    stringsAsFactors = FALSE)
  rownames(edge_attr_map) <- c("Mother", "Father", "Child", "Mate")
  edge_attr_loci <- data.frame(
    color = c(rep("green", 3), "red", "gray"),
    stringsAsFactors = FALSE)
  rownames(edge_attr_loci) <- c("Allele1", "Allele2", "Either", "Neither",
                                 "Unknown")
  igraph::E(g)$arrow.mode  <- edge_attr_map[d$Relationship, "arrow.mode"]
  igraph::E(g)$color       <- ifelse(is.na(d$OtherAllele),
                                     edge_attr_loci["Unknown", "color"],
                                     edge_attr_loci[as.character(d$OtherAllele),
                                                    "color"])
  igraph::E(g)$width       <- edge_attr_map[d$Relationship, "width"]
  igraph::E(g)$arrow.width <- edge_attr_map[d$Relationship, "arrow.width"]
  igraph::E(g)$arrow.size <- 0.3
  # Match nodes to attributes
  node_attrs <- chimp_attrs[match(names(V(g)), chimp_attrs$Name), ]
  node_shape = c(Male="square", Female="circle", Unknown="none")
  igraph::V(g)$shape <- node_shape[ifelse(is.na(node_attrs$Sex),
                                          "Unknown",
                                          node_attrs$Sex)]
  # Plot
  igraph::plot.igraph(g)
  return(g)
}

# Detailed per-locus graph, showing alleles, color-coding by allele
# match/mismatch.
example_graph_locus <- function(d, locus) {
  d <- subset(d, Locus == locus | is.na(Locus))
  g <- igraph::graph_from_data_frame(d)
  ### Match edge attributes
  # Mapping of relationship to display settings.
  edge_attr_map <- data.frame(
    arrow.mode = c(">", ">", "-", "-"),
    color      = c("red", "blue", "green", "black"),
    width      = c(1, 1, 0, 0),
    arrow.width = c(1, 1, 0, 0),
    stringsAsFactors = FALSE)
  rownames(edge_attr_map) <- c("Mother", "Father", "Child", "Mate")
  edge_attr_loci <- data.frame(
    color = c(rep("green", 3), "red", "gray"),
    stringsAsFactors = FALSE)
  rownames(edge_attr_loci) <- c("Allele1", "Allele2", "Either", "Neither",
                                "Unknown")
  igraph::E(g)$arrow.mode  <- edge_attr_map[d$Relationship, "arrow.mode"]
  igraph::E(g)$color       <- ifelse(is.na(d$OtherAllele),
                                     edge_attr_loci["Unknown", "color"],
                                     edge_attr_loci[as.character(d$OtherAllele),
                                                    "color"])
  igraph::E(g)$width       <- edge_attr_map[d$Relationship, "width"]
  igraph::E(g)$arrow.width <- edge_attr_map[d$Relationship, "arrow.width"]
  igraph::E(g)$arrow.size <- 0.3
  # TODO pass in kg as argument?
  id  <- paste(names(V(g)), rep(locus, length(V(g))))
  id2 <- paste(kg$Name,     kg$Locus)
  a1 <- make_allele_name(kg[match(id, id2), "Allele1Seq"])
  a2 <- make_allele_name(kg[match(id, id2), "Allele2Seq"])
  igraph::V(g)$label <- paste(names(V(g)), a1, a2, sep = "\n")
  igraph::V(g)$shape <- "rectangle"
  igraph::V(g)$size <- 20
  igraph::V(g)$size2 <- 15
  igraph::V(g)$label.cex = 0.7
  # Plot
  igraph::plot.igraph(g)
  return(g)
}

kinship_worked_example <- function(kg, individual_attrs) {
  kg <- cbind(kg, match_kin(kg, individual_attrs))
  d <- tabulate_kinship(kg)
  d_sub <- kinship_subset(d, "Aris", dist = 3)
  example_graph_locus(d_sub, "NMS2")
}

# There are other graph structures that may be useful besides kinship.  For
# example, a representation of the all-to-all distance matrix.

# does this work like I think it does?  Not sure.
distmat_subset <- function(dm, name, cutoff=8, dist=1) {
  nms <- name
  for (i in 1:dist) {
    nms2 <- nms
    for (nm in nms) {
      x <- dm[nm, ]
      nms2 <- c(nms2, names(x[x<cutoff]))
    }
    nms <- unique(nms2)
  }
  dm[nms, nms]
}

distmat_worked_example <- function(kg) {
  d <- make_dist_mat(subset(kg,
                            Locus %in% c("A", "B", "C", "D", 1, 2, 3, 4) &
                              Location == "Gombe"))

  #d_sub <- distmat_subset(d, "Ch-095", 8, 3)

  nms <- colnames(d)[(d["Ch-095", ] < 12)]
  d_sub <- d[nms, nms]

  # Translate the distance value integer into a closeness value decimal, from an
  # arbitrary 1 to 1000.  (This seems to work nicely with the layout functions.)
  adj_sub <- 10^((16-d_sub)/16*3)

  g <- igraph::graph_from_adjacency_matrix(adj_sub,
                                           mode = "undirected",
                                           weighted = TRUE,
                                           diag = FALSE)
  # Use a range of 1 to 10 for the displayed edge thickness, but still on that
  # log scale.
  width_min <- 1
  width_max <- 20
  width_from_weights <- function(w) {
    (w - min(w))/(max(w) - min(w))*(width_max - width_min) + width_min
  }
  E(g)$width <- width_from_weights(E(g)$weight)

  # Use a range of (0.5 to 1)*255 for the displayed edge alpha transparency, but
  # still on that log scale.
  alpha_min <- 15
  alpha_max <- 255
  color_from_weights <- function(w) {
    alpha <- (w - min(w))/(max(w) - min(w))*(alpha_max - alpha_min) + alpha_min
    rgb(0, 0, 0, alpha, maxColorValue = 255)
  }

  # Fade the edges
  E(g)$color <- color_from_weights(E(g)$weight)

  g <- add_layout_(g, with_fr())

  plot(g)

  g
}

# !!!! We can define clustering with this!  Perfect for seeing where a given
# individual gets assigned.  Does it make more sense to define graphs with
# number of edges as the number of matches instead of using weights?  Not sure.
# https://stackoverflow.com/questions/9471906/what-are-the-differences-between-community-detection-algorithms-in-igraph/
community_detection <- function(kg) {
  LOCI <-  c("A", "B", "C", "D",
             1, 2, 3, 4,
             "FP1", "FP2", "FP3", "FP4",
             "NMS2", "NMS5", "NMS11", "AMEL")
  d <- make_dist_mat(subset(kg,
                            Locus %in% LOCI))
  d_sub <- d

  # ##############################
  # All this as above...
  weighted <- TRUE
  if (weighted) {
    adj_sub <- 10^((length(LOCI)*2-d_sub)/(length(LOCI)*2)*3)
    g <- igraph::graph_from_adjacency_matrix(adj_sub,
                                             mode = "undirected",
                                             weighted = TRUE,
                                             diag = FALSE)
  } else {
    adj_sub <- max(d_sub) - d_sub
    g <- igraph::graph_from_adjacency_matrix(adj_sub,
                                             mode = "undirected",
                                             weighted = NULL,
                                             diag = FALSE)
  }

  width_min <- 1
  width_max <- 20
  width_from_weights <- function(w) {
    (w - min(w))/(max(w) - min(w))*(width_max - width_min) + width_min
  }
  E(g)$width <- width_from_weights(E(g)$weight)
  alpha_min <- 15
  alpha_max <- 255
  color_from_weights <- function(w) {
    alpha <- (w - min(w))/(max(w) - min(w))*(alpha_max - alpha_min) + alpha_min
    rgb(0, 0, 0, alpha, maxColorValue = 255)
  }
  E(g)$color <- color_from_weights(E(g)$weight)
  g <- add_layout_(g, with_fr())


  # Since it's huge, shrink things down a bit.
  V(g)$size <- 0.1
  E(g)$width <- E(g)$width/5
  V(g)$label.cex <- 0.7


  x <- cluster_fast_greedy(g)

  cols <- c("black", "red", "green", "blue", "yellow", "purple")
  V(g)$label.color <- cols[membership(x)]

  plot(g)
  g
}


# kinship2 ----------------------------------------------------------------

# Helper functions for working with the pedigree objects from the kinship2
# package.

pedigree_setup <- function(d) {
  idx <- is.na(d$Mother) | is.na(d$Father)
  d$Mother[idx] <- NA
  d$Father[idx] <- NA
  kinship2::pedigree(id = d$Name,
                     dadid = d$Father,
                     momid = d$Mother,
                     sex = d$Sex)
}

# subset a kinship2 pedigree object from a named individual out to a given
# distance.
pedigree_subset <- function(p, name, dist=1) {
  nms <- name
  d <- as.data.frame(p)
  d$dadid[d$dadid == 0] <- NA
  d$momid[d$momid == 0] <- NA
  for (i in 1:dist) {
    d_sub <- subset(d, id %in% nms | dadid %in% nms | momid %in% nms)
    nms <- list(d_sub$id, d_sub$momid, d_sub$dadid)
    nms <- unique(unlist(lapply(nms, function(nm) as.character(nm[!is.na(nm)]))))
  }
  d_sub <- subset(d, id %in% nms |
                    id %in% momid & ! is.na(momid) |
                    id %in% dadid & ! is.na(dadid))
  idx_parents_outside <-
    (! d_sub$momid %in% d_sub$id) |
    (! d_sub$dadid %in% d_sub$id)
  d_sub$momid[idx_parents_outside] <- NA
  d_sub$dadid[idx_parents_outside] <- NA
  do.call(pedigree, d_sub)
}

# proof of concept showing an approximation of the family tree shown in the
# paper for Locus 3.  With some more work we could make filled polygons and
# match it even better.
figure_from_paper <- function(p, kg, chimp_attrs) {
  cols <- c(
    "234-a" = "yellow",
    "234-b" = "blue",
    "234-c" = "red",
    "234-d" = "green",
    "226-a" = "gray"
    )
  p2 <- pedigree_subset(p, "Glama", 2)

  kg_kin <- cbind(kg, match_kin(kg, chimp_attrs))
  idx <- match(paste(p2$id, "3"), paste(kg_kin$Name, kg_kin$Locus))
  x <- normalize_alleles(kg_kin[idx, c("Allele1Name", "Allele2Name")])

  items <- plot(p2)
  points(items$x-items$boxw/4, items$y+items$boxh/2, pch=19, cex=2, col=cols[x$V1])
  points(items$x+items$boxw/4, items$y+items$boxh/2, pch=19, cex=2, col=cols[x$V2])
}

# plot an example pedigree tree for a given locus and mark those individuals who
# have mismatches with a parent.
# TODO finish this; draw filled half squares/circles with labels like Hannah
# did.
figure_with_mismatches <- function(locus="NMS2") {
  kg <- load_genotypes("~/analysis/chimp-microsat-combined/from-google/known_genotypes.csv")
  chimp_attrs <- load_genotypes("~/analysis/chimp-microsat-combined/from-google/individual_attrs.csv")
  kg_kin <- cbind(kg, match_kin(kg, chimp_attrs))

  chimp_attrs$Mother[chimp_attrs$Name %in% c("Baroza", "Bahati")] <- NA
  chimp_attrs$Father[chimp_attrs$Name %in% c("Baroza", "Bahati")] <- NA
  chimp_attrs$Sex[chimp_attrs$Name == "Baroza"] <- "F"
  p <- pedigree_setup(chimp_attrs)

  p2 <- pedigree_subset(p, "Nyota", 3)


  idx <- match(paste(p2$id, locus), paste(kg_kin$Name, kg_kin$Locus))

  idx2 <- kg_kin$MotherAllele[idx] %in% "Neither" | kg_kin$FatherAllele[idx] %in% "Neither"
  mismatched <- rep(0, length(idx))
  mismatched[idx2] <- 1

  x <- normalize_alleles(kg_kin[idx, c("Allele1Name", "Allele2Name")])
  lvls <- unique(c(x$V1, x$V2))
  x$V1 <- factor(x$V1, levels = lvls)
  x$V2 <- factor(x$V2, levels = lvls)
  cols <- rainbow(length(levels(x$V1)))
  cols[lvls == ""] <- "#ffffff"

  items <- plot(p2, col=mismatched+1, lwd=2)
  # Draw color-coded dots for all known alleles.  Note that we need to use
  # xpd=TRUE to prevent clipping of things drawn at the edges.
  #
  # Still not sure how best to manage the duplicated shapes. items$plist$nid
  # shows them but not items$x, y.
  points(items$x-items$boxw/4, items$y+items$boxh/2, pch=19, cex=2, col=cols[x$V1], xpd=T)
  points(items$x+items$boxw/4, items$y+items$boxh/2, pch=19, cex=2, col=cols[x$V2], xpd=T)
  points(items$x[mismatched==1], items$y[mismatched==1], pch=1, cex=2, col="red")
  items
}

# Old Stuff ---------------------------------------------------------------



# create a graph object.  If no locus is specified this will aggregate
# connections for all loci available.  If a locus is given, the graph will be
# specific to that locus.
kinship_graph_old <- function(kinship_data, locus=NA, layoutType="dot") {
  # TODO don't do this!
  kinship_data <- subset(kinship_data, ! is.na(Name))
  kinship_data <- kinship_data[match(unique(kinship_data$Name),
                                     kinship_data$Name), ]

  gnodes <- kinship_data$Name
  names(gnodes) <- kinship_data$Name

  gnodes <- sort(unique(c(kinship_data$Name,
                          kinship_data$Mother,
                          kinship_data$Father)))
  gnodes[gnodes == ""] <- NA
  gnodes <- gnodes[! is.na(gnodes)]
  names(gnodes) <- gnodes

  gedges <- lapply(gnodes, function(nm) {
    idx <- match(nm, kinship_data$Name)
    x <- as.character(kinship_data[idx, c("Mother", "Father")])
    names(x) <- c("Mother", "Father")
    x[x == ""] <- NA
    list(edges=x[!is.na(x)])
  })

  # Find every other node that has this node as a parent, and make an edge
  gedges_children <- lapply(gnodes, function(nm) {
    children <- sapply(gedges, function(node) { nm %in% node$edges })
    children <- names(children[children])
    if (length(children) > 0)
      names(children) <- paste0(rep("Child", length(children)),
                                1:length(children))
    list(edges=children)
  })

  gedges_both <- gedges
  for (nm in names(gedges)) {
    gedges_both[[nm]]$edges <- c(gedges[[nm]]$edges, gedges_children[[nm]]$edges)
  }

  # if using gedges_both, "undirected" works:
  g <- graph::graphNEL(nodes=gnodes, edgeL=gedges_both, edgemode = "directed")

  # remove orphaned nodes
  # Directed version:
  idxl <- graph::degree(g)[[1]] == 0 & graph::degree(g)[[2]] == 0
  # Undirected version:
  #idxl <- graph::degree(g) == 0
  g <- graph::removeNode(gnodes[idxl], g)

  # Graph attributes

  node_attrs <- list(color = "black")
  #edge_attrs <- list(arrowhead=c("Apple~Aphro"="box"))
  edge_attrs <- list(color = "black")
  graph_attrs <- list(node = node_attrs, edge=edge_attrs) #list(node=list(fontsize = 4))
  g <- Rgraphviz::layoutGraph(g, attrs=graph_attrs, layoutType=layoutType)
  graph::nodeRenderInfo(g) <- list(fontsize=60, shape="box")

  return(g)
}

# Create a sub-graph of the given name and  all those nodes directly connected.
kinship_subset_old <- function(g, name) {
  g_sub <- subGraph(c(name, unlist(adj(g, name))), g)
  # Need to re-layout the subgraph so it's ready to render.
  Rgraphviz::layoutGraph(g_sub)
}

# Util --------------------------------------------------------------------


# six-col character matrix: two cols of alleles, two cols for mother, two for
# father.  NA for the 2nd/4th/6th implies homozygous for that individual.
match_kin_alleles <- function(alleles) {
  alleles[is.na(alleles[, 2]), 2] <- alleles[is.na(alleles[, 2]), 1]
  alleles[is.na(alleles[, 4]), 4] <- alleles[is.na(alleles[, 4]), 3]
  alleles[is.na(alleles[, 6]), 6] <- alleles[is.na(alleles[, 6]), 5]
  # For each row, define each allele's matching to either parent.
  case <- c("Neither", "Father", "Mother", "Either")
  result <- apply(alleles, 1, function(r) {
    if (is.na(r[1]))
      return(c(NA, NA))
    in_m <- r[1:2] %in% r[3:4]
    in_f <- r[1:2] %in% r[5:6]
    a1 <- c(in_m[1], in_f[1])
    a2 <- c(in_m[2], in_f[2])
    a1 <- case[1 + sum(a1 * 2:1)]
    a2 <- case[1 + sum(a2 * 2:1)]
    # If we're missing any parent data, that could explain the lack of a match.
    # Backpedal for that case to a "no comment."
    a <- c(a1, a2)
    if (any(is.na(r[3:6])))
      a[a == "Neither"] <- NA
    #if (sum(in_m, in_f) > 0)
    #  browser()
    return(a)
  })
  result <- as.data.frame(t(result))
  colnames(result) <- c("Allele1Source", "Allele2Source")
  # TODO define level order?
  result[[1]] <- factor(result[[1]], levels = case)
  result[[2]] <- factor(result[[2]], levels = case)
  return(result)
}

# translate match_kin_alleles output (columns identifying parents for each
# allele) into columns identifying alleles for each parent.
# four-col data frame: two cols for results of match_kin_alleles, one col each
# for mother and father index in known genoyptes (i.e., do we have genotypes for
# them?)
match_alleles_kin <- function(a2) {
  case <- c("Neither", "Allele2", "Allele1", "Either")
  kin_to_alleles <- apply(a2, 1, function(r) {
    if (all(is.na(r[1:2])))
      return(c(NA, NA))
    m_in <- r[1:2] %in% c("Mother", "Either")
    f_in <- r[1:2] %in% c("Father", "Either")
    m <- case[1 + sum(m_in * 2:1)]
    f <- case[1 + sum(f_in * 2:1)]
    if (is.na(r[3]))
      m <- NA
    if (is.na(r[4]))
      f <- NA
    p <- c(m, f)
    return(p)
  })
  kin_to_alleles <- as.data.frame(t(kin_to_alleles))
  kin_to_alleles[[1]] <- factor(kin_to_alleles[[1]], levels = case)
  kin_to_alleles[[2]] <- factor(kin_to_alleles[[2]], levels = case)
  colnames(kin_to_alleles) <- c("MotherAllele", "FatherAllele")
  return(kin_to_alleles)
}
