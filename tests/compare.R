library(vegan)
library(rhdf5)

import_biom <- function(file) {
  bm <- read_biom(file)
  bmotu <- as.data.frame(as.matrix(biom_data(bm)), check.names=FALSE)
  bmtax <- as.data.frame(as.matrix(observation_metadata(bm)), check.names=FALSE)
  return(list(otu=bmotu, tax=bmtax, biom=bm))
}

fix_biom <- function(biom) {
  biom$rows <- lapply(biom$rows, function(x) { mm <- x$metadata; x$metadata <- NULL; x$metadata$taxonomy <- mm; x })
  if (is.null(names(biom$columns[[1]]$metadata))) return(biom)
  biom$columns <- lapply(biom$columns, function(x) { x$metadata <- as.list(x$metadata); return(x) })
  return(biom)
}

write_biom <- function(biom, file, pretty=FALSE) {
  jsonlite::write_json(biom, pretty=pretty, auto_unbox = TRUE, na = "string", path = file)
}

write_biom_hdf5 <- function(biom, file) {
  observation <- list()
  sample <- list()

#  observation$`group-metadata` <- list()
#  sample$`group-metadata` <- list()
  observation$ids <- as.array(unlist(lapply(biom$rows, function(x) x$id)))
  y <- attributes(as(as.matrix(biom_data(biom)), "dgRMatrix"))
  observation$matrix <- list(data=y$x, indices=y$j, indptr=y$p )
  shape <- dim(biom_data(biom))
  observation$metadata$taxonomy <- as.matrix(t(observation_metadata(biom)))
  rownames(observation$metadata$taxonomy) <- NULL
  colnames(observation$metadata$taxonomy) <- NULL
  
  sample$ids <- as.array(unlist(lapply(biom$columns, function(x) x$id)))
  y <- attributes(as(as.matrix(biom_data(biom)), "dgCMatrix"))
  sample$matrix <- list(data=y$x, indices=y$i, indptr=y$p )
  sample$metadata <- as.list(sample_metadata(biom))

  if (file.exists(file)) file.remove(file)
  loc <- H5Fcreate(file)
  h5save(observation, sample, file=file, createnewfile = FALSE)
  
  atts <- list(biom$date, biom$format_url, c(2,1), biom$generated_by, "No Table ID", length(sample$matrix$data), shape, x3$type)
  att_names <- c("creation-date", "format-url", "format-version", "generated-by", "id", "nnz", "shape", "type")
  out <- list(observation=observation, sample=sample)
  names(atts) <- att_names
  attributes(out) <- atts
  
 
  lapply(1:8, function (x) h5writeAttribute(atts[[x]], loc, att_names[x]))
  H5close()
  invisible(out)


}

biomv2 <- "testdata/rich_dense_otu_table.v2.biom"

x <- h5read("testdata/rich_dense_otu_table.v2.biom", "/", read.attributes = TRUE)
x3 <- read_biom("testdata/rich_dense_otu_table.v2.biom")

