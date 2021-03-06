################################################################################
################################################################################
#' Read a biom-format file, returning a \code{biom-class}.
#'
#' Import the data from a biom-format file into R, represented as an instance
#' of the \code{\link{biom-class}}; essentially a \code{\link{list}} with 
#' special constraints that map to \href{http://biom-format.org/documentation/biom_format.html}{the biom-format definition}.
#' 
#' The BIOM file format (canonically pronounced biome) is designed to be a general-use format for representing biological sample by observation contingency tables. BIOM is a recognized standard for the \href{http://www.earthmicrobiome.org/}{Earth Microbiome Project} and is a \href{http://gensc.org/}{Genomics Standards Consortium} candidate project. Please see \href{http://biom-format.org/}{the biom-format home page} for more details.
#' 
#' It is tempting to include an argument identifying the 
#' biom-format version number of the data file being imported.
#' However, the biom-format version number is a required
#' field in the biom-format definition. 
#' Rather than duplicate this formal specification
#' and allow the possibility of a conflict, the version 
#' number of the biom format will be referred to only by
#' the "format" field in the biom formatted data,
#' or its representation in R.
#'
#' @usage read_biom(biom_file)
#'
#' @param biom_file (Required). A character string indicating the 
#'  file location of the biom formatted file. This is a HDF5 or JSON formatted file
#'  specific to biological datasets. 
#'  The format is formally defined at \href{http://biom-format.org/documentation/biom_format.html}{the biom-format definition}
#'  and depends on the versioning.
#'
#' @return An instance of the \code{biom-class}.
#'
#' @seealso 
#' 
#' Function to create a biom object from R data,
#' \code{\link{make_biom}}.
#' 
#' Definition of the
#' \code{\link{biom-class}}. 
#' 
#' Function to write a biom format file from a biom object,
#' \code{\link{write_biom}}
#'
#' Accessor functions like \code{\link{header}}.
#'
#' @references \url{http://biom-format.org/}
#'
#' @importFrom jsonlite fromJSON
#' @export
#' @examples
#' # # # import with default parameters, specify a file
#' biom_file <- system.file("extdata", "rich_sparse_otu_table.biom", package = "biomformat")
#' biom_file
#' read_biom(biom_file)
#' biom_file <- system.file("extdata", "min_sparse_otu_table.biom", package = "biomformat")
#' biom_file
#' read_biom(biom_file)
#' ## The previous examples use system.file() because of constraints in specifying a fixed
#' ##   path within a reproducible example in a package. 
#' ## In practice, however, you can simply provide "hard-link"
#' ## character string path to your file:
#' # mybiomfile <- "path/to/my/biomfile.biom"
#' # read_biom(mybiomfile)
read_biom <- function(biom_file){	
	# Read the file, storing as a list 
	# generated by jsonlite w/ default JSON parsing params
	trash = try(silent=TRUE,
	            expr = {
	              x <- fromJSON(biom_file, simplifyDataFrame = FALSE, simplifyMatrix = FALSE)
	            })
	if(inherits(trash, "try-error")){
	  # If JSON interpretation attempt failed, try HDF5
	  trash = try(silent=TRUE,
	              expr = {x <- read_hdf5_biom(biom_file)})
	}
	if(inherits(trash, "try-error")){
	  # If still bad, throw helpful error.
	  stop("Both attempts to read input file:\n", 
	       biom_file, "\n", 
	       "either as JSON (BIOM-v1) or HDF5 (BIOM-v2).\n",
	       "Check file path, file name, file itself, then try again.")
	}
	# Use the biom() constructor function to 
	# instantiate a biom-class, perform validity checks. Return.
	return( biom(x) )
}
################################################################################
#' Write a biom-format v1 file, returning a \code{biom-class}.
#'
#' @param biom (Required). A biom object that is going to be written to file
#'  as a proper biom formatted file, adhering to 
#'  \href{http://biom-format.org/documentation/biom_format.html}{the biom-format definition}.
#'  
#' @param biom_file (Required). A character string indicating the 
#'  file location of the biom formatted file. This is a JSON formatted file
#'  specific to biological datasets. 
#'  The format is formally defined at 
#'  \href{http://biom-format.org/documentation/biom_format.html}{the biom-format definition}
#'  
#' @param pretty logical; Should biom output be pretty printed?
#'
#' @return Nothing. The first argument, \code{x}, is written to a file.
#'
#' @seealso 
#' 
#' Function to create a biom object from R data,
#' \code{\link{make_biom}}.
#' 
#' Definition of the
#' \code{\link{biom-class}}. 
#' 
#' The \code{\link{read_biom}} import function.
#'
#' Accessor functions like \code{\link{header}}.
#'
#' @references \url{http://biom-format.org/}
#'
#' @export
#' @importFrom jsonlite toJSON write_json
#' @examples
#' biom_file <- system.file("extdata", "rich_sparse_otu_table.biom", package = "biomformat")
#' x = read_biom(biom_file)
#' outfile = tempfile()
#' write_biom(x, outfile)
#' y = read_biom(outfile)
#' identical(x, y)
write_biom <- function(biom, biom_file, pretty=FALSE) {
  write_json(biom, path=biom_file, pretty=pretty, auto_unbox = TRUE, na = "null")
}

################################################################################
#' Read in a biom-format vs 2 file, returning a \code{list}.
#'
#' This function is meant only to be used if the user knows the file is
#' a particular version / hdf5 format. Otherwise, the `read_biom` file should be used.
#'
#' @param biom_file (Required). A biom object that is going to be written to file
#'  as a proper biom formatted file, adhering to 
#'  \href{http://biom-format.org/documentation/biom_format.html}{the biom-format definition}.
#'  
#' @return Nothing. The first argument, \code{x}, is written to a file.
#'
#' @seealso 
#' 
#' Function to create a biom object from R data,
#' \code{\link{make_biom}}.
#' 
#' Definition of the
#' \code{\link{biom-class}}. 
#' 
#' The \code{\link{read_hdf5_biom}} import function.
#'
#' Accessor functions like \code{\link{header}}.
#'
#' @references \url{http://biom-format.org/}
#'
#' @export
#' @importFrom rhdf5 h5read
#' @examples
#' biom_file <- system.file("extdata", "rich_sparse_otu_table_hdf5.biom", package = "biomformat")
#' x = read_hdf5_biom(biom_file)
#' x = biom(x)
#' outfile = tempfile()
#' write_biom(x, outfile)
#' y = read_biom(outfile)
#' identical(observation_metadata(x),observation_metadata(y))
#' identical(sample_metadata(x),sample_metadata(y))
#' identical(biom_data(x), biom_data(y))
read_hdf5_biom<-function(biom_file){
	x = h5read(biom_file,"/",read.attributes = TRUE)
	data = generate_matrix(x)
	rows = generate_metadata(x$observation)
	columns = generate_metadata(x$sample)
	shape = c(length(data),length(data[[1]]))
 
	id = attr(x,"id")
	vs = attr(x,"format-version")
	format = sprintf("Biological Observation Matrix %s.%s",vs[1],vs[2])
	format_url = attr(x,"format-url")
	type = "OTU table"
	generated_by = attr(x,"generated-by")
	date = attr(x,"creation-date")
	matrix_type = "dense"
	matrix_element_type = "int"
	
	namedList(id,format,format_url,type,generated_by,date,matrix_type,matrix_element_type,
		rows,columns,shape,data)
}
################################################################################

#' write string variables to scalar hdf5 space
#'
#' @param attr Character string.  Attribute to be saved.
#' @param h5obj An object of class H5IdComponent representing a H5 object identifier (file, group, or dataset). See H5Fcreate, H5Fopen, H5Gcreate, H5Gopen, H5Dcreate, or H5Dopen to create an object of this kind.
#' @param name Character string. Name of attribute.
#'
#' @return
#' @export
#' @importFrom rhdf5 H5Screate H5Tcopy H5Tset_size H5Acreate H5Awrite H5Aclose H5Sclose
#'
#' @examples
h5writeAttribute.scalar <- function(attr, h5obj, name) {
  h5space4 = H5Screate("H5S_SCALAR")
  typeid <- H5Tcopy("H5T_C_S1")
  H5Tset_size(typeid, size=nchar(attr)+1)
  ga <- H5Acreate(h5obj, name, typeid, h5space4)
  ccc <- H5Awrite(ga, attr)
  H5Aclose(ga)
  H5Sclose(h5space4)
  
}



#' Write biom file in HDF5 format
#'
#' @param biom A biom object.
#' @param file A character string indicating the 
#'  file location of the biom formatted file. This is a HDF5 formatted file
#'  specific to biological datasets. 
#'  The format is formally defined at 
#'  \href{http://biom-format.org/documentation/biom_format.html}{the biom-format definition}
#'
#' @return
#' @export
#' @import rhdf5
#' @examples
write_hdf5_biom <- function(biom, file) {
  observation <- list()
  sample <- list()
  
  observation$`group-metadata` <- list()
  sample$`group-metadata` <- list()
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
  print(loc)
  h5save(observation, sample, file=file, createnewfile = FALSE)
  
  atts <- list(biom$date, biom$format_url, c(2,1), biom$generated_by, "No Table ID", length(sample$matrix$data), shape, biom$type)
  names(atts) <- c("creation-date", "format-url", "format-version", "generated-by", "id", "nnz", "shape", "type")
  
  lapply(c(1,2,4,5,8),function (x) h5writeAttribute.scalar(atts[[x]], loc, names(atts)[x]) )
  lapply(c(3,7),function (x) h5writeAttribute(atts[[x]], loc, names(atts)[x]) )
  h5writeAttribute.integer(atts[[6]], loc, names(atts)[6])
  H5close()
  
  out <- list(observation=observation, sample=sample)
  attributes(out) <- atts
  invisible(out)
  
}
