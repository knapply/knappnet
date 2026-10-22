\name{read_network}
\alias{read_network}
\alias{write_network}
\title{Read and write network files}
\usage{
read_network(file, format = NULL, directed = NULL, nodes = NULL)

write_network(x, file, format = NULL, weight = NULL, nodes_file = NULL)
}
\arguments{
\item{file}{Path to the file.}
\item{format}{Format name (\code{"graphml"}, \code{"gexf"}, \code{"gml"},
\code{"pajek"}, \code{"dl"}, \code{"json"}, \code{"csv"}, \code{"tsv"}, or
\code{"dot"}), or \code{NULL} to use the file extension.}
\item{directed}{For formats that do not record direction (UCINET DL and edge
lists), whether edges are directed. Defaults to \code{TRUE}.}
\item{nodes}{For edge lists only: a data frame of nodes, or the path to a CSV
file of nodes whose first column is the node identifier.}
\item{x}{A \code{\link{netframe}}.}
\item{weight}{Edge column to write as the link value in formats with a single
positional value (Pajek, UCINET DL) and as the GEXF edge weight. Defaults to
the network's weight column.}
\item{nodes_file}{For edge lists only: optional path for a CSV file of nodes,
so isolates and node attributes are kept.}
}
\value{
\code{read_network()} returns a \code{\link{netframe}}.
\code{write_network()} returns \code{file}, invisibly.
}
\description{
\code{read_network()} reads a file into a \code{\link{netframe}};
\code{write_network()} writes one out. The format is taken from the file
extension unless \code{format} is given.
}
\details{
\tabular{lllll}{
\strong{Format} \tab \strong{Extensions} \tab \strong{Read} \tab \strong{Write} \tab \strong{Common source} \cr
GraphML \tab .graphml, .xml \tab yes \tab yes \tab igraph, NetworkX, Gephi, yEd, Cytoscape \cr
GEXF \tab .gexf \tab yes \tab yes \tab Gephi \cr
GML \tab .gml \tab yes \tab yes \tab igraph, NetworkX, Cytoscape \cr
Pajek \tab .net, .paj \tab yes \tab yes \tab Pajek, igraph, NetworkX \cr
UCINET DL \tab .dl \tab yes \tab yes \tab UCINET, igraph \cr
JSON node-link \tab .json \tab yes \tab yes \tab NetworkX, D3 \cr
Edge list \tab .csv, .tsv \tab yes \tab yes \tab spreadsheets, any tool \cr
DOT \tab .dot, .gv \tab no \tab yes \tab Graphviz \cr
}
Attributes keep their types where the format records them (GraphML, GEXF,
JSON) and are inferred otherwise. Node identifiers that are all whole numbers
are returned as integers. Formats differ in what they can hold: Pajek and
UCINET DL store only structure, node labels, and a single edge value, and DOT
is for drawing.
}
\examples{
g <- netframe(data.frame(from = c("a", "b"), to = c("b", "c"), w = c(1, 2)),
              weight = "w", directed = FALSE)
path <- tempfile(fileext = ".graphml")
write_network(g, path)
read_network(path)
}
