\name{net_read}
\alias{net_read}
\alias{net_write}
\title{Read and write network files}
\description{
\code{net_read()} reads a file into a \code{\link{net_frame}};
\code{net_write()} writes one out. The format is taken from the file
extension unless \code{format} is given.
}
\usage{
net_read(file, format = NULL, directed = NULL, nodes = NULL, source = "source",
         target = "target", id = "id", weight = NULL, weight_type = NULL,
         quiet = FALSE)

net_write(x, file, format = NULL, weight = NULL, nodes_file = NULL)
}
\arguments{
\item{file}{Path to the file.}
\item{format}{Format name (\code{"graphml"}, \code{"gexf"}, \code{"gml"},
\code{"pajek"}, \code{"dl"}, \code{"json"}, \code{"dynetml"}, \code{"csv"},
\code{"tsv"}, \code{"gpkg"}, or \code{"dot"}), or \code{NULL} to use the
file extension.}
\item{directed}{For formats that do not record direction (UCINET DL and edge
lists), whether edges are directed. Defaults to \code{TRUE}.}
\item{nodes}{For edge lists only: a data frame of nodes, or the path to a CSV
file of nodes.}
\item{source, target, id}{For edge lists (and GeoPackages without a
\code{network} layer): the names of the endpoint columns of the edge file
and of the identifier column of the node file, as in
\code{\link{net_frame}}. Other formats name these columns themselves:
\code{source}, \code{target}, and \code{id}.}
\item{x}{A \code{\link[=net_frame]{net_frame}}.}
\item{weight}{When reading, the edge column holding weights; by default a
numeric column called \code{weight} is used if there is one, and
\code{NA} ignores it. When writing, the edge column written as the link
value in formats with a single positional value (Pajek, UCINET DL, DyNetML)
and as the GEXF edge weight; defaults to the network's weight column.}
\item{weight_type}{When reading, what the weights mean (\code{"distance"} or
\code{"strength"}); see \code{\link{net_frame}}. Required when the file
has a weight column, since no format records the meaning (a GeoPackage
written by \code{net_write()} does, and needs no argument).}
\item{quiet}{Whether to suppress the advisory notes of validation.}
\item{nodes_file}{For edge lists only: optional path for a CSV file of nodes,
so isolates and node attributes are kept.}
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
DyNetML \tab .xml \tab yes \tab yes \tab ORA (CASOS, Carnegie Mellon) \cr
Edge list \tab .csv, .tsv \tab yes \tab yes \tab spreadsheets, any tool \cr
DOT \tab .dot, .gv \tab no \tab yes \tab Graphviz \cr
GeoPackage \tab .gpkg \tab yes \tab yes \tab QGIS, GDAL, sf (spatial networks; needs sf) \cr
}
Attributes keep their types where the format records them (GraphML, GEXF,
JSON) and are inferred otherwise. Node identifiers that are all whole numbers
are returned as integers. Formats differ in what they can hold: Pajek and
UCINET DL store only structure, node labels, and a single edge value, and DOT
is for drawing.

Files are read and written as UTF-8, and file paths may contain any
character, so identifiers and attributes in any script round-trip on every
platform. The exception is writing edge lists, which uses
\code{\link[utils]{write.table}}: in R before 4.2 on Windows, and in other
locales that are not UTF-8, it can only write characters of the native
encoding.

GraphML and DyNetML files both use \code{.xml}. When reading a \code{.xml}
file without \code{format}, the root element decides (\code{<graphml>} or a
DyNetML \code{<DynamicNetwork>}, \code{<DynamicMetaNetwork>},
\code{<MetaNetwork>}, or \code{<MetaMatrix>}). Writing a \code{.xml} file
without \code{format} produces GraphML; use \code{format = "dynetml"} (or
\code{"ora"}) to write DyNetML.
}
\value{
\code{net_read()} returns a \code{\link[=net_frame]{net_frame}}.
\code{net_write()} returns \code{file}, invisibly.
}
\section{GeoPackage}{
A GeoPackage (\code{.gpkg}) holds a spatial network (see
\code{\link{net_to_sf}}) as two layers, \code{nodes} and \code{edges}, each
written with \code{\link[sf]{st_write}} when it has geometry and as a plain
table otherwise, plus a small table layer \code{network} with the
direction, the identifier, endpoint, and weight columns, and network-level
attributes (as text). Reading needs only the \code{nodes} and \code{edges}
layers: without \code{network}, the \code{id}, \code{source}, and
\code{target} arguments name the columns, the network is directed unless
\code{directed} says otherwise, and the weight is found as for other
formats. Geometries and their CRS round-trip; attribute types follow GDAL's
GeoPackage driver (logical columns come back as integers). Requires the sf
package. \code{net_write()} replaces an existing file; on Windows a
GeoPackage that GDAL still holds open (for example one read earlier in the
session) cannot be replaced, so write to a new file name instead.
}
\section{DyNetML}{
DyNetML (Tsvetovat, Reminga, and Carley 2004) is the XML format of ORA for
multimode meta-networks: node classes such as agents, organizations,
locations, tasks, and resources, with several networks among them, and
optionally one meta-network per time period. It is read into one network:
\itemize{
\item Nodes get \code{class} (the node class id) and \code{type} (its ORA
type, such as \code{Agent}), \code{title}, and one column per node property,
typed from the file's property identities (\code{number} as double;
\code{text}, \code{categoryText}, \code{URI}, and \code{date} as text).
\item The same node id may appear in two classes. Only when it does are
identifiers written \code{class::id}, with the original in \code{node_id}.
\item Edges get \code{network} (the network id) and \code{weight} (the link
value; links without a value get 1).
\item Each network says whether it is directed. A file mixing directed and
undirected networks is read as directed, each undirected link stored as two
opposite edges, with a logical \code{directed} column marking them.
\item A file with several meta-networks (time periods) gets a \code{period}
column on nodes and edges, and identifiers \code{period::id}; each period's
nodes are separate rows.
}
Ports, survey (AnthroPac) data, and properties on networks are not supported
and give an error naming the line. Writing produces one meta-network per
\code{period} value, node classes from the \code{class} and \code{type}
columns (a single \code{Agent} class without them), and one network per
\code{network} value (or per pair of node classes); other node columns become
node properties. Logical node columns are written as ORA category text and
read back as text, and integer columns read back as double.
}
\references{
Tsvetovat, M., Reminga, J., and Carley, K. M. (2004). DyNetML: Interchange
format for rich social network data. Technical report CMU-ISRI-04-105,
Carnegie Mellon University. \doi{10.1184/R1/6622022.v1}

Altman, N., Carley, K. M., & Reminga, J. (2020). ORA User's Guide 2020.
Technical report CMU-ISR-20-110, Carnegie Mellon University.
}
\examples{
g <- net_frame(data.frame(source = c("a", "b"), target = c("b", "c"), w = c(1, 2)),
              weight = "w", weight_type = "distance", directed = FALSE)
path <- tempfile(fileext = ".graphml")
net_write(g, path)
net_read(path, weight_type = "distance")
}
