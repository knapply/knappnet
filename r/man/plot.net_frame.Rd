\name{plot.net_frame}
\alias{plot.net_frame}
\title{Plot a network}
\description{
Draws a network with base R graphics. With no arguments beyond the
network, nodes are placed by \code{\link{net_layout}}, sized by degree,
and labelled when the network is small; edges are semi-transparent, with
arrows when the network is directed. Node and edge appearance can be set
from columns of \code{x$nodes} and \code{x$edges}.
}
\usage{
\method{plot}{net_frame}(x, layout = NULL, node_size = NULL, node_color = NULL,
     node_label = NULL, edge_width = NULL, edge_color = NULL,
     edge_linetype = NULL, curved = NULL, ..., legend = TRUE, seed = NULL)
}
\arguments{
\item{x}{A \code{\link{net_frame}}.}
\item{layout}{Node positions: \code{NULL} for the default of
\code{\link{net_layout}} (or the node locations of a spatial network;
see Details), \code{"geo"} for the node locations, the name of a
\code{net_layout()} method (such as \code{"fr"} or \code{"circle"}),
the names of two numeric columns of \code{x$nodes} (for example
\code{c("lon", "lat")}), a matrix or data frame with two columns and one
row per node (as returned by \code{net_layout()}), or a function that
takes the network and returns such a matrix.}
\item{node_size, node_color, node_label, edge_width, edge_color}{Appearance
of nodes and edges. Each can be the name of a column of \code{x$nodes} (for
node arguments) or \code{x$edges} (for edge arguments), a vector with one
value per node or edge, or a single value for all. Columns and vectors are
mapped to sizes, widths, or colors (see Details) and get a legend; wrap a
vector in \code{\link{I}()} to use its values as they are (sizes as
\code{cex} multipliers, widths as \code{lwd}, colors as R colors).
\code{NULL} uses the defaults (see Details). \code{node_label = TRUE}
labels nodes with their identifiers and \code{FALSE} draws no labels.}
\item{edge_linetype}{Line types of edges: the name of a column of
\code{x$edges} or a vector with one value per edge (at most six
categories, mapped to solid, dashed, dotted, dot-dash, long dash, and
two-dash, with a legend), a single line type for all, or a vector wrapped
in \code{\link{I}()} of R line types. \code{NULL} draws negative ties of
a signed network dashed and everything else solid.}
\item{curved}{Whether edges are drawn as arcs: \code{NULL} (the default)
curves only edges that share a pair of nodes (parallel edges, and the two
directions of a reciprocated tie), so that each stays visible;
\code{TRUE} curves every edge; \code{FALSE} draws straight lines.}
\item{...}{Further arguments passed to \code{\link[graphics]{title}}, such
as \code{main}.}
\item{legend}{Whether to draw legends for mapped columns and vectors.}
\item{seed}{Seed for the layout, as in \code{\link{net_layout}}, when
\code{layout} is \code{NULL} or a method name.}
}
\details{
Defaults:
\itemize{
  \item Node size grows with degree (as area), and nodes are smaller in
    larger networks.
  \item Nodes are a single blue with a white border, except in networks
    with modes, where each mode has its own color; two-mode networks
    (\code{\link{net_twomode}}) are laid out in two rows, one per mode.
  \item Nodes are labelled with their identifiers in networks of at most
    40 nodes.
  \item Edges are gray and more transparent the more edges there are,
    except that edges are coloured by relation in a multiplex network and
    by layer in a multilayer network (with edges between layers in gray).
    In directed networks, arrows end at the border of the target node.
  \item Negative ties of a signed network are dashed, and unobserved
    pairs (\code{edge_missing} in \code{\link{net_frame}}) are dotted and
    light gray, whatever else is mapped.
  \item Edge width is constant, except in networks whose weights are
    strengths (\code{weight_type = "strength"} in \code{\link{net_frame}}),
    where it grows with the weight. Distance weights do not change the
    width, since a thick line would suggest a strong tie.
}

Mapping: numeric values give sizes and widths that grow with the value
(node area grows linearly) and colors from the sequential viridis palette
(\code{grDevices::hcl.colors(, "viridis")}; light for small and dark for
large values). Factors, character vectors, and logical vectors give colors
from the qualitative palette of Okabe and Ito (2008), as
\code{grDevices::palette.colors(, "Okabe-Ito")} ships it with black left
out, in the order of the factor levels or of the sorted values; beyond
eight categories the HCL palette \code{"Dark 3"} of Zeileis et al. (2020)
is used instead, which stays distinguishable in sequence but is not
guaranteed to be safe for every kind of color vision deficiency. Every
plot in the package uses these two palettes, so a category keeps its color
from one plot to the next. Missing values are drawn in light gray (colors)
or at the smallest size or width. Sizes and widths must be numeric.

Legends are drawn in the right margin. Coordinates keep their aspect
ratio, so distances in the plot are proportional to distances in the
layout. Self-loops are drawn as small circles beside their node, and
parallel edges as arcs (see \code{curved}). For one panel per layer,
time window, or ego, see \code{\link{net_plot_layers}},
\code{\link{net_plot_snapshots}}, and \code{\link{net_plot_egos}}.

The plot is reproducible for a given \code{seed} (or after
\code{\link{set.seed}}).

Spatial networks (those with sf node geometry, or coordinate columns
recorded by \code{\link{net_spatial}} or
\code{\link{net_proximity}}) are drawn at their node locations by
\code{\link{net_plot_spatial}} when \code{layout} is \code{NULL} or
\code{"geo"}, with the same appearance arguments; any other \code{layout}
draws them like any network. \code{"geo"} is an error for a network without
locations.
}
\value{The layout used, invisibly: a matrix with columns \code{x} and
\code{y} and one row per node, which can be passed back as \code{layout}.
For a spatial network drawn at its locations, \code{x} invisibly, as
\code{\link{net_plot_spatial}} returns.}
\references{
Okabe, M., & Ito, K. (2008). Color Universal Design (CUD): How to make
figures and presentations that are friendly to colorblind people.
\url{https://jfly.uni-koeln.de/color/}

Zeileis, A., Fisher, J. C., Hornik, K., Ihaka, R., McWhite, C. D., Murrell,
P., Stauffer, R., & Wilke, C. O. (2020). colorspace: A toolbox for
manipulating and assessing colors and palettes. \emph{Journal of
Statistical Software}, 96(1), 1-49. \doi{10.18637/jss.v096.i01}
}
\seealso{\code{\link{net_layout}}, \code{\link{net_plot_spatial}}, and the
other plotting functions: \code{\link{net_plot_matrix}},
\code{\link{net_plot_centrality}}, \code{\link{net_plot_partition}},
\code{\link{net_plot_degree}}.}
\examples{
g <- net_frame(data.frame(source  = c("a", "a", "b", "c", "c", "d"),
                         target    = c("b", "c", "c", "d", "e", "e"),
                         calls = c(1, 4, 2, 8, 3, 5)),
              directed = FALSE)
g$nodes$group <- c("x", "x", "y", "y", "y")
plot(g, seed = 1)
plot(g, node_color = "group", edge_width = "calls", seed = 1, main = "Calls")
xy <- plot(g, layout = "circle", node_label = FALSE)
plot(g, layout = xy, node_size = I(c(3, 1, 1, 1, 1)))

attended <- data.frame(person = c("ann", "ann", "bob", "cat"),
                       event  = c("gala", "fair", "gala", "fair"))
# Two modes in two rows.
plot(net_twomode(attended, source = "person", target = "event"))

# Line types from a column; parallel edges as arcs and a self-loop.
m <- net_frame(data.frame(source = c("a", "a", "a", "b", "c"),
                          target = c("b", "b", "a", "c", "a"),
                          kind = c("call", "visit", "call", "call", "visit")))
plot(m, edge_linetype = "kind", seed = 1)

# Edges coloured by relation in a multiplex network.
plot(net_example("noordin", relation = c("classmates", "kinship")), node_label = FALSE,
     seed = 1)
}
