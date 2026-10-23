\name{net_split_layers}
\alias{net_split_layers}
\alias{net_edge_layers}
\alias{net_is_multilayer}
\alias{net_as_multilayer}
\title{Multilayer networks}
\description{
A multilayer network has one set of actors, several layers, and edges
that name the layer of each endpoint (Kivela et al. 2014). Actors stay in
the node table, one row each; the \code{source_layer} and
\code{target_layer} columns of \code{\link{net_frame}} give the layers.
\code{net_edge_layers()} classifies every edge, \code{net_split_layers()}
returns one network per layer, \code{net_is_multilayer()} tests for the
class, and \code{net_as_multilayer()} lifts a multiplex network (one with
a \code{relation} role) to layer form, adding coupling edges.
}
\usage{
net_split_layers(x)

net_edge_layers(x)

net_is_multilayer(x)

net_as_multilayer(x, couplings = c("categorical", "ordinal", "none"),
                  relation = NULL, coupling_weight = 1)
}
\arguments{
\item{x}{A multilayer net_frame (class \code{multilayer_net_frame}); for
\code{net_as_multilayer()}, a net_frame with a relation column; for
\code{net_is_multilayer()}, any object.}
\item{couplings}{Which coupling edges join each actor to itself across
layers: \code{"categorical"} (the default) between every pair of layers,
\code{"ordinal"} between consecutive layers in order of first appearance
(for layers that are ordered, such as time periods), or \code{"none"}.
In a directed network each coupling is added in both directions.}
\item{relation}{The name of the column of \code{x$edges} whose values
become the layers; \code{NULL} uses the relation role.}
\item{coupling_weight}{The weight of the coupling edges of a weighted
network (positive when the weights are distances).}
}
\details{
An edge with endpoints \eqn{(a, b)} in layers \eqn{(l_a, l_b)} is
\emph{intralayer} when \eqn{l_a = l_b}, a \emph{coupling} when
\eqn{l_a \ne l_b} and \eqn{a = b}, and \emph{interlayer} otherwise
(Kivela et al. 2014; De Domenico et al. 2013). A coupling is not a
self-loop, and two edges in different layers are never parallel. Every
actor is taken to be present in every layer (the network is
node-aligned). Mucha et al. (2010) call layers slices.

A multiplex network (several relations on one node set, the
\code{relation} role) is the special case whose only edges between layers
are couplings; it needs no layer columns, and most measures take
\code{relation =} directly. \code{net_as_multilayer()} copies the relation
into both layer columns, adds the couplings, records the layer labels in
a layer table, and declares \code{couplings_only = TRUE}. Temporal
networks are refused, since coupling edges would need times of their own.

\code{net_split_layers()} keeps every actor in every layer, with that
layer's intralayer edges and the attributes of the original rows; couplings
and interlayer edges are left out. \code{\link{net_plot_layers}} draws the
layers side by side with a shared layout.
}
\value{
\code{net_split_layers()}: a named list of net_frames, one per layer, in
the order of the layer labels. \code{net_edge_layers()}: a data frame
aligned with the rows of \code{x$edges}, with columns \code{source_layer},
\code{target_layer}, and \code{class} (\code{"intralayer"},
\code{"coupling"}, or \code{"interlayer"}). \code{net_is_multilayer()}:
\code{TRUE} or \code{FALSE}. \code{net_as_multilayer()}: a multilayer
net_frame.
}
\references{
De Domenico, M., Sole-Ribalta, A., Cozzo, E., Kivela, M., Moreno, Y.,
Porter, M. A., Gomez, S., & Arenas, A. (2013). Mathematical formulation
of multilayer networks. \emph{Physical Review X}, 3(4), 041022.
\doi{10.1103/PhysRevX.3.041022}

Kivela, M., Arenas, A., Barthelemy, M., Gleeson, J. P., Moreno, Y., &
Porter, M. A. (2014). Multilayer networks. \emph{Journal of Complex
Networks}, 2(3), 203-271. \doi{10.1093/comnet/cnu016}

Mucha, P. J., Richardson, T., Macon, K., Porter, M. A., & Onnela, J.-P.
(2010). Community structure in time-dependent, multiscale, and multiplex
networks. \emph{Science}, 328(5980), 876-878.
\doi{10.1126/science.1184819}
}
\seealso{
\code{\link{net_frame}} for the layer roles, the layer table, and the
\code{couplings_only} declaration; \code{\link{net_plot_layers}}.
}
\examples{
# Two actors in two layers, with a coupling and a layer table.
edges <- data.frame(source = c("ann", "ann", "ann"), target = c("bob", "bob", "ann"),
                    from_layer = c("work", "online", "work"),
                    to_layer = c("work", "online", "online"))
layers <- data.frame(layer = c("work", "online"), medium = c("office", "chat"))
g <- net_frame(edges, directed = FALSE, source_layer = "from_layer",
               target_layer = "to_layer", layers = layers)
g
class(g)
net_edge_layers(g)
net_is_simple(g)   # the coupling is not a loop

# Krackhardt's three relations as layers coupled across every pair.
krack <- net_example("krackhardt")
ml <- net_as_multilayer(krack)
ml
table(net_edge_layers(ml)$class)
lapply(net_split_layers(ml), net_density)
net_plot_layers(ml, node_label = FALSE, seed = 1)
}
