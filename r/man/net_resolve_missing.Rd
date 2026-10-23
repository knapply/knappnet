\name{net_resolve_missing}
\alias{net_resolve_missing}
\title{Resolve missing ties and non-respondents}
\description{
A network may record what was not observed: unobserved pairs, flagged by
the \code{edge_missing} column of \code{\link{net_frame}}, and
non-respondents, flagged by its \code{node_missing} column.
\code{net_resolve_missing()} returns the network with that missingness
resolved by an explicit rule, so that measures can use it.
}
\usage{
net_resolve_missing(x, na_action)
}
\arguments{
\item{x}{A \code{\link{net_frame}}.}
\item{na_action}{How to resolve the missing data, with no default:
\code{"absent"} treats every unobserved pair as no tie (its rows are
dropped from \code{x$edges}) and the unreported ties of non-respondents
as absent; \code{"present"} treats every unobserved pair as a tie (its
rows are kept as ordinary edges).}
}
\details{
An absent tie was observed not to exist and has no row. A missing tie is
a pair whose status was not observed: a row of the edge table with
\code{edge_missing = TRUE}, kept in the table (with its relation, layers,
and spell, which say where it is unobserved) but left out of the graph
that measures use. A non-respondent is a node whose own reports were not
collected; ties to it that others reported are known. In a directed
network the out-ties of a non-respondent are unobserved; in an undirected
one, the pairs between two non-respondents are. These implied pairs are
never written out (Kossinets 2006; Rawlings et al. 2023, ch. 4).

Treating missing ties as absent is what most software does silently, and
it biases density, degree, and path lengths downwards when missingness is
substantial (Kossinets 2006). It is offered here as an explicit choice so
that the analysis states it. \code{"present"} is refused when there are
non-respondents, since it would tie each of them to every node whose tie
it did not report; impute those ties as rows of \code{x$edges} instead.

The result has neither missing-data role. The flag columns stay in the
tables as ordinary attributes.
}
\value{A \code{\link{net_frame}} without missing-data roles.}
\references{
Kossinets, G. (2006). Effects of missing data in social networks.
\emph{Social Networks}, 28(3), 247-268.
\doi{10.1016/j.socnet.2005.07.002}

Rawlings, C. M., Smith, J. A., Moody, J., & McFarland, D. A. (2023).
\emph{Network Analysis: Integrating Social Network Theory, Method, and
Application with R}. Cambridge University Press.
\url{https://inarwhal.github.io/NetworkAnalysisR-book/}
}
\seealso{
\code{\link{net_frame}} for the \code{edge_missing} and
\code{node_missing} roles.
}
\examples{
edges <- data.frame(source = c("ann", "ann", "bob", "cy"),
                    target = c("bob", "cy", "cy", "dee"),
                    unknown = c(FALSE, TRUE, FALSE, FALSE))
g <- net_frame(edges, directed = FALSE, edge_missing = "unknown")
g

# The unobserved pair is dotted; resolving it either way changes the graph.
op <- par(mfrow = c(1, 3))
xy <- plot(g, main = "unobserved pair", seed = 1)
plot(net_resolve_missing(g, "absent"), layout = xy, main = "as absent")
plot(net_resolve_missing(g, "present"), layout = xy, main = "as present")
par(op)

net_density(net_resolve_missing(g, "absent"))
net_density(net_resolve_missing(g, "present"))

# A non-respondent.
people <- data.frame(id = c("ann", "bob", "cy", "dee"),
                     no_reply = c(FALSE, FALSE, TRUE, FALSE))
h <- net_frame(edges[-2, ], people, directed = TRUE, node_missing = "no_reply")
h
net_resolve_missing(h, "absent")
}
