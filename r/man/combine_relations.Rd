\name{combine_relations}
\alias{combine_relations}
\title{Combine the relations of a multiplex network}
\usage{
combine_relations(x, relation = NULL,
                  method = c("union", "intersection", "sum", "difference", "select"),
                  which = NULL)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{relation}{The edge column that names each edge's relation, by name or
position; \code{NULL} looks for a column named \code{relation} or
\code{network}.}
\item{method}{How to combine; see Description.}
\item{which}{Which relations (values of the relation column) to combine;
\code{NULL} means all of them. \code{"difference"} needs exactly two, in
order, and \code{"select"} exactly one.}
}
\value{A netframe.}
\description{
A multiplex network stores several relations in one edge table, told apart
by a column such as \code{relation} (as \code{example_network("noordin")}
and files read from DyNetML have). \code{combine_relations()} builds a
one-relation network from them:
\describe{
  \item{\code{"union"}}{A pair is tied if it is tied in any of the relations
    (UCINET: Transform > Matrix Operations > Between Datasets > Statistical
    Summaries, Maximum; Boolean OR).}
  \item{\code{"intersection"}}{A pair is tied if it is tied in every one of
    the relations (Statistical Summaries, Minimum; Boolean AND).}
  \item{\code{"sum"}}{A pair is tied with weight equal to the number of
    relations in which it is tied, or the sum of its weights over the
    relations when the network is weighted (Statistical Summaries, Sum).}
  \item{\code{"difference"}}{A pair is tied if it is tied in the first of
    two relations and not in the second (Transform > Matrix Operations >
    Between Datasets > Boolean Combinations).}
  \item{\code{"select"}}{The edges of one relation, with all their columns
    (Data > Unpack).}
}
}
\details{
Within one relation parallel edges count once (a pair is tied or not),
except that \code{"sum"} adds the weights of parallel edges when the network
is weighted. Self-loops are kept as ties of a node with itself. Pairs are
ordered in directed networks and unordered in undirected ones. The result
keeps the node table and the network's direction; its edge table has
\code{from}, \code{to}, and \code{relations} (the number of selected
relations in which the pair is tied), plus \code{weight} for \code{"sum"},
which becomes the weight column (with the network's \code{weight_type}, or
\code{"strength"} for counts). Other edge columns are dropped, since they
cannot be merged across relations; \code{"select"} keeps every column.

The UCINET menu paths follow Hanneman and Riddle (2005, chapter 16), which
also describes Transform > Multiplex, the inverse operation that stacks
relations into one nominal-valued matrix.
}
\references{
Hanneman, R. A., & Riddle, M. (2005). \emph{Introduction to Social Network
Methods}. University of California, Riverside.

Borgatti, S. P., Everett, M. G., & Freeman, L. C. (2002). \emph{UCINET for
Windows: Software for Social Network Analysis}. Analytic Technologies.
}
\seealso{\code{\link{example_network}}, \code{\link{dichotomize}},
\code{\link{symmetrize}}}
\examples{
noordin <- example_network("noordin")
table(noordin$edges$relation)

# Every relation in which a pair is tied: the "alive" network of Everton
# (2012), and the number of relations per pair as a weight.
combine_relations(noordin, method = "union")
combined <- combine_relations(noordin, method = "sum")
table(combined$edges$weight)

# Pairs tied by both friendship and an operation, friends who never shared
# an operation, and a single relation.
combine_relations(noordin, method = "intersection", which = c("friendship", "operations"))
combine_relations(noordin, method = "difference", which = c("friendship", "operations"))
combine_relations(noordin, method = "select", which = "kinship")
}
