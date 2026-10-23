\name{net_auto_equiv}
\alias{net_auto_equiv}
\title{Automorphic equivalence}
\description{
Nodes are automorphically equivalent when some automorphism of the
network, a relabelling of the nodes that maps the network onto itself,
maps one onto the other (Everett 1985; Borgatti and Everett 1992): they
are interchangeable, and have the same value on every measure that depends
on the structure alone. The classes are the orbits of the automorphism
group, computed exactly by an individualization-refinement search (McKay
and Piperno 2014).
}
\usage{
net_auto_equiv(x, max_steps = 1e6)
}
\arguments{
\item{x}{A \code{\link{net_frame}}.}
\item{max_steps}{A cap on the number of refined colourings the search
examines. The search is exponential in the worst case (highly symmetric
networks such as strongly regular graphs), so it stops with an error,
rather than run indefinitely, once it has examined this many; results are
never approximate.}
}
\details{
The equivalence is defined on binary ties: weights are ignored and parallel
edges count once. Direction is respected: an automorphism must preserve
every arc. A self-loop is a property of its node, so a node with a loop
can only map to a node with a loop. Classes are numbered from 1 in order
of their first node.

Automorphic equivalence implies regular equivalence and is implied by
structural equivalence. Networks of a few hundred nodes are usually quick.
igraph and NetworkX have no function for the orbits, and UCINET's routine
(MAXSIM) computes an approximation based on distance profiles rather than
the exact orbits; the validation suite checks the result against the
orbits found by trying every permutation of small networks.
}
\value{An integer vector of class numbers (starting at 1) aligned with the
rows of \code{x$nodes}.}
\references{
Everett, M. G. (1985). Role similarity and complexity in social networks.
\emph{Social Networks}, 7(4), 353-359. \doi{10.1016/0378-8733(85)90013-9}

Borgatti, S. P., & Everett, M. G. (1992). Notions of position in social
network analysis. \emph{Sociological Methodology}, 22, 1-35.
\doi{10.2307/270991}

McKay, B. D., & Piperno, A. (2014). Practical graph isomorphism, II.
\emph{Journal of Symbolic Computation}, 60, 94-112.
\doi{10.1016/j.jsc.2013.09.003}
}
\seealso{\code{\link{net_regular_equiv}},
\code{\link{net_struct_equiv}}.}
\examples{
g <- net_frame(data.frame(source = c(1, 1, 2, 2, 3, 3), target = c(2, 3, 4, 5, 6, 7)))
net_auto_equiv(g)  # {1}, {2, 3}, {4, 5, 6, 7}

# Two triangles joined through a middle node: the outer corners are all
# interchangeable, as are the two inner ones.
h <- net_frame(data.frame(source = c(1, 1, 2, 3, 4, 5, 5, 6),
                         target   = c(2, 3, 3, 4, 5, 6, 7, 7)), directed = FALSE)
net_auto_equiv(h)  # {1, 2, 6, 7}, {3, 5}, {4}
}
