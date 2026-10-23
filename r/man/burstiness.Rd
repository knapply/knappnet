\name{burstiness}
\alias{burstiness}
\alias{inter_event_times}
\title{Inter-event times and burstiness}
\description{
\code{inter_event_times()} lists the times between consecutive events of
each node or pair of nodes; \code{burstiness()} summarizes them with Goh
and Barabási's (2008) burstiness and memory coefficients.
}
\usage{
burstiness(x, by = c("node", "pair"), mode = c("all", "out", "in"))

inter_event_times(x, by = c("node", "pair"), mode = c("all", "out", "in"))
}
\arguments{
\item{x}{A temporal network (\code{\link{temporal_network}}).}
\item{by}{\code{"node"} or \code{"pair"}.}
\item{mode}{For nodes of directed networks: \code{"all"}, \code{"out"},
or \code{"in"}. Ignored for pairs and for undirected networks.}
}
\details{
The events of a node are the start times of its contacts: all of them
(\code{mode = "all"}), or in directed networks those it sends
(\code{"out"}) or receives (\code{"in"}); a self-loop is one event. The
events of a pair are the start times of the contacts between the two
nodes (ordered pairs in directed networks, unordered otherwise).
Simultaneous events give inter-event times of 0. Durations, weights, and
node presence intervals are not used.

With inter-event times \eqn{\tau_1, \dots, \tau_k} of mean \eqn{\mu} and
standard deviation \eqn{\sigma}, the burstiness is
\deqn{B = \frac{\sigma - \mu}{\sigma + \mu},}{B = (sigma - mu) / (sigma + mu),}
which is -1 for perfectly regular events, about 0 for a Poisson process,
and approaches 1 for very bursty activity (Goh and Barabási 2008, Eq. 1).
The memory coefficient is the correlation of consecutive inter-event
times (their Eq. 4),
\deqn{M = \frac{1}{k - 1} \sum_{i=1}^{k-1} \frac{(\tau_i - \mu_1)(\tau_{i+1} - \mu_2)}{\sigma_1 \sigma_2},}{M = (1 / (k - 1)) sum_i (tau_i - mu_1)(tau_{i+1} - mu_2) / (sigma_1 sigma_2),}
where \eqn{\mu_1, \sigma_1} (\eqn{\mu_2, \sigma_2}) are the mean and
standard deviation of \eqn{\tau_1, \dots, \tau_{k-1}}
(\eqn{\tau_2, \dots, \tau_k}). Standard deviations divide by the number
of values (the standard deviation of the observed distribution), so
\eqn{M} is exactly the Pearson correlation; \code{sd()} divides by one
less, which gives a slightly different \eqn{B}. \eqn{B} is \code{NaN}
with fewer than two inter-event times or when all are 0, and \eqn{M} with
fewer than three or when either sequence is constant. Neither coefficient
corrects for the number of events (Kim and Jo 2016 discuss finite-size
effects), so compare nodes with similar activity.
}
\value{
\code{inter_event_times()} returns a data frame with the node (or the
pair, as \code{from} and \code{to}) and \code{gap}, the inter-event time
in the unit of the times, ordered by node (or pair) and then time.
\code{burstiness()} returns a data frame with one row per node (aligned
with \code{x$nodes}) or per pair (in order of first contact), with
\code{events}, \code{burstiness}, and \code{memory}.
}
\references{
Goh, K.-I., & Barabási, A.-L. (2008). Burstiness and memory in complex
systems. \emph{EPL (Europhysics Letters)}, 81(4), 48002.
\doi{10.1209/0295-5075/81/48002}

Kim, E.-K., & Jo, H.-H. (2016). Measuring burstiness for finite event
sequences. \emph{Physical Review E}, 94(3), 032311.
\doi{10.1103/PhysRevE.94.032311}
}
\seealso{
\code{\link{temporal_network}}, \code{\link{aggregate_network}} for
contact counts and durations per pair.
}
\examples{
ev <- data.frame(from = c("a", "a", "a", "a", "b", "b", "b"),
                 to = c("b", "b", "c", "c", "c", "c", "c"),
                 time = c(0, 1, 2, 10, 0, 5, 10))
g <- temporal_network(ev, directed = FALSE)
inter_event_times(g)
burstiness(g)
burstiness(g, by = "pair")
# Events a node sends or receives in a directed network.
burstiness(temporal_network(ev), mode = "out")
}
