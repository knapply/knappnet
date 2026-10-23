\name{net_as_signed}
\alias{net_as_signed}
\title{Split signed values into sign and magnitude}
\description{
Signed ties often arrive as one numeric column whose negative values are
negative ties (as in Pajek and UCINET files). A signed network here keeps
the sign and the magnitude apart, so that a negative tie is never read as
a negative length or strength. \code{net_as_signed()} splits such a
column of an edge table into a sign column (\code{-1} or \code{1}) and a
non-negative weight column, ready for \code{net_frame(sign = , weight = )}.
}
\usage{
net_as_signed(edges, value, sign = "sign", weight = "weight")
}
\arguments{
\item{edges}{A data frame of edges.}
\item{value}{The name of the numeric column of \code{edges} holding the
signed values. Zero and missing values are an error, since a zero is
neither positive nor negative; drop those rows first.}
\item{sign, weight}{The names of the columns to write (replacing columns of
those names): the sign as an integer \code{-1} or \code{1}, and the
absolute value.}
}
\details{
A signed graph gives every edge a positive or negative sign (Harary 1953;
Cartwright and Harary 1956). The \code{sign} role of
\code{\link{net_frame}} holds it; a weight, if any, is the magnitude and
must not be negative. Plots draw negative ties dashed. The sign role is
provisional: it is validated and stored, and no measure uses it yet.
}
\value{\code{edges} with the two columns added or replaced.}
\references{
Cartwright, D., & Harary, F. (1956). Structural balance: A generalization
of Heider's theory. \emph{Psychological Review}, 63(5), 277-293.
\doi{10.1037/h0046049}

Harary, F. (1953). On the notion of balance of a signed graph.
\emph{Michigan Mathematical Journal}, 2(2), 143-146.
\doi{10.1307/mmj/1028989917}
}
\seealso{
\code{\link{net_frame}} for the \code{sign} role.
}
\examples{
ratings <- data.frame(source = c("ann", "ann", "bob", "cy"),
                      target = c("bob", "cy", "cy", "dee"),
                      value = c(3, -2, 1, -1))
signed <- net_as_signed(ratings, "value")
signed
g <- net_frame(signed, directed = FALSE, sign = "sign", weight = "weight",
               weight_type = "strength")
g
plot(g, seed = 1)

# Alliances (positive) and enmities (negative) among the Gahuku-Gama.
tribes <- net_example("highland_tribes")
tribes
plot(tribes, seed = 1)
}
