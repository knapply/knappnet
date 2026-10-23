\name{net_recode}
\alias{net_recode}
\alias{net_join}
\title{Recode and join node attributes}
\description{
\code{net_recode()} replaces values of a column of \code{x$nodes} (or
\code{x$edges}) according to a mapping, as UCINET's Transform > Recode does
for data values (Hanneman and Riddle 2005, chapter 6).
\code{net_join()} adds the columns of a data frame to
\code{x$nodes}, matched by node identifier (UCINET: Data > Join > Columns,
or loading an attribute dataset).
}
\usage{
net_recode(x, column, map, default = NULL, table = c("nodes", "edges"))

net_join(x, data, by = NULL)
}
\arguments{
\item{x}{A \code{\link{net_frame}}.}
\item{column}{The column to recode, by name or position.}
\item{map}{A named vector (old value = new value) or a function.}
\item{default}{Value for entries not named in \code{map} (ignored when
\code{map} is a function); \code{NULL} leaves them unchanged.}
\item{table}{\code{"nodes"} or \code{"edges"}.}
\item{data}{A data frame with one row per node to join.}
\item{by}{The column of \code{data} holding node identifiers, by name or
position; \code{NULL} uses a column with the same name as the network's
identifier column, or else the first column.}
}
\details{
\code{map} is a named vector, old value to new (\code{c(male = "M", female =
"F")}), or a function applied to the whole column. Values not named in the
vector are left unchanged unless \code{default} is given, in which case they
become \code{default}; a factor column is recoded through its labels and
returned as character. Recoding the identifier column rewrites the edge
endpoints to match, so new identifiers must stay unique.

\code{net_join()} matches \code{data[[by]]} to \code{x$nodes[[id]]};
nodes without a matching row get \code{NA}, rows of \code{data} without a
matching node are dropped, and duplicate keys in \code{data} are an error.
Columns of \code{data} that already exist in \code{x$nodes} replace them.
}
\value{A net_frame.}
\references{
Hanneman, R. A., & Riddle, M. (2005). \emph{Introduction to Social Network
Methods}. University of California, Riverside.
}
\seealso{\code{\link{net_filter_nodes}}, \code{\link{net_frame}}}
\examples{
ji <- net_example("jemaah_islamiyah")
table(ji$nodes$role)
# Collapse roles into two categories.
ji2 <- net_recode(ji, "role", c(`command team` = "command"), default = "operative")
table(ji2$nodes$role)

# Join an attribute table by name.
extra <- data.frame(name = ji$nodes$name[1:3], arrested = c(TRUE, FALSE, TRUE))
head(net_join(ji, extra)$nodes)
}
