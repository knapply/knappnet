\name{saom_effects}
\alias{saom_effects}
\alias{saom_effect_names}
\title{Effects of a stochastic actor-oriented model}
\usage{
saom_effects(formula, data = NULL)

saom_effect_names()
}
\arguments{
\item{formula}{A one-sided formula whose terms are effect names, for
example \code{~ density + recip + transTrip + simX("smoke")}. Structural
effects are bare names; covariate effects are calls with the covariate name
as a string or symbol; the GWESP effects accept an internal parameter,
\code{gwespFF(parameter = 69)} (the default, giving
\eqn{\alpha = 0.69 = \log 2}{alpha = 0.69 = log 2}). A character vector of
structural effect names, or an effects data frame returned by this
function, is also accepted.}
\item{data}{Optional \code{\link{saom_data}} whose covariates the formula
may use; covariate names are then checked.}
}
\value{
\code{saom_effects()}: a data frame with one row per effect: \code{effect}
(the short name), \code{covariate}, \code{parameter}, \code{label}, and
\code{name} (the coefficient name). \code{saom_effect_names()}: the short
names of all available effects.
}
\description{
Turns a model formula into the table of effects \code{\link{saom_fit}}
estimates, and lists the effects available. Effects carry the short names
of RSiena (Ripley et al. 2026, chapter 13), and their definitions match
RSiena's exactly; see Details.
}
\details{
With \eqn{x_{ij} = 1}{x_ij = 1} for a tie from \eqn{i} to \eqn{j},
\eqn{x_{i+}}{x_i+} the out-degree, \eqn{x_{+i}}{x_+i} the in-degree,
\eqn{v_i} a centered actor covariate and \eqn{w_{ij}}{w_ij} a centered
dyadic covariate, the evaluation function of actor \eqn{i} is
\eqn{f_i(x) = \sum_k \beta_k s_{ik}(x)}{f_i(x) = sum_k beta_k s_ik(x)} with
these statistics \eqn{s_{ik}(x)}{s_ik(x)}:
\tabular{lll}{
  \code{density} \tab out-degree \tab \eqn{x_{i+}}{x_i+} \cr
  \code{recip} \tab reciprocity \tab \eqn{\sum_j x_{ij} x_{ji}}{sum_j x_ij x_ji} \cr
  \code{transTrip} \tab transitive triplets \tab
    \eqn{\sum_{j,h} x_{ij} x_{ih} x_{hj}}{sum_jh x_ij x_ih x_hj} \cr
  \code{transTies} \tab transitive ties \tab
    \eqn{\sum_j x_{ij} \max_h x_{ih} x_{hj}}{sum_j x_ij max_h x_ih x_hj} \cr
  \code{transRecTrip} \tab transitive reciprocated triplets \tab
    \eqn{\sum_{j,h} x_{ij} x_{ji} x_{ih} x_{hj}}{sum_jh x_ij x_ji x_ih x_hj} \cr
  \code{cycle3} \tab 3-cycles \tab
    \eqn{\sum_{j,h} x_{ij} x_{jh} x_{hi}}{sum_jh x_ij x_jh x_hi} \cr
  \code{gwespFF}, \code{gwespFB}, \code{gwespBF}, \code{gwespBB} \tab
    geometrically weighted edgewise shared partners \tab
    \eqn{\sum_j x_{ij} e^\alpha \{1 - (1 - e^{-\alpha})^{esp_{ij}}\}}{sum_j x_ij e^a (1 - (1 - e^-a)^esp_ij)} \cr
  \code{inPop}, \code{inPopSqrt} \tab in-degree popularity \tab
    \eqn{\sum_j x_{ij} x_{+j}}{sum_j x_ij x_+j},
    \eqn{\sum_j x_{ij} \sqrt{x_{+j}}}{sum_j x_ij sqrt(x_+j)} \cr
  \code{outPop}, \code{outPopSqrt} \tab out-degree popularity \tab
    \eqn{\sum_j x_{ij} x_{j+}}{sum_j x_ij x_j+},
    \eqn{\sum_j x_{ij} \sqrt{x_{j+}}}{sum_j x_ij sqrt(x_j+)} \cr
  \code{inAct}, \code{inActSqrt} \tab in-degree activity \tab
    \eqn{x_{i+} x_{+i}}{x_i+ x_+i}, \eqn{x_{i+} \sqrt{x_{+i}}}{x_i+ sqrt(x_+i)} \cr
  \code{outAct}, \code{outActSqrt} \tab out-degree activity \tab
    \eqn{x_{i+}^2}{x_i+^2}, \eqn{x_{i+}^{1.5}}{x_i+^1.5} \cr
  \code{egoX} \tab covariate ego \tab \eqn{v_i x_{i+}}{v_i x_i+} \cr
  \code{altX} \tab covariate alter \tab \eqn{\sum_j x_{ij} v_j}{sum_j x_ij v_j} \cr
  \code{simX} \tab covariate similarity \tab
    \eqn{\sum_j x_{ij} (sim_{ij} - \overline{sim})}{sum_j x_ij (sim_ij - mean sim)} \cr
  \code{sameX} \tab same covariate \tab
    \eqn{\sum_j x_{ij} I\{v_i = v_j\}}{sum_j x_ij I(v_i = v_j)} \cr
  \code{egoXaltX} \tab covariate ego x alter \tab
    \eqn{v_i \sum_j x_{ij} v_j}{v_i sum_j x_ij v_j} \cr
  \code{X} \tab dyadic covariate \tab \eqn{\sum_j x_{ij} w_{ij}}{sum_j x_ij w_ij}
}
For the GWESP effects \eqn{esp_{ij}}{esp_ij} counts the actors \eqn{h} with
\eqn{i \to h \to j}{i -> h -> j} (\code{gwespFF}),
\eqn{i \to h \leftarrow j}{i -> h <- j} (\code{gwespFB}),
\eqn{i \leftarrow h \to j}{i <- h -> j} (\code{gwespBF}), or
\eqn{i \leftarrow h \leftarrow j}{i <- h <- j} (\code{gwespBB}), and
\eqn{\alpha}{alpha} is \code{parameter / 100} (Hunter 2007; Ripley et al.
2026). As in RSiena these are elementary effects: the change statistic of
a tie is that tie's own weight, not the change in the sum over all of the
actor's ties. The similarity is
\eqn{sim_{ij} = 1 - |v_i - v_j| / range}{sim_ij = 1 - |v_i - v_j| / range}
with the range taken over all actors and waves, centered by its mean over
all ordered pairs of distinct actors (and waves).

Covariates are centered as described in \code{\link{saom_data}}. The
statistics of every effect, and the change statistics
\eqn{s_{ik}(x^{+ij}) - s_{ik}(x^{-ij})}{s_ik(x+ij) - s_ik(x-ij)} that
ministeps use, are checked against a direct evaluation of these formulas
and against RSiena's own target statistics in the package's validation
suite.
}
\references{
Ripley, R. M., Snijders, T. A. B., Boda, Z., and Voros, A. (2026).
\emph{Manual for RSiena}. University of Oxford, Department of Statistics and
Nuffield College; University of Groningen, Department of Sociology.
\url{https://www.stats.ox.ac.uk/~snijders/siena/RSiena_Manual.pdf}

Hunter, D. R. (2007). Curved exponential family models for social networks.
\emph{Social Networks}, 29(2), 216-230. \doi{10.1016/j.socnet.2006.08.005}
}
\seealso{
\code{\link{saom_fit}}
}
\examples{
saom_effects(~ density + recip + transTrip + gwespFF(parameter = 69) + simX("smoke"))
saom_effect_names()
}
