# Temporary: runs the R code of the articles (not part of the package build).

test_that("article knappnet.Rmd runs", {
  skip_on_cran()
  skip_if_not_installed("igraph")
  skip_if_not_installed("network")
  on_null_device(local({
has_igraph <- requireNamespace("igraph", quietly = TRUE)
has_network <- requireNamespace("network", quietly = TRUE)
library(knappnet)
ties <- data.frame(
  source = c("Acciaiuoli", "Medici", "Medici", "Medici", "Medici", "Medici",
           "Castellani", "Castellani", "Castellani", "Peruzzi", "Peruzzi",
           "Strozzi", "Strozzi", "Ridolfi", "Tornabuoni", "Albizzi", "Albizzi",
           "Salviati", "Bischeri", "Guadagni"),
  target   = c("Medici", "Barbadori", "Ridolfi", "Tornabuoni", "Albizzi", "Salviati",
           "Peruzzi", "Strozzi", "Barbadori", "Strozzi", "Bischeri",
           "Ridolfi", "Bischeri", "Tornabuoni", "Guadagni", "Ginori", "Guadagni",
           "Pazzi", "Guadagni", "Lamberteschi")
)
families <- data.frame(name = sort(unique(c(ties$source, ties$target, "Pucci"))))
flo <- net_frame(ties, nodes = families, directed = FALSE, id = "name")
flo
try({
bad <- data.frame(source = c("Medici", NA, "Strozzi"), target = c("Pazzi", "Ginori", "Medcii"))
net_frame(bad, nodes = families, directed = FALSE, id = "name")
}, silent = TRUE)
plot(flo, seed = 1)
head(flo$nodes)
head(flo$edges)
flo$nodes$degree <- net_degree(flo)
flo$nodes$betweenness <- net_betweenness(flo)
flo$nodes$harmonic <- round(net_harmonic(flo, normalized = TRUE), 3)
flo$nodes$core <- net_coreness(flo)
flo$nodes[order(-flo$nodes$betweenness), ]
plot(flo, node_size = "betweenness", node_color = "core", seed = 1)
net_closeness(flo)[flo$nodes$name == "Pucci"]
net_density(flo)
net_centralization(flo, "degree")
comp <- net_components(flo)
comp$no
comp$csize
net_max_cliques(flo, min_size = 3)
groups <- net_louvain(flo, seed = 1)
table(groups)
roads <- data.frame(
  source    = c("a", "b", "a"),
  target      = c("b", "c", "c"),
  minutes = c(5, 5, 20)
)
trips <- net_frame(roads, directed = FALSE, weight = "minutes", weight_type = "distance")
trips
net_betweenness(trips)                    # a to c is quicker through b
net_betweenness(trips, weighted = FALSE)  # ignoring travel time
calls <- data.frame(
  source  = c("a", "b", "a"),
  target    = c("b", "c", "c"),
  calls = c(1, 1, 20)
)
phone <- net_frame(calls, directed = FALSE, weight = "calls", weight_type = "strength")
net_betweenness(phone)  # a and c talk directly, so b brokers nothing
net_strength(phone)
op <- par(mfrow = c(1, 2))
xy <- plot(trips, main = "minutes (distance)", edge_color = I(c("black", "black", "gray")),
           seed = 2)
plot(phone, layout = xy, main = "calls (strength)",
     edge_color = I(c("gray", "gray", "black")))
par(op)
women <- net_example("southern_women")
women
op <- par(mfrow = c(1, 2))
plot(women, node_label = FALSE, main = "two-mode: women and events")
plot(net_example("krackhardt"), node_label = FALSE, main = "multiplex: three relations",
     seed = 1)
par(op)
m <- matrix(c(0, 2, 0,
              2, 0, 1,
              0, 1, 0), 3, byrow = TRUE,
            dimnames = list(c("x", "y", "z"), c("x", "y", "z")))
small <- as_net_frame(m)
small
small$edges
as.matrix(small, attr = "weight")
if (has_igraph) {
ig <- net_to_igraph(flo)
igraph::vcount(ig)
head(igraph::vertex_attr(ig, "betweenness"))
karate <- as_net_frame(igraph::make_graph("Zachary"))
karate
}
if (has_network) {
net <- net_to_network(flo)
network::network.size(net)
}
path <- tempfile(fileext = ".graphml")
net_write(flo, path)
back <- net_read(path)
back
head(back$nodes)
edges_csv <- tempfile(fileext = ".csv")
nodes_csv <- tempfile(fileext = ".csv")
net_write(flo, edges_csv, nodes_file = nodes_csv)
net_read(edges_csv, nodes = nodes_csv, directed = FALSE)
  }))
  succeed()
})

test_that("article network-types.Rmd runs", {
  skip_on_cran()
  on_null_device(local({

library(knappnet)
net_examples()[, c("name", "type", "source")]
friends <- net_example("krackhardt", relation = "friendship")
head(friends$edges)
friends
undirected <- net_frame(friends$edges, friends$nodes, directed = FALSE, source = "from",
                        target = "to", id = "name", quiet = TRUE)
undirected
op <- par(mfrow = c(1, 2))
xy <- plot(friends, main = "directed", node_label = FALSE, seed = 1)
plot(undirected, layout = xy, main = "undirected", node_label = FALSE)
par(op)
women <- net_example("southern_women")
head(women$nodes)
head(women$edges)
women
op <- par(mfrow = c(1, 2))
plot(women, node_label = FALSE, main = "women and events")
plot(net_project(women), node_label = FALSE, main = "women, by shared events", seed = 1)
par(op)
cb <- net_example("country_borders")
africa <- cb$nodes[cb$nodes$continent == "Africa", ]
regions <- sort(unique(africa$subregion))
nodes <- data.frame(id = c(africa$code, regions, "Africa"),
                    mode = c(rep("country", nrow(africa)), rep("subregion", length(regions)),
                             "continent"))
edges <- data.frame(source = c(africa$code, regions),
                    target = c(africa$subregion, rep("Africa", length(regions))))
three <- net_frame(edges, nodes, directed = FALSE, mode = "mode", between_modes_only = TRUE)
three
plot(three, layout = "modes", node_label = FALSE, seed = 1)
w <- data.frame(source = c("a", "b", "a"), target = c("b", "c", "c"), w = c(1, 1, 5))
as_strength <- net_frame(w, directed = FALSE, weight = "w", weight_type = "strength")
as_distance <- net_frame(w, directed = FALSE, weight = "w", weight_type = "distance")
as_strength
net_distances(as_strength)
net_distances(as_distance)
op <- par(mfrow = c(1, 2))
xy <- plot(as_strength, main = "strength: a-c is strong", seed = 3)
plot(as_distance, layout = xy, main = "distance: a-c is far",
     edge_color = I(c("black", "black", "gray80")))
par(op)
mis <- net_example("les_miserables", weighted = TRUE)
mis
plot(mis, node_label = FALSE, seed = 1)
tribes <- net_example("highland_tribes")
head(tribes$edges)
tribes
plot(tribes, seed = 1)
krack <- net_example("krackhardt")
table(krack$edges$relation)
krack
plot(krack, node_label = FALSE, seed = 1)
xy <- net_layout(krack, seed = 1)
relations <- net_split_relations(krack)
op <- par(mfrow = c(1, 3))
for (r in names(relations)) plot(relations[[r]], layout = xy, main = r, node_label = FALSE)
par(op)
edges <- data.frame(source = c("ann", "ann", "bob", "ann", "ann"),
                    target = c("bob", "bob", "cy", "ann", "cy"),
                    from_layer = c("work", "online", "online", "work", "work"),
                    to_layer = c("work", "online", "online", "online", "online"))
layers <- data.frame(layer = c("work", "online"), medium = c("office", "chat"))
ml <- net_frame(edges, directed = FALSE, source_layer = "from_layer",
                target_layer = "to_layer", layers = layers)
ml
cbind(edges, net_edge_layers(ml)["class"])
krack_layers <- net_as_multilayer(krack)
krack_layers
net_plot_layers(krack_layers, node_label = FALSE, seed = 1)
calls <- data.frame(source = c("a", "b", "a", "c", "b", "d", "a"),
                    target = c("b", "c", "c", "d", "d", "a", "d"),
                    time = c(1, 2, 2, 4, 5, 7, 8))
contacts <- net_frame(calls, directed = FALSE, start = "time")
contacts
net_plot_timeline(contacts)
net_plot_snapshots(contacts, width = 3, seed = 1)
spells <- data.frame(source = c("ann", "ann", "bob", "cy"),
                     target = c("bob", "bob", "cy", "dee"),
                     from = c(0, 12, 3, 6), to = c(10, 15, 9, 14))
people <- data.frame(id = c("ann", "bob", "cy", "dee"), joined = c(0, 0, 2, 5),
                     left = c(NA, NA, NA, 16))
stream <- net_frame(spells, people, directed = FALSE, start = "from", end = "to",
                    node_start = "joined", node_end = "left")
stream
net_plot_snapshots(stream, breaks = c(0, 4, 8, 12, 16), seed = 1)
liking <- net_example("sampson", relation = "liking")
liking
net_plot_snapshots(liking, node_label = FALSE, seed = 1)
waves <- net_snapshots(contacts, breaks = c(1, 5, 9))
panel <- net_panel(waves)
panel
sapply(net_waves(panel), function(g) nrow(g$edges))
miles <- net_example("miles")
ne <- miles$nodes[miles$nodes$state %in% c("ME", "NH", "VT", "MA", "RI", "CT", "NY", "NJ",
                                           "PA"), ]
near <- miles$edges[miles$edges$from %in% ne$name & miles$edges$to %in% ne$name &
                      miles$edges$miles <= 120, ]
roads <- net_frame(near, ne, directed = FALSE, source = "from", target = "to", id = "name",
                   weight = "miles", weight_type = "distance",
                   coords = c("longitude", "latitude"), quiet = TRUE)
roads
plot(roads, node_label = FALSE)
cb
plot(cb, node_label = FALSE)
banded <- near
banded$band <- banded$miles
st <- net_frame(banded, ne, directed = FALSE, source = "from", target = "to", id = "name",
                start = "band", coords = c("longitude", "latitude"), quiet = TRUE)
class(st)
st
net_plot_snapshots(st, breaks = c(0, 60, 90, 121), layout = c("longitude", "latitude"),
                   node_label = FALSE)
advice <- net_example("krackhardt", relation = "advice")
egos <- net_ego_frame(advice, egos = c(2, 7, 21), threads = 1)
egos
head(egos$nodes)
net_plot_egos(egos, node_label = FALSE, seed = 1)
fr_nodes <- friends$nodes
fr_nodes$no_reply <- fr_nodes$name %in% c(7, 15)
fr_edges <- friends$edges[!friends$edges$from %in% c(7, 15), ]
fr_edges$unknown <- FALSE
unasked <- data.frame(from = c(1, 3), to = c(5, 9), relation = "friendship",
                      unknown = TRUE)
fr_edges <- rbind(fr_edges, unasked)
partial <- net_frame(fr_edges, fr_nodes, source = "from", target = "to", id = "name",
                     node_missing = "no_reply", edge_missing = "unknown", quiet = TRUE)
partial
plot(partial, node_color = "no_reply", node_label = FALSE, seed = 1)
net_density(net_resolve_missing(partial, "absent"))
multi <- data.frame(source = c("a", "a", "a", "b", "c", "c"),
                    target = c("b", "b", "a", "c", "a", "a"),
                    calls = c(1, 2, 1, 4, 1, 1))
g <- net_frame(multi, weight = "calls", weight_type = "strength")
g
net_multiplicity(g)
simple <- net_simplify(g, remove_loops = TRUE)
simple
op <- par(mfrow = c(1, 2))
xy <- plot(g, main = "multigraph", seed = 2)
plot(simple, layout = xy, main = "simplified (weights summed)")
par(op)
try({
net_frame(multi, loops = FALSE, multiple = FALSE)
}, silent = TRUE)
net_example("sampson")
  }))
  succeed()
})

test_that("article definitions.Rmd runs", {
  skip_on_cran()
  on_null_device(local({

library(knappnet)
nodes <- data.frame(name = c("ana", "ben", "cy", "dee"),
                    age  = c(31, 45, 28, 52))
edges <- data.frame(source = c("ana", "ana", "ben", "cy"),
                    target = c("ben", "cy",  "cy",  "dee"),
                    years  = c(4, 2, 7, 1))
edges
g <- net_frame(edges, nodes, directed = FALSE,
               source = "source", target = "target", id = "name",
               weight = "years", weight_type = "strength")
g
plot(g, seed = 1)
try({
net_frame(data.frame(source = c("ana", NA, "ben"), target = c("ben", "cy", "zed")), nodes,
          id = "name")
}, silent = TRUE)
multiplex <- data.frame(source   = c("ana", "ana", "ana", "ben"),
                        target   = c("ben", "ben", "ben", "cy"),
                        relation = c("friendship", "advice", "advice", "advice"))
multiplex
m <- net_frame(multiplex, relation = "relation")
net_multiplicity(m)
net_is_simple(m)
op <- par(mfrow = c(1, 2))
xy <- plot(m, main = "multigraph", seed = 1)
plot(net_simplify(m), layout = xy, main = "simplified")
par(op)
try({
net_frame(edges, weight = "years")                           # error: say what the weights mean
net_frame(edges, weight = "years", weight_type = "strength") # years of acquaintance
}, silent = TRUE)
as_distance <- net_frame(edges, weight = "years", weight_type = "distance", directed = FALSE)
net_distances(g)["ana", ]
net_distances(as_distance)["ana", ]
affiliations <- data.frame(source = c("ana", "ana", "ben", "cy"),
                           target = c("club A", "club B", "club A", "club B"))
members <- data.frame(name = c("ana", "ben", "cy", "club A", "club B"),
                      mode = c("person", "person", "person", "club", "club"))
members
clubs <- net_twomode(affiliations, members, mode = "mode", id = "name")
clubs
c(declared_two_mode = net_is_twomode(clubs), bipartite = net_is_bipartite(clubs))
plot(clubs)
path <- net_frame(data.frame(source = c("a", "b"), target = c("b", "c")), directed = FALSE)
c(declared_two_mode = net_is_twomode(path), bipartite = net_is_bipartite(path))
ties <- data.frame(source = c("ana", "ana", "ana"), target = c("ben", "ben", "ana"),
                   from_layer = c("work", "online", "work"),
                   to_layer = c("work", "online", "online"))
ml <- net_frame(ties, directed = FALSE, source_layer = "from_layer",
                target_layer = "to_layer")
ml
net_edge_layers(ml)
ratings <- data.frame(source = c("ana", "ana", "ben"), target = c("ben", "cy", "cy"),
                      value = c(2, -1, 3))
signed <- net_frame(net_as_signed(ratings, "value"), directed = FALSE, sign = "sign",
                    weight = "weight", weight_type = "strength")
signed
plot(signed, seed = 1)
contacts <- data.frame(source = c("ana", "ana", "ben"),
                       target = c("ben", "cy",  "cy"),
                       start  = c(1, 1, 3),
                       end    = c(3, 2, 6))
contacts   # active on [1, 3), [1, 2), [3, 6)
stream <- net_frame(contacts, directed = FALSE, start = "start", end = "end")
stream
net_plot_timeline(stream)
waves <- list(net_frame(data.frame(source = c("ana", "ben"), target = c("ben", "cy"))),
              net_frame(data.frame(source = c("ana", "cy"), target = c("cy", "ben"))))
panel <- net_panel(waves)
panel
net_plot_snapshots(panel, seed = 1)
places <- data.frame(id = c("Baltimore", "Washington", "Philadelphia"),
                     lon = c(-76.61, -77.04, -75.17), lat = c(39.29, 38.91, 39.95))
trains <- data.frame(source = c("Washington", "Baltimore"),
                     target = c("Baltimore", "Philadelphia"), minutes = c(35, 75))
rail <- net_frame(trains, places, directed = FALSE, weight = "minutes",
                  weight_type = "distance", coords = c("lon", "lat"))
rail
plot(rail)
whole <- net_frame(data.frame(source = c("ana", "ana", "ben", "cy", "dee"),
                              target = c("ben", "cy", "cy", "dee", "eli")),
                   directed = FALSE)
egos <- net_ego_frame(whole, egos = c("ana", "dee"), threads = 1)
egos
net_plot_egos(egos, seed = 1)
survey <- data.frame(source = c("ana", "ana", "ben"), target = c("ben", "cy", "cy"),
                     unknown = c(FALSE, TRUE, FALSE))
partial <- net_frame(survey, directed = FALSE, edge_missing = "unknown")
partial
op <- par(mfrow = c(1, 2))
xy <- plot(partial, main = "unobserved pair dotted", seed = 1)
plot(net_resolve_missing(partial, "absent"), layout = xy, main = "resolved as absent")
par(op)
hyperedges <- data.frame(source = c("ana", "ben", "cy", "ana", "dee"),
                         target = c("h1",  "h1",  "h1", "h2",  "h2"))
hyperedges   # h1 = {ana, ben, cy}, h2 = {ana, dee}
hyper <- net_twomode(hyperedges)
net_project(hyper)$edges   # the clique expansion, weighted by shared hyperedges
  }))
  succeed()
})

test_that("article coming-from-igraph.Rmd runs", {
  skip_on_cran()
  skip_if_not_installed("igraph")
  on_null_device(local({

library(knappnet)
ig <- igraph::make_graph("Zachary")
karate <- as_net_frame(ig)
karate
head(karate$edges)
ig_scores <- data.frame(
  degree      = igraph::degree(ig),
  betweenness = igraph::betweenness(ig),
  closeness   = igraph::closeness(ig),
  eigenvector = igraph::eigen_centrality(ig)$vector,
  core        = igraph::coreness(ig),
  constraint  = igraph::constraint(ig)
)
karate$nodes$degree      <- net_degree(karate)
karate$nodes$betweenness <- net_betweenness(karate)
karate$nodes$closeness   <- net_closeness(karate)
karate$nodes$eigenvector <- net_eigen(karate)
karate$nodes$core        <- net_coreness(karate)
karate$nodes$constraint  <- net_constraint(karate)
head(karate$nodes)
sapply(names(ig_scores), function(m) {
  isTRUE(all.equal(karate$nodes[[m]], ig_scores[[m]], check.attributes = FALSE,
                   tolerance = 1e-7))
})
c(igraph = igraph::edge_density(ig), knappnet = net_density(karate))
c(igraph = igraph::centr_betw(ig)$centralization,
  knappnet = net_centralization(karate, "betweenness"))
c(igraph = igraph::clique_num(ig), knappnet = net_clique_number(karate))
set.seed(1)
ig_groups <- igraph::cluster_louvain(ig)
kn_groups <- net_louvain(karate, seed = 1)
c(igraph = igraph::modularity(ig_groups), knappnet = attr(kn_groups, "modularity"))
igraph::E(ig)$weight <- rep(c(1, 2, 3), length.out = igraph::ecount(ig))
weighted_karate <- as_net_frame(ig)
weighted_karate
all.equal(net_betweenness(weighted_karate), igraph::betweenness(ig))
back <- net_to_igraph(karate)
head(igraph::vertex_attr(back, "betweenness"))
  }))
  succeed()
})

test_that("article coming-from-networkx.Rmd runs", {
  skip_on_cran()
  on_null_device(local({

library(knappnet)
edges <- data.frame(source = c("ana", "ana", "ben", "cy", "cy"),
                    target   = c("ben", "cy", "cy", "dee", "eli"))
g <- net_frame(edges, directed = FALSE)
net_betweenness(g)
split <- net_frame(data.frame(source = c("a", "b", "d"), target = c("b", "c", "e")),
                  directed = FALSE)
data.frame(node     = split$nodes$id,
           networkx = round(net_closeness(split, convention = "networkx"), 3),
           igraph   = round(net_closeness(split, normalized = TRUE), 3))
multi <- net_frame(data.frame(source = c("a", "a", "b", "a", "d"),
                             target   = c("b", "b", "c", "d", "c")), directed = FALSE)
data.frame(node     = multi$nodes$id,
           igraph   = round(net_betweenness(multi), 3),
           networkx = round(net_betweenness(multi, convention = "networkx"), 3))
kite <- net_frame(data.frame(source = c(1, 1, 1, 2, 2, 3, 4),
                            target   = c(2, 3, 4, 3, 4, 4, 5)), directed = FALSE)
data.frame(node     = kite$nodes$id,
           igraph   = round(net_eigen(kite), 3),
           networkx = round(net_eigen(kite, convention = "networkx"), 3))
  }))
  succeed()
})

test_that("article coming-from-statnet.Rmd runs", {
  skip_on_cran()
  skip_if_not_installed("network")
  skip_if_not_installed("sna")
  on_null_device(local({

library(knappnet)
# Padgett's Florentine marriage ties, with the Pucci family added as an
# isolate (it married outside this group); see ?net_example.
ties <- net_example("florentine")$edges
families <- data.frame(name = sort(unique(c(ties$source, ties$target, "Pucci"))))
flo <- net_frame(ties, nodes = families, directed = FALSE, id = "name")

net <- net_to_network(flo)
network::network.size(net)
network::network.edgecount(net)
as_net_frame(net)

m <- as.matrix(flo)   # a symmetric 0/1 sociomatrix with node names
m[1:4, 1:4]
multi <- net_frame(data.frame(source = c("a", "a", "b", "c"),
                             target   = c("b", "b", "c", "c")), directed = FALSE)
data.frame(node = multi$nodes$id,
           default = net_degree(multi),
           sna = net_degree(multi, convention = "sna"))
sna_closeness <- sna::closeness(m, gmode = "graph", cmode = "undirected")
data.frame(family = flo$nodes$name,
           sna = sna_closeness,
           knappnet_sna = net_closeness(flo, convention = "sna"),
           knappnet_default = round(net_closeness(flo, normalized = TRUE), 3))
star <- net_frame(data.frame(source = "hub", target = c("a", "b", "c", "d", "e")),
                 directed = FALSE)
c(igraph = net_centralization(star, "degree"),
  sna = net_centralization(star, "degree", convention = "sna"),
  sna_itself = sna::centralization(as.matrix(star), sna::degree, mode = "graph",
                                   cmode = "indegree"))
sna_scores <- data.frame(
  degree      = sna::degree(m, gmode = "graph", cmode = "indegree"),
  betweenness = sna::betweenness(m, gmode = "graph", cmode = "undirected"),
  harmonic    = sna::closeness(m, gmode = "graph", cmode = "suminvundir")
)
kn_scores <- data.frame(
  degree      = net_degree(flo, convention = "sna"),
  betweenness = net_betweenness(flo, convention = "sna"),
  harmonic    = net_harmonic(flo, convention = "sna")
)
sapply(names(sna_scores), function(v) {
  isTRUE(all.equal(kn_scores[[v]], sna_scores[[v]], check.attributes = FALSE))
})
side <- ifelse(flo$nodes$name < "M", 1, 2)   # A-L and M-Z
kn_brokerage <- net_brokerage(flo, side)
sna_brokerage <- sna::brokerage(m, side)$raw.nli[, c("w_I", "w_O", "b_IO", "b_OI", "b_O", "t")]
all.equal(unname(as.matrix(kn_brokerage[, -1])), unname(sna_brokerage))
all.equal(unname(net_struct_equiv(flo, "hamming")),
          unname(as.matrix(sna::sedist(m, method = "hamming"))))
staff <- data.frame(
  id = c("ann", "bob", "cat", "dan", "eve", "fay"),
  dept = c("sales", "sales", "eng", "eng", "eng", "legal")
)
advice <- net_frame(
  data.frame(source = c("ann", "bob", "cat", "dan", "eve", "bob", "fay", "ann"),
             target   = c("bob", "cat", "dan", "eve", "cat", "fay", "dan", "cat")),
  nodes = staff
)
friends <- net_frame(
  data.frame(source = c("ann", "ann", "bob", "cat", "dan"),
             target   = c("bob", "eve", "eve", "fay", "eve")),
  nodes = staff, directed = FALSE
)
same_dept <- outer(staff$dept, staff$dept, "==") * 1

Y <- as.matrix(advice)
X <- as.matrix(friends)
test <- net_qap_cor(advice, friends, exact = TRUE)
c(knappnet = test$statistic, sna = sna::gcor(Y, X))

fit <- net_mrqap(advice, list(friends = friends, same_dept = same_dept), exact = TRUE)
set.seed(1)
sna_fit <- sna::netlm(Y, list(X, same_dept), reps = 1000)
rbind(knappnet = fit$coefficients, sna = sna_fit$coefficients)
set.seed(1)
sna_test <- sna::qaptest(list(Y, X), sna::gcor, g1 = 1, g2 = 2, reps = 1000)
rbind(knappnet_exact = test$p.value[c("Pr(<=r)", "Pr(>=r)")],
      sna_1000_draws = c(sna_test$pleeq, sna_test$pgreq))

rbind(knappnet_exact = fit$p.value[, "Pr(>=|b|)"],
      sna_1000_draws = sna_fit$pgreqabs)
  }))
  succeed()
})

test_that("article coming-from-ucinet-ora-pajek-gephi.Rmd runs", {
  skip_on_cran()
  skip_if_not_installed("igraph")
  on_null_device(local({

library(knappnet)
dl <- tempfile(fileext = ".dl")
writeLines(c(
  "dl n = 5",
  "format = fullmatrix",
  "labels:",
  "ann,bob,cat,dan,eve",
  "data:",
  "0 1 1 0 0",
  "1 0 1 0 0",
  "1 1 0 1 0",
  "0 0 1 0 1",
  "0 0 0 1 0"
), dl)
team <- net_read(dl, directed = FALSE)
team
team$edges
out <- tempfile(fileext = ".dl")
net_write(team, out)
cat(readLines(out), sep = "\n")
org <- net_frame(
  data.frame(source    = c("ann", "bob", "ann", "bob", "cat"),
             target      = c("bob", "cat", "sales", "sales", "eng"),
             network = c("advice", "advice", "member", "member", "member")),
  nodes = data.frame(id    = c("ann", "bob", "cat", "sales", "eng"),
                     class = c("people", "people", "people", "units", "units"),
                     type  = c("Agent", "Agent", "Agent", "Organization", "Organization"))
)
xml <- tempfile(fileext = ".xml")
net_write(org, xml, format = "dynetml")
meta <- net_read(xml)
meta$nodes
meta$edges
advice <- net_frame(meta$edges[meta$edges$network == "advice", ],
                   nodes = meta$nodes[meta$nodes$class == "people", ],
                   directed = meta$directed)
advice
net_betweenness(advice)
member <- net_frame(meta$edges[meta$edges$network == "member", ], nodes = meta$nodes,
                    directed = meta$directed)
member <- net_as_twomode(member, mode = "class")
member
net_incidence(member)
net <- tempfile(fileext = ".net")
writeLines(c(
  "*Vertices 4",
  '1 "ann"',
  '2 "bob"',
  '3 "cat"',
  '4 "dan"',
  "*Edges",
  "1 2",
  "1 3",
  "2 3",
  "3 4"
), net)
net_read(net)
team$nodes$betweenness <- net_betweenness(team)
team$nodes$group <- as.vector(net_louvain(team, seed = 1))
gexf <- tempfile(fileext = ".gexf")
net_write(team, gexf)
net_read(gexf)$nodes
split <- net_frame(data.frame(source = c("a", "b", "d"), target = c("b", "c", "e")),
                  directed = FALSE)
data.frame(node   = split$nodes$id,
           ucinet = round(net_closeness(split, convention = "ucinet"), 3),
           igraph = round(net_closeness(split, normalized = TRUE), 3))
karate <- as_net_frame(igraph::make_graph("Zachary"))
n <- nrow(karate$nodes)
ucinet_style <- data.frame(
  node          = karate$nodes$id,
  degree        = net_degree(karate, convention = "ucinet"),
  nDegree       = 100 * net_degree(karate, convention = "ucinet") / (n - 1),
  nBetweenness  = 100 * net_betweenness(karate, normalized = TRUE, convention = "ucinet"),
  nCloseness    = 100 * net_closeness(karate, convention = "ucinet"),
  eigenvector   = net_eigen(karate, convention = "ucinet")
)
head(round(ucinet_style, 3))
100 * net_centralization(karate, "degree", convention = "ucinet")
100 * net_centralization(karate, "betweenness", convention = "ucinet")
star <- net_frame(data.frame(source = "hub", target = c("a", "b", "c", "d", "e")),
                 directed = FALSE)
c(igraph = net_centralization(star, "degree"),
  ucinet = net_centralization(star, "degree", convention = "ucinet"))
net_fragmentation(split)
  }))
  succeed()
})

test_that("article conventions.Rmd runs", {
  skip_on_cran()
  on_null_device(local({

library(knappnet)
multi <- net_frame(data.frame(source = c("a", "a", "b", "c"),
                             target   = c("b", "b", "c", "c")), directed = FALSE)
multi$edges
data.frame(node   = multi$nodes$id,
           igraph = net_degree(multi),
           sna    = net_degree(multi, convention = "sna"))
paths <- net_frame(data.frame(source = c("a", "a", "b", "a", "d"),
                             target   = c("b", "b", "c", "d", "c")), directed = FALSE)
data.frame(node     = paths$nodes$id,
           igraph   = round(net_betweenness(paths), 3),
           networkx = round(net_betweenness(paths, convention = "networkx"), 3),
           sna      = round(net_betweenness(paths, convention = "sna"), 3))
g <- net_frame(data.frame(source = c("a", "a", "a", "b"), target = c("b", "c", "d", "c")),
              directed = FALSE)
data.frame(node             = g$nodes$id,
           betweenness      = net_betweenness(g),
           normalized       = round(net_betweenness(g, normalized = TRUE), 3),
           ucinet_percent   = round(100 * net_betweenness(g, normalized = TRUE,
                                                      convention = "ucinet"), 1),
           sna_rescale      = round(net_betweenness(g) / sum(net_betweenness(g)), 3))
split <- net_frame(data.frame(source = c("a", "b", "d"), target = c("b", "c", "e")),
                  directed = FALSE)
data.frame(
  node     = split$nodes$id,
  igraph   = round(net_closeness(split), 3),
  networkx = round(net_closeness(split, convention = "networkx"), 3),
  sna      = net_closeness(split, convention = "sna"),
  ucinet   = round(net_closeness(split, convention = "ucinet"), 3)
)
data.frame(node       = split$nodes$id,
           igraph     = net_harmonic(split),
           networkx   = net_harmonic(split, convention = "networkx"),
           sna_ucinet = net_harmonic(split, convention = "sna"))
d <- net_frame(data.frame(source = c("a", "b", "c", "b", "c"),
                         target   = c("b", "c", "a", "a", "b")))
data.frame(node     = d$nodes$id,
           igraph   = round(net_eigen(d), 3),
           networkx = round(net_eigen(d, convention = "networkx"), 3),
           sna      = round(net_eigen(d, convention = "sna"), 3),
           ucinet   = round(net_eigen(d, convention = "ucinet"), 3))
try({
net_eigen(split, convention = "networkx")
}, silent = TRUE)
star <- net_frame(data.frame(source = "hub", target = c("a", "b", "c", "d", "e")),
                 directed = FALSE)
c(igraph          = net_centralization(star, "degree"),
  igraph_no_loops = net_centralization(star, "degree", loops = FALSE),
  sna             = net_centralization(star, "degree", convention = "sna"),
  ucinet          = net_centralization(star, "degree", convention = "ucinet"))
c(igraph = net_density(multi), sna = net_density(multi, convention = "sna"))
  }))
  succeed()
})

test_that("article plotting.Rmd runs", {
  skip_on_cran()
  on_null_device(local({

library(knappnet)
karate <- net_example("karate")
florentine <- net_example("florentine")
plot(karate, node_color = "faction", seed = 1, main = "Zachary's karate club")
net_plot_centrality(karate, net_betweenness(karate), title = "betweenness", seed = 1)
net_plot_partition(karate, net_louvain(karate, seed = 1), title = "community", seed = 1)
net_plot_scores(net_betweenness(florentine), florentine, main = "Betweenness, Florentine families")
net_plot_matrix(karate, order = net_louvain(karate, seed = 1), main = "Sorted by community")
bm <- net_blockmodel(florentine, net_concor(florentine, depth = 2))
plot(bm, network = florentine)
net_plot_mixing(karate, "faction", main = "Ties within and between factions")
net_assort_nominal(karate, "faction")
net_plot_degree(karate)
net_plot_degree(karate, type = "ccdf", log = TRUE)
friends <- net_example("noordin", relation = "friendship")
net_plot_components(friends)
advice <- net_frame(data.frame(source = c(1, 2, 3, 3, 4, 5, 5, 6, 1, 2),
                              target   = c(2, 3, 1, 4, 5, 3, 6, 4, 6, 1)))
net_plot_triads(advice, threads = 2)
net_plot_coreness(karate, seed = 1, threads = 2, main = "k-core shells")
net_plot_ego(florentine, "Medici", seed = 1, main = "The Medici and their neighbours")
ji <- net_example("jemaah_islamiyah")
kp <- net_key_players(ji, 3, "fragmentation", seed = 1, threads = 2)
kp$set
net_plot_robustness(ji, nodes = kp$set, label = "key players", seed = 1, threads = 2)
d <- net_mds_distances(florentine)
fit <- net_mds(d, k = 2)
fit
plot(fit, main = "Geodesic distances, metric scaling")
plot(fit, which = "shepard", d = d)
set.seed(1)
m <- as.matrix(florentine)
noise <- matrix(rbinom(length(m), 1, 0.15), nrow(m))
noisy <- pmax(m, (noise + t(noise)) > 0)
q <- net_qap_cor(m, noisy, mode = "graph", reps = 1000, seed = 1, threads = 2)
plot(q)
set.seed(2)
x1 <- matrix(rnorm(144), 12)
x2 <- matrix(rnorm(144), 12)
y <- 0.8 * x1 + matrix(rnorm(144), 12)
fit <- net_mrqap(y, list(same = x1, other = x2), reps = 1000, seed = 1, threads = 2)
plot(fit)
plot(fit, which = "distributions")
logins <- data.frame(user = c("u1", "u2", "u1", "u3", "u2", "u3", "u1", "u2"),
                     host = c("h1", "h1", "h2", "h2", "h3", "h1", "h1", "h2"),
                     time = c(0, 1, 5, 6, 12, 13, 14, 15))
g <- net_temporal(logins, source = "user", target = "host")
net_plot_timeline(g, main = "Logins")
net_plot_rolling(net_rolling(g, net_degree, width = 6, step = 3), main = "Degree in sliding windows")
  }))
  succeed()
})

test_that("article spatial-networks.Rmd runs", {
  skip_on_cran()
  skip_if_not_installed("sf")
  on_null_device(local({
has_sf <- requireNamespace("sf", quietly = TRUE)
has_sfnetworks <- has_sf && requireNamespace("sfnetworks", quietly = TRUE)

library(knappnet)
stations <- sf::st_as_sf(
  data.frame(id = c("north", "east", "south", "west"),
             x = c(500000, 503000, 500000, 497000),
             y = c(6003000, 6000000, 5997000, 6000000),
             elevation = c(120, 85, 90, 140)),
  coords = c("x", "y"), crs = 32633)
links <- data.frame(source = c("north", "east", "south"),
                    target = c("east", "south", "west"))
g <- net_frame(links, stations, directed = FALSE)
g
g$roles[c("id", "node_geometry", "edge_geometry")]
g$roles$crs$Name
g$nodes
edge_lines <- net_to_sf(g, "edges")
edge_lines
h <- net_frame(edge_lines, stations, directed = FALSE)
h
net_degree(h)
net_betweenness(h)
net_edge_lengths(g)
net_edge_lengths(h)
bend <- sf::st_linestring(rbind(c(500000, 6003000), c(503000, 6003000), c(503000, 6000000)))
sf::st_geometry(edge_lines)[1] <- sf::st_sfc(bend, crs = 32633)
sf::st_geometry(edge_lines)[3] <- sf::st_sfc(sf::st_linestring(), crs = 32633)
mixed <- net_frame(edge_lines, stations, directed = FALSE)
mixed
net_edge_lengths(mixed)
mixed$edges$metres <- as.vector(net_edge_lengths(mixed))
w <- net_frame(mixed$edges, stations, directed = FALSE, weight = "metres", weight_type = "distance")
net_distances(w)
square <- function(x, y) {
  sf::st_polygon(list(rbind(c(x, y), c(x + 1, y), c(x + 1, y + 1), c(x, y + 1), c(x, y))))
}
regions <- sf::st_sf(region = c("sw", "se", "nw", "ne"),
                     row = c("bottom", "bottom", "top", "top"),
                     geometry = sf::st_sfc(square(0, 0), square(1, 0), square(0, 1), square(1, 1),
                                           crs = 32633))
r <- net_from_polygons(regions, id = "region")
net_group_geometry(r, "row")
net_from_polygons(regions, id = "region")$edges
net_from_polygons(regions, contiguity = "rook", id = "region")$edges
segment <- function(x0, y0, x1, y1) sf::st_linestring(rbind(c(x0, y0), c(x1, y1)))
streets <- sf::st_sf(
  name = c("1st Ave", "2nd Ave", "3rd Ave", "A St", "B St", "C St"),
  geometry = sf::st_sfc(segment(0, 0, 200, 0), segment(0, 100, 200, 100), segment(0, 200, 200, 200),
                        segment(0, 0, 0, 200), segment(100, 0, 100, 200), segment(200, 0, 200, 200),
                        crs = 32633))
roads <- net_from_lines(streets)
roads
roads$edges
net_distances(roads)[1, ]
net_proximity(regions, within = 1, id = "region")$edges
net_proximity(regions, within = 0, id = "region", exact = TRUE)$edges
inner <- sf::st_bbox(c(xmin = 498000, ymin = 5998000, xmax = 504000, ymax = 6004000),
                     crs = sf::st_crs(32633))
net_spatial_filter(h, inner)$nodes$id

zones <- sf::st_sf(zone = c("uptown", "downtown"),
                   geometry = sf::st_sfc(
                     sf::st_polygon(list(rbind(c(495000, 6001500), c(505000, 6001500),
                                               c(505000, 6005000), c(495000, 6005000),
                                               c(495000, 6001500)))),
                     sf::st_polygon(list(rbind(c(495000, 5995000), c(505000, 5995000),
                                               c(505000, 6001500), c(495000, 6001500),
                                               c(495000, 5995000))))),
                   crs = 32633)
net_spatial_join(h, zones)$nodes[, c("id", "zone")]

ll <- net_transform(h, 4326)
ll$roles$crs$Name
net_edge_lengths(ll)
net_check_geometry(h)
net_check_geometry(mixed)
net_plot_spatial(r, node_color = "row", main = "Queen contiguity of four regions")
net_plot_spatial(roads, node_label = FALSE, edge_color = "length", main = "Six streets")
net_plot_spatial(mixed, edge_width = 2, main = "Stations, with one bent and one implicit edge")
net_plot_edge_lengths(mixed, scale = 3000)
path <- tempfile(fileext = ".gpkg")
net_write(roads, path)
sf::st_layers(path)$name
back <- net_read(path)
identical(net_distances(back), net_distances(roads))
unlink(path)
if (has_sfnetworks) {
net <- net_to_sfnetwork(roads)
net
as_net_frame(net)
}
  }))
  succeed()
})
