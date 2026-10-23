# Contact timelines and measures over windows (R/plot-temporal.R).

test_that("temporal plots draw timelines and measures over time", {
  calls <- data.frame(source = c("a", "a", "b", "c", "a", "b"),
                      target   = c("b", "c", "c", "d", "b", "c"),
                      time = c(1, 2, 2, 5, 7, 9), end = c(3, 2, 4, 6, 7, 10))
  g <- net_temporal(calls, end = "end", directed = FALSE)
  on_null_device({
    expect_silent(out <- net_plot_timeline(g))
    expect_identical(out$row[1:2], c("a -- b", "b -- c"))
    expect_identical(out$contacts[1:2], c(2L, 2L))
    expect_identical(sum(out$contacts), 6L)
    expect_silent(out <- net_plot_timeline(g, by = "node", max_rows = 2, main = "Nodes"))
    expect_identical(nrow(out), 2L)
    expect_identical(out$row[1], "b")
    expect_restores_par(net_plot_timeline(g))
    directed <- net_temporal(calls, end = "end")
    expect_silent(out <- net_plot_timeline(directed))
    expect_identical(out$row[1], "a -> b")
    instant <- net_temporal(data.frame(source = "a", target = "a", time = 3))
    expect_silent(net_plot_timeline(instant))
    dated <- net_temporal(data.frame(source = c("a", "b"), target = c("b", "c"),
                                         time = as.Date(c("2024-01-01", "2024-03-01"))))
    expect_silent(net_plot_timeline(dated))
    stamped <- net_temporal(data.frame(source = c("a", "b"), target = c("b", "c"),
                                           time = as.POSIXct(c("2024-01-01 10:00", "2024-01-01 12:30"),
                                                             tz = "UTC")))
    expect_silent(net_plot_timeline(stamped))
    empty <- net_temporal(calls[0, ], nodes = data.frame(id = c("a", "b")), end = "end")
    expect_error(net_plot_timeline(empty), "no contacts")
    expect_error(net_plot_timeline(net_frame(calls)), "temporal network")
  })
  logins <- data.frame(user = c("u1", "u2", "u1", "u3", "u2", "u3", "u1"),
                       host = c("h1", "h1", "h2", "h2", "h3", "h1", "h1"),
                       time = c(0, 1, 5, 6, 12, 13, 14))
  h <- net_temporal(logins, source = "user", target = "host")
  r <- net_rolling(h, net_degree, width = 5, step = 2)
  on_null_device({
    expect_silent(out <- net_plot_rolling(r))
    expect_identical(out, r)
    expect_silent(net_plot_rolling(r, nodes = c("u1", "h1"), main = "Two nodes"))
    expect_silent(net_plot_rolling(r, top = 2))
    expect_silent(net_plot_rolling(net_rolling(h, function(w) nrow(w$edges), width = 5), ylab = "contacts"))
    expect_restores_par(net_plot_rolling(r))
    expect_error(net_plot_rolling(r, nodes = "zz"), "not in the results")
    expect_error(net_plot_rolling(data.frame(a = 1)), "net_rolling()")
    expect_error(net_plot_rolling(r[0, ]), "no windows")
    dated <- net_temporal(data.frame(source = c("a", "b", "a"), target = c("b", "c", "c"),
                                         time = as.Date(c("2024-01-01", "2024-01-05", "2024-01-09"))))
    expect_silent(net_plot_rolling(net_rolling(dated, net_degree, width = 4)))
  })
})
