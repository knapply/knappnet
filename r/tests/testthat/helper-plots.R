# Plots are checked by behaviour on a null device, not by pixels: they must
# draw without error or warning, return what they document, and restore
# the graphical parameters they change.

on_null_device <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  force(code)
}

expect_restores_par <- function(code) {
  before <- graphics::par(c("mar", "mfrow", "xpd"))
  force(code)
  expect_identical(graphics::par(c("mar", "mfrow", "xpd")), before)
}
