## -----------------------------------------------------------------------------
#| label: setup
#| include: false
if (!requireNamespace("fiber", quietly = TRUE)) {
  # R CMD check re-builds vignette outputs in a fresh R process where the
  # just-installed package is not always on the library path (seen on Windows).
  # Skip the executable chunks rather than fail the build in that case.
  knitr::opts_chunk$set(eval = FALSE)
} else {
  library(fiber)
}

# The dti package (a Suggests) pulls in rgl on load; in headless environments
# rgl cannot open a display and emits a warning ('rgl.init' failed / X11 error).
# That warning concerns rgl's graphics, not fiber, so we route rgl to its null
# device to keep the coercion example below clean.
Sys.setenv(RGL_USE_NULL = "true")


## -----------------------------------------------------------------------------
# A helix in 3-D, sampled at 50 points
t <- seq(0, 2 * pi, length.out = 50)
sl <- streamline(
  points = cbind(X = cos(t), Y = sin(t), Z = t / (2 * pi))
)
sl


## -----------------------------------------------------------------------------
sl <- streamline(
  points = cbind(X = cos(t), Y = sin(t), Z = t / (2 * pi)),
  point_data = list(FA = runif(50)),
  streamline_data = list(mean_FA = mean(runif(50)), label = "helix")
)
sl@n_points # 50
sl@point_attributes # "FA"
sl@streamline_attributes # c("mean_FA", "label")


## -----------------------------------------------------------------------------
get_curvilinear_length(sl)
get_sinuosity(sl)


## -----------------------------------------------------------------------------
head(get_curvature(sl))


## -----------------------------------------------------------------------------
sl2 <- add_shape_descriptors(sl)
sl2@streamline_data$euclidean_length
head(sl2@point_data$curvature)


## -----------------------------------------------------------------------------
sl_a <- streamline(points = cbind(X = cos(t), Y = sin(t), Z = t / (2 * pi)))
sl_b <- streamline(
  points = cbind(X = 1.1 * cos(t), Y = 1.1 * sin(t), Z = t / (2 * pi))
)
b <- bundle(streamlines = list(sl_a, sl_b))
b
b@n_streamlines # 2
b@streamline_data$mean_FA # lifted automatically, if present


## -----------------------------------------------------------------------------
b2 <- bind_bundles(sl_a, b)
b2@n_streamlines # 3


## -----------------------------------------------------------------------------
sl_20 <- reparametrize(sl_a, n_points = 20L)
sl_20@n_points # 20

b_20 <- reparametrize(b, n_points = 20L)
b_20[[1]]@n_points # 20
b_20[[2]]@n_points # 20


## -----------------------------------------------------------------------------
# streamlines -> scalar
compute_hausdorff_distance(sl_a, sl_b)

# bundle (with no second argument) -> pairwise dist object
d <- compute_hausdorff_distance(b)
as.matrix(d)

# bundle vs streamline -> vector
compute_hausdorff_distance(b, sl_a)


## -----------------------------------------------------------------------------
# Single subject: several white-matter bundles of one subject
bs_single <- bundle_set(list(
  callosum = bundle(list(sl_a)),
  arcuate = bundle(list(sl_b))
))
bs_single
bs_single@n_bundles # 2
bs_single[["callosum"]] # extract one bundle (push-down)

# Multi-subject: the same bundle across two subjects / sessions
b_sub01 <- bundle(
  streamlines = list(sl_a),
  bundle_data = list(subject = "sub-01")
)
b_sub02 <- bundle(
  streamlines = list(sl_b),
  bundle_data = list(subject = "sub-02")
)

bs <- bundle_set(list("sub-01" = b_sub01, "sub-02" = b_sub02))
bs
bs@bundle_data$subject # c("sub-01", "sub-02")
bs[["sub-01"]] # extract one subject (push-down)


## -----------------------------------------------------------------------------
cc <- read.csv(
  system.file("extdata", "CC_SingleTensor.csv", package = "fiber"),
  header = TRUE
)

make_streamline <- function(d) {
  d <- d[order(d$PointId), ]
  streamline(
    points = as.matrix(d[c("X", "Y", "Z")]),
    point_data = list(FA = d$FA)
  )
}
cc_bundle <- bundle(
  streamlines = lapply(split(cc, cc$StreamlineId), make_streamline)
)

cc_bundle@n_streamlines # 100
cc_bundle[[1]]@n_points # 150
cc_bundle[[1]]@point_attributes # "FA"
range(cc_bundle[[1]]@point_data$FA) # per-point FA on the first stream


## -----------------------------------------------------------------------------
cc_bundle2 <- bundle(
  streamlines = lapply(
    split(cc, cc$StreamlineId),
    function(d) {
      d <- d[order(d$PointId), ]
      streamline(
        points = as.matrix(d[c("X", "Y", "Z")]),
        point_data = list(FA = d$FA),
        streamline_data = list(mean_FA = mean(d$FA))
      )
    }
  )
)
cc_bundle2@streamline_data$mean_FA
mean(cc_bundle2@streamline_data$mean_FA) # mean across the 100 streamlines


## -----------------------------------------------------------------------------
if (requireNamespace("dti", quietly = TRUE)) {
  dfi <- as_dwifiber(b)
  class(dfi) # "dwiFiber"
  b_back <- as_bundle(dfi)
  b_back@n_streamlines
}

