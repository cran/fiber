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
# device to keep the coercion examples below clean.
Sys.setenv(RGL_USE_NULL = "true")


## -----------------------------------------------------------------------------
sl <- streamline(
  points = cbind(X = 0:4, Y = 0:4, Z = 0:4),
  point_data = list(FA = runif(5)),
  streamline_data = list(label = "CST")
)
sl
sl@n_points
sl@point_attributes


## -----------------------------------------------------------------------------
sl1 <- streamline(
  points = cbind(X = 0:4, Y = 0:4, Z = 0:4),
  streamline_data = list(mean_FA = 0.6, label = "A")
)
sl2 <- streamline(
  points = cbind(X = 1:3, Y = 1:3, Z = 1:3),
  streamline_data = list(mean_FA = 0.7, label = "A")
)
b <- bundle(streamlines = list(sl1, sl2))
b@streamline_data$mean_FA # c(0.6, 0.7) — lifted automatically
b@streamline_attributes # c("mean_FA", "label")


## -----------------------------------------------------------------------------
b[[1]]@streamline_data$mean_FA # 0.6


## -----------------------------------------------------------------------------
b[1]@streamline_data$mean_FA # c(0.6)
b[1]@streamline_data$label # "A"


## -----------------------------------------------------------------------------
b_sub01 <- bundle(
  streamlines = list(sl1),
  bundle_data = list(subject = "sub-01")
)
b_sub02 <- bundle(
  streamlines = list(sl2),
  bundle_data = list(subject = "sub-02")
)
bs <- bundle_set(list("sub-01" = b_sub01, "sub-02" = b_sub02))
bs@bundle_data$id_from_input_names # c("sub-01", "sub-02")


## -----------------------------------------------------------------------------
bs[["sub-01"]]@bundle_data$subject # "sub-01" (pushed down)
bs[1]@bundle_data$subject # "sub-01"


## -----------------------------------------------------------------------------
is_streamline(sl) # TRUE
is_bundle(sl) # FALSE
is_bundle_set(bs) # TRUE


## -----------------------------------------------------------------------------
b_all <- bind_bundles(sl1, b)
b_all@n_streamlines # 3


## -----------------------------------------------------------------------------
bs_merged <- bind_bundle_sets("sub-A" = b_sub01, bs)
bs_merged@n_bundles # 3


## -----------------------------------------------------------------------------
sl_20 <- reparametrize(sl1, n_points = 20L)
sl_20@n_points # 20
b_20 <- reparametrize(b, n_points = 20L)
b_20[[1]]@n_points # 20


## -----------------------------------------------------------------------------
sl_curve <- streamline(
  points = cbind(
    X = cos(seq(0, 2 * pi, length.out = 21)),
    Y = sin(seq(0, 2 * pi, length.out = 21)),
    Z = seq(0, 2 * pi, length.out = 21)
  )
)
# Scalar descriptors are fine on any streamline
get_euclidean_length(sl1)
get_curvilinear_length(sl1)
get_sinuosity(sl1)

# Per-point profiles are meaningful on a truly curved tract
head(get_curvature(sl_curve))
head(get_torsion(sl_curve))


## -----------------------------------------------------------------------------
b_long <- bundle(
  streamlines = list(
    streamline(points = cbind(X = 0:9, Y = 0:9, Z = 0:9)),
    streamline(points = cbind(X = 1:10, Y = 1:10, Z = 1:10))
  )
)
b_desc <- add_shape_descriptors(b_long)
names(b_desc@streamline_data)
names(b_desc[[1]]@point_data)


## -----------------------------------------------------------------------------
sb <- streamline(points = cbind(X = 0:4, Y = 0:4, Z = 0:4))
sc <- streamline(points = cbind(X = 1:5, Y = 1:5, Z = 1:5))

# streamline x streamline -> scalar
compute_hausdorff_distance(sb, sc)

# bundle x missing -> dist object
bd <- bundle(streamlines = list(sb, sc))
as.matrix(compute_hausdorff_distance(bd))

# bundle x streamline -> vector
compute_hausdorff_distance(bd, sb)

# bundle x bundle -> combined dist object
bd2 <- bundle(streamlines = list(sc))
as.matrix(compute_hausdorff_distance(bd, bd2))


## -----------------------------------------------------------------------------
if (requireNamespace("dti", quietly = TRUE)) {
  dfi <- as_dwifiber(bd)
  class(dfi)
  # round-trip
  as_bundle(dfi)@n_streamlines
}

