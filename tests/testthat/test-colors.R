test_that("scale_color_paletteer_d works", {
  p <- scale_fill_pnw()
  expect_s3_class(p, c("ScaleDiscrete", "Scale", "ggproto", "gg"))
})
