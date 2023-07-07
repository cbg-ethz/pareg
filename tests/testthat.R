library(testthat)
library(pareg)

cl <- basiliskStart(
  pareg_env,
  testload = c("tensorflow", "tensorflow_probability")
)
basiliskRun(
  proc = cl,
  fun = test_check,
  "pareg"
)
basiliskStop(cl)
