library(testthat)
library(pareg)

cl <- basilisk::basiliskStart(
  pareg::pareg_env,
  testload = c("tensorflow", "tensorflow_probability")
)
basilisk::basiliskRun(
  proc = cl,
  fun = test_check,
  "pareg"
)
basilisk::basiliskStop(cl)
