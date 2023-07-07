pareg_env <- basilisk::BasiliskEnvironment(
    envname = "pareg",
    pkgname = "pareg",
    packages = c(
        "tensorflow==2.10.0",
        "tensorflow-probability==0.19.0"
    ),
    channels = c("anaconda")
)
