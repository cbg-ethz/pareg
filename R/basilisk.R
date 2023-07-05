my_env <- basilisk::BasiliskEnvironment(
    envname = "pareg_env",
    pkgname = "pareg",
    packages = c(
        "python=3.8",
        "tensorflow==2.4.0",
        "tensorflow_probability==0.12.0"
    )
)
