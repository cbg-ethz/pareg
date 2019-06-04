library(purrr)
library(ggplot2)

test_that("idea works", {
  # setup terms
  df.terms <- rbind(
    data.frame(
      name="foo",
      gene=paste("g", 100:200, sep="")
    )
  )

  # compute enrichments
  set.seed(13)

  window.size <- 10
  genes <- paste("g", 90:109, sep="")

  df.enr <- map_df(0:10, function (x) {
    df.genes <- data.frame(
      gene=genes,
      pvalue=c(
        rbeta(x, 1, 1), # uniform
        rbeta(window.size, .1, 1), # "significant" p-values
        rbeta(20-x-window.size, 1, 1) # uniform
      )
    )

    df.res <- pareg(df.genes, df.terms)

    cbind(
      df.res,
      data.frame(x=x)
    )
  })

  # validation plot
  ggplot(df.enr, aes(x=x, y=enrichment)) +
    geom_line() +
    theme_minimal()

  # assertions
  expect_true(df.enr[df.enr$x == 0,]$enrichment > 0)
  expect_true(df.enr[df.enr$x == 10,]$enrichment < 0)
})
