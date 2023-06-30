# convaQ

devtools::install_github("SimonLarsen/convaq")

x <- high_risk_losses
y <- low_risk_losses

library(convaq)
results <- convaq(
        x, y,
        model = "statistical",
        name1 = "high", name2="low",
        p.cutoff = 0.0001
)


head(regions(results))

result0.001 <- regions(results)
