library(ggplot2)
library(tibble)
library(hexSticker)

n <- 1000
x <- seq(0, 1, length.out = n)
m <- function(x) (2*(x- 0.5))^3
Y_no_error <-  m(x)
varepsilon <- 1-2*rbinom(n,1, p=0.5)
Y <- (Y_no_error + varepsilon) %>% sample(n)
dat <- tibble(x=x, Y=Y, Y_no_error = Y_no_error)

p <- ggplot(dat, aes(x=x, y=Y_no_error)) +
  geom_line(col="red") +
  scale_y_continuous("") +
  geom_point(aes(x=x, y=Y), size = .1) +
  theme_void()

sticker(p, package="UncoupledIsoReg", p_size=15, s_x=1, s_y=.75, s_width=1.3, s_height=1,
        h_color = "#3234a8",
        filename="inst/figures/sticker.png")
