library(tidyverse)  # We are working in the tidyverse
library(pracma)     # Numerical/mathmatical functions, 'practical math'
library(magrittr)   # Some more usefull tidyverse things, e.g. '%<>%'
library(tikzDevice) # Output graphics to TikZ for inclusion in LaTeX
library(ggthemes)   # Additional themes for ggplot, has solarized colour scheme
                    # for discrete plots
library(viridis)    # Good continous colour scheme
library(parallel)   # Parallel mapply and lapply=


# The relative size of sample under mask as a function of the incident and
# detection directions
mask <- function(theta_d, theta_i) tan(theta_i*pi/180) + tan(theta_d*pi/180)

# The number of points in each direction
N <- 200

# The range of incidence and detection directions
theta_is <- seq(0, 60, length.out = N)
theta_ds <- seq(0, 60, length.out = N)

# Calculate the mask size for these and plot it as a heat map
df_masks <- tibble(theta_i = rep(theta_is, times = N),
                   theta_d = rep(theta_ds, each = N)) %>%
    mutate(rel_size = mask(theta_d, theta_i))

# Plot of the relative size of the masks
plt_masks <- df_masks %>%
    ggplot(aes(x = theta_i, y = theta_d)) +
    geom_raster(aes(fill = rel_size), interpolate = TRUE) +
    geom_contour(aes(z = rel_size), colour = "white", size = 0.5, binwidth = 0.5) +
    geom_hline(aes(yintercept = 45), colour = "firebrick1") +
    geom_vline(aes(xintercept = 45), colour = "firebrick1") +
    scale_fill_viridis() +
    coord_fixed(expand = FALSE) + 
    xlab(expression("Incidence angle/"*degree)) + 
    ylab(expression("Detection angle/"*degree)) +
    labs(fill = "Relative\nmask size")
