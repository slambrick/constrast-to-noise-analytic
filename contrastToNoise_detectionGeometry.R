library(tidyverse)  # We are working in the tidyverse
library(pracma)     # Numerical/mathmatical functions, 'practical math'
library(magrittr)   # Some more usefull tidyverse things, e.g. '%<>%'
library(tikzDevice) # Output graphics to TikZ for inclusion in LaTeX
library(ggthemes)   # Additional themes for ggplot, has solarized colour scheme
                    # for discrete plots
library(viridis)    # Good continous colour scheme
library(parallel)   # Parallel mapply and lapply
library(akima)      # Bicubic interpolation

#------------------------------- Parameters -----------------------------------#

# How much faster is pumping by helium than nitrogen
helium_factor <- 3.5

# The pumping from the stagnation region to the ioniser. l/s, for helium
pump_ion <- 0.4*2

# Working distance, specify a working distance, for now we have three variables
# to loop over
working_dist <- 1

# The number of points to split each variable into
N <- 200

# The range of angles around the normal to the lab floor we want to optimise for
optimise_range <- c(-30, 30)

#-------------------------------- Functions -----------------------------------#

# Contrast value from  a range of intensity values
grad_contrsat <- function(ints, psi_range) {
    grads <- abs(diff(ints)/diff(psi_range))
    grads[is.na(grads)] <- 0
    return(mean(grads))
}

#-------------------------------- Calculate -----------------------------------#


# We calculate the intensities as a function of beta and psi, then we have a
# large range of intensity data to perform our statistics on. Intensities are
# normalised such that a surface at 45deg to a detector of 10deg half cone,
# which is a flat surface in the current microscope.
N_int <- intensity(45, 10)

# The range of psis
psis <- as.vector(as.numeric(seq(0, 179)))

# The range of betas
betas <- as.vector(as.numeric(seq(0, 90)))

df_intensities <- tibble(beta = rep(betas, times = length(psis)),
                         psi = rep(psis, each = length(betas)))

df_intensities %<>% mutate(
    Is = mcmapply(intensity, psi = psi, beta = beta, mc.cores = 4)) %>%
    mutate(Is = Is/N_int)

#--------------------------- Analysis functions -------------------------------#

# Matrix of intensities
Is_mat <- matrix(df_intensities[["Is"]], nrow=length(psis), ncol=length(betas))

mtx <- outer(psis, betas, Vectorize(intensity))

# For a range of psis and fixed beta calculate the contrast to noise ratio
calc_CNR <- function(beta_value, psi_range) {
    psis_range <- seq(psi_range[1], psi_range[2], length.out = 200)
    beta_range <- vector(mode = "numeric", length = 200) + beta_value
    ints <- bicubic(x = psis, y = betas, z = mtx, 
                    x0 = psis_range, y0 = beta_range)
    return(c(grad_contrsat(ints$z, psi_range), mean(ints$z)))
}

# Calculate the CNR for a series of values of beta over a range of psis
calc_CNR_overBeta <- function(the_betas, psi_range) {
    CNR <- lapply(the_betas, calc_CNR, psi_range = psi_range)
    CNR <- matrix(unlist(CNR), nrow = 2)
    df_betas <- tibble(beta = the_betas, sig = CNR[1,], av = CNR[2,])
    return(df_betas)
}

# Convert angles to the lab frame to psi given the detection angle
psi_of_ang <- function(det_ang, ang_range) ang_range + det_ang

#-------------------------- Calculation functions -----------------------------#

# For a fixed detection angle calculate the contrast varying the aperture size
fixedDetection_beta <- function(det_angle, surf_angs, working_dist) {
    if (!working_dist) {
        # The working distance provided is 0, ignore pumping
    }
}

# For a fixed beta calculate the contrast to noise varying the detection angle
fixedBeta_detection <- function(beta_size, surf_angs, working_dist) {
    # Detection angles to look through
    det_angs <- seq(0, 60, length.out = 250)

    # Convert the surface angles to a range of psi
    psi_ranges <- lapply(det_angs, psi_of_ang, ang_range = surf_angs)

    # Calculate the contrast to noise
    CNR <- mclapply(psi_ranges, calc_CNR, beta_value = beta_size, mc.cores = 4)
    CNR <- matrix(unlist(CNR), nrow = 2)

    # Data frame the results
    df_detectionAngle <- tibble(det_ang = det_angs, sig = CNR[1,], av = CNR[2,]) %>%
        mutate(CNR = sig/sqrt(av))

    # Normalise the contrast to noise such that 45deg is 1
    tmp <- calc_CNR(10, psi_of_ang(45, surf_angs))
    CNR_norm <- tmp[1]/sqrt(tmp[2])
    df_detectionAngle %<>% mutate(CNR = CNR/CNR_norm)

    plt_sig_av_det <- df_detectionAngle %>% 
        gather(key = "variable", value = "value", -det_ang, -CNR) %>%
        ggplot(aes(x = det_ang, y = value, colour = variable)) +
        geom_line() + 
        scale_colour_manual(values = colorRampPalette(solarized_pal()(8))(2)) +
        xlab(expression("Detection angle/"*degree)) +
        theme(legend.title = element_blank())
    
    plt_CNR_det <- df_detectionAngle %>%
        ggplot(aes(x = det_ang, y = CNR)) + 
        geom_line(colour = colorRampPalette(solarized_pal()(8))(1)) + 
        xlab(expression("Detection angle/"*degree)) +
        ylab("Contrast to noise")

    # display plots
    dev.new()
    print(plt_sig_av_det)
    
    dev.new()
    print(plt_CNR_det)
    
    return(df_detectionAngle)
}

# For a series of beta calculate the contrast to noise varying the detection
# angle
fixedD_betaDetection <- function(range_betas, surf_angs, working_dist) {
    # Detection angles to look through
    det_angs <- seq(0, 60, length.out = 50)

    # Aperture sizes
    the_betas <- seq(range_betas[1], range_betas[2], length.out = 50)
    
    # Convert the surface angles to a range of psi
    psi_ranges <- lapply(det_angs, psi_of_ang, ang_range = surf_angs)

    # Calculate the contrast to noise
    CNR <- mclapply(psi_ranges, calc_CNR_overBeta, the_betas = the_betas, mc.cores = 4)
    
    for (i in seq(50)) {
        CNR[[i]] %<>% mutate(det_ang = det_angs[i])
    }
    
    CNR <- bind_rows(CNR) %>% 
        mutate(pumping = conductance_ap(beta, working_dist)) %>%
        mutate(transmission = pump_ion/(pump_ion + pumping)) %>%
        mutate(av = av*transmission, sig = sig*transmission) %>%
        mutate(CNR = sig/sqrt(av))
    
    CNR %<>% mutate(CNR = replace(CNR, av == 0, 0))
    
    # Normalise the contrast to noise such that 45deg is 1
    tmp <- calc_CNR(10, psi_of_ang(45, surf_angs))
    conduct_now <- conductance_ap(10, 3)
    transmission_now <- 0.4/(0.4 + conduct_now)
    CNR_norm <- tmp[1]*sqrt(transmission_now)/sqrt(tmp[2])
    CNR %<>% mutate(CNR = CNR/CNR_norm)
    
    # Some portions of the space we disregard as it brings the aperture too
    # close to the surface
    CNR %<>% mutate(CNR = 
        replace(
            CNR, 
            working_dist*(cos(det_ang*pi/180) - 
                tan(beta*pi/180)*sin(det_ang*pi/180)) < 0.5, 
            NaN
        )
    )
    
    # A string for the title
    title_str <- paste("Working distance = ", as.character(working_dist), "mm", 
                       ", Surface angles = (", as.character(surf_angs[1]),
                       "\u00b0,", as.character(surf_angs[2]), 
                       "\u00b0)", sep = "")
    
    # Plot of the CNR as a contour plot and heatmap
    plt_contour <- CNR %>%
        ggplot(aes(x = det_ang, y = beta)) +
        geom_raster(aes(fill = CNR), interpolate = TRUE) +
        geom_contour(aes(z = CNR), colour = "white", size = 0.5) +
        scale_fill_viridis() +
        coord_fixed(expand = FALSE) + 
        xlab(expression("Detection angle/"*degree)) + 
        ylab(expression("Half cone angle/"*degree)) +
        labs(fill = "CNR") +
        ggtitle(title_str)
    
    dev.new()
    print(plt_contour)
    
    return_list <- list(CNR, plt_contour)
    
    return(return_list)
}
