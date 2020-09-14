library(tidyverse)  # We are working in the tidyverse
library(pracma)     # Numerical/mathmatical functions, 'practical math'
library(magrittr)   # Some more usefull tidyverse things, e.g. '%<>%'
library(tikzDevice) # Output graphics to TikZ for inclusion in LaTeX
library(ggthemes)   # Additional themes for ggplot, has solarized colour scheme
                    # for discrete plots
library(viridis)    # Good continous colour scheme

#------------------------------- Parameters -----------------------------------#

# Range of betas (half cone angle of apeture)
psi_range <- c(0, 90 + 89)

# Range of psis (angle of the surface to the apeture)
beta_range <- c(0, 89)

# Should tex files be made of the plots?
save_tex <- TRUE

# The pumping from the stagnation region to the ioniser. l/s, for helium
pump_ion <- 0.4

# Pumping factor for helium, how much faster does helium travel through a pipe
# than nitrogen.
helium_factor <- 3.5

# The distance from the sample to the aperture, defines the area, mm^2, of the
# aperture and thus the pumping from the stagnation region back to the sample
# chamber. In mm.
working_dists <- seq(0.1, 5, by = 0.05)

# The present half cone angle
present_half_cone <- 10

# The present working distance
present_working_dist <- 3

#-------------------------------- Functions -----------------------------------#

# Calculate the signal detected as a function of psi and beta, the function only
# applies where beta + psi < 90deg
intensityAnl <- function(psi, beta) 0.25*pi*cos(psi)*(1 - cos(2*beta) )

# Cosine of the angle to the surface normal
cosChi <- function(theta, phi, psi, beta) {
    dDotn <- cos(theta)*(sin(psi)*tan(theta)*cos(phi) + cos(psi))
    result <- if_else(dDotn < 0, 0, dDotn)
    return(result)
}

# Calculate the intensity where the analytic formula cannot be applied.
intensityInt <- function(psi, beta) {
    integrand <- function(theta, phi) sin(theta)*cosChi(theta, phi, psi, beta)
    result <- quad2d(integrand, 0, beta, 0, pi)
    return(result)
}

# Calculate the intensity using the analytic if possible and the integration if
# not
intensity <- function(psi, beta) {
    if (psi + beta > pi/2)
        result <- intensityInt(psi, beta)
    else
        result <- intensityAnl(psi, beta)
    return(result)
}
intensity <- Vectorize(intensity, "psi")

# Calculates the Michelson contrast to noise between two values
Michelson <- function(I1, I2) {
    difference <- if_else(I1 > I2, I1 - I2, I2 - I1)
    CNR <- difference/(I1 + I2)^(3/2)
    return(CNR)
}

# Calculates the difference to noise between two values
diffNoise <- function(I1, I2) {
    difference <- if_else(I1 > I2, I1 - I2, I2 - I1)
    DNR <- difference/sqrt(I1 + I2)
    return(DNR)
}

# Calculates the pumping through the aperture
conductance_ap <- function(beta, working_dist) {
    # Area of the aperture, cm^2
    d <- working_dist/10
    r <- d*sin(beta*pi/180)
    
    # Pumping through that aperture
    C <- helium_factor*9.3*(2*r)^2
    return(C)
}

# Get the standard deviation from the other data
sigma_calc <- function(signals, beta_v) (signals %>% filter(beta == beta_v))[["sigma"]]
sigma_calc <- Vectorize(sigma_calc, vectorize.args = "beta_v")

# Get the mean of the absolute value of 
diff_calc <- function(signals, beta_v) (signals %>% filter(beta == beta_v))[["diff"]]
diff_calc <- Vectorize(diff_calc, vectorize.args = "beta_v")

# Get the average from the other data
av_calc <- function(signals, beta_v) (signals %>% filter(beta == beta_v))[["av"]]
av_calc <- Vectorize(av_calc, vectorize.args = "beta_v")

#-------------------------------- Calculate -----------------------------------#

# All the different betas
betas <- seq(beta_range[1], beta_range[2], by = 1)

# For each value of beta calculate the intensities for a range of psis
psis <- seq(psi_range[1], psi_range[2], by = 1)

# Overall data frame
df_signal <- tibble(beta = NA, I = NA, psi = NA)

for (i in seq(1, length(betas))) {
    dfPsi <- tibble(psi = psis, beta = betas[i]) %>%
        mutate(I = intensity(psis*pi/180, betas[i]*pi/180))
    df_signal <- bind_rows(df_signal, dfPsi)
}

# Remove the first row that are NA's
df_signal <- df_signal[-1,]

#---------------------------- Analyse contrast --------------------------------#

# Normalise the signals such that 1 corresponds to a 45deg surface angle and a
# 5deg half cone angle, i.e. the current level of signal
N <- df_signal %>% filter(beta == 5, psi == 45)
df_signal %<>% mutate(I = I/N[["I"]])

# Calculate the signal to noise, then normalise that as above
df_signal %<>% mutate(SNR = I/sqrt(I))
df_signal %<>% mutate(SNR = ifelse(is.nan(SNR), 0, SNR))

# Calculate the standard deviation contrast and noise for each value of beta
df_sigma <- tibble(beta = NA, contrast = NA, av = NA)
df_diff <- tibble(beta = NA, contrast = NA, av = NA)
for (i in seq(1, length(betas))) {
    df_oneBeta <- df_signal %>% filter(beta == betas[i], psi <= beta + 90)
    sig <- sd(df_oneBeta[["I"]])
    dif <- mean(abs(diff(df_oneBeta[["I"]])))
    av <- mean(df_oneBeta[["I"]])
    df_sigma %<>% add_row(beta = betas[i], contrast = sig, av = av)
    df_diff %<>% add_row(beta = betas[i], contrast = dif, av = av)
}
df_sigma <- df_sigma[-1,]
df_diff <- df_diff[-1,]

# RMS and differential contrast to noise
df_sigma %<>% mutate(CNR = contrast/sqrt(av))
df_diff %<>% mutate(CNR = contrast/sqrt(av))

# Normalise such that a value of 1 is the present aperture size
N1 <- (df_sigma %>% filter(beta == present_half_cone))[["CNR"]]
N2 <- (df_diff %>% filter(beta == present_half_cone))[["CNR"]]
df_sigma %<>% mutate(CNR = CNR/N1)
df_diff %<>% mutate(CNR = CNR/N2)

# Normalise the contrast too
N1 <- (df_sigma %>% filter(beta == present_half_cone))[["contrast"]]
N2 <- (df_diff %>% filter(beta == present_half_cone))[["contrast"]]
df_sigma %<>% mutate(contrast = contrast/N1)
df_diff %<>% mutate(contrast = contrast/N2)

# Combine the RMs and differential contrast to noise
df_contrast <- bind_rows("sigma" = df_sigma, "diff" = df_diff, .id = "method")

# Add the equivalent signal increase
df_contrast %<>% mutate(signal_increase = CNR^2)

# Calculate the Contrast to noise and difference to noise between 0/45, 0/90, 
# 45/90, for each beta
df_difference <- tibble(beta = NA, CNR = NA, I_0_45 = NA, I_0_90 = NA, 
                        I_45_90 = NA)
for (i in seq(1, length(betas))) {
    df_oneBeta <- df_signal %>% filter(beta == betas[i])
    I_0 <- (df_oneBeta %>% filter(psi == 0))[["I"]]
    I_45 <- (df_oneBeta %>% filter(psi == 45))[["I"]]
    I_90 <- (df_oneBeta %>% filter(psi == 90))[["I"]]
    df_difference %<>% add_row(beta = betas[i], CNR = "Difference", 
                               I_0_45 = I_0 - I_45,
                               I_0_90 = I_0 - I_90,
                               I_45_90 = I_45 - I_90)
    df_difference %<>% add_row(beta = betas[i], CNR = "DNR", 
                               I_0_45 = diffNoise(I_0, I_45),
                               I_0_90 = diffNoise(I_0, I_90),
                               I_45_90 = diffNoise(I_45, I_90))
}
df_difference <- df_difference[-1,]

# Normalise so that the current set up at a surface of 45deg has a value of 1
# for the difference between 0/45
N1 <- (df_difference %>% filter(beta == present_half_cone, 
                                CNR == "Difference"))[["I_0_45"]]
N2 <- (df_difference %>% filter(beta == present_half_cone, 
                                CNR == "DNR"))[["I_0_45"]]
df_difference %<>% mutate(
    I_0_45 = if_else(CNR == "Difference", I_0_45/N1, I_0_45/N2),
    I_0_90 = if_else(CNR == "Difference", I_0_90/N1, I_0_90/N2),
    I_45_90 = if_else(CNR == "Difference", I_45_90/N1, I_45_90/N2)
)

#--------------------------------- Pumping ------------------------------------#
# 
# # Data frame
# df_pumping <- tibble(beta = rep(betas, times = length(working_dists)),
#                      working_dist = rep(working_dists, each = length(betas)))
# 
# # Calculate the pumping for each 
# df_pumping %<>% mutate(pumping = conductance_ap(beta, working_dist))
#                      
# # Calculate the standard deviation and the RMS contrast to noise
# df_pumping %<>% mutate(
#     sigma = df_contrast %>% sigma_calc(beta),
#     av = df_contrast %>% av_calc(beta)
# )
# df_pumping %<>% mutate(
#     sigma = sigma*pump_ion/(pump_ion + pumping),
#     av = av*pump_ion/(pump_ion + pumping),
#     CNR_sigma = sigma/sqrt(av)
# )
# 
# # Normalise such that a value of 1 is the present aperture size
# N <- (df_pumping %>% filter(beta == present_half_cone, 
#     abs(working_dist - present_working_dist) < 0.0001))[["CNR_sigma"]]
# df_pumping %<>% mutate(CNR_sigma = CNR_sigma/N)
# 
# # Add the equivalent signal increase
# df_pumping %<>% mutate(signal_increase = CNR_sigma^2)

#----------------------------------- Plot -------------------------------------#

# Plot of the signal as a contour plot and heatmap
plt_contour <- df_signal %>%
    ggplot(aes(x = beta, y = psi)) +
    geom_raster(aes(fill = I), interpolate = TRUE) +
    geom_contour(aes(z = I), colour = "white", size = 0.5, binwidth = 10) +
    geom_vline(aes(xintercept = present_half_cone), colour = "firebrick1") +
    scale_fill_viridis() +
    scale_x_continuous(breaks = c(0,20,40,60,80)) +
    scale_y_continuous(breaks = c(0,20,40,60,80,100,120,140,160,180)) +
    coord_fixed(expand = FALSE) + 
    xlab(expression("Half cone angle/"*degree)) + 
    ylab(expression("Surface angle/"*degree)) +
    labs(fill = "Signal")

# Plot of the signal to noise ratio as a contour plot and heatmap
plt_contour_SNR <- df_signal %>%
    ggplot(aes(x = beta, y = psi)) +
    geom_raster(aes(fill = SNR), interpolate = TRUE) +
    geom_contour(aes(z = SNR), colour = "white", size = 0.5, binwidth = 1) +
    geom_vline(aes(xintercept = present_half_cone), colour = "firebrick1") +
    scale_fill_viridis() +
    scale_x_continuous(breaks = c(0,20,40,60,80)) +
    scale_y_continuous(breaks = c(0,20,40,60,80,100,120,140,160,180)) +
    coord_fixed(expand = FALSE) + 
    xlab(expression("Half cone angle/"*degree)) + 
    ylab(expression("Surface angle/"*degree)) +
    labs(fill = "SNR")
    
# Plot of the RMS contrst to noise and the standard deviation
facet_names <- c(
    'CNR'="Contrast to noise",
    'contrast'="Contrast"
)
leg_labs <- list(
    'sigma'="RMS",
    'diff'="Differential"
)
plt_CNR_sig <- df_contrast %>%
    gather(key = "variable", value = "value", -beta, -av, -method, -signal_increase) %>%
    ggplot(aes(x = beta, y = value, colour = method)) +
    geom_line() + 
    geom_vline(aes(xintercept = present_half_cone)) +
    scale_x_continuous(breaks = c(0,20,40,60,80)) +
    scale_colour_manual(values = colorRampPalette(solarized_pal()(8))(2),
                        labels = leg_labs)  +
    facet_wrap( ~ variable, ncol = 2, scale = "free_y", 
               labeller = as_labeller(facet_names)) +
    xlab(expression("Half cone angle/"*degree)) + 
    ylab(NULL) +
    theme(aspect.ratio = 2/(1+sqrt(5)), legend.title = element_blank())

# Plot the required signal increase
plt_signal_increase <- df_contrast %>%
    ggplot(aes(x = beta, y = signal_increase, colour = method)) +
    geom_line() + 
    geom_vline(aes(xintercept = present_half_cone)) +
    scale_x_continuous(breaks = c(0,20,40,60,80)) +
    scale_colour_manual(values = colorRampPalette(solarized_pal()(8))(2),
                        labels = leg_labs)  +
    xlab(expression("Half cone angle/"*degree)) + 
    ylab("Equivalent overall signal increase") +
    theme(aspect.ratio = 2/(1+sqrt(5)), legend.title = element_blank())

# Plot the difference to noise and the difference
leg_labs <- list(
    bquote(0^degree * "," ~ 45^degree), 
    bquote(0^degree * "," ~ 90^degree),
    bquote(45^degree * "," ~ 90^degree)
)
plt_CNR_diff <- df_difference %>%
    gather(key = "pair", value = "value", -beta, -CNR) %>%
    ggplot(aes(x = beta, y = value, colour = pair)) + 
    geom_line() +
    geom_vline(aes(xintercept = present_half_cone)) +
    scale_x_continuous(breaks = c(0,20,40,60,80)) +
#    scale_colour_manual(values = colorRampPalette(solarized_pal()(8))(3), 
#                        labels = leg_labs)  +
    facet_wrap( ~ CNR, ncol = 2, scale = "free_y") +
    xlab(expression("Half cone angle/"*degree)) + 
    ylab(NULL) +
    theme(legend.title = element_blank()) +
    theme(aspect.ratio = 2/(1+sqrt(5)))

# Plot of the contrast to noise with the pumping taken into account
#facet_names <- c(
#    'CNR_sigma'="RMS contrast to noise",
#    'sigma'="Standard deviation"
#)
#plt_CNR_pumping <- df_pumping %>%
#    filter(working_dist == 3) %>%
#    gather(key = "variable", value = "value", -beta, -av, -signal_increase) %>%
#    ggplot(aes(x = beta, y = value)) +
#    geom_line(colour = colorRampPalette(solarized_pal()(8))(1)) + 
#    geom_vline(aes(xintercept = present_half_cone)) +
#    scale_x_continuous(breaks = c(0,20,40,60,80)) +
#    facet_wrap( ~ variable, ncol = 2, scale = "free_y", labeller = as_labeller(facet_names)) +
#    xlab(expression("Half cone angle/"*degree)) + 
#    ylab(NULL) +
#    theme(aspect.ratio = 2/(1+sqrt(5)))

# # Plot the equivalent signal increase
# plt_signal_increase_pumping <- df_pumping %>%
#     filter(working_dist == 3) %>%
#     ggplot(aes(x = beta, y = signal_increase)) +
#     geom_line(colour = colorRampPalette(solarized_pal()(8))(1)) + 
#     geom_vline(aes(xintercept = present_half_cone)) +
#     scale_x_continuous(breaks = c(0,20,40,60,80)) +
#     xlab(expression("Half cone angle/"*degree)) + 
#     ylab("Equivalent overall signal increase") +
#     theme(aspect.ratio = 2/(1+sqrt(5)))
# 
#     
# # Remove the NaNs for plotting
# df_pumping[is.na(df_pumping)] <- 0
# 
# # Plot of the signal improvment as a function of both half cone angle and the
# # working distance.
# plt_contour_pumping <- df_pumping %>%
#     ggplot(aes(x = beta, y = working_dist)) +
#     geom_raster(aes(fill = signal_increase), interpolate = TRUE) +
# #    geom_contour(aes(z = signal_increase), colour = "white", size = 0.5, binwidth=1) +
#     geom_vline(aes(xintercept = present_half_cone), colour = "firebrick1") +
#     geom_hline(aes(yintercept = present_working_dist), colour = "firebrick1") +
#     scale_fill_viridis() +
#     coord_cartesian(expand = FALSE) + 
#     xlab(expression("Half cone angle/"*degree)) + 
#     ylab(expression("Working distance/mm")) +
#     labs(fill = "Equivalent signal\nincrease")
# 
# # Plot of the RMS CNR as a function of both half cone angle and the
# # working distance.
# plt_contour_pumpingCNR <- df_pumping %>%
#     ggplot(aes(x = beta, y = working_dist)) +
#     geom_raster(aes(fill = CNR_sigma), interpolate = TRUE) +
#     geom_contour(aes(z = CNR_sigma), colour = "white", size = 0.5, binwidth=1) +
#     geom_vline(aes(xintercept = present_half_cone), colour = "firebrick1") +
#     geom_hline(aes(yintercept = present_working_dist), colour = "firebrick1") +
#     scale_fill_viridis() +
#     coord_cartesian(expand = FALSE) + 
#     xlab(expression("Half cone angle/"*degree)) + 
#     ylab(expression("Working distance/mm")) +
#     labs(fill = "CNR")


#----------------------------------- Save -------------------------------------#

# TODO: save to text file

# TODO: save the plots to eps

# Save the plots as TikZ pictures for inclusion in latex. for 160mm wide page,
# and 10pt font with the figures/figure captions \small
if (save_tex) {
    # The 2D plots
    plt_contour <- plt_contour + 
        xlab("Half cone angle/$^\\circ$") +
        ylab("Surface angle/$^\\circ$") +
        theme(text = element_text(size = 9))
    tikz("2Dplot.tex", standAlone = FALSE, width = 3.15, height = 4.5)
    print(plt_contour)
    dev.off()
    
    plt_contour_SNR <- plt_contour_SNR + 
        xlab("Half cone angle/$^\\circ$") +
        ylab("Surface angle/$^\\circ$") +
        theme(text = element_text(size = 9))
    tikz("2Dplot_SNR.tex", standAlone = FALSE, width = 3.15, height = 4.5)
    print(plt_contour_SNR)
    dev.off()
    
    # Plot of the RMS contrast to noise
    plt_CNR_sig <- plt_CNR_sig +
        xlab("Half cone angle/$^\\circ$") +
        theme(text = element_text(size = 9))
    tikz("RMS_contrastNoise.tex", standAlone = FALSE, width = 6.3, 
         height = 2.15)
    print(plt_CNR_sig)
    dev.off()
    
    # Plot of the difference and the DNR
    tex_labs <- list("0$^\\circ$, 45$^\\circ$", "0$^\\circ$, 90$^\\circ$", "45$^\\circ$, 90$^\\circ$")
    plt_CNR_diff <- plt_CNR_diff +
        scale_colour_manual(values = colorRampPalette(solarized_pal()(8))(3), 
                            labels = tex_labs) +
        xlab("Half cone angle/$^\\circ$") +
        theme(text = element_text(size = 9))
    tikz("Differnece_contrastNoise.tex", standAlone = FALSE, width = 6.3, 
         height = 2.15)
    print(plt_CNR_diff)
    dev.off()
    
    # Plot of the RMS contrast to noise
    plt_signal_increase <- plt_signal_increase +
        xlab("Half cone angle/$^\\circ$") +
        theme(text = element_text(size = 9))
    tikz("signal_increase.tex", standAlone = FALSE, width = 3.15, height = 2.1)
    print(plt_signal_increase)
    dev.off()
}
