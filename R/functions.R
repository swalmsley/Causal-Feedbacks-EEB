
# Custom functions to be called in '_targets.R'


# Save plot ---------------------------------------------------------------
save_figure <- function(path, w, h, call) {
  png(path, width=w, height=h, units='in', res=2001)
  print(call)
  dev.off()
}


# Plot densities ----------------------------------------------------------
density_plot <- function(fit, true_effect, var_name, var_label, xloc, my_cols, colNum) {
  
  g <- ggplot(data = as_draws_df(fit), aes(x=get(var_name))) + 
    geom_density(color=my_cols[colNum], fill=my_cols[colNum], alpha=0.75) +
    xlim(-2.5,2.5)+
    geom_vline(xintercept=true_effect,linetype='dashed',linewidth=0.5)+
    labs(x=paste(var_label),y='Posterior density')+
    annotate('text', x=xloc, y=2.5,label='True effect',size=2.5)+
    theme_classic()
  g
  
}
# fit <- iv_model
# true_effect <- -1
# var_name <- 'b_violenceobserved_development_applied'
# var_label <- 'Estimated effect'
# xloc <- -1.4
# colNum <- 1



# Simulate x causes y, and "correct' regression of y on x -----------------
regression_sim_x_to_y <- function() {
  
  n <- 50
  beta <- seq(-2,2,by=0.1)
  output <- data.table(index=(1:length(beta)))
  
  for (i in 1:length(beta)){
    
    x <- rnorm(n,0,1) # simulate x values
    y <- beta[i]*x + rnorm(n) # simulate y values
    
    m <- lm(y ~ x) # fit simple regression
    
    # save results to output 
    output[i, true_xy:=beta[i],]
    output[i, coefficient:=m$coefficients[2],]
    output[i, lowCI:=confint(m)[2],]
    output[i, highCI:=confint(m)[4],]
    
  }
  return(output)
}



# Simulate y causes x, and "incorrect" regression of y on x ---------------
regression_sim_reverse <- function() {
  
  n <- 50
  beta <- seq(-2,2,by=0.1)
  output <- data.table(index=(1:length(beta)))
  
  for (i in 1:length(beta)){
    
    # Note that this is the reverse of typical causal expection of x -> y
    y <- rnorm(n,0,1) # simulate y values
    x <- beta[i]*y + rnorm(n) # simulate x values
    
    m <- lm(y ~ x) # fit simple regression
    
    # save results to output 
    output[i, true_yx:=beta[i],]
    output[i, coefficient:=m$coefficients[2],]
    output[i, lowCI:=confint(m)[2],]
    output[i, highCI:=confint(m)[4],]
    
  }
  return(output)
}



# Simulate time series data for cross-lagged effect examples --------------
sim_temporal_data <- function() {
  
  # Initialize variables
  n_time_steps <- 20
  n_individuals <- 1
  
  # auto effects (dampening)
  a_to_a <- 0.25
  b_to_b <- 0.25
  
  # cross effects
  a_to_b <- 0.3
  b_to_a <- -0.3
  
  # Create an empty data frame to store the results
  df <- data.frame(
    individual = integer(),
    time = integer(),
    a = numeric(),
    b = numeric())
  
  # Loop over individuals
  for (i in 1:n_individuals) {
    
    # Initialize starting values
    a <- numeric(n_time_steps)
    b <- numeric(n_time_steps)
    a[1] <- 1
    b[1] <- 1
    
    # Loop over time steps
    for (t in 2:n_time_steps) {
      
      # A changes as a function of previous values of A (auto) and B (cross) plus random variation
      a[t] <- a[t-1] + b_to_a * b[t-1] + rnorm(1, mean = 0, sd = 0.25)
      
      # A changes as a function of previous values of A (auto) and B (cross) plus random variation
      b[t] <- b[t-1] + a_to_b * a[t-1] + rnorm(1, mean = 0, sd = 0.25)
      
    }
    
    # Append the results for this individual to the dataframe
    df <- rbind(df, data.frame(
      individual = rep(i, n_time_steps),
      time = 1:n_time_steps,
      a = a,
      b = b
    ))
  }
  
  # Convert from wide to long format for plotting with ggplot
  df_long <- df %>%
    tidyr::pivot_longer(c(a, b), names_to = "trait", values_to = "value")
  
  df <-
    df %>%
    mutate(a_1 = lag(a),
           b_1 = lag(b))
  
  return(df)
  
}



# Fit model with cross-lagged effects -------------------------------------
fit_panel_model <- function(df) {
  
  f1 <- bf(b ~ 0 + Intercept + b_1 + a_1)
  f2 <- bf(a ~ 0 + Intercept + a_1 + b_1)
  
  m <-
    brm(data = df,
        family = 'Gaussian',
        f1 + f2 + set_rescor(FALSE),
        prior = c(prior(normal(0, 1), class = b, resp = b, coef = Intercept),
                  prior(normal(0, 1), class = b, resp = a, coef = Intercept),
                  prior(normal(0, 1), class = b, resp = b),
                  prior(normal(0, 1), class = b, resp = a),
                  prior(exponential(1), class = sigma, resp = b),
                  prior(exponential(1), class = sigma, resp = a)),
        control = list(max_treedepth = 15),
        iter = 4000, warmup = 2000, chains = 4, cores = 4)
  
  return(m)
  
}


# Simulate x and y cause each other, and fit both possible regress --------
regression_sim_reciprocal <- function() {
  
  n <- 50
  beta_xy <- seq(-2,2,by=0.2)
  beta_yx <- seq(-2,2,by=0.2)
  
  output <- data.table(expand.grid(beta_xy=beta_xy, beta_yx=beta_yx))
  
  for (i in 1:nrow(output)){
    
    x_initial <- rnorm(n,0,1) # starter
    y_initial <- rnorm(n,0,1) # starter
    
    x <- x_initial + (output[i,beta_yx,])*y_initial + rnorm(n,0,1) # simulate x values
    y <- y_initial + (output[i,beta_xy,])*x_initial + rnorm(n,0,1) # simulate y values
    
    m_y <- lm(y ~ x) # fit simple regression
    m_x <- lm(x ~ y) # fit simple regression
    
    # save results to output 
    output[i, m_y_coefficient:=m_y$coefficients[2],]
    output[i, m_y_lowCI:=confint(m_y)[2],]
    output[i, m_y_highCI:=confint(m_y)[4],]
    
    output[i, m_x_coefficient:=m_x$coefficients[2],]
    output[i, m_x_lowCI:=confint(m_x)[2],]
    output[i, m_x_highCI:=confint(m_x)[4],]
    
  }
  return(output)
}



# Draw together above simulations to create panels for Figure 2 -----------
figure_2 <- function(){
  
  # SETUP 
  my_cols <- c('steelblue','red')
  
  # PANEL A
  output <- regression_sim_x_to_y()
  
  A <- ggplot(data=output, aes(x=true_xy, y=coefficient))+
    geom_abline(intercept=0, slope=1)+
    geom_errorbar(aes(ymin=lowCI, ymax=highCI), width=0, linewidth=0.75, color='grey')+
    geom_point(color=my_cols[1])+
    labs(x = expression("Causal effect of x on y"),
         y = bquote("Effect of x on y from regression (" * beta[1] * ")")) + 
    xlim(-2,2)+ylim(-2,2)+
    
    annotate('text', x=-0.8, y=1.3,label='(Causal simulation)',size=3.25,col='black')+
    annotate('text', x=-0.8, y=0.8,label='(Regression)',size=3.25,col=my_cols[1])+
    
    annotate('text', x=-0.8, y=1,label=expression(y == beta[0] + beta[1]*x + epsilon),size=4.5,col=my_cols[1])+
    annotate('text', x=-0.8, y=1.5,label=expression(y %<-% x),size=4.5,col='black')+
    annotate('text', x=-0.6, y=2,label='Correct causal direction',size=4,col='black', fontface='bold.italic')+
    
    theme_classic()
  
  A
  
  
  # PANEL B
  output <- regression_sim_reverse()
  
  B <- ggplot(output, aes(x=true_yx, y=coefficient))+
    geom_abline(intercept=0, slope=0)+
    geom_errorbar(aes(ymin=lowCI, ymax=highCI), width=0, linewidth=0.75, color='grey')+
    geom_point(color=my_cols[1])+
    labs(x = expression("(Reverse) causal effect of y on x"),
         y = bquote("Effect of x on y from regression (" * beta[1] * ")")) + 
    xlim(-2,2)+ylim(-2,2)+
    annotate('text', x=1, y=-0.6,label='True effect of x on y is 0',size=3.25,col='black')+
    annotate('text', x=-0.8, y=1.2,label=expression(y == beta[0] + beta[1]*x + epsilon),size=4.5,col=my_cols[1])+
    annotate('text', x=-0.8, y=1.5,label=expression(y %->% x),size=4.5,col='black')+
    annotate('text', x=-0.8, y=2,label='Reverse causation',size=4,col='black', fontface='bold.italic')+
    annotate(
      "curve",
      x = 1, y = -0.5,      # arrow tail near your text
      xend = 1, yend = -0.05,   # arrow head touching the zero line
      curvature = 0.25,
      arrow = arrow(type='closed', length = unit(0.04, "inches")),
      linewidth = 0.2
    )+
    theme_classic()
  B
  
  # PANEL C
  output <- regression_sim_reciprocal()
  
  C <- ggplot(output, aes(x=beta_xy, y=m_y_coefficient))+
    
    geom_abline(intercept=0, slope=1)+
    #geom_smooth(method='lm',se=F, color='black', linewidth=0.75)+
    #geom_errorbar(aes(ymin=m_y_lowCI, ymax=m_y_highCI), width=0, linewidth=0.75, color='grey')+
    geom_point(color=my_cols[1],alpha=0.5)+
    
    #geom_smooth(method='lm',se=F, color='black', linewidth=0.75)+
    #geom_errorbar(inherit.aes=FALSE, aes(x=beta_yx, y=m_x_coefficient, ymin=m_y_lowCI, ymax=m_y_highCI), width=0, linewidth=0.75, color='grey')+
    geom_point(inherit.aes=FALSE, aes(x=beta_yx, y=m_x_coefficient),color=my_cols[2],alpha=0.5)+
    
    labs(x = bquote("Causal effect (" * x ~ "on" ~ y * ") or (" * y ~ "on" ~ x * ")"),
         y=bquote(paste('Effects estimated from regression (', beta[1], ')')))+
    xlim(-2,2)+ylim(-2,2)+
    annotate('text', x=1, y=-1,label='Neither regression recovers \n direct causal effects',size=3.25,col='black')+
    annotate('text', x=-0.8, y=1.1+0.05,label=expression(y == beta[0] + beta[1]*x + epsilon),size=4.5,col=my_cols[1])+
    annotate('text', x=-0.8, y=0.875+0.05,label=expression(x == beta[0] + beta[1]*y + epsilon),size=4.5,col=my_cols[2])+
    annotate('text', x=-0.8, y=1.5,label=expression(y %<->% x),size=4.5,col='black')+
    annotate('text', x=-0.8, y=1.5,label=expression(y %<->% x),size=4.5,col='black')+
    
    annotate('text', x=-0.65, y=2,label='Bidirectional causation',size=4,col='black', fontface='bold.italic')+
    
    theme_classic()
  
  # COMBINE PANELS
  
  figure <- (A | B | C) + plot_annotation(tag_levels = c('A'))
  
  return(figure)
  
}



# Portion of panel model figure -------------------------------------------
panel_model_figure <- function(fit, pf1, pf2) {
  
  my_cols <- c('steelblue','red')
  
  p2 <- density_plot(fit, 0.3, 'b_b_a_1', 'Effect of boldness on size', 1.5, my_cols, 2) + xlim(-0.6,0.6)+
    geom_density(inherit.aes=FALSE,data=as_draws_df(pf2),aes(x=get('b_a')),color='grey50',fill='grey50',alpha=0.35)+
    theme(axis.title = element_text(size=8))
  
  p3 <- density_plot(fit, -0.3, 'b_a_b_1', 'Effect of size on boldness', 1.5, my_cols, 1) + xlim(-0.6,0.6)+
    geom_density(inherit.aes=FALSE,data=as_draws_df(pf1),aes(x=get('b_b')),color='grey50',fill='grey50',alpha=0.35)+
    theme(axis.title = element_text(size=8))
  
  plot <- (p2/p3)
  
  return(plot)
  
}
# fit <- tar_read(pfit)
# pf1 <- tar_read(panel_regular_1)
# pf2 <- tar_read(panel_regular_2)


# Simulate and model eco-evo feedback (continuous time example) -----------
eco_evo_sim <- function() {
  
  # Set up time points
  times <- seq(from=0, to=200, by=1)
  Nobs <- length(times)
  
  # System parameters for perfect oscillations
  A_pop <- -0.01         # Modest self-regulation for population 
  A_trait <- -0.01      # Modest self-regulation for trait
  
  # Exactly balanced cross-effect parameters
  Eco_to_Evo <- -0.1    # Effect of population on trait
  Evo_to_Eco <- 0.1   # Effect of trait on population (negative feedback)
  
  # Equal continuous intercepts
  B_pop <- 0
  B_trait <- 0
  
  # System noise
  G_pop <- 0.2
  G_trait <- 0.2
  G_cross <- 0.2
  
  # Create data frame
  data <- data.frame(
    Time = times,
    PopDensity = rep(NA, Nobs),
    TraitValue = rep(NA, Nobs)
  )
  
  # Specify initial values
  initialPopDensity <- 0
  initialTraitValue <- 0 
  
  # Initialize states
  PopDensityState <- initialPopDensity
  TraitValueState <- initialTraitValue
  
  # Generate data with eco-evolutionary feedback
  for(obsi in 1:Nobs) {
    
    # if first observation, just use initial values
    if(obsi == 1) {
      data$PopDensity[obsi] <- PopDensityState
      data$TraitValue[obsi] <- TraitValueState
      
    } else {
      
      # Calculate deterministic changes:
      dPopDensity <- A_pop * PopDensityState + Evo_to_Eco * TraitValueState + B_pop
      dTraitValue <- A_trait * TraitValueState + Eco_to_Evo * PopDensityState + B_trait
      
      # Generate system noise
      systemNoisePopDensity <- rnorm(n=1, mean=0, sd=1)
      systemNoiseTraitValue <- rnorm(n=1, mean=0, sd=1)
      systemNoiseCross <- rnorm(n=1, mean=0, sd=1)
      
      # Update states
      PopDensityState <- PopDensityState + dPopDensity + 
        G_pop * systemNoisePopDensity + G_cross * systemNoiseCross
      
      TraitValueState <- TraitValueState + dTraitValue + 
        G_trait * systemNoiseTraitValue + G_cross * systemNoiseCross
    }
    
    data$PopDensity[obsi] <- PopDensityState
    data$TraitValue[obsi] <- TraitValueState
    
  }
  
  # Add measurement error
  data$PopDensity <- data$PopDensity + rnorm(n=Nobs, mean=0, sd=0.5)
  data$TraitValue <- data$TraitValue + rnorm(n=Nobs, mean=0, sd=0.5)
  
  # Add ID variable for ctsem format
  data$ID <- 1
  
  # Specify continuous time model for eco-evolutionary feedback
  ct_model <- ctModel(
    # Variable naming
    manifestNames = c("PopDensity", "TraitValue"),  # Observed variables
    latentNames = c("PopDensity", "TraitValue"),    # Latent processes (same as manifest)
    
    # Data structure parameters
    time = "Time",                                  # Time column name
    id = "ID",                                      # ID column name
    type = "stanct",                                # Continuous time model
    
    # Measurement model parameters
    LAMBDA = diag(1, 2),                            # Identity mapping from latent to observed
    MANIFESTMEANS = 0,                              # No measurement intercept
    MANIFESTVAR = c(
      "measError_Pop", 0,                          # Measurement error for population
      0, "measError_Trait"                         # Measurement error for trait
    ),
    
    # Dynamics matrix - key eco-evolutionary parameters
    DRIFT = c(
      "A_pop", "Evo_to_Eco",                       # Population dynamics row
      "Eco_to_Evo", "A_trait"                      # Trait dynamics row
    ),
    
    # Continuous intercepts (baseline rates)
    CINT = c("B_pop", "B_trait"),                  
    
    # Initial state parameters
    T0MEANS = c("init_PopDensity", "init_TraitValue"),
    
    # System noise parameters (stochastic components)
    DIFFUSION = c(
      "G_pop", 0,                                  # Pop noise and upper triangle 0
      "G_cross", "G_trait"                         # Correlation & trait noise
    )
  )
  
  # Fit the model to data
  ct_fit <- ctStanFit(datalong = data, ctstanmodel = ct_model, optimize=FALSE, priors=TRUE)
  
  return(ct_fit)
  
}



# Plot eco-evo dynamics (continuous time example) -------------------------
plot_dynamics <- function(ct_fit) {
  
  # Summarize results
  summary(ct_fit, parmatrices = FALSE)
  
  # Extract predictions without plotting
  kalman_data <- ctKalman(fit=ct_fit, 
                          kalmanvec=c('y', 'yprior'),  # Original data and predictions
                          errorvec='yprior',          # Include uncertainty estimates 
                          plot=FALSE)                 # Don't plot yet
  
    # Create a custom plot with correct column names
    eco_evo_plot <- ggplot() +
      # Add raw data points (observational data)
      geom_point(data=kalman_data[kalman_data$Element == 'y' & !is.na(kalman_data$value), ],
                 aes(x=Time, y=value, color=Row),
                 alpha=0.5, size=1.5) +
      
      # Add model prediction lines
      geom_line(data=kalman_data[kalman_data$Element == 'yprior', ],
                aes(x=Time, y=value, color=Row),
                size=1) +
      
      # Optional: Add uncertainty ribbons with less opacity
      geom_ribbon(data=kalman_data[kalman_data$Element == 'yprior', ],
                  aes(x=Time, y=value, ymin=value-sd, ymax=value+sd, fill=Row),
                  alpha=0.1, linetype=0) +
      
      # Set custom colors
      scale_color_manual(values=c("PopDensity"="steelblue", "TraitValue"="red")) +
      scale_fill_manual(values=c("PopDensity"="steelblue", "TraitValue"="red")) +
      
      # Labels and theming
      labs(title="Eco-evolutionary dynamics",
           x="Time", y="Value") +
      annotate('text', 28, -4.5, label='Trait value', color='red')+
      annotate('text', 20, 3, label='Population density', color='steelblue')+
      theme_classic() +
      ylim(-6, 6) +
      theme(legend.position="none",
            legend.title=element_blank())
    
    eco_evo_plot
    
    # Display the plot
    return(eco_evo_plot)
    
  }



# Plot eco-evo effects (continuous time example) --------------------------
plot_effects <- function(ct_fit) {
  
  # extract original data from model
  data <- ct_fit$standata$Y
  
  ## Extract posteriors for ctsem model.
  post <- ctExtract(ct_fit)
  drift_samples <- post$pop_DRIFT   # posterior samples: draws × 2 × 2
  
  ct_Eco_to_Evo <- drift_samples[, 2, 1]  # Pop → Trait
  ct_Evo_to_Eco <- drift_samples[, 1, 2]  # Trait → Pop
  
  # Weakly informative priors
  reg_prior <- prior(normal(0, 1), class = "b") +
    prior(exponential(1), class = "sigma")
  
  # Trait ~ Pop
  m_trait_on_pop <- brm(
    TraitValue ~ PopDensity,
    data   = data,
    family = gaussian(),
    prior  = reg_prior,
    chains = 2, iter = 1000, cores = 2,
    refresh = 0
  )
  
  # Pop ~ Trait
  m_pop_on_trait <- brm(
    PopDensity ~ TraitValue,
    data   = data,
    family = gaussian(),
    prior  = reg_prior,
    chains = 2, iter = 1000, cores = 2,
    refresh = 0
  )
  
  # Posterior draws for slopes
  draws_trait_on_pop  <- as.data.frame(m_trait_on_pop)
  draws_pop_on_trait  <- as.data.frame(m_pop_on_trait)
  
  reg_Eco_to_Evo <- draws_trait_on_pop$`b_PopDensity`   # Pop -> Trait
  reg_Evo_to_Eco <- draws_pop_on_trait$`b_TraitValue`   # Trait -> Pop
  
  post_df <- rbind(
    data.frame(Value = ct_Eco_to_Evo,  Effect = "Eco → Evo", Method = "ctsem"),
    data.frame(Value = reg_Eco_to_Evo, Effect = "Eco → Evo", Method = "regression"),
    data.frame(Value = ct_Evo_to_Eco,  Effect = "Evo → Eco", Method = "ctsem"),
    data.frame(Value = reg_Evo_to_Eco, Effect = "Evo → Eco", Method = "regression")
  )
  
  Eco_to_Evo <- -0.1
  Evo_to_Eco <- 0.1
  
  truth_df <- data.table(
    Effect = c("Eco → Evo", "Evo → Eco"),
    True   = c(Eco_to_Evo,   Evo_to_Eco)
  )
  
  post_df$Method <- factor(post_df$Method, levels = c('regression', 'ctsem'))
  
  post_df <- data.table(post_df)
  
  g1 <- ggplot(post_df[Effect=="Eco → Evo",,], aes(x = Value, colour = Method, fill = Method)) +
    geom_density(alpha = 0.3, linewidth = 0.8, adjust=2) +
    geom_vline(data = truth_df[Effect=="Eco → Evo",,], aes(xintercept = True), 
               linetype = "dashed", colour = "black") +
    # facet_wrap(~ Effect, nrow = 2) +
    labs(x = "Effect", y = "Density of posterior")+
    theme_classic(base_size = 12) +
    scale_fill_manual(values=c('grey90', 'black'))+
    scale_color_manual(values=c('grey50','black'))+
    xlim(-0.25, 0.25)+
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"))
  
  g2 <- ggplot(post_df[Effect=="Evo → Eco",,], aes(x = Value, colour = Method, fill = Method)) +
    geom_density(alpha = 0.3, linewidth = 0.8, adjust=2) +
    geom_vline(data = truth_df[Effect=="Evo → Eco",,], aes(xintercept = True), 
               linetype = "dashed", colour = "black") +
    # facet_wrap(~ Effect, nrow = 2) +
    labs(x = "Effect", y = "Density of posterior")+
    theme_classic(base_size = 12) +
    scale_fill_manual(values=c('grey90', 'black'))+
    scale_color_manual(values=c('grey60','black'))+
    xlim(-0.25, 0.25)+
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"))
  
  return(g1 / g2)
  
}



# Simulate data for toy instrumental variable examples --------------------
run_IV_sim <- function(){
  
  n <- 100
  
  Instrument_1 <- rnorm(n)
  Instrument_2 <- rnorm(n)
  
  A_initial <- rnorm(n) + 1*Instrument_2 # initial values of A and B are influenced by their respective instruments
  B_initial <- rnorm(n) + 1*Instrument_1
  
  A <- A_initial + (-1)*B_initial + rnorm(n) # A is influenced by instrument_2 and B 
  B <- B_initial + (1)*A_initial  + rnorm(n) # B is influenced by instrument_1 and A
  
  # Create a dataframe
  data <- data.frame(A, B, Instrument_1, Instrument_2)
  
  return(data)
  
}



# Fit instrumental variable models ----------------------------------------
fit_IV_model_brms <- function(data, option){
  
  # Option 1: direct effect of B -> A
  
  if (option==1) {
    
    f1 <- bf(B ~ 1 + Instrument_2)
    f2 <- bf(A ~ 1 + B)
    model <- brm(f1 + f2 + set_rescor(TRUE),
                 data=data, 
                 prior = c(# First model
                   prior(normal(0, 0.5), class = Intercept, resp = B),
                   prior(normal(0, 0.5), class = b, resp = B),
                   prior(exponential(0.5), class = sigma, resp = B),
                   
                   # Second model
                   prior(normal(0, 0.5), class = Intercept, resp = A),
                   prior(normal(0, 0.5), class = b, resp = A),
                   prior(exponential(0.5), class = sigma, resp = A),
                   
                   # rho
                   prior(lkj(2), class = rescor)),
                 iter = 10000, warmup = 5000, chains = 4, cores = 4,
                 family = 'Gaussian')
    
  }
  
  # Option 2: direct effect of A -> B
  
  if (option==2) {
    
    f1 <- bf(A ~ 1 + Instrument_1)
    f2 <- bf(B ~ 1 + A)
    model <- brm(f1 + f2 + set_rescor(TRUE),
                 data=data, 
                 prior = c(# First model
                   prior(normal(0, 0.5), class = Intercept, resp = B),
                   prior(normal(0, 0.5), class = b, resp = B),
                   prior(exponential(0.5), class = sigma, resp = B),
                   
                   # Second model
                   prior(normal(0, 0.5), class = Intercept, resp = A),
                   prior(normal(0, 0.5), class = b, resp = A),
                   prior(exponential(0.5), class = sigma, resp = A),
                   
                   # rho
                   prior(lkj(2), class = rescor)),
                 iter = 10000, warmup = 5000, chains = 4, cores = 4,
                 family = 'Gaussian')
  }
  
  # Option 3: standard regression B -> A
  
  if (option==3) {
    
    model <- brm(A ~ 1 + B,
                 data=data,
                 chains = 4, cores = 4,
                 family = 'Gaussian')
  }
  
  
  # Option 4: standard regression A -> B
  
  if (option==4) {
    
    model <- brm(B ~ 1 + A,
                 data=data,
                 chains = 4, cores = 4,
                 family = 'Gaussian')
  }
  
  return(model)
  
}



# Create code portion of figure for instrumental variable plot ------------
instrumental_variable_plot <- function(fit1, fit2, r1, r2) {
  
  # Set up colors
  my_cols <- c('steelblue','red')
  
  pm1 <- density_plot(fit1, 1, 'b_A_B', 'Estimated effect', 1.5, my_cols, 1) +
    geom_density(inherit.aes=FALSE,data=as_draws_df(r1),aes(x=get('b_B')),
                 color='white',fill='grey50',alpha=0.6) +
    annotate('text', x=-0.75,color='grey50', y=2,label='Naive regression',size=2.5)+
    annotate('text', x=1.6,color=my_cols[1], y=1,label='IV analysis',size=2.5)
  pm1
  
  pm2 <- density_plot(fit2, -1, 'b_B_A', 'Estimated effect', -1.5, my_cols, 2) +
    geom_density(inherit.aes=FALSE,data=as_draws_df(r2),aes(x=get('b_A')),
                 color='white',fill='grey50',alpha=0.6) +
    annotate('text', x=0.75,color='grey50', y=2,label='Naive regression',size=2.5)+
    annotate('text', x=-1.5,color=my_cols[2], y=1,label='IV analysis',size=2.5)
  
  
  pr1 <- density_plot(r1, -1, 'b_B', 'Estimated effect', -1.6, my_cols, 2)
  pr2 <- density_plot(r2, 1, 'b_A', 'Estimated effect', -1.6, my_cols, 2)
  
  plot <- ((pm1 / pm2)) 
  return(plot)
  
}
# fit1 <- tar_read(IV_m1)
# fit2 <- tar_read(IV_m2)
# r1 <- tar_read(reg_m1)
# r2 <- tar_read(reg_m2)





# Create supplemental figure 1 --------------------------------------------
framework_figure <- function() {
  
  my_cols <- c('steelblue','red')
  
  # Oversimplified DAG
  dag_text <- 'dag{
  A -> B;
  A -> C;
  C -> D;
  C -> E;
  E -> F;
  E -> G;
  F -> H;
  F -> I;
  G -> J;
  G -> K;
  }'
  
  dag <- dagitty(dag_text)
  
  # Manually set the coordinates of the nodes
  coordinates(dag) <- list(
    x = c(A=0, B=-1, C=1, D=0, E=2, F=0, G=4, H=-1, I=1, J=3, K=5),
    y = c(A=0, B=-1, C=-1, D=-2, E=-2, F=-4, G=-4, H=-5, I=-5, J=-5, K=-5)
  )
  
  p <- ggdag(dag)+
    geom_dag_edges(edge_color='white')+
    
    geom_dag_node(colour='white') +
    
    
    geom_dag_edges(edge_width=0.25, 
                   start_cap=ggraph::circle(6, 'mm'),
                   end_cap=ggraph::circle(6, 'mm'))+
    
    annotate('text', x=0, y=0.3,label="Assume a researcher wants to estimate \n the causal effect of X on Y",size=2.25,fontface='plain',col=my_cols[1])+
    #annotate('text', x=0, y=0.25,label=expression(y == beta[0] + beta[1]*x + epsilon),size=5,fontface='bold',col=my_cols[1])+
    annotate('text',0 ,0, label='Could Y be a cause of X?',size=2,fontface='bold',col='black')+
    annotate('text',-1 ,-1, label='Proceed to infer causal effect, \n assuming given DAG.',size=2,col='black')+
    annotate('text',1 ,-1, label='Are the directions of causation \n between X and Y relevant \n to the scientific question?',size=2,col='black')+
    annotate('text',0 ,-2, label='Interpret effect of Y on X \n as non-causal correlation',size=2,col='black')+
    annotate('text',2 ,-2, label='Are there repeated measurements \n of Y and X across time?',size=2,col='black')+
    annotate('text',4 ,-4, label='Does the hypothesized effect \n of X on Y occur at discrete, \n separable time periods?',size=2,col='black')+
    annotate('text',0 ,-4, label= 'Is there an instrument \n available for Y?',size=2,col='black')+
    annotate('text',-1 ,-5, label='Not possible to identify \n direct causal effects',size=2,col='black')+
    annotate('text',1 ,-5, label= 'Instrumental \n variable methods',size=2,col='black')+
    annotate('text',3 ,-5, label= 'Continuous time \n methods',size=2,col='black')+
    annotate('text',5 ,-5, label= 'Time-indexed \n acyclic diagrams',size=2,col='black')+
    
    annotate('text',0.7 ,-0.4, label= 'Yes',size=1.5,col='black')+
    annotate('text',1.75,-1.4, label= 'Yes',size=1.5,col='black')+
    annotate('text',3.25 ,-2.9, label= 'Yes',size=1.5,col='black')+
    annotate('text',0.7 ,-4.4, label= 'Yes',size=1.5,col='black')+
    annotate('text',4.75 ,-4.4, label= 'Yes',size=1.5,col='black')+
    
    annotate('text',-0.7 ,-0.4, label= 'No',size=1.5,col=my_cols[2])+
    annotate('text',0.25 ,-1.4, label= 'No',size=1.5,col=my_cols[2])+
    annotate('text',0.75 ,-2.9, label= 'No',size=1.5,col=my_cols[2])+
    annotate('text',-0.7 ,-4.4, label= 'No',size=1.5,col=my_cols[2])+
    annotate('text',3.25 ,-4.4, label= 'No',size=1.5,col=my_cols[2])+
    
    annotate('text',1 ,-5.25, label= 'Method 3',fontface='bold',size=2,col=my_cols[1])+
    annotate('text',3 ,-5.25, label= 'Method 2',fontface='bold',size=2,col=my_cols[1])+
    annotate('text',5 ,-5.25, label= 'Method 1',fontface='bold',size=2,col=my_cols[1])+
    
    annotate('rect',xmin=0.5, xmax=1.5, ymin=-5.5, ymax=-4.75, alpha=0.2, fill=my_cols[1])+
    annotate('rect',xmin=2.5, xmax=3.5, ymin=-5.5, ymax=-4.75, alpha=0.2, fill=my_cols[1])+
    annotate('rect',xmin=4.5, xmax=5.5, ymin=-5.5, ymax=-4.75, alpha=0.2, fill=my_cols[1])+
    
    xlim(-1.5,5.5)+
    
    theme_dag()
  
  p
  
}
# framework_figure()



print('Functions successfully loaded')
