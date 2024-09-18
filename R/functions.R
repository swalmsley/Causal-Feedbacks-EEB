
# Custom functions to be called in '_targets.R'



# Save plot ---------------------------------------------------------------

save_figure <- function(path, w, h, call) {
  png(path, width=w, height=h, units='in', res=2001)
  print(call)
  dev.off()
}



# Make a basic plot -------------------------------------------------------

plot_xy <- function(df, xcol, ycol) {
  
  stopifnot("X variable not found" = xcol %in% names(df))
  stopifnot("Y variable not found" = ycol %in% names(df))
  
  g <- ggplot(df) + geom_point(aes(x=.data[[xcol]], y = .data[[ycol]]))
  
  stopifnot(class(g)[1]=="gg")
  
  return(g)
  
}


# # Lynx-hare function ------------------------------------------------------
# sim_lynx_hare <- function( n_steps , init , theta , dt) {
#   L <- rep(NA,n_steps)
#   H <- rep(NA,n_steps)
#   L[1] <- init[1]
#   H[1] <- init[2]
#   for ( i in 2:n_steps ) {
#     H[i] <- H[i-1] + dt*H[i-1]*( theta[1] - theta[2]*L[i-1] )
#     L[i] <- L[i-1] + dt*L[i-1]*( theta[3]*H[i-1] - theta[4] )
#   }
#   return( cbind(L,H) )
# }
# 
# library(rethinking)
# data(Lynx_Hare)
# theta <- c(0.5, 0.05, 0.025, 0.5)
# z <- sim_lynx_hare(1e4, init=c(4,30), theta, dt=0.002)
# 
# plot(z[,2], type='l', ylim=c(0, max(z[,2])), lwd=2, xaxt='n',
#      ylab='number (thousands)', xlab='')
# lines(z[,1], col=rangi2, lwd=2)
# mtext('time',1)

# 
# # Simulate time series of trait evolution ---------------------------------
# temporal_simulation <- function(n_timesteps,
#                                 n_individuals, 
#                                 a_to_b, 
#                                 b_to_a) {
#   df <- data.frame(
#     individual = integer(),
#     time = integer(),
#     a = numeric(),
#     b = numeric())
#   
#   for (i in 1:n_individuals) {
#     
#     # Initialize starting values of 0 for simplicity
#     a <- numeric(n_timesteps)
#     b <- numeric(n_timesteps)
#     a[1] <- 0
#     b[1] <- 0
#     
#     for (t in 2:n_timesteps) {
#       a[t] <- a[t-1] + b_to_a*b[t-1] + rnorm(1, mean = 0, sd = 0.1) 
#       b[t] <- b[t-1] + a_to_b*a[t-1] + rnorm(1, mean = 0, sd = 0.1) 
#     }
#     
#     # Append the results for this individual to the dataframe
#     df <- rbind(df, data.frame(
#       individual = rep(i, n_timesteps),
#       time = 1:n_timesteps,
#       a = a,
#       b = b))
#     
#   }
#   
#   df_long <- df %>% tidyr::pivot_longer(c(a, b),
#                                         names_to = "trait", values_to = "value")
#   return(data.table(df_long))
# }
# 


# Plot results of temporal trait simulation -------------------------------
plot_temporal <- function(simulation_results, title) {

  d <- simulation_results[trait=='a',,] # only plotting 'a'
  
  g <- ggplot(d, aes(x = time, y = value, colour = trait)) +
    geom_line(aes(group = interaction(individual, trait)), alpha = 1, linewidth = 0.25) +
    labs(title = paste(title),
         x = "Time Step", y = "A") +
    theme_classic() +
    ylim(-10,10)+
    theme(legend.position='none', title=element_text(size=10), plot.tag = element_text(size=12))
  
  return(g)

}





library(data.table)
library(ggplot2)
# Simultaneity  -----------------------------------------------------------
Example1_simultaneity <- function() {
  
  # set up parameters for simulation
  n_betas <- 50
  n <- 200 # number of observations
  list_bxy <- seq(from=-3,to=3,length.out=n_betas)
  list_byx <- seq(from=-3,to=3,length.out=n_betas)
  
  combinations <- data.table(as.data.frame(t(combn(c(list_bxy,list_byx),2))))
  colnames(combinations) <- c('true_xy','true_yx')
  
  # NOTE currently effects are always equivalent, need to create situations with different combinations
  # this is interesting though: coefficients seem to flatten around 0.6...
  
  # set up datatable to hold results
  d <- combinations[,iteration:=.I,]
  
  for (i in 1:nrow(d)) {
    
    # pull out 'true' effects
    effect_xy <- d[i,true_xy,]
    effect_yx <- d[i,true_yx,]
    
    # generate initial and final values of x and y
    x1 <- rnorm(n, 0, 1)
    y1 <- rnorm(n, 0, 1)
    
    iv_1 <- rnorm(n, 0, 1)
    
    x <- rnorm(n, x1 + effect_yx*y1, 1)
    y <- rnorm(n, y1 + effect_xy*x1, 1)
    
    # fill in output for plotting
    d[i,B_xy:=summary(lm(y~x))$coefficients[2],]
    #d[i,B_yx:=summary(lm(x~y))$coefficients[2],]
    
    #d[i,B_adj:=summary(lm(y~x + x1 + y1))$coefficients[2],]
    #d[i,B_adj:=summary(lm(y~ x + x1 + y1))]
    
    # Could be nice to show both "true direct effect" and "true total effect", if/when they differ?
    
    # Not sure I'm doing this properly
    # key question is what do we mean when we ask about effect of X on Y? do we mean Xnow and Ynow?
    # no of course not, we mean historicall: how have these traits influenced one another
    # like should I be setting the exposure to y1 -> x and x1 -> y? are these actually what we're after? 
    
    # plot is coming along but should make it so that true x->y and true y->x are visible so we can see effects (e.g., lattice plot)
    # especially when one of the effects is 0, in which case this might work, right? 
    
    # maybe I'm getting at making this into a basic panel model?
    
    # have I shown endogeneity bias though? Maybe if we pretend that x1 and y1 don't actually exist in the DAG, it's just X and Y
    
    # what is reverse causation bias, then?
    
    
  }

  # subset to look at x on y effects and coefficients only
  d <- d[,c('iteration','true_xy','B_xy','true_yx'),]

  g <- ggplot(d, aes(x=true_xy, y=B_xy)) +
    geom_point(alpha=0.1) + 
    #facet_wrap(~true_yx)+
    geom_abline(slope=1, intercept=0)+
    ylim(-2,2)
  g
  
}


# 
#   d <- simulation_results[trait=='a',,] # only plotting 'a'
#   
#   g <- ggplot(d, aes(x = time, y = value, colour = trait)) +
#     geom_line(aes(group = interaction(individual, trait)), alpha = 1, linewidth = 0.25) +
#     labs(title = paste(title),
#          x = "Time Step", y = "A") +
#     theme_classic() +
#     ylim(-10,10)+
#     theme(legend.position='none', title=element_text(size=10), plot.tag = element_text(size=12))
#   
#   return(g)
#   
# }


# n <- 1000
# 
# effect_xy <- 5
# effect_yx <- 0
# 
# # generate initial and final values of x and y
# x1 <- rnorm(n, 0, 1)
# y1 <- rnorm(n, 0, 1)
# 
# iv_1 <- rnorm(n, 0, 1)
# 
# # x <- rnorm(n, x1 + effect_yx*y1, 1)
# # y <- rnorm(n, y1 + effect_xy*x1, 1)
# 
# # x <- x1 + effect_yx * y1 + rnorm(n, 0, 1)
# # y <- y1 + effect_xy * x1 + rnorm(n, 0, 1)
# 
# x <- x1 + rnorm(n, 0, 2)   # Generating x with random noise
# y <- y1 + effect_xy * x1 + rnorm(n, 0, 1)  # Generating y with the effect of x1 and random noise
# 
# 
# summary(lm(y~x))
# summary(lm(x~y))
# 
# plot(x,y)
# 
# # first, why are we not getting the expected effect out of this?
# # getting exactly half the effect...
# 
# # works if we regress on x1
# 
# # maybe the tricky thing is simulating simultaneity?






# 
# # Number of observations
# n <- 100
# 
# # Generate variables
# z <- rnorm(n, 0, 1)  # Variable affecting x directly
# y <- rnorm(n,0,1) 
# 
# # Simultaneity between x and y
# x <- rnorm(n, 0, 1) + (-1 * y) + rnorm(n, 0, 1)  # x affects y positively
# y <- rnorm(n, 0, 1) + (0.3 * x) + rnorm(n, 0, 1)  # y affects x positively
# 
# # Adding direct effect of z on x
# x <- x + 0.7 * z + rnorm(n, 0, 1)
# 
# # Check correlations (optional)
# cor(x, y)  # Should show a positive correlation due to simultaneity
# cor(x, z)  # Should show a positive correlation due to direct effect
# cor(y, z)  # Should show little to no correlation
# 
# # Regression models
# summary(lm(y ~ x + z))  # y regressed on x and z
# summary(lm(x ~ y + z))  # x regressed on y and z
# 
# d <- data.frame(x=x,y=y,z=z)
# 
# # Trying an IV model
# library(brms)
# f1 <- bf(x ~ z)
# f2 <- bf(y ~ x)
# iv_model <- brm(f1 + f2 + set_rescor(TRUE),
#                 data=d, 
#                 family = 'Gaussian')
# 
# # examine model output
# iv_model
# 
# 
# 
# 
# # Generate variables
# set.seed(123)  # for reproducibility
# n <- 100
# z <- rnorm(n)
# x <- rnorm(n) + 0.7 * z + rnorm(n)
# y <- rnorm(n) + 0.3 * x + rnorm(n)
# 
# # IV model using `ivreg` from `AER` package
# library(AER)
# iv_model <- ivreg(y ~ x | z)
# summary(iv_model)
# 
# 
# 
# 
# library(tidyr)
# library(dplyr)
# 
# # make a standardizing function
# standardize <- function(x) {
#   (x - mean(x)) / sd(x)
# }
# 
# # simulate
# set.seed(73) 
# 
# n <- 500
# 
# dat_sim <-
#   tibble(u_sim = rnorm(n, mean = 0, sd = 1),
#          q_sim = sample(1:4, size = n, replace = T)) %>% 
#   mutate(e_sim = rnorm(n, mean = u_sim + q_sim, sd = 1)) %>% 
#   mutate(w_sim = rnorm(n, mean = u_sim + 0 * e_sim, sd = 1)) %>% 
#   mutate(w = standardize(w_sim),
#          e = standardize(e_sim),
#          q = standardize(q_sim))
# 
# dat_sim
# 
# 
# e_model <- bf(e ~ 1 + q)
# w_model <- bf(w ~ 1 + e)
# 
# b14.6 <-
#   brm(data = dat_sim, 
#       family = gaussian,
#       e_model + w_model + set_rescor(TRUE),
#       prior = c(# E model
#         prior(normal(0, 0.2), class = Intercept, resp = e),
#         prior(normal(0, 0.5), class = b, resp = e),
#         prior(exponential(1), class = sigma, resp = e),
#         
#         # W model
#         prior(normal(0, 0.2), class = Intercept, resp = w),
#         prior(normal(0, 0.5), class = b, resp = w),
#         prior(exponential(1), class = sigma, resp = w),
#         
#         # rho
#         prior(lkj(2), class = rescor)),
#       iter = 2000, warmup = 1000, chains = 4, cores = 4,
#       seed = 14)




# Plot simple DAG ---------------------------------------------------------


# plot_simple_dag <- function(dag_text) {
# 
#   dag <- dagitty(dag_text)
# 
#   p1 <- ggdag(dag,
#               layout='circle',
#               edge_type =,
#               node = FALSE,
#               text = FALSE) + theme_dag() +
#     geom_dag_text(color=c(c1,c2),aes(hjust='inward', vjust='inward'))
# 
# }
# 
# 
# 
# plot_simple_dag <- function(dag_text) {
# 
#   dag <- dagitty(dag_text)
# 
#   # Manually set the coordinates of the nodes
#   coordinates(dag) <- list(
#     x = c(A = 1, B = 2),
#     y = c(A = 1, B = 1)
#   )
# 
# 
#   p <- ggdag(dag)+
#     geom_dag_node(colour='white')+
#     geom_dag_text(colour='black')+
#     theme_dag()
# 
#   return(p)
# 
# }
# 
# 
# dag_text <- 'dag{A -> B}'
# 
# plot_simple_dag(dag_text)



# Plot reciprocal dag -----------------------------------------------------
plot_reciprocal_dag <- function(sign1, sign2) {
  
  dag_text <- 'dag{
  A1 -> A2;
  A1 -> B2;
  B1 -> B2;
  B1 -> A2;
  A2 -> A3;
  A2 -> B3;
  B2 -> B3;
  B2 -> A3;
  }'
  
  dag <- dagitty(dag_text)
  
  # Manually set the coordinates of the nodes
  coordinates(dag) <- list(
    x = c(A1 = 1, A2 = 2, A3 = 3, B1 = 1, B2 = 2, B3 = 3),
    y = c(A1 = 2, A2 = 2, A3 = 2, B1 = 1, B2 = 1, B3 = 1)
  )
  
  p <- ggdag(dag) +
    geom_dag_node(colour='white') +
    geom_dag_text(
      label = expression(A[t], A[t+1], A[t + ...], B[t], B[t+1], B[t + ...]),
      parse = TRUE,
      size=4,
      colour = 'black'
    ) +
    annotate('text', x=1.55, y=1.7,label=paste('(',sign1,')',sep=''),size=3,col='red')+
    annotate('text', x=1.55, y=1.3,label=paste('(',sign2,')',sep=''),size=3,col='red')+
    annotate('text', x=2.55, y=1.7,label=paste('(',sign1,')',sep=''),size=3,col='red')+
    annotate('text', x=2.55, y=1.3,label=paste('(',sign2,')',sep=''),size=3,col='red')+
    ylim(0.5,2.5)+
    theme_dag()
  
  return(p)
  
  
}



# Plot autocorrelated dag -------------------------------------------------
plot_auto_dag <- function() {
  
  dag_text <- 'dag{
  A1 -> A2;
  B1 -> B2;
  A2 -> A3;
  B2 -> B3;
  }'
  
  dag <- dagitty(dag_text)
  
  # Manually set the coordinates of the nodes
  coordinates(dag) <- list(
    x = c(A1 = 1, A2 = 2, A3 = 3, B1 = 1, B2 = 2, B3 = 3),
    y = c(A1 = 2, A2 = 2, A3 = 2, B1 = 1, B2 = 1, B3 = 1)
  )
  
  p <- ggdag(dag) +
    geom_dag_node(colour='white') +
    geom_dag_text(
      label = expression(A[t], A[t+1], A[t + ...], B[t], B[t+1], B[t + ...]),
      parse = TRUE,
      size=4,
      colour = 'black'
    ) +
    ylim(0.5,2.5)+
    theme_dag()
  
  return(p)
  
  
}



# Plot single auto --------------------------------------------------------
plot_single_auto_dag <- function() {
  
  dag_text <- 'dag{
  A1 -> A2;
  A2 -> A3;
  }'
  
  dag <- dagitty(dag_text)
  
  # Manually set the coordinates of the nodes
  coordinates(dag) <- list(
    x = c(A1 = 1, A2 = 2, A3 = 3),
    y = c(A1 = 1.5, A2 = 1.5, A3 = 1.5))
  
  p <- ggdag(dag) +
    geom_dag_node(colour='white') +
    geom_dag_text(
      label = expression(A[t], A[t+1], A[t + ...]),
      parse = TRUE,
      size=4,
      colour = 'black'
    ) +
    ylim(0.5,2.5)+
    theme_dag()
  
  return(p)
  
  
}









instrumental_variable_plot <- function(fit1, fit2, r1, r2) {
  
  # Set up colors
  # my_cols <- viridis_pal(option='A')(12)[c(4,8)]
  # my_cols <- viridis_pal(option='A')(12)[c(4,8)]
  my_cols <- c('steelblue','red')
  
  
  # # Time-dependent DAG
  # 
  # dag_text <- 'dag{
  # I1 -> A0;
  # I2 -> B0;
  # A0 -> A;
  # A0 -> B;
  # B0 -> B;
  # B0 -> A;
  # }'
  # 
  # dag <- dagitty(dag_text)
  # 
  # # Manually set the coordinates of the nodes
  # coordinates(dag) <- list(
  #   x = c(I1 = 0.5, I2 = 2.5, A0 = 1, A = 1, B0 = 2, B = 2),
  #   y = c(I1 = 2, I2 = 2, A0 = 2, A = 1, B0 = 2, B = 1)
  # )
  # 
  # p1 <- ggdag(dag) +
  #   geom_dag_node(colour='white') +
  #   #geom_dag_edges(edge_colour=c('black','blue','blue','black','black','black'))+
  #   geom_dag_text(
  #     label = expression(A[observed], A[ t-1], B[observed], B[ t-1], Instrument[1], Instrument[2]),
  #     parse = TRUE,
  #     size=c(5,5,5,5,3.5,3.5),
  #     colour = c('black','grey50','black','grey50','black','black')
  #   ) +
  #   # annotate('text', x=1.55, y=1.7,label=paste('(',sign1,')',sep=''),size=3,col='red')+
  #   # annotate('text', x=1.55, y=1.3,label=paste('(',sign2,')',sep=''),size=3,col='red')+
  #   #ylim(0.5,2.5)+
  #   theme_dag()
  # 
  # #p1
  
  
  
  # Oversimplified DAG
  dag_text <- 'dag{
  A -> B;
  B -> A;
  I1 -> A;
  I2 -> B;
  }'
  
  dag <- dagitty(dag_text)
  
  ggdag(dag)
  
  # Manually set the coordinates of the nodes
  coordinates(dag) <- list(
    x = c(A = 1, B = 2, I1=0, I2=3),
    y = c(A = 1, B = 1, I1=1, I2=1)
  )
  
  # 
  # pFull <- ggdag(dag) +
  #   geom_dag_node(colour='white') +
  #   geom_dag_edges(edge_width=rep(0.6,400),edge_color=c(rep('white',100),rep('white',100),rep('black',100),rep('black',100)))+
  #   geom_dag_edges_arc(curvature = 0.00000001,edge_width=rep(0.6,400),edge_color=c(rep(my_cols[1],100),rep(my_cols[2],100),rep('white',100),rep('white',100)))+
  #   geom_dag_text(
  #     parse=FALSE,
  #     size=c(7,7,7,6),
  #     colour = c('black','black','black','black'),
  #     fontface=c('bold','bold','bold','bold')
  #   ) +
  #   theme_dag()+
  #   theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
  # 
  # pFull
  
  
  # p0 <- ggdag(dag) +
  #   geom_dag_node(colour='white') +
  #   geom_dag_edges(edge_width=c(rep(0.6,400)),edge_color=c(rep('white',400)))+
  #   geom_dag_edges_arc(curvature = 0.00000001,edge_width=rep(0.6,400),
  #                      #edge_linetype=c(rep('dashed',400)),
  #                      edge_color=c(rep(my_cols[1],100),rep(my_cols[2],100),rep('white',100),rep('white',100)))+
  #   geom_dag_text(
  #     parse=FALSE,
  #     size=c(7,7,7,6),
  #     colour = c('black','black','white','white'),
  #     fontface=c('bold','bold','bold','plain')
  #   ) +
  #   theme_dag()
  # #theme(plot.margin = unit(c(2,2,2,2), "cm"))
  # 
  # p0
  
  
  p1 <- ggdag(dag) +
    geom_dag_node(colour='white') +
    geom_dag_edges(edge_width=c(rep(0.6,200),rep(1,100),rep(0.6,100)),edge_color=c(rep('white',100),rep('white',100),rep('black',100),rep('white',100)))+
    geom_dag_edges_arc(curvature = 0.00000001,edge_width=c(rep(1,100),rep(0.6,300)),edge_color=c(rep(my_cols[1],100),rep('black',100),rep('white',100),rep('white',100)))+
    geom_dag_text(
      label=c('Y','X','I.1','I.2'),
      parse=FALSE,
      size=c(5,5,5,5),
      colour = c('black','black','black','white'),
      fontface=c('bold','bold','bold','plain')
    ) +
    theme_dag()+
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
  
  p1

  p2 <- ggdag(dag) +
    geom_dag_node(colour='white') +
    geom_dag_edges(edge_width=c(rep(0.6,200),rep(0.6,100),rep(1,100)),edge_color=c(rep('white',100),rep('white',100),rep('white',100),rep('black',100)))+
    geom_dag_edges_arc(curvature = 0.00000001,edge_width=c(rep(0.6,100),rep(1,100),rep(0.6,200)),edge_color=c(rep('black',100),rep(my_cols[2],100),rep('white',100),rep('white',100)))+
    geom_dag_text(
      label=c('Y','X','I.1','I.2'),
      parse=FALSE,
      size=c(5,5,5,5),
      colour = c('black','black','white','black'),
      fontface=c('bold','bold','plain','bold')
    ) +
    theme_dag()+
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
    
  p2
  # 
  # p2 <- ggdag(dag) +
  #   geom_dag_node(colour='white') +
  #   geom_dag_edges(edge_width=c(rep(0.6,200),rep(0.6,100),rep(1.5,100)),edge_color=c(rep('white',100),rep('white',100),rep('black',100),rep('black',100)))+
  #   geom_dag_edges_arc(curvature = 0.1,edge_width=c(rep(0.6,100),rep(1.5,100),rep(0.6,200)),edge_color=c(rep('black',100),rep(my_cols[2],100),rep('white',100),rep('white',100)))+
  #   geom_dag_text(
  #     label=c('MPA','Biodiversity','Political Will','Nutrients'),
  #     parse=FALSE,
  #     size=c(3,3,3,3),
  #     colour = c('black','black','black','black'),
  #     fontface=c('bold','bold','plain','bold')
  #   ) +
  #   theme_dag()+
  #   theme(plot.margin = unit(c(3,0.5,3,0.5), "cm"))
  # 
  # p2
  
  pm1 <- density_plot(fit1, 1, 'b_A_B', 'Estimated effect', 1.5, my_cols, 1) +
    geom_density(inherit.aes=FALSE,data=as_draws_df(r1),aes(x=get('b_B')),
                 color='white',fill='grey',alpha=0.35) +
    annotate('text', x=-0.75,color='grey', y=2,label='Naive regression',size=4)+
    annotate('text', x=1.6,color=my_cols[1], y=1,label='IV analysis',size=4)
  pm1
  
  pm2 <- density_plot(fit2, -1, 'b_B_A', 'Estimated effect', -1.5, my_cols, 2) +
    geom_density(inherit.aes=FALSE,data=as_draws_df(r2),aes(x=get('b_A')),
                 color='white',fill='grey',alpha=0.35) +
    annotate('text', x=1,color='grey', y=2,label='Naive regression',size=4)+
    annotate('text', x=-1.5,color=my_cols[2], y=1,label='IV analysis',size=4)
  
  # fit <- iv_model
  # true_effect <- -1
  # var_name <- 'b_violenceobserved_development_applied'
  # var_label <- 'Estimated effect'
  # xloc <- -1.4
  # colNum <- 1
  
  # standard regression for comparison

  pr1 <- density_plot(r1, -1, 'b_B', 'Estimated effect', -1.6, my_cols, 2)
  pr2 <- density_plot(r2, 1, 'b_A', 'Estimated effect', -1.6, my_cols, 2)
  
  
  pVoid <- ggplot() + theme_void()
  
  # Density plot
  p3 <- ggplot(data.frame(v = rnorm(100)),aes(x=v))+
    geom_density()
  
  plot <- ((p1 / p2) | (pm1 / pm2)) + plot_annotation(tag_levels='A')
  return(plot + plot_layout(widths=c(1,2)))
  #plot <- (pVoid | pVoid / pFull / pVoid | pVoid) | (p3 / p3) 
  # plot + plot_layout(widths=c(1,3,1,10))
  #return(plot)

  # plot <- (p1 / p2) | (pm1 / pm2)
  # return(plot + plot_layout(widths=c(1,2)))
  
}
# fit1 <- tar_read(IV_m1)
# fit2 <- tar_read(IV_m2)
# r1 <- tar_read(reg_m1)
# r2 <- tar_read(reg_m2)


density_plot <- function(fit, true_effect, var_name, var_label, xloc, my_cols, colNum) {

  g <- ggplot(data = as_draws_df(fit), aes(x=get(var_name))) + 
    geom_density(color='white', fill=my_cols[colNum], alpha=0.5) +
    xlim(-2.5,2.5)+
    geom_vline(xintercept=true_effect,linetype='dashed',linewidth=0.5)+
    labs(x=paste(var_label),y='Posterior density')+
    annotate('text', x=xloc, y=2.5,label='True effect',size=4)+
    theme_classic()
  g
  
}
# fit <- iv_model
# true_effect <- -1
# var_name <- 'b_violenceobserved_development_applied'
# var_label <- 'Estimated effect'
# xloc <- -1.4
# colNum <- 1







create_ODE_example <- function() {
  
  
  # Define the ODE system
  ode_system <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
      dx <- a*X - b*X*Y
      dy <- c*Y*X - d*Y
      return(list(c(dx, dy)))
    })
  }
  
  # Set initial guess for parameters
  initial_params <- c(a = 0.1, b = 0.1, c = 0.1, d = 0.1)
  
  # Simulated data (with added noise)
  set.seed(123)
  times <- seq(0, 50, by = 0.1)
  true_parameters <- c(a = 1.1, b = 0.4, c = 0.1, d = 0.4)
  state <- c(X = 1, Y = 1)
  simulated_data <- ode(y = state, times = times, func = ode_system, parms = true_parameters)
  simulated_data <- as.data.frame(simulated_data)  # Ensure data is a data.frame
  simulated_data$X <- simulated_data$X + rnorm(nrow(simulated_data), sd = 1)
  simulated_data$Y <- simulated_data$Y + rnorm(nrow(simulated_data), sd = 1)
  
  
  
  
  lotka_model <- odemodel(
    name="Lotka Volterra model",
    model=list(
      u ~ a * u - b * u * v,
      v ~ c * u * v - d * v
    ),
    observation=list(
      X ~ dpois(lambda=u),
      Y ~ dpois(lambda=v)
    ),
    initial=list(
      u ~ u0,
      v ~ v0
    ),
    par=c("a", "b", "c", "d", "u0", "v0")
  )
  
  harestart <- c(a = 1.1,
                 b = 0.4, 
                 c = 0.1, 
                 d = 0.4, 
                 u0 = 1,
                 v0 = 1)
  
  # harefit <- fitode(lotka_model, data=simulated_data,
  #                   start=harestart,
  #                   tcol="time")
 
  proposal.vcov <- matrix(0, 6, 6)
  diag(proposal.vcov) <- c(1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4)
  proposal.vcov
  
  harefitMCMC <- fitodeMCMC(lotka_model, data=simulated_data,
                            start=harestart,
                            proposal.vcov = proposal.vcov,
                            tcol="time")
  
  return(harefitMCMC)
  
}










ODEs_plot <- function(fitMCMC) {
  
  
  # Set up colors

  my_cols <- c('steelblue','red')
  

  # Oversimplified DAG
  dag_text <- 'dag{
  A -> B;
  B -> A;
  }'
  
  dag <- dagitty(dag_text)
  
  ggdag(dag)
  
  # Manually set the coordinates of the nodes
  coordinates(dag) <- list(
    x = c(A = 1, B = 2),
    y = c(A = 1, B = 1)
  )

  
  p1 <- ggdag(dag) +
    geom_dag_node(colour='white') +
    geom_dag_edges(edge_color='white')+
    geom_dag_edges_arc()+
    geom_dag_text(
      label=c('Prey','Predator'),
      parse=FALSE,
      size=c(2,2),
      colour = c('white','white')
    ) +
    annotate('text', x=1, y=1,label='Prey',color=my_cols[1],size=3)+
    annotate('text', x=1.9, y=1,label='Predator',color=my_cols[2],size=3)+
    theme_dag()
  p1
  # density plots
  
  # # extract draws
  # draws <- data.frame(fitMCMC@mcmc[[1]])
  # draws$draw <- 1:1000
  # 
  # g1 <- ggplot(data = draws, aes(x=b)) + 
  #   geom_density(color='white', fill=my_cols[1], alpha=0.5) +
  #   xlim(-1,1)+
  #   geom_vline(xintercept=0.4,linetype='dashed',linewidth=0.5)+
  #   labs(x=paste('label'),y='Posterior density')+
  #   annotate('text', x=0.2, y=20,label='True effect',size=4)+
  #   theme_classic()
  # g1
  
  # g1 <- plot(fitMCMC)

  # plot <- ((p1 / p2) | (pm1 / pm2)) + plot_annotation(tag_levels='A')
  # return(plot + plot_layout(widths=c(1,2)))
  #plot <- (pVoid | pVoid / pFull / pVoid | pVoid) | (p3 / p3) 
  # plot + plot_layout(widths=c(1,3,1,10))
  #return(plot)
  
  
  pred <- predict(fitMCMC,0.95,simplify=TRUE)
  data <- data.table(fitMCMC@data)
  dataLong <- melt(data, id.vars='times')
    
  
  g1 <- ggplot(dataLong, aes(x=times, y=value, group=variable, color=variable))+
    geom_point(alpha=0.3, size=1.25)+
    geom_line(inherit.aes = FALSE, data=pred$X, aes(x=times, y=estimate), alpha=1, linewidth=1.1, color=my_cols[1])+
    geom_line(inherit.aes = FALSE, data=pred$Y, aes(x=times, y=estimate), alpha=1, linewidth=1.1, color=my_cols[2])+
    labs(x='Time', y='Population abundance')+
    annotate('text', x=54, y=14,label='Prey',color=my_cols[1],size=3)+
    annotate('text', x=56, y=3,label='Predator',color=my_cols[2],size=3)+
    theme_classic()+
    scale_color_manual(values=c(my_cols))+
    xlim(0,60)+
    theme(legend.position='none', axis.title=element_text(size=8))

  
  g2 <- ggplot() +
    annotate('text', x=0, y=0.4, label=expression(frac(dy, dt) == Predator(r[Predator]*Prey - beta[Mortality])),size=1.75, col=my_cols[2])+
    annotate('text', x=0, y=-0.4, label=expression(frac(dx, dt) == Prey(r[Prey] - beta[Consumption]*Predator)), size=1.75, col=my_cols[1])+
    ylim(-1,1)+
    xlim(-10,10)+
    theme_no_axes()+ 
    theme(panel.border = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"))
  g2
  
  
  plot <- ((p1 | g2) / g1) + plot_annotation(tag_levels='A')
  # plot <- (p1 | g1) + plot_layout(widths=c(1,4))
  
  return(plot)

  # return(plot + plot_layout(widths=c(1,2)))
  
}
# fitMCMC <- tar_read(fitMCMC)




density_plot <- function(fit, true_effect, var_name, var_label, xloc, my_cols, colNum) {
  
  g <- ggplot(data = as_draws_df(fit), aes(x=get(var_name))) + 
    geom_density(color='white', fill=my_cols[colNum], alpha=0.5) +
    xlim(-2.5,2.5)+
    geom_vline(xintercept=true_effect,linetype='dashed',linewidth=0.5)+
    labs(x=paste(var_label),y='Posterior density')+
    annotate('text', x=xloc, y=2.5,label='True effect',size=4)+
    theme_classic()
  g
  
}
# fit <- iv_model
# true_effect <- -1
# var_name <- 'b_violenceobserved_development_applied'
# var_label <- 'Estimated effect'
# xloc <- -1.4
# colNum <- 1














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




regression_sim_reverse <- function() {
  
  n <- 50
  beta <- seq(-2,2,by=0.1)
  output <- data.table(index=(1:length(beta)))
  
  for (i in 1:length(beta)){
    
    x <- rnorm(n,0,1) # simulate x values
    y <- beta[i]*x + rnorm(n) # simulate y values
    
    m <- lm(x ~ y) # fit simple regression
    
    # save results to output 
    output[i, true_xy:=beta[i],]
    output[i, coefficient:=m$coefficients[2],]
    output[i, lowCI:=confint(m)[2],]
    output[i, highCI:=confint(m)[4],]
    
  }
  return(output)
}




sim_temporal_data <- function() {
  
  # Initialize variables
  n_time_steps <- 25
  n_individuals <- 25
  a_to_b <- -0.2
  b_to_a <- 0.2
  
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
  
  df <- data.table(df)
  df <- df[individual==1]
  
  df <-
    df %>%
    mutate(a_1 = lag(a),
           b_1 = lag(b))
  
  return(df)
  
}


fit_panel_model <- function(df) {
  
  f1 <- bf(b ~ 0 + Intercept + b_1 + a_1)
  f2 <- bf(a ~ 0 + Intercept + a_1 + b_1)
  
  m <-
    brm(data = df,
        family = 'Gaussian',
        f1 + f2 + set_rescor(FALSE),
        prior = c(prior(normal(log(10), 1), class = b, resp = b, coef = Intercept),
                  prior(normal(log(10), 1), class = b, resp = a, coef = Intercept),
                  prior(normal(0, 0.5), class = b, resp = b),
                  prior(normal(0, 0.5), class = b, resp = a),
                  prior(exponential(1), class = sigma, resp = b),
                  prior(exponential(1), class = sigma, resp = a)),
        iter = 2000, warmup = 1000, chains = 4, cores = 4)
  
  return(m)
  
}




panel_model_figure <- function(fit, pf1, pf2) {
  
  my_cols <- c('steelblue','red')

  # Oversimplified DAG
  dag_text <- 'dag{
  A1 -> A2;
  A1 -> B2;
  B1 -> B2;
  B1 -> A2;
  A2 -> A3;
  A2 -> B3;
  B2 -> B3;
  B2 -> A3;
  A3 -> A4;
  A3 -> B4;
  B3 -> B4;
  B3 -> A4;
  }'
  
  dag <- dagitty(dag_text)
  
  # Manually set the coordinates of the nodes
  coordinates(dag) <- list(
    x = c(A1=1, B1=1, A2=2, B2=2, A3=3, B3=3, A4=4, B4=4),
    y = c(A1=1, B1=2, A2=1, B2=2, A3=1, B3=2, A4=1, B4=2)
  )
  
  p <- ggdag(dag)+
    geom_dag_edges(edge_color='white')+
    geom_dag_node(colour='white') +
    geom_dag_edges(edge_color=c(rep('black',100),
                                rep(my_cols[1],100),
                                rep('black',100),
                                rep(my_cols[1],100),
                                rep('black',100),
                                rep(my_cols[1],100),
                                rep(my_cols[2],100),
                                rep('black',100),
                                rep(my_cols[2],100),
                                rep('black',100),
                                rep(my_cols[2],100),
                                rep('black',100)),
                   edge_width=c(rep(0.5,100),
                                rep(1,100),
                                rep(0.5,100),
                                rep(1,100),
                                rep(0.5,100),
                                rep(1,100),
                                rep(1,100),
                                rep(0.5,100),
                                rep(1,100),
                                rep(0.5,100),
                                rep(1,100),
                                rep(0.5,100)))+
    geom_dag_text(colour='black', size=5,
                  label=c(expression(Y[t]),
                          expression(Y[t+1]),
                          expression(Y[t+2]),
                          expression(Y[t+3]),
                          expression(X[t]),
                          expression(X[t+1]),
                          expression(X[t+2]),
                          expression(X[t+3])))+
    theme_dag()
  p
  
  
  p2 <- density_plot(fit, -0.2, 'b_b_a_1', 'Estimated effect', 1.5, my_cols, 1) + xlim(-0.5,0.5)+
     geom_density(inherit.aes=FALSE,data=as_draws_df(pf1),aes(x=get('b_b')),color='white',fill='grey',alpha=0.35)+
    annotate('text', x=1,color='grey', y=2,label='Naive regression',size=4)+
    annotate('text', x=-1.5,color=my_cols[2], y=1,label='IV analysis',size=4) +
    theme(axis.title = element_text(size=8))
    
    
    
    # geom_density(inherit.aes=FALSE, data=as_draws_df(pf1),aes(x=get('b_A')),
  
  p3 <- density_plot(fit, 0.2, 'b_a_b_1', 'Estimated effect', 1.5, my_cols, 2) + xlim(-0.5,0.5)+
    geom_density(inherit.aes=FALSE,data=as_draws_df(pf2),aes(x=get('b_a')),color='white',fill='grey',alpha=0.35)+
    annotate('text', x=1,color='grey', y=2,label='Naive regression',size=4)+
    annotate('text', x=-1.5,color=my_cols[2], y=1,label='IV analysis',size=4)+
    theme(axis.title = element_text(size=8))
  
  plot <- (p/p2/p3) + plot_annotation(tag_levels='A') 
  
  return(plot)
  
}
# fit <- tar_read(pfit)
# pf1 <- tar_read(panel_regular_1)
# pf2 <- tar_read(panel_regular_2)

# two types of simple models to consider: multiple sites, take a pinpoint sample of each, look at relationship (cross-sectional)
# take longitudinal data and fit standard regression 
# key question then is if IV would work for multiple-site cross-sectional situation




panel_model_confound_figure <- function(fit, pf1, pf2) {
  
  my_cols <- c('steelblue','red')
  
  # Oversimplified DAG
  dag_text <- 'dag{
  A1 -> A2;
  A1 -> B2;
  B1 -> B2;
  B1 -> A2;
  A2 -> A3;
  A2 -> B3;
  B2 -> B3;
  B2 -> A3;
  A3 -> A4;
  A3 -> B4;
  B3 -> B4;
  B3 -> A4;
  Z1 -> A2;
  Z1 -> B2;
  Z2 -> A3;
  Z2 -> B3;
  Z3 -> A4;
  Z3 -> B4;
  }'
  
  dag <- dagitty(dag_text)
  
  # Manually set the coordinates of the nodes
  coordinates(dag) <- list(
    x = c(A1=1, B1=1, A2=2, B2=2, A3=3, B3=3, A4=4, B4=4, Z1=1.5, Z2=2.5, Z3=3.5),
    y = c(A1=1, B1=2, A2=1, B2=2, A3=1, B3=2, A4=1, B4=2, Z1=1.5, Z2=1.5, Z3=1.5)
  )
  
  p <- ggdag(dag)+
    geom_dag_edges(edge_color='white')+
    geom_dag_node(colour='white') +
    geom_dag_edges(edge_color=c(rep('black',100),
                                rep(my_cols[1],100),
                                rep('black',100),
                                rep(my_cols[1],100),
                                rep('black',100),
                                rep(my_cols[1],100),
                                rep(my_cols[2],100),
                                rep('black',100),
                                rep(my_cols[2],100),
                                rep('black',100),
                                rep(my_cols[2],100),
                                rep('black',100),
                                rep('black',600)),
                   edge_width=c(rep(0.5,100),
                                rep(1,100),
                                rep(0.5,100),
                                rep(1,100),
                                rep(0.5,100),
                                rep(1,100),
                                rep(1,100),
                                rep(0.5,100),
                                rep(1,100),
                                rep(0.5,100),
                                rep(1,100),
                                rep(0.5,100),
                                rep(0.5,600)))+
    geom_dag_text(colour='black', size=5,
                  label=c(expression(Y[t]),
                          expression(Y[t+1]),
                          expression(Y[t+2]),
                          expression(Y[t+3]),
                          expression(X[t]),
                          expression(X[t+1]),
                          expression(X[t+2]),
                          expression(X[t+3]),
                          expression(Z[1]),
                          expression(Z[2]),
                          expression(Z[3])))+
    theme_dag()
  p

  return(p)
  
}
# fit <- tar_read(pfit)
# pf1 <- tar_read(panel_regular_1)
# pf2 <- tar_read(panel_regular_2)

# two types of simple models to consider: multiple sites, take a pinpoint sample of each, look at relationship (cross-sectional)
# take longitudinal data and fit standard regression 
# key question then is if IV would work for multiple-site cross-sectional situation










# output <- regression_sim_x_to_y()
# output <- regression_sim_reverse()
# 
# my_cols <- 'steelblue'
# 
# ggplot(output, aes(x=true_xy, y=coefficient))+
#   #geom_smooth(method='lm',se=F, color='black', linewidth=0.75)+
#   geom_errorbar(aes(ymin=lowCI, ymax=highCI), width=0, linewidth=0.75, color='grey')+
#   geom_point(color=my_cols[1])+
#   geom_abline(intercept=0, slope=1)+
#   labs(x='True causal effect (simulated)', y=bquote(paste('Regression coefficient (', beta[1], ')')))+
#   xlim(-2,2)+ylim(-2,2)+
#   annotate('text', x=0.75, y=2,label=expression(x == beta[0] + beta[1]*y + epsilon),size=8,col=my_cols[1])+
#   theme_classic()
# 



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
# 
# output <- regression_sim_reciprocal()
# 
# my_cols[2] <- 'red'
# g <- ggplot(output, aes(x=beta_xy, y=m_y_coefficient))+
#   
#   #geom_smooth(method='lm',se=F, color='black', linewidth=0.75)+
#   #geom_errorbar(aes(ymin=m_y_lowCI, ymax=m_y_highCI), width=0, linewidth=0.75, color='grey')+
#   geom_point(color=my_cols[1],alpha=0.5)+
# 
#   #geom_smooth(method='lm',se=F, color='black', linewidth=0.75)+
#   #geom_errorbar(inherit.aes=FALSE, aes(x=beta_yx, y=m_x_coefficient, ymin=m_y_lowCI, ymax=m_y_highCI), width=0, linewidth=0.75, color='grey')+
#   geom_point(inherit.aes=FALSE, aes(x=beta_yx, y=m_x_coefficient),color=my_cols[2],alpha=0.5)+
#   
#   geom_abline(intercept=0, slope=1)+
#   labs(x='True causal effect (simulated)', y=bquote(paste('Regression coefficient (', beta[1], ')')))+
#   xlim(-2,2)+ylim(-2,2)+
#   annotate('text', x=-1, y=1.25,label=expression(y == beta[0] + beta[1]*x + epsilon),size=6,col=my_cols[1])+
#   annotate('text', x=-1, y=1,label=expression(x == beta[0] + beta[1]*y + epsilon),size=6,col=my_cols[2])+
#   
#   theme_classic()
# 
# 
# g


figure_2 <- function(){
  
  # SETUP 
  my_cols <- c('steelblue','red')
  
  # PANEL A
  output <- regression_sim_x_to_y()
  
  A <- ggplot(data=output, aes(x=true_xy, y=coefficient))+
    #geom_smooth(method='lm',se=F, color='black', linewidth=0.75)+
    geom_errorbar(aes(ymin=lowCI, ymax=highCI), width=0, linewidth=0.75, color='grey')+
    geom_point(color=my_cols[1])+
    geom_abline(intercept=0, slope=1)+
    labs(x='True causal effect (simulated)', y=bquote(paste('Regression coefficient (', beta[1], ')')))+
    xlim(-2,2)+ylim(-2,2)+
    annotate('text', x=-0.8, y=1.25,label=expression(y == beta[0] + beta[1]*x + epsilon),size=5,col=my_cols[1])+
    annotate('text', x=-0.8, y=1.75,label=expression(y %<-% x),size=5,col='black')+
    theme_classic()
  
  
  
  # PANEL B
  output <- regression_sim_reverse()
  
  B <- ggplot(output, aes(x=true_xy, y=coefficient))+
    #geom_smooth(method='lm',se=F, color='black', linewidth=0.75)+
    geom_errorbar(aes(ymin=lowCI, ymax=highCI), width=0, linewidth=0.75, color='grey')+
    geom_point(color=my_cols[2])+
    geom_abline(intercept=0, slope=1)+
    labs(x='True causal effect (simulated)', y=bquote(paste('Regression coefficient (', beta[1], ')')))+
    xlim(-2,2)+ylim(-2,2)+
    annotate('text', x=-0.8, y=1.25,label=expression(x == beta[0] + beta[1]*y + epsilon),size=5,col=my_cols[2])+
    annotate('text', x=-0.8, y=1.75,label=expression(y %<-% x),size=5,col='black')+
    theme_classic()
  
  
  
  
  # PANEL B
  output <- regression_sim_reverse()
  
  B <- ggplot(output, aes(x=true_xy, y=coefficient))+
    #geom_smooth(method='lm',se=F, color='black', linewidth=0.75)+
    geom_errorbar(aes(ymin=lowCI, ymax=highCI), width=0, linewidth=0.75, color='grey')+
    geom_point(color=my_cols[2])+
    geom_abline(intercept=0, slope=1)+
    labs(x='True causal effect (simulated)', y=bquote(paste('Regression coefficient (', beta[1], ')')))+
    xlim(-2,2)+ylim(-2,2)+
    annotate('text', x=-0.8, y=1.25,label=expression(x == beta[0] + beta[1]*y + epsilon),size=5,col=my_cols[2])+
    annotate('text', x=-0.8, y=1.75,label=expression(y %<-% x),size=5,col='black')+
    theme_classic()
  
  
  
  
  
  # PANEL C
  output <- regression_sim_reciprocal()
  
  C <- ggplot(output, aes(x=beta_xy, y=m_y_coefficient))+
    
    #geom_smooth(method='lm',se=F, color='black', linewidth=0.75)+
    #geom_errorbar(aes(ymin=m_y_lowCI, ymax=m_y_highCI), width=0, linewidth=0.75, color='grey')+
    geom_point(color=my_cols[1],alpha=0.5)+
    
    #geom_smooth(method='lm',se=F, color='black', linewidth=0.75)+
    #geom_errorbar(inherit.aes=FALSE, aes(x=beta_yx, y=m_x_coefficient, ymin=m_y_lowCI, ymax=m_y_highCI), width=0, linewidth=0.75, color='grey')+
    geom_point(inherit.aes=FALSE, aes(x=beta_yx, y=m_x_coefficient),color=my_cols[2],alpha=0.5)+
    
    geom_abline(intercept=0, slope=1)+
    labs(x='True causal effect (simulated)', y=bquote(paste('Regression coefficient (', beta[1], ')')))+
    xlim(-2,2)+ylim(-2,2)+
    annotate('text', x=-0.8, y=1.25,label=expression(y == beta[0] + beta[1]*x + epsilon),size=5,col=my_cols[1])+
    annotate('text', x=-0.8, y=1,label=expression(x == beta[0] + beta[1]*y + epsilon),size=5,col=my_cols[2])+
    
    annotate('text', x=-0.8, y=1.75,label=expression(y %<->% x),size=5,col='black')+
    
    theme_classic()
  
  
  # COMBINE PANELS
  
  figure <- (A | B | C) + plot_annotation(tag_levels = c('A'))
  return(figure)
  
}







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
  #ggdag(dag)


  p <- ggdag(dag)+
    geom_dag_edges(edge_color='white')+

    geom_dag_node(colour='white') +
    # geom_dag_edges_arc(curvature=0.4,
    #                    edge_color=c(rep('white',1000),
    #                                 rep('black',100)))+
    # geom_dag_edges(edge_color=c(rep('black',1000),
    #                             rep('white',100)))+
    # geom_dag_text(color='black',
    #               size=2,
    #               label=c('Could Y be a cause of X?',
    #                       'Proceed to \n infer causal effect, \n assuming given DAG.',
    #                       'Are the \n directions of \n causation between \n X and Y relevant \n to the scientific question?',
    #                       'Interpret effect \n of Y on X \n as non-causal \n correlation',
    #                       'Are there \n repeated measurements \n of Y and X \n across time?s',
    #                       'Are there suitable \n instruments for X?',
    #                       'Does the \n hypothesized effect \n of X on Y \n occur at discrete, \n separable time \n periods?',
    #                       'Not possible \n to identify \n direct causal \n effects',
    #                       'Does the hypothesized \n effect of X on \n Y occur at discrete, \n separable time periods',
    #                       'Continuous time \n methods',
    #                       'Time-indexed \n acyclic diagrams'))+
    #xlim(-4,6)+
  
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
  #return(p)
  
  # p + geom_dag_label_repel(aes(label = 'test'), colour = "white", show.legend = FALSE) 

  
}
# framework_figure()












print('Functions successfully loaded')
