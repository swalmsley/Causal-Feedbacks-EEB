
library(targets)

# Source
tar_source('R') # will do all in 'R' folder

# Options
tar_option_set()

# Seed
tar_option_set(seed = 12345)

# Variables
c1 <- col.desat( "blue" , 0.8 )
c2 <- col.desat( "red" , 0.8 )



# Targets
list(
  
  # CREATE FIGURE 2 - Bidirectionality and simple regression
  tar_target(Figure2, save_figure('./Manuscript/Figures/Figure2.png',w=11,h=5,
                                  figure_2())),
  
  # METHOD 1 - Temporally explicit model
  tar_target(panel_sim, sim_temporal_data()),
  tar_target(pfit, fit_panel_model(panel_sim)),
  tar_target(panel_regular_1, brm(a~b,data=panel_sim)),
  tar_target(panel_regular_2, brm(b~a,data=panel_sim)),
  tar_target(Figure4, save_figure('./Manuscript/Figures/Figure4.png',w=5,h=4,
                                  panel_model_figure(pfit, panel_regular_1, panel_regular_2))),

  # METHOD 2 - ODE analysis
  tar_target(ct_fit, eco_evo_sim()),
  tar_target(ct_plot1, save_figure('./Manuscript/Figures/Figure5-part1.png',w=8,h=4.5,
                                   plot_dynamics(ct_fit))),
  tar_target(ct_plot2, save_figure('./Manuscript/Figures/Figure5-part2.png',w=8*0.8,h=4*0.8,
                                   plot_effects(ct_fit))),
  
  # METHOD 3 - Instrumental variable analyses
  tar_target(IV_sim, run_IV_sim()),
  tar_target(IV_m1, fit_IV_model_brms(IV_sim, 1)),
  tar_target(IV_m2, fit_IV_model_brms(IV_sim, 2)),
  tar_target(reg_m1, fit_IV_model_brms(IV_sim, 3)), # regular regression for comparison
  tar_target(reg_m2, fit_IV_model_brms(IV_sim, 4)), # regular regression for comparison
  tar_target(Figure3, save_figure('./Manuscript/Figures/Figure3.png',w=5.5,h=5.5,
                                  instrumental_variable_plot(IV_m1, IV_m2, reg_m1, reg_m2))),

  ## Write vignettes - Github
  tar_quarto(v1g, file.path('Vignettes', 'Vignette-DiscreteTime_github.qmd')),
  tar_quarto(v2g, file.path('Vignettes', 'Vignette-ContinuousTime_github.qmd')),
  tar_quarto(v3g, file.path('Vignettes', 'Vignette-InstrumentalVariables_github.qmd')),
  tar_quarto(v4g, file.path('Vignettes', 'Vignette-ReverseEffects_github.qmd')),
  
  ## Write vignettes - HTML
  tar_quarto(v1g_html, file.path('Vignettes', 'Vignette-DiscreteTime_HTML.qmd')),
  tar_quarto(v2g_html, file.path('Vignettes', 'Vignette-ContinuousTime_HTML.qmd')),
  tar_quarto(v3g_html, file.path('Vignettes', 'Vignette-InstrumentalVariables_HTML.qmd')),
  tar_quarto(v0g_html, file.path('Vignettes', 'Vignette-ReverseEffects_HTML.qmd')),
  
  # Supplemental figures
  tar_target(FigureS1, save_figure('./Manuscript/Figures/FigureS1.png',w=5,h=5,
                                   framework_figure())),

  # Write supplement
  tar_quarto(supplement, file.path('Manuscript', 'Supplement.qmd')) # full version
  

)
