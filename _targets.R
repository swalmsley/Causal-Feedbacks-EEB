
library(targets)

# Source
tar_source('R') # will do all in 'R' folder

# Options
tar_option_set()

# Seed
tar_option_set(seed = 1234)

# Variables
c1 <- col.desat( "blue" , 0.8 )
c2 <- col.desat( "red" , 0.8 )



# Targets
list(
  
  
  # General figures
  tar_target(Figure1, save_figure('./Manuscript/Figures/Figure1.png',w=5,h=5,
                                  framework_figure())),
  tar_target(Figure2, save_figure('./Manuscript/Figures/Figure2.png',w=11,h=5,
                                  figure_2())),
  tar_target(FigureCoral, save_figure('./Manuscript/Figures/FigureCoral.png',w=9,h=3,
                                  coral_dag())),

  
  # METHOD 1 - Temporally explicit model
  tar_target(panel_sim, sim_temporal_data()),
  tar_target(pfit, fit_panel_model(panel_sim)),
  tar_target(panel_regular_1, brm(a~b,data=panel_sim)),
  tar_target(panel_regular_2, brm(b~a,data=panel_sim)),
  tar_target(Figure4, save_figure('./Manuscript/Figures/Figure4.png',w=4.5,h=5,
                                  panel_model_figure(pfit, panel_regular_1, panel_regular_2))),

  # METHOD 2 - ODE analysis
  tar_target(fitMCMC, create_ODE_example()),
  tar_target(Figure5, save_figure('./Manuscript/Figures/Figure5.png',w=4.5,h=4.5,
                                  ODEs_plot(fitMCMC))),
  
  # METHOD 3 - Instrumental variable analyses
  tar_target(IV_sim, run_IV_sim()),
  tar_target(IV_m1, fit_IV_model_brms(IV_sim, 1)),
  tar_target(IV_m2, fit_IV_model_brms(IV_sim, 2)),
  tar_target(reg_m1, fit_IV_model_brms(IV_sim, 3)), # regular regression for comparison
  tar_target(reg_m2, fit_IV_model_brms(IV_sim, 4)), # regular regression for comparison
  tar_target(Figure3, save_figure('./Manuscript/Figures/Figure3.png',w=10,h=5.5,
                                  instrumental_variable_plot(IV_m1, IV_m2, reg_m1, reg_m2))),
  
  ## Write vignettes
  tar_quarto(v1g, file.path('Vignettes', 'Vignette-PanelModel_github.qmd')),
  tar_quarto(v2g, file.path('Vignettes', 'Vignette-ODEs_github.qmd')),
  tar_quarto(v3g, file.path('Vignettes', 'Vignette-InstrumentalVariables_github.qmd')),
  
  # tar_quarto(v1, file.path('Vignettes', 'Vignette-PanelModel.qmd')),
  # tar_quarto(v2, file.path('Vignettes', 'Vignette-InstrumentalVariables.qmd')),
  # tar_quarto(v3, file.path('Vignettes', 'Vignette-ODEs.qmd')),
  
  # Write supplement
  # tar_quarto(supplement_1, file.path('Manuscript', 'Supplement_CausalFeedbacks.qmd')),
  tar_quarto(supplement, file.path('Manuscript', 'Supplement.qmd')),
  
  # Write manuscript
  tar_quarto(
    paper,
    file.path('Manuscript','MS_CausalFeedbacks.qmd'))
  
)
