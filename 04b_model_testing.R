library(dplyr)
library(tidyr)
library(brms)
library(rstan)

sim_func <- function(nsquare = 50, nplotpersquare = 5){
  nplot <- nsquare*nplotpersquare
  
  sim_data <- data.frame(SQUARE = gl(nplot/nplotpersquare,nplotpersquare)) %>%
    mutate(REP_ID = paste0(SQUARE,"X",1:nplot),    
           Improved = rbinom(nplot,1,0.007*as.numeric(SQUARE))) %>%
    mutate(pH = rnorm(nplot,mean = 4.5 + Improved + rep(rnorm(nplot/nplotpersquare,0,0.5),each = nplotpersquare), 1)) %>%
    mutate(Ell = rnorm(nplot, mean = pH + rep(rnorm(nplot/nplotpersquare,0,1),each = nplotpersquare), 1)) %>%
    mutate(C_year12 = rnorm(nplot, rep(rnorm(nplot/nplotpersquare,0,0.2),each = nplotpersquare), 1)) %>%
    mutate(N_year12 = rnorm(nplot, C_year12, 0.1)) %>%
    mutate(Sdep_year12 = rep(rnorm(nplot/nplotpersquare,0,1),each = nplotpersquare)) %>%
    mutate(Ndep_year12 = rep(rnorm(nplot/nplotpersquare,-Sdep_year12,1),each = nplotpersquare)) %>%
    mutate(rain_year12 = rep(rnorm(nplot/nplotpersquare,0,1),each = nplotpersquare),
           PH_SE_year12 = 0.2,
           ELL_SE_year12 = 0.15 + Improved*0.05,
           PH_SE_year23 = 0.15,
           ELL_SE_year23 = 0.2 + Improved*0.05,
           PH_SE_year34 = 0.15,
           ELL_SE_year34 = 0.15 + Improved*0.05) %>%
    mutate(pH_diffYear12 = rstudent_t(nplot, 4,
                                      0.5*Improved+(1-Improved)*0.5*Sdep_year12 + 
                                        0.3*rain_year12,PH_SE_year12))%>%
    mutate(Ell_diffYear12 = rstudent_t(nplot, 4,
                                       pH_diffYear12*(1-Improved) + 
                                         Improved*0.5*pH_diffYear12 + 
                                         (pH>5.5)*-0.5*N_year12,ELL_SE_year12)) %>%
    mutate(pH_year2 = pH + pH_diffYear12,
           Ell_year2 = Ell + Ell_diffYear12,
           C_year23 = C_year12 + rnorm(nplot,-0.15*C_year12,0.2),
           rain_year23 = rep(rnorm(nplot/nplotpersquare,0,1),each = nplotpersquare),
           Sdep_year23 = rep(rnorm(nsquare, Sdep_year12[seq(1,nplot,nplotpersquare)], 0.1),each=5),
           Ndep_year23 = rep(rnorm(nsquare, Ndep_year12[seq(1,nplot,nplotpersquare)], 0.1), each = 5)) %>%
    mutate(N_year23 = N_year12 + rnorm(nplot, -0.1*N_year12 +
                                         0.5*Ndep_year23 + (C_year23 - C_year12),
                                       0.1)) %>%
    mutate(pH_diffYear23 = rstudent_t(nplot, 4,
                                      0.5*Improved+(1-Improved)*0.5*Sdep_year23 + 
                                        0.3*rain_year23, PH_SE_year23)) %>%
    mutate(Ell_diffYear23 = rstudent_t(nplot, 4,
                                       pH_diffYear23*(1-Improved) + 
                                         Improved*0.5*pH_diffYear23 + 
                                         (pH_year2>5.5)*-0.5*N_year23,ELL_SE_year23)) %>%
    mutate(pH_year3 = pH_year2 + pH_diffYear23,
           Ell_year3 = Ell_year2 + Ell_diffYear23,
           C_year34 = C_year12 + rnorm(nplot,-0.15*C_year12,0.2),
           rain_year34 = rep(rnorm(nplot/nplotpersquare,-0.5,1),each = nplotpersquare),
           Sdep_year34 = rep(rnorm(nsquare, Sdep_year23[seq(1,nplot,nplotpersquare)], 0.1),each=5),
           Ndep_year34 = rep(rnorm(nsquare, Ndep_year23[seq(1,nplot,nplotpersquare)], 0.1), each = 5)) %>%
    mutate(N_year34 = N_year23 + rnorm(nplot, -0.1*N_year23 +
                                         0.5*Ndep_year34 + (C_year34 - C_year23),
                                       0.1)) %>%
    mutate(pH_diffYear34 = rstudent_t(nplot, 4,
                                      0.5*Improved+(1-Improved)*0.5*Sdep_year34 + 
                                        0.3*rain_year34, PH_SE_year34)) %>%
    mutate(Ell_diffYear34 = rstudent_t(nplot, 4,
                                       pH_diffYear34*(1-Improved) + 
                                         Improved*0.5*pH_diffYear34 + 
                                         (pH_year3>5.5)*-0.5*N_year34,ELL_SE_year34)) %>%
    select(SQUARE, REP_ID, Improved,
           Sdep_year12, Sdep_year23, Sdep_year34, Ndep_year12, Ndep_year23, Ndep_year34, 
           rain_year12, rain_year23, rain_year34,
           C_year12, C_year23, C_year34, N_year12, N_year23, N_year34,
           PH_SE_year12, PH_SE_year23, PH_SE_year34,
           ELL_SE_year12, ELL_SE_year23, ELL_SE_year34,
           PH_year12 = pH_diffYear12, PH_year23 = pH_diffYear23, PH_year34 = pH_diffYear34,
           Year1_pH_year12 = pH, Year1_pH_year23 = pH_year2, Year1_pH_year34 = pH_year3,
           Ell_year12 = Ell_diffYear12, Ell_year23 = Ell_diffYear23, Ell_year34 = Ell_diffYear34) %>%
    pivot_longer(contains("_year"), names_to = c("Variable","Time_period"),
                 names_sep = "_year") %>%
    pivot_wider(names_from = Variable, values_from = value) %>%
    mutate(YRnm = as.numeric(as.factor(Time_period)))
  
  
  mod_pr <- c(prior(normal(0,0.5), class = "b", resp = "Ell"),
              prior(normal(0,0.5), class = "b", resp = "PH"),
              prior(normal(0,0.5), class = "b", resp = "N"),
              prior(normal(0.8,0.2), class = "b", coef = "C", resp = "N"),
              prior(student_t(3, 0, 2.5), class = "sds", resp = "Ell"),
              prior(normal(0,0.25), class = "Intercept", resp = "Ell"),
              prior(normal(0,0.25), class = "Intercept", resp = "PH"),
              prior(student_t(3,0,1), class = "Intercept", resp = "N"),
              prior(gamma(4,1), class = "nu", resp = "Ell"),
              prior(gamma(4,1), class = "nu", resp = "PH"),
              prior(normal(0,0.1), class = "ar", resp = "Ell"),
              prior(normal(0,0.1), class = "ar", resp = "PH"),
              prior(normal(0,0.1), class = "ar", resp = "N"),
              prior(student_t(3,0,0.5), class = "sd", resp = "Ell"),
              prior(student_t(3,0,0.5), class = "sd", resp = "PH"),
              prior(student_t(3,0,0.5), class = "sd", resp = "N"),
              prior(student_t(3,0,0.5), class = "sigma", resp = "Ell"),
              prior(student_t(3,0,0.5), class = "sigma", resp = "PH"),
              prior(student_t(3,0,0.5), class = "sigma", resp = "N"))
  
  sim_mod <- brm(bf(Ell | mi(ELL_SE) ~ Improved*mi(PH) + s(Year1_pH, N, Improved) + 
                      (1|SQUARE) +
                      ar(time = YRnm, gr = REP_ID),
                    family = "student") +
                   bf(PH | mi(PH_SE)  ~ Improved*Sdep + rain + (1|SQUARE) +
                        ar(time = YRnm, gr = REP_ID), family = "student") + 
                   bf(N ~ Ndep + C + (1|SQUARE) +
                        ar(time = YRnm, gr = REP_ID)) +
                   set_rescor(FALSE), data = sim_data, prior = mod_pr,
                 save_pars = save_pars(all = TRUE, latent = TRUE), 
                 cores = 4, iter = 4000)
  
  sampler_params <- rstan::get_sampler_params(sim_mod$fit, inc_warmup = FALSE)
  divergent <- c(sampler_params[[1]][,'divergent__'],
                 sampler_params[[2]][,'divergent__'],
                 sampler_params[[3]][,'divergent__'],
                 sampler_params[[4]][,'divergent__'])
  
  # rerun if there are divergent transitions 
  # happens once, if setting adapt delta to 0.99 doesn't fix it then discard model results later
  if(sum(divergent)>0){
    sim_mod <- update(sim_mod, cores = 4, iter = 4000, 
                      save_pars = save_pars(all = TRUE, latent = TRUE), 
                      control = list(adapt_delta = 0.99))
    
    divergent <- rstan::get_sampler_params(sim_mod$fit, inc_warmup=FALSE)[[1]][,'divergent__']
    
  }
  
  # rerun with more iterations if neff too low - only do this if Rhat acceptable
  if(min(neff_ratio(sim_mod))<0.05 & min(neff_ratio(sim_mod))>0.02 & max(rhat(sim_mod))<1.05){
    sim_mod <- update(sim_mod, cores = 4, iter = 2000, 
                      save_pars = save_pars(all = TRUE, latent = TRUE))
    
  }
  
  # if divergent okay, rhat okay and neff okay then calculate if estimates
  # within the 50 and 90% CI
  sim_pars <- c(0,0,0,0.5,0.5,0.5,0.3,-0.5,0.5,1,0,0,0,0,0,0,0,0,0,1,-0.5)
  if(sum(divergent) == 0 & max(rhat(sim_mod))<1.05 & min(neff_ratio(sim_mod))>0.05){
    par_vals <- fixef(sim_mod, probs = c(0.05,0.25,0.75,0.95))
    
    CI90 <- sim_pars > par_vals[,"Q5"] & sim_pars < par_vals[,"Q95"]
    CI50 <- sim_pars > par_vals[,"Q25"] & sim_pars < par_vals[,"Q75"]
  } else{
    CI90 <- rep(NA, 21)
    CI50 <- rep(NA, 21)
  }
  
  return(cbind(CI50,CI90))
}

set.seed(7)
brm_model_test <- replicate(10, sim_func(), simplify = FALSE)


test <- make_conditions(sim_data, vars = c("Improved","Year1_pH"))[1:6,]
test[,"Improved"] <- c(0,0,0,1,1,1)
test[,"cond__"] <- gsub("-0.2","0",test[,"cond__"])
test[,"cond__"] <- gsub("0.2","1",test[,"cond__"])

plot(conditional_effects(sim_mod, conditions = test, effects = "N", 
                         select_points = 0.5), 
     points = TRUE, plot = FALSE)[[1]]
