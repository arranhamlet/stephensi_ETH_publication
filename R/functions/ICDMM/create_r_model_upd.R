create_r_model_upd <- function(odin_model_path = system.file("extdata/odin_model.R",package="hanojoel"),
         het_brackets = 5,
         age = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,7.5,10,15,20,30,40,50,60,70,80),
         init_EIR = 10,
         init_ft = 0.4,
         country = NULL,
         admin2 = NULL,
         ...){
  
  ## create model param list using necessary variables
  mpl <- model_param_list_create(...)
  
  # generate initial state variables from equilibrium solution
  state <- equilibrium_init_create(age_vector=age, EIR=init_EIR,ft=init_ft,
                                   model_param_list = mpl, het_brackets=het_brackets,
                                   country = country,
                                   admin_unit = admin2)
  
  # create odin generator
  gen <- odin::odin(odin_model_path, verbose=FALSE)
  state <- state[names(state) %in% coef(gen)$name]
  
  # return mod
  return(list("state" = state, "generator" = gen, "mpl" = mpl))
}
