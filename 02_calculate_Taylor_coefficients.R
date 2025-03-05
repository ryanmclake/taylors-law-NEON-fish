gc()


# setup function to manually predict fish abundance using the three pass method
predict_function <- function(x, beta, alpha){
  est <- (alpha*x)/(beta+x)
  return(est)
}



rows <- unique(joined_fish$row)

# creat an output list
out <- list()

# loop through rows  and calculate the fish abundance for each independent 3-pass

for(i in 1:length(rows)){
  
  tryCatch({
    
  # create a vector of the 3-pass EF data
  for_vector <- joined_fish %>%
    filter(row == rows[i])
  vector = c(for_vector$`1`, for_vector$`2`, for_vector$`3`)
  
  # use the FSA package to calculate abundance for the Carle-Strube,
  # the Moran, Schnute, Seber, and Burnham logistic model methods
  carlestrube_predict <- removal(vector)
  
  carlestrube_estimate <- data.frame(estimate = carlestrube_predict$est[1], 
                                     lowerCI = carlestrube_predict$est[3],
                                     upperCI = carlestrube_predict$est[4]) %>%
    tibble() %>%
    mutate(row = rows[i]) %>%
    mutate(method = "Carle-Strube")
  
  moran_predict <- removal(vector, method="Moran")
  
  moran_estimate <- data.frame(estimate = moran_predict$est[1], 
                               lowerCI = moran_predict$est[2],
                               upperCI = moran_predict$est[3]) %>%
    tibble() %>%
    mutate(row = rows[i]) %>%
    mutate(method = "Moran")
  
  schnute_predict <- removal(vector,method="Schnute")
  
  schnute_estimate <- data.frame(estimate = schnute_predict$est[1], 
                                 lowerCI = schnute_predict$est[2],
                                 upperCI = schnute_predict$est[3]) %>%
    tibble() %>%
    mutate(row = rows[i]) %>%
    mutate(method = "Schnute")
  
  seber_predict <- removal(vector,method="Seber3")
  
  seber_estimate <- data.frame(estimate = seber_predict$est[1], 
                               lowerCI = seber_predict$est[3],
                               upperCI = seber_predict$est[4]) %>%
    tibble() %>%
    mutate(row = rows[i]) %>%
    mutate(method = "Seber")
  
  burnham_predict <- removal(vector,method="Burnham")
  
  burnham_estimate <- data.frame(estimate = burnham_predict$est[1], 
                                 lowerCI = burnham_predict$est[3],
                                 upperCI = burnham_predict$est[4]) %>%
    tibble() %>%
    mutate(row = rows[i]) %>%
    mutate(method = "Burnham")
  
  # bind each of the calculated methods to the list
  out[[i]] <- rbind(carlestrube_estimate, moran_estimate, 
                    schnute_estimate, seber_estimate, burnham_estimate)
  
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

# unpack the list and output the high, low, and mean estimated abundance of fish
predicted_abundance <- do.call("rbind", out) %>%
  group_by(row) %>%
  summarise_at(vars(estimate),
               list(Q1=~quantile(., probs = 0.5,na.rm = T),
                    mean=~mean(.,na.rm = T), 
                    Q3=~quantile(., probs = 0.95,na.rm = T),
                    variance=~var(estimate, na.rm = T))) %>%

estimated_abundance <- left_join(joined_fish, predicted_abundance, by = "row") %>%
  filter(!is.na(variance))

ggplot(estimated_abundance, aes(mean,variance, color = domainID))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  geom_smooth(method = "lm", se = F)


sites <- unique(estimated_abundance$siteID)

# creat an output list
out2 <- list()

for(i in 1:length(sites)){

  for_TPL <- estimated_abundance %>%
    filter(siteID == sites[i])
# calibrate the dmodel with the pass data with the nlsLM package
TPL = nlsLM(variance ~ alpha*mean^beta,
                             start = list(alpha = 3, 
                                          beta = 1.5),
                             data = for_TPL,
                             control = nls.lm.control(maxiter=1024))

 alpha_estimate <- data.frame(parameter = "ALPHA",
                              estimate = TPL$m$getPars()[1], 
                              estimate_error = summary(TPL)$coefficients[,2][1],
                              siteID = sites[i])

 beta_estimate <- data.frame(parameter = "BETA",
                             estimate = TPL$m$getPars()[2], 
                             estimate_error = summary(TPL)$coefficients[,2][2],
                             siteID = sites[i])

 out2[[i]] <- rbind(alpha_estimate, beta_estimate)

 }

predicted_parameters <- do.call("rbind", out2) %>%
  filter(parameter == "BETA")


predicted_parameters <- do.call("rbind", out2) %>%
  filter(parameter == "ALPHA")

estimated_abundance <- estimated_abundance %>%
  mutate(years = lubridate::year(boutEndDate))
  
year <- unique(estimated_abundance$years)
  
# creat an output list
out3 <- list()

for(f in 1:length(year)){
  
  for_TPL <- estimated_abundance %>%
    filter(years == year[f])
  # calibrate the dmodel with the pass data with the nlsLM package
  TPL = nlsLM(variance ~ alpha*mean^beta,
              start = list(alpha = 3, 
                           beta = 1.5),
              data = for_TPL,
              control = nls.lm.control(maxiter=1024))
  
  alpha_estimate <- data.frame(parameter = "ALPHA",
                               estimate = TPL$m$getPars()[1], 
                               estimate_error = summary(TPL)$coefficients[,2][1],
                               YEAR = year[f])
  
  beta_estimate <- data.frame(parameter = "BETA",
                              estimate = TPL$m$getPars()[2], 
                              estimate_error = summary(TPL)$coefficients[,2][2],
                              YEAR = year[f])
  
  out3[[f]] <- rbind(alpha_estimate, beta_estimate)
  
}

predicted_temporal_parameters <- do.call("rbind", out3) %>%
  filter(parameter == "ALPHA") %>%
  arrange(YEAR)
# Development #
# this is a Mechalis-Menton model that is used to estimate abundance
# in principle, we start with a zero pass with zero fish, the following 
# passes are an accumulated abundance curve, with each pass
# having fewer and fewer total fish. This creates and asymptotic curve that
# can be presented with a MM model. the model function is pulled into the
# # R environment in the first script. 
# 
# # setup a DF with an additonal pass that represents "0" pass and "0" fish
# three_pass_w_zero <- joined_fish %>%
#   tidyr::gather(., passNumber, count, `1`:`3`, factor_key=TRUE) %>%
#   group_by(taxonID, siteID, year, month, row) %>%
#   mutate(pass_sum = cumsum(count)) %>%
#   mutate(pass = as.numeric(passNumber)) %>%
#   group_modify(~ add_row(.x,.before=0)) %>%
#   select(taxonID, siteID, year, month, row, pass, pass_sum) %>%
#   mutate(pass = ifelse(is.na(pass), 0, pass)) %>%
#   mutate(pass_sum = ifelse(is.na(pass_sum), 0, pass_sum))
# 
# # like above, run a for loop that applies the model to each seperate 
# # fishing event that included a 3-pass
# rows <- c(unique(three_pass_w_zero$row))
# out_hollings <- list()
# three_prediction_raw <- list()
# 
# for(i in 1:length(rows)){
#   
#   tryCatch({
#     
#     data_for_model <- three_pass_w_zero %>%
#       filter(row == rows[i])
#     
#     # calibrate the dmodel with the pass data with the nlsLM package
#     three_pass_predictor = nlsLM(pass_sum ~ (alpha*pass)/(beta+pass),
#                                  start = list(alpha = data_for_model$pass_sum[4], 
#                                               beta = 3),
#                                  data = data_for_model,
#                                  control = nls.lm.control(maxiter=1024))
#     
#     # Generate the estimate from the function and pull "JUST THE PARAMERTERS"
#     # Here we are pulling 500 random parameter estimates
#     
#     # alpha is the atymptote (i.e. V Max)
#     alpha_estimate <- data.frame(estimate = three_pass_predictor$m$getPars()[1], 
#                                  estimate_error = summary(three_pass_predictor)$coefficients[,2][1]) %>%
#       tibble() %>%
#       mutate(row = rows[i]) %>%
#       mutate(method = "Type 2 Vmax") %>%
#       mutate(estimate_low = estimate - estimate_error,
#              estimate_high = estimate + estimate_error) %>%
#       dplyr::select(-estimate_error) %>%
#       dplyr::select(estimate, estimate_low, estimate_high, row, method)
#     
#     # Beta is K or the half saturation constant. 
#     beta_estimate <- data.frame(estimate = three_pass_predictor$m$getAllPars()[2], 
#                                 estimate_error = summary(three_pass_predictor)$coefficients[,2][2]) %>%
#       tibble() %>%
#       mutate(row = rows[i]) %>%
#       mutate(method = "Beta 3-Pass") %>%
#       mutate(estimate_low = estimate - estimate_error,
#              estimate_high = estimate + estimate_error) %>%
#       dplyr::select(-estimate_error) %>%
#       dplyr::select(estimate, estimate_low, estimate_high, row, method)
#     
#     # Randomly select 500 parameter values for the Vmax and K
#     three_pass_alpha <- rnorm(500, three_pass_predictor$m$getAllPars()[1], summary(three_pass_predictor)$coefficients[,2][1])
#     three_pass_beta <- rnorm(500, three_pass_predictor$m$getAllPars()[2], summary(three_pass_predictor)$coefficients[,2][2])
#     
#     # create a data frame of these outputs
#     three_pass_parms <- data.frame(alpha = three_pass_alpha, 
#                                    beta = three_pass_beta) %>%
#       tibble(.)
#     
#     # setup passes from 3 to 50 (these represent the projections into the future passes)
#     passes <- seq(from = 3, to = 50, by = 1)
#     # create an output list 
#     prediction <- list()
#     
#     for(s in 1:length(passes)){
#       # run through each pass and predict the fish abundance for each successive pass from 
#       # 3 through 50. You cn go as high a you want, but the values will never exceed 
#       # the Vmax
#       prediction_raw <- predict_function(x = passes[s],
#                                          beta = three_pass_parms$beta,
#                                          alpha = three_pass_parms$alpha) %>%
#         data.frame(abundance = ., pass = passes[s]) %>%
#         tibble(.) %>%     
#         mutate(abundance = round(abundance, digits = 0)) %>%
#         mutate(abundance = ifelse(abundance < 0, 0, abundance)) %>%
#         mutate(iteration = rep(row_number()))
#       
#       prediction[[s]] <- prediction_raw
#       
#     }
#     
#     # output the predictions and bind them together, creates a big DF
#     three_prediction_raw[[i]] <- do.call("rbind", prediction) %>%
#       mutate(row = rows[i])
#     
#     # Filter the 50th prediction and pull the summary stats to get the abundance
#     hollings_predict <- do.call("rbind", prediction) %>%
#       filter(pass == 50) %>%
#       summarise(estimate = mean(abundance, na.rm = TRUE),
#                 estimate_low = quantile(abundance, probs = 0.05),
#                 estimate_high = quantile(abundance, probs = 0.95)) %>%
#       mutate(row = rows[i]) %>%
#       mutate(method = "Type 2 Posterior Prediction 3-Pass")
#     
#     out_hollings[[i]] <-  rbind(alpha_estimate, hollings_predict, beta_estimate)
#     
#   }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
#   
# }
# 
# # bind and compile all of the comparisons from FSA package and this new modelling approach
# three_pass_prediction <- do.call("rbind", out_hollings) %>%
#   dplyr::select(estimate, estimate_low, estimate_high, row, method) %>%
#   arrange(row)
# 
# readr::write_csv(three_pass_prediction, "./fish/input/three_pass_predictions_raw_taxon.csv")
# 
# prediction_comparisons <- bind_rows(three_pass_prediction, three_pass_predict_methods_compare) %>%
#   arrange(row)
# 
# readr::write_csv(prediction_comparisons, "./fish/input/predictions_comparisons_file_taxon.csv")
# 
# 
