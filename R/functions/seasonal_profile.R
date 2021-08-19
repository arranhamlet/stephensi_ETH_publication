seasonal_profile <- function(fourier_transform_data, scalar = 1){
  
  sapply(1:365, function(t) pmax((fourier_transform_data$a0 + 
                                   scalar * fourier_transform_data$a1 * cos(2 * pi * t/365) + 
                                   scalar * fourier_transform_data$a2 * cos(2 * 2 * pi * t/365) + 
                                   scalar * fourier_transform_data$a3 * cos(3 * 2 * pi * t/365) + 
                                   scalar * fourier_transform_data$b1 * sin(2 * pi * t/365) + 
                                   scalar * fourier_transform_data$b2 * sin(2 * 2 * pi * t/365) + 
                                   scalar * fourier_transform_data$b3 * sin(3 * 2 * pi * t/365) )/(fourier_transform_data$theta_c), 
                                0.001))
}