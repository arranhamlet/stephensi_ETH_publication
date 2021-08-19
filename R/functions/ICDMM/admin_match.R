#------------------------------------------------
#' match clean
#'
#' @param a First string to compare. Default = NULL
#' @param b Second string to compare. Default = NULL
#'
#' @importFrom RecordLinkage levenshteinSim
#'
#' @export

match_clean <- function(a, b){
  
  a <- gsub("[[:punct:][:space:]]", "", tolower(stringi::stri_trans_general(a, "latin-ascii")))
  b <- gsub("[[:punct:][:space:]]", "", tolower(stringi::stri_trans_general(b, "latin-ascii")))
  
  ret <- which(b %in% a)
  
  if(length(ret) == 0){
    distance <- levenshteinSim(a, b)
    ret <- which.max(distance)
  }
  return(ret)
}

#------------------------------------------------
#' match admin region
#'
#' \code{admin_match} Matches the user input admin unit and country with data
#'
#' @param country Character for country within which admin unit is in.
#'   Default = NULL
#' @param admin_unit Character for admin region. Some fuzzy logic will be used to
#'   match. If not provided then no seasonality is introduced. Default = NULL
#' @param admin_units_seasonal Dataframe of seasonality data for country and admin unit
#'
#' @export

admin_match <- function(admin_unit = NULL, country = NULL,
                        admin_units_seasonal){
  
  # intialise admin match as no match
  admin_matches <- 0
  
  if (!is.null(admin_unit)) {
    
    # if there is no country given then search for the admin unit
    if (is.null(country)) {
      
      # find exact match
      admin_matches <- which(tolower(admin_units_seasonal$admin1) %in% tolower(admin_unit))
      
      # if exact does not match try closest match
      if (length(admin_matches) < 1) {
        admin_matches <- match_clean(admin_unit, admin_units_seasonal$admin1)
      } else if (length(admin_matches) > 1){
        stop('Please specify the country of admin unit.  There are multiple with same name.')
      }
      
      # if we do have a country though find that match first and then find admin
    } else {
      
      # first find an exact match
      country_matches <- which(tolower(admin_units_seasonal$country) %in% tolower(country))
      
      if (length(country_matches) < 1) {
        country_name <- admin_units_seasonal$country[match_clean(country, admin_units_seasonal$country)]
        country_matches <- which(tolower(admin_units_seasonal$country) %in% tolower(country_name))
      }
      
      sub_admin_units_seasonal <- admin_units_seasonal[country_matches,]
      
      # find exact match
      admin_sub_matches <- which(tolower(sub_admin_units_seasonal$admin1) %in% tolower(admin_unit))
      
      # if exact does not match try closest match
      if (length(admin_sub_matches) != 1) {
        admin_sub_matches <- match_clean(admin_unit,
                                         sub_admin_units_seasonal$admin1)
      } else if (length(admin_sub_matches) > 1){
        stop('There are multiple admins with same name - check the data file!')
      }
      
      admin_matches <- country_matches[admin_sub_matches]
    }
    
    message("Requested: ", admin_unit,
            "\nReturned: ", admin_units_seasonal$admin1[admin_matches], ", ",
            admin_units_seasonal$country[admin_matches])
  }
  
  return(admin_matches)
}