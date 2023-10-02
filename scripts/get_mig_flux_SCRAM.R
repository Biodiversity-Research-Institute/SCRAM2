#' Migration Flux factor
#'
#' Returns the migratory flux factor, expressing the total number of
#' bird flights through the rotors of the wind farm per month, if all flights
#' occur within the rotor's circle area of all turbines, and assuming birds take no
#' avoiding action.
#'
#' @details The flux factor is used for other model calculations.
#'    Methodology and assumptions underpinning
#'    \code{get_mig_flux_factor_SCRAM} are described in "Annex6" of
#'    \href{https://www.bto.org/sites/default/files/u28/downloads/Projects/Final_Report_SOSS02_Band1ModelGuidance.pdf}{Band (2012)}
#'
#' @param n_turbines An integer. The number of turbines on the wind farm (\eqn{T}).
#' @param rotor_radius An integer. The radius of the rotor (\eqn{R}), in meters
#' @param density_est A numeric. The density estimate (birds/km) derived from modeling occupancy and the cell width for the modeled cell or population moving across broad front
#'
#' @return
#' The number of bird flights potentially transiting through rotors
#' at each time period (assuming no avoidance), if all flights occur within the
#' rotor's circular area
#'
#' @examples
#' get_mig_flux_SCRAM(
#'      n_turbines = 100,
#'      rotor_radius = 120,
#'      density_est = 0.51
#' )
#' @export
get_mig_flux_SCRAM <- function(n_turbines, rotor_radius, density_est){
  
  tot_frontal_area <- n_turbines * pi * rotor_radius^2
  
  # Get an estimate of birds/km and then convert that to birds/m
  bird_dens_per_m <- density_est / 1000
  
  #active_secs <- (daynight_hrs$Day + noct_activity * daynight_hrs$Night) * 3600
  
  (bird_dens_per_m/(2*rotor_radius)) * tot_frontal_area #* active_secs
  
}
