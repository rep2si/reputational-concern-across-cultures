# Returns a vector of names that can be used in ggplot
get_site_labels<- \() {
  return(
    c(
      total = "All communities",
      te = "TE - India",
      wl = "WL - Colombia",
      wh = "WH - Colombia",
      tz = "TZ - Morocco",
      az = "AZ - India",
      tm = "TM - Mexico",
      np = "NP - Nepal",
      up = "UP - India"
    )
  )
}
