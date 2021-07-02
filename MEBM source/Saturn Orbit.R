orbit <- read.csv("Saturn_Orbital_Elements.csv", as.is = TRUE, 
                  stringsAsFactors = FALSE)

colnames(orbit) <- c("ecc", "inc", "w_p", "capom", "L_p", "obl", "w", "f", "M",
                     "a", "L_s", "t")
tail(orbit)

plot(rev(orbit$t)/1000, rev(orbit$ecc), type = "l", lwd = "0.25", xlab = "kYr", 
     ylab = "Eccentricity", ylim = c(0, 0.1))

# KEY:

# Column 1; ecc: eccentricty
# Column 2; inc: incidence angle [radians]
# Column 3; w_p: longitude of perihelion [radians]
# Column 4; capom: longitude of ascending node [radians]
# Column 5; L_p: Ls of perihelion [radians]
# Column 6; obl: Obliquity of Saturn [radians] (calculated from secular theory)
# Column 7; w: argument of periapse [radians]
# Column 8; f: true anomaly [radians]
# Column 9; M: mean anomaly [radians]
# Column 10; a: semi-major axis (AU)
# Column 11; L_s: planetocentric solar longitude [radians]
# Column 12; t: input time vector [years since J2000]

