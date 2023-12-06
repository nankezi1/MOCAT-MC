get the latest DISCOS data for physical parameters.
script: discosweb.m

get TLEs for the period you care about (modify the URL to choose date you want), download as CSV.
TLE data part in script: discos_plus_TLE.m 

merge the two via discos_plus_TLE  to create initialized.mat or such (mat_sats variable)
script: discos_plus_TLE.m
