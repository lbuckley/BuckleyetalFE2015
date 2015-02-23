


#-------------------

biophys=function(Taira, Ta, WIND, SVL, MASS, psi, rho_S, elevation, J){
#Taira: air temperature, C
#Ta: surface temperature, C
#Wind speed, m/s
#SVL, mm
#MASS= g 
#psi = radian
#rho_S: albedo percent
#Elevation (m)
#J is Julian Day

# Biophysical models from Campbell & Norman 1998
# constants
sigma=5.67*10^-8 # stefan-boltzman constant, W m^-2 K^-4
c_p=29.3 # specific heat of air, J/mol degrees K or C

# absorptivity
alpha_S=0.9  # solar absorptivity (Gates 1980, Table 11.4)
alpha_L=0.965 # thermal absoptivity, Bartlett & Gates 1967 
epsilon_s=0.965 # surface emisivity (p163), Bartlett & Gates 1967

F_d=0.8  # diffuse view angle, Bartlett & Gates 1967
F_r=0.5  # reflected solar radiation
F_a=0.5  # atmospheric radiation
F_g=0.5  # ground thermal radation

tau=0.65 # atmospheric transmisivity
S_p0=1360 # extraterrestrial flux density, W/m^2 (p159)

rd=180/pi  # factor to convert radians into degrees

# Calculate radiation
# view angles, parameterize for animal suspended above ground (p181), on ground- adjust F_e, F_r, and F_g
h=SVL/1000 # length of cylinder in m
theta = psi # angle between solar beam and a normal to the plane in radians, = psi for horizontal surfaces

# F_p=(cos (theta)+(4*h*sin (theta))/(pi*d))/(2+4*h/d)  # beam view angle, Fig 11.6
A=0.121*MASS^0.688   # total lizard area, roughgarden 1981 from Norris (1965) and Porter and James (1979)
A_p= (-1.1756810^-4*psi^2-9.2594*10^-2*psi+26.2409)*A/100      # projected area
F_p=A_p/A
                
# radiation
p_a=101.3* exp (-elevation/8200)  # atmospheric pressure
m_a=p_a/(101.3*cos (psi))  # (11.12) optical air mass
m_a[which(psi>(80*pi/180))]=5.66
                
# Flux densities
#without climate change
epsilon_ac= 9.2*10^-6*(Taira+273)^2 # (10.11) clear sky emissivity
L_a=epsilon_ac*sigma*(Taira+273)^4  # (10.7) long wave flux densities from atmosphere 
L_g=epsilon_s*sigma*(Ta+273)^4  # (10.7) long wave flux densities from ground

#S_p=S_p0*tau^m_a # (11.11) direct irradience , W/m^2
dd2= 1+2*0.1675*cos(2*pi*J/365)
S_p=S_p0*tau^m_a*dd2 *cos(psi)  #Sears and Angilletta 2012 #dd is correction factor accounting for orbit

S_d=0.3*(1-tau^m_a)* S_p  # (11.13) diffuse radiation
#S_t=S_p*cos (psi)+S_d # solar irradience 
S_r= rho_S*S_p # (11.10) reflected radiation


			   
#__________________________________________________
# conductance

dim=SVL/1000 # characteristic dimension in meters (Table 9.5)
g_r= 4*sigma*(Taira+273)^3/c_p # (12.7) radiative conductance

g_Ha=1.4*0.135*sqrt(WIND/dim) # boundary conductance, factor of 1.4 to account for increased convection (Mitchell 1976)
                
#__________________________________________________
# operative environmental temperature

#calculate with both surface and air temp (on ground and in tree)

sprop=1 #proportion of radiation that is direct, Sears and Angilletta 2012
R_abs= sprop*alpha_S*(F_p*S_p+ F_d*S_d + F_r*S_r)+alpha_L*(F_a*L_a+F_g*L_g) # (11.14) Absorbed radiation
Te=Taira+(R_abs-epsilon_s*sigma*(Taira+273)^4)/(c_p*(g_r+g_Ha))                       
Te_surf= Ta+(R_abs-epsilon_s*sigma*(Ta+273)^4)/(c_p*(g_r+g_Ha))        

# calculate in shade, no direct radiation
sprop=0 #proportion of radiation that is direct, Sears and Angilletta 2012
R_abs= sprop*alpha_S*(F_p*S_p+ F_d*S_d + F_r*S_r)+alpha_L*(F_a*L_a+F_g*L_g) # (11.14) Absorbed radiation
TeS=Taira+(R_abs-epsilon_s*sigma*(Taira+273)^4)/(c_p*(g_r+g_Ha))                       
TeS_surf=Ta+(R_abs-epsilon_s*sigma*(Ta+273)^4)/(c_p*(g_r+g_Ha))  

return(c(Te, TeS,Te_surf,TeS_surf))
}
