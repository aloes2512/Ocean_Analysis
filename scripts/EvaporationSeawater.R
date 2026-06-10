#Water evaporation:
## Clausius- Clapeyron) e saturation evaporation pressure
d_ln(e)= (-dh/(RT^2)+C )*dT# h latent heat; R gas constant
h=2.5*10^6 #J/kg for water
R=4615 #J/kg·K
# integrating with constant h
e = e0*exp(-h0/R(1/T-1/T0))
# evaporation E links to e :
E ~ (e- ea)*u # ea = air vapor pressure; u = wind speed
A #surface area m^2
E*A=E*(e- ea)*u*A # 𝐴⋅(𝑥𝑠−𝑥)
# approximation (standard bulk formula)
v #'windspeed'm/s
gh=(25+19*v)*(Xs-X)*A # evaporation rate (kg/h)
Xs=e # [kg/kg] saturation huminity level of seawater at Temperature T
Xs # calculated using Clausius Clapeyron
X #actual surrounding air huminity measured (kg H₂O/kg dry air)

Θ = 25+19*v #= evaporation coefficient (kg/(m²·h)),
Θ=25+19*v #(with v = wind speed in m/s)
A #surface area of water (m²)
Xs=e #humidity ratio of saturated air at the water surface temperature (kg H₂O/kg dry air)
X # humidity ratio of the surrounding air (kg H₂O/kg dry air)
# evaporation rate
Θ*(Xs-X)*A
