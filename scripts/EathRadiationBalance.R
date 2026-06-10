# Heat balance
𝑑ℎ=(1−α)∗𝑆−4∗σ∗ε∗𝑇4
# h is the earth heat-content,
#S solar radiation 𝑊/𝑚2,
α=0.3 # reflection constant (albedo),
#T temperature in K,
σ = 5.67/10^8 #𝑊/𝑚^2 Stephan-Boltzmann-Constant,
A #Surface area
M = σ*T^4 # Stephan-Boltzmann Law for black body
ε=0.6 # estimate of average emissivity of earth surface.
Pw=M*ε # emitted power of real body
#$ ~ 1.36 kW/m^2$ #S at the top of the atmosphere
S= 1.36

# @ radiative equilibrium dh=0
α=0.3
T=((1- α)*𝑆/4*σ*ε)^0.25
0.7 = (1- α) # absorbed radiation; if reflection α=0.3
S= 1368  # W/m^2 im Mittel. 0.25*1368 Wm^-2
ε=0.6    #an estimate of emissivity
σ=5.67*10 -8 #𝑊*𝑚^−2𝐾^−4 Stephan-Boltzman-Constant
4*5.67*10^-8 #/10^8*0.6
T=sqrt(sqrt(0.7*S/(4*5.67*10^-8*0.6)))-273.2 # [°C]
T # 16.43261
