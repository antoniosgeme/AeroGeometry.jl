
"""
Convert feet and inches to meters.
"""
ft2m(feet, inches = 0) = 0.3048 * feet + 0.0254 * inches

"""
Convert pounds (lb) to kilograms (kg).
"""
lb2kg(pounds) = pounds * 0.453592

"""
Convert ounces (oz) to kilograms (kg).
"""
oz2kg(ounces) = ounces * 0.0283495

"""
Convert gallons (US) to liters (L).
"""
gal2l(gallons) = gallons * 3.78541

"""
Convert pounds per square inch (psi) to Pascals (Pa).
"""
psi(pressure) = pressure * 6894.76

"""
Convert atmospheres (atm) to Pascals (Pa).
"""
atm2pa(atm) = atm * 101325

"""
Convert inches to meters.
"""
in2m(inches) = inches * 0.0254

"""
Convert miles to meters.
"""
mi2m(miles) = miles * 1609.34

"""
Convert miles per hour (mph) to meters per second (m/s).
"""
mph2mps(speed) = speed * 0.44704

"""
Convert knots (kts) to meters per second (m/s).
"""
kts2mps(knots) = knots * 0.514444

"""
Convert horsepower (hp) to watts (W).
"""
hp2w(horsepower) = horsepower * 745.7

"""
Convert yards (yd) to meters (m).
"""
yd2m(yards) = yards * 0.9144

"""
Convert pound-feet (lb-ft) to Newton-meters (Nm).
"""
lbft2nm(torque) = torque * 1.35582

"""
Convert foot-pounds force (ft-lbf) to joules (J).
"""
ftlbf2j(energy) = energy * 1.35582

"""
Convert Rankine (°R) to Kelvin (K).
"""
rankine2kelvin(temp_r) = temp_r * (5/9)
