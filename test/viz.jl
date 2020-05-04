using Controlz

# https://lpsa.swarthmore.edu/Root_Locus/RLocusExamples.html
g = (s + 1) / (s^3 + 4 * s^2 + 6 * s + 4)
g = 1 / (s * (s + 3))
g = (s + 3) / (s^2 - s - 2)
g = (s^2 + 2 * s + 2) / (s * (s^4 + 9 * s^3 + 33 * s^2 + 51 * s + 26))

root_locus(g, max_mag_Kc=3.0, nb_pts=5000)
