from Method import epslon, Cmin, Thin, Tcin, Tcout, Thout, cp1c, cp1h, mctotal, est

q = epslon * Cmin * (Thin - Tcin)
qc = mctotal * cp1c * (Tcin - Tcout)
qh = mctotal * cp1h * (Thin - Thout)
    
print("est, cp1h, q, qc, qh", est, cp1h, q, qc, qh)