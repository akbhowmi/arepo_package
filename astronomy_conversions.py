def bolometric_correction_hardXray(Lbol):
            Lsun=3.8e33
            c1,k1,c2,k2= 4.073,-0.026, 12.60, 0.278
            Lbol_by_Lband = c1*(Lbol/1e10/Lsun)**k1 + c2*(Lbol/1e10/Lsun)**k2
            #print("Bol",Lbol_by_Lband)
            return Lbol_by_Lband
