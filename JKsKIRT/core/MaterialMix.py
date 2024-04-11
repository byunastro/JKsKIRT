import numpy as np
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Dust Material Mix
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

def MatMix(wfile,
           Dust_Type,
           Dust_Env,
           N_Si,
           N_C,
           N_PAH,
           indent,
           indent_base
           ):
    
    N_idt = indent_base

    print((N_idt)*indent+'<materialMix type="MaterialMix">',file=wfile)
    
    if Dust_Type == 'mrn':
        print((N_idt+1)*indent+'<MRNDustMix environment="MilkyWay" numSilicateSizes="%s" numGraphiteSizes="%s"/>'%(N_Si,N_C),file=wfile) 
    elif Dust_Type == 'zubko':
        print((N_idt+1)*indent+'<ZubkoDustMix numSilicateSizes="%s" numGraphiteSizes="%s" numPAHSizes="%s"/>'%(N_Si,N_C,N_PAH),file=wfile) 
    elif Dust_Type == 'draine_li':
        print((N_idt+1)*indent+'<DraineLiDustMix numSilicateSizes="%s" numGraphiteSizes="%s" numPAHSizes="%s"/>'%(N_Si,N_C,N_PAH),file=wfile)
    elif Dust_Type == 'weingartner_draine':
        if Dust_Env in ['MilkyWay','MW',',mw','Milkyway','milkyway','milkyWay', 'Milky Way', 'milky way', 'Milky-Way', 'milky-way']:
            print((N_idt+1)*indent+'<WeingartnerDraineDustMix environment="MilkyWay" numSilicateSizes="%s" numGraphiteSizes="%s" numPAHSizes="%s"/>'%(N_Si,N_C,N_PAH),file=wfile)
        elif Dust_Env in ['LMC','lmc','Large Magellanic Cloud', 'Large-Magellanic-Cloud','large magellanic cloud', 'large-magellanic-cloud']:
            print((N_idt+1)*indent+'<WeingartnerDraineDustMix environment="LMC" numSilicateSizes="%s" numGraphiteSizes="%s" numPAHSizes="%s"/>'%(N_Si,N_C,N_PAH),file=wfile)
        elif Dust_Env in ['SMC','smc','Small Magellanic Cloud', 'Small-Magellanic-Cloud','small magellanic cloud', 'small-magellanic-cloud']:
            print((N_idt+1)*indent+'<WeingartnerDraineDustMix environment="LMC" numSilicateSizes="%s" numGraphiteSizes="%s" numPAHSizes="%s"/>'%(N_Si,N_C,N_PAH),file=wfile)
    elif Dust_Type == 'trust':
        print((N_idt+1)*indent+'<TrustBenchmarkDustMix numSilicateSizes="%s" numGraphiteSizes="%s" numPAHSizes="%s"/>'%(N_Si,N_C,N_PAH),file=wfile)
    elif Dust_Type == 'themis':
        print((N_idt+1)*indent+'<ThemisDustMix numSilicateSizes="%s" numHydrocarbonSizes="%s"/>'%(N_Si,N_C),file=wfile)   

    print((N_idt)*indent+'</materialMix>',file=wfile)
    
    return N_idt


    








