library(IsoSpecR)
D %>% filter(mz_L > 1079.5, mz_R < 1080)

Iso = function(mol)
    IsoSpecify(molecule = mol, stopCondition = .99) %>%
    as_data_frame() %>%
    mutate( prob = exp(logProb) ) %>%
    select( mass, prob)

z11 = c(C=63, H=96, N=17, O=13, S=1) 
z9  = c(C=52, H=77, N=12, O=11, S=1)


precursor = c(C=63,H=98,N=18,O=13,S=1)

Q = 1.0

Iso(z9) %>%
    mutate( mz = (mass + Q)/Q ) %>%
    data.frame

Iso(precursor) %>%
    mutate( mz = (mass + Q)/Q ) %>%
    data.frame


Iso(precursor) %>%
    mutate( mz = (mass + 3.)/3. ) %>%
    data.frame
