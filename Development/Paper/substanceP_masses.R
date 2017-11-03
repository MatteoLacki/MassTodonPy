library(IsoSpecR)
library(tidyverse)

ENVELOPE = IsoSpecify(molecule = c(C=63, H=97, N=17, O=14, S=1), .99, showCounts = T)
ENVELOPE = ENVELOPE %>% 
           tbl_df() %>%
           mutate(prob = exp(logProb))

monoisotopic = ENVELOPE %>% filter(prob == max(prob)) 
m_mono = monoisotopic$mass

sprintf("%.10f", m_mono) 


data(isotopicData)
isotopicData$IsoSpec %>% filter()

m_H = IsoSpecify(molecule = c(H=1), 2.0, showCounts = T)[1,1]

round((m_mono + 3*m_H), 3)
round((m_mono + 3*m_H)/2, 3)
round((m_mono + 3*m_H)/3, 3)

m_C13_peak = ENVELOPE[2,1]

round((m_C13_peak + 2*m_H)/2, 3)
round((m_C13_peak + 2*m_H)/2, 3)