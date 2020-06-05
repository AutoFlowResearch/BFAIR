import pandas as pd
data = pd.read_csv (r'C:\Users\Teeradon Phlairaharn\Desktop\m_MTBLS1455_LC-MS_positive_reverse-phase_metabolite_profiling_v2_maf.tsv',sep='\t')   
df = pd.DataFrame(data)
print(df)