import urllib.request
import sys
import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile


filename = sys.argv[-1]

#The method read_excel() reads the data into a Pandas Data Frame, where the first parameter is the filename and the second parameter is the sheet.
df = pd.read_excel(filename, sheetname='Sheet1')

#To iterate over the list we can use a loop:
for i in df.index:
    #print(df['Gene name'][i])
    x = df['Gene name'][i]
    print(x)
    with urllib.request.urlopen("https://hive.biochemistry.gwu.edu/biomuta/api/genesearch?gene=x") as response:
   		biomuta_result = response.read()
   		print(biomuta_result)