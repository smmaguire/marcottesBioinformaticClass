import pandas as pd
files =!ls CpG*.txt
fullDat=pd.DataFrame()
for nextFile in files:
    data = pd.read_table(nextFile, header=None, skiprows=1, sep= r'\s*')
    data=data.ix[:,1:]
    data=data.rename(columns={'X.2':'meth1','X.3':'chrom','X.4':'site','X.5':'meth2'})
    grouped = data.groupby('site')
    test2=grouped.count()
    grouped2 =data.groupby(['meth2','site'])
    test3= grouped2.count()
    test4=test3.xs('Z',level='meth2')
    test5=pd.DataFrame(test4['meth1'])
    test6=pd.DataFrame(test2['chrom'])
    test7=test5.join(test6)
    test7['per_methyl']=test7['meth1'].astype(float)/(test7['chrom'])
    test7['sample']=nextFile
    test7=test7.rename(columns={'chrom':'coverage'})
    test7=test7.reset_index()
    fullDat=pd.concat([fullDat,test7])
print fullDat.head()