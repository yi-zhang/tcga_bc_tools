
""" RNA table class
This class works for TCGA RNA table as downloaded from gdc. The header line is "gene" followed by barcodes for each RNA-seq sample.

Example:
    # create a rnatable class from a pandas dataframe.
    myrna = rnatable(dataframe, NAME)
    # take log for provided columns, overwritten.
    myrna.log(COLUMN_LIST)

Example of reading BRCA expressions: 

import sys

MyPackageDIR=<DIR to rnatable_tcga>
if(MyPackageDIR not in sys.path):
        sys.path.append(MyPackageDIR)

        import rnatable_tcga.rnatable_func
        genecol = 'gene'
        # The following files on biocluster as examples
        #expFile='/home/groups/bd2k/TCGA/BRCA/RNA-seq/processed/Level_3/FPKM_hg38/tcga_brca_fpkm_hg38_er_pos_tumor_table.txt'
        expFile = '/home/groups/BD2K/TCGA/BRCA/RNA-seq/processed/Level_3/FPKM_UQ/TCGA_BRCA_FPKMUQ.txt'

        rna = pd.read_csv(expFile, header=0, sep='\t')
        rnapd = rnatable_tcga.rnatable_func.rnatable(rna, 'rnalog')
        rnapd.log(rnapd.pd.columns.values[1:])
        duppats = rnatable_tcga.rnatable_func.check_patdupsamp(rnapd.pd.columns.values)
        rnapd.duppat_takemean(duppats)
        rnapd.barcode_TP_keep12()
        rnapd.rmlowexp()
"""

def pcorr(x1,x2,x3):
    """ partial correlation"""
    import scipy.stats as ss
    [r12,pal] = ss.pearsonr(x1,x2)
    [r13,pal] = ss.pearsonr(x1,x3)
    [r23,pal] = ss.pearsonr(x2,x3)
    pcorr=(r12-r23*r13)/(1-r13*r13)
    return(pcorr)

def findgene( list1, genelist):
    found = [t for t in genelist if t in list1]
    notfound = list(set(list1)-set(genelist))
    return([found,notfound])

def subselectpats(dat, pats, keepcol = ['gene']):
    """ given the pats list, keep the keepcol and the given pats. can assign the returned subpd object as a new rnatable class.
     subname example: RRpd
    """
    found = [t for t in pats if t in dat.columns.values[1:]] # ordered
    notfound = set(pats) - set(found)
    if(len(notfound)!=0):
        print("Ignored pats returned.")
    dat = dat.loc[:,list(keepcol)+list(found)]
    return({"sub":dat,"notfound":notfound})
        
def subselectgenes(dat, genecol, genesets):
    """ given the gene list. can assign the returned subpd object as a new rnatable class. # unordered"""
    res = dat[dat[genecol].isin(list(genesets))]
                
    notfound = [t for t in genesets if t not in res[genecol].values]
    if(len(notfound)!=0):
        print("Ignored genes returned.")
    return({"sub":res,"notfound":notfound})

def gettfgenes(tfgenefiles):
    import pandas as pd
    # example:
    #tf_db_file = '/home/groups/song/songlab2/yzhan201/DataTools/AnimalTFDB/TFHumanGenes.txt'
    #tf_mkup_file = '/home/groups/song/songlab2/yzhan201/DataTools/AnimalTFDB/TfHumanMissed.mkup.txt'
    tfhuman=pd.read_csv(tfgenefiles[0],sep='\t',header=0,index_col=None)
    GeneList= list(set([t for t in tfhuman['Symbol'].values if str(t)!='nan']))
    AllTf_GeneList = GeneList
    if(len(tfgenefiles)==2):
        tfhumanmkup=pd.read_csv(tfgenefiles[1], sep='\t',header=0,index_col=None)
        GeneListMkup=list(tfhumanmkup['TcgaName'].values)
        AllTf_GeneList = list(set(GeneList+GeneListMkup))
    return(AllTf_GeneList) 

def getcorr_1target(dat, tf_A_genelist, tf_B_genelist, genecol = 'gene'):
    """ calling : subselectgenes """
    targetgenenow = tf_B_genelist[0]
    rnaAB = subselectgenes(dat, genecol, tf_A_genelist+tf_B_genelist) # new pandas containing list A and B. 
    rnaABcorr = rnaAB['sub'].set_index(genecol).transpose().corr()
    rnaABcorr = rnaABcorr.loc[tf_A_genelist,tf_B_genelist]
    #rnaABcorr.reset_index(inplace=True)
    rnaABcorr[genecol] = rnaABcorr.index
    rnaABcorr.rename(columns={targetgenenow:'target'}, inplace=True)
    return(rnaABcorr)

def check_dupgenes(all_genelist):
    """Not completed.# pick : genes mapping to multiple positions in gencode"""
    # all_genelist = rnapd.pd[genecol].values
    import collections
    cnt = collections.Counter(all_genelist)
    dupgenes = [t[0] for t in cnt.most_common() if t[1]>1]
    return(dupgenes)
        
def check_nolowexp(rna, cutoff=1):
    """ calculate mean expression and return mask """
    mask = rna.mean(axis=1,numeric_only=True, skipna=True) >1
    return(list(rna.reset_index().loc[mask,:]['gene'].values))

def check_patdupsamp(allbarcode_list):
    """check patients with duplicate samples; return a dict with patients and their duplicated samples.
    # An example for input barcodelist: 
    # tcgaTPcol=[re.match('TCGA.+\-01[AB]\-.+',term).group() for term in allbarcode_list if re.match('TCGA.+\-01[AB]\-.+',term)] 
    """
    seenpats= []
    duppats = {}
    for barcode in allbarcode_list:
        if(barcode[:12] not in seenpats):
            seenpats.append(barcode[:12])
        else:
            duppat = barcode[:12]
            duppat_barcode = [t for t in allbarcode_list if t[:12] == duppat ] 
            duppats.update({duppat:duppat_barcode})
    return(duppats)

def check_01RNA(allbarcode_list):
    """check tumor (01) RNA-seq samples; return a list of tumor RNA-seq barcodes.
    """
    rna01_list = [t for t in allbarcode_list if ("-" in t and t.split("-")[3][:-1] == "01" ) ] 
    ## example: TCGA-AA-AAAA-01A-XX...
    return(rna01_list)
    

class rnatable:

    def log(self, colnames):
        import numpy as np
        """ take the table, transform all numerical to log2(value+1)"""
        self.pd[colnames]= self.pd[colnames].copy().fillna(0).applymap(lambda x: np.log2(x+1))
        #return(self)

    def get_gene_mask(self, genelist):
        """ take the table, get the row masks for input genelist. unordered."""
        mask = self.pd['gene'].isin(genelist)
        return(mask)

    def barcode_keep12(self):
        """ take the table, cut it into TCGA-XX-XXXX format."""
        cols = self.pd.columns.values
        cols_new = {}
        for i in range(len(cols)):
            if(cols[i][0:4]=='TCGA'):
                cols_new.update( {cols[i] : cols[i][0:12] })
        # be careful of the order
        self.pd.rename(columns = cols_new, inplace=True )
        #return(self)
        
    def barcode_TP_keep12(self,keepcol=['gene']):
        import re
        """ take the table, cut it into TCGA-XX-XXXX format, only for TP samples"""
        cols = self.pd.columns.values
        cols_new = {}
        for i in range(len(cols)):
        	#correct the regex which does not match barcode length 16.
        	#if(re.match('TCGA.+\-01[AB]\-.+',cols[i])):
        	#    new = re.match('TCGA.+\-01[AB]\-.+',cols[i]).group()[0:12]
            if(re.match('TCGA.+\-01[AB]',cols[i])):
                new = re.match('TCGA.+\-01[AB]',cols[i]).group()[0:12]
                cols_new.update( {cols[i] : new })
        # be careful of the order
        self.pd = self.pd[list(keepcol)+list(cols_new.keys())].copy()
        self.pd.rename(columns = cols_new, inplace=True ) 
    
    def rmlowexp(self, cutoff=1):
        """ calculate mean expression and filter """
        mask = self.pd.mean(axis=1,numeric_only=True, skipna=True) >1
        self.pd = self.pd.loc[mask,:]
        #return(self)

    
    def duppat_takemean(self,duppats={}, keepcol=['gene']):
        """combine patients with duplicates mean. duppats: a dict.
        # example of pats: duppats = check_patdupsamp(TPcol)
        """
        for pat in duppats.keys():
            self.pd[pat] = self.pd[duppats[pat]].mean(axis=1)
            for t in duppats[pat]:
                del self.pd[t] 
            self.pd.rename(columns = {pat:duppats[pat][0]}, inplace = True) # still using the first name

    def keep_barcodes(self, keep_list, keepcol=['gene']):
        """# For example: before combining patients with duplicates, ensure only tumor RNA-seq (01) is kept.
        # example of keep_list: keep_list = check_tcga01(rnapd.pd.columns.values)
        """
        #allbarcode_list = list(self.pd.columns.values)
        #keep_list = check_01RNA(allbarcode_list)
        self.pd = self.pd.loc[:, keepcol+ keep_list]

    def __init__(self, pd, name):
        self.pd = pd.copy()
        self.name = name
        
        
def __init__():
    import pandas as pd
    import numpy as np
    import re
