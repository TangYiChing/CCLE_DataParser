"""
Data parser for CCLE
     # filename # source/link
     baseline expression: CCLE_expression.csv, source link: https://depmap.org/portal/download/ (19Q3: log2(RPKM+1))
     baseline expression: CCLE_expression.csv, source link: https://depmap.org/portal/download/ (20Q1: log2(TPM+1))
     mutation: CCLE_mutations.csv, source link: https://depmap.org/portal/download/
     copy number variation: CCLE_gene_cn.csv, source link: https://depmap.org/portal/download/
     cell line: sample_info.csv, source link: https://depmap.org/portal/download/
     compound: , source link:
     drug sensitivity: CCLE_NP24.2009_Drug_data_2015.02.24.csv, source link:
     model list: https://cellmodelpassports.sanger.ac.uk/downloads

"""


#import built-in pkgs
import os
import sys
import numpy as np
import pandas as pd
import itertools as itrs
import pubchempy as pcp
from pubchempy import Compound, get_compounds
#import customized pkgs
import util as ut
import plot_util as pt
import externaldb as exdb

DB_PATH = '/data/DR/db/CCLE/'
DB_FILE = {'MODEL':'model_list_20200204.csv',
           'EXP_tpm':'CCLE_expression.csv', # log2(TPM+1)
           'EXP_rpkm': 'CCLE_expression_log2RPKM.csv', # log2(RPKM+1)
           'CNV':'CCLE_gene_cn.csv',
           'MUT':'CCLE_mutations.csv',
           'RESP':'CCLE_NP24.2009_Drug_data_2015.02.24.csv',
           'CELL': 'sample_info.csv',
           'DRUG':''}

class UseCCLE:
    # initialize 
    def __init__(self, dbPathStr=DB_PATH, dbFileDict=DB_FILE):
        """
        Data parser for CCLE dataset: https://depmap.org/portal/download/

        :param dbPathStr: string representing path to the database
        :param dbFileDict: dict representing file location to each data
        """
        self.folder_str = dbPathStr
        self.file_dict = dbFileDict

    # parse raw data
    def parseRESP(self):
        """
        Read file: CCLE_NP24.2009_Drug_data_2015.02.24.csv

        :return df: dataframe with headers
        """
        f_str = self.folder_str + self.file_dict['RESP']
        df = pd.read_csv(f_str, header=0,  sep=",")
        # clean data
        df.columns = ['CCLE Cell Line Name', 'Primary Cell Line Name', 'Compound', 'Target', 'Doses (uM)',
       'Activity Data (median)', 'Activity SD', 'Num Data', 'FitType',
       'EC50', 'IC50', 'Amax', 'ActArea']
        return df

    def parseCELL(self):
        """
        Read file: sample_info.csv

        :return df: dataframe with headers        
        """
        f_str = self.folder_str + self.file_dict['CELL']
        df = pd.read_csv(f_str, header=0, sep=",")
        return df

    def parseEXP(self, use='TPM'):
        """
        Read file:
            EXP_tpm: CCLE_expression.csv,           # log2(TPM+1)
            EXP_rpkm: CCLE_expression_log2RPKM.csv, # log2(RPKM+1)

        :param use: string representing data type, options=[ TPM | RPKM ], default:TPM
        :return df: dataframe with headers
        """
        if use == 'TPM':
            f_str = self.folder_str + self.file_dict['EXP_tpm']
        elif use == 'RPKM':
            f_str = self.folder_str + self.file_dict['EXP_rpkm']
        else:
            print( 'ERROR: use=[TPM | RPKM], got={:}'.format(use) )
            sys.exit(1) 
        df = pd.read_csv(f_str, header=0, index_col=0, sep=",")
        return df

    def parseCNV(self):
        """
        Read file: CCLE_gene_cn.csv

        :return df: dataframe with headers
        """
        f_str = self.folder_str + self.file_dict['CNV']
        df = pd.read_csv(f_str, header=0, index_col=0, sep=",")
        return df

    def parseMUT(self):
        """
        Read file: CCLE_mutations.csv

        :return df: dataframe with headers
        """
        f_str = self.folder_str + self.file_dict['MUT']
        df = pd.read_csv(f_str, header=0, index_col=0, sep=",", low_memory=False)
        return df

    # retrieve processed data
    def getRESP(self):
        """
        Retrieve processed response data
        :return df: dataframe with headers

        Note:
        =====
        Due to cname:cid replacement, output of getRESP() may have fewer cells than output of  parseRESP()
        """
        # load data
        resp = self.parseRESP()
        # load data dictionary to replace cname with cid
        cid_dict = self._replace_cell_name(use=['resp'])
        # clean data
        cname_list = sorted(list(set(resp['Primary Cell Line Name']) & set(cid_dict['resp'].keys())))
        resp = resp.loc[resp['Primary Cell Line Name'].isin(cname_list)]
        # replace cname with cid
        df = resp.copy()
        df['Primary Cell Line Name'] = df['Primary Cell Line Name'].replace(to_replace=cid_dict['resp'])
        # add LN_IC50
        df['LN_IC50'] = np.log( df['IC50'] )
        # return
        return df 
   
    def getMODEL(self):
        """
        Retrieve model information.

        :return df: dataframe with headers

        Note:
        =====
        Return dataframe was merged from two files:
            1. CCLE's sample_info.csv
            2. Cell Model Passport's model_list_20200204.csv
        This method calls the retrieveMergedMODEL() method of externaldb.py
        """
        df = exdb.retrieveMergedMODEL(db='CCLE')
        return df

    def getDRUG(self, use=['target', 'smile', 'resp']):
        """
        Retrieve drug that has data in given data types

        :param use: list containing data types, options=[target | smile], default: [target | smile]
        :return inter_drugname_list: list containing drug name that has data in given data types
        """
        # load data
        target = self.getCompoundTARGET(use_exdb=True)
        smile = self.getCompoundSMILE(use_exdb=True)
        resp = self.getRESP()
        drugname_dict = { 'target': list(target['drug'].unique()), 'smile': list(smile['drug'].unique()),
                          'resp': list(resp['Compound'].unique()) }
        # get intersection drug names
        use_drugname_list = []
        for data_str in use:
            use_drugname_list.append( set(drugname_dict[data_str]) )
            #print(data_str, len(drugname_dict[data_str]))
        inter_drugname_list = sorted(list( set.intersection( *use_drugname_list ) ))
        #print( 'data used={:}, shared drugs={:}'.format(use, len(inter_drugname_list)) )
        return inter_drugname_list

    def getCELL(self, use=['exp', 'cnv', 'mut', 'resp']):
        """
        Retrieve cell that has data in given data types

        :param use: list containing data types, options=[exp | cnv | mut | resp], default: [exp, cnv, mut, resp]
        :return inter_modelID_list: list containing cell modelIDs that has data in given data types
        """
        # load data 
        exp = self.getEXP(use='TPM')
        cnv = self.getCNV()
        mut = self.parseMUT()
        resp = self.getRESP() # Primary Cell Line Name is ACH-ID
        modelID_dict = { 'exp': exp.index.tolist(), 'cnv': cnv.index.tolist(), 
                         'mut': list(mut['DepMap_ID'].unique()), 'resp': list(resp['Primary Cell Line Name'].unique()) }
        # get intersection cell lines
        use_modelID_list = []
        for data_str in use:
            use_modelID_list.append( set(modelID_dict[data_str]) )
            #print(data_str, len(modelID_dict[data_str]))
        inter_modelID_list = sorted(list( set.intersection( *use_modelID_list ) ))
        #print( 'data used={:}, shared cells={:}'.format(use, len(inter_modelID_list)) )
        return inter_modelID_list

    def getCELLMap(self, cellList, colStr):
        """
        Retrieve cell mapping for the given cells

        :param cellList: list containing a list of cells (e.g., model_id, model_name)
        :param colStr: string representing column name of the model where itemList can be retrieved.
                       options=[model_id, model_name, BROAD_ID, RRID, COSMIC_ID]
        :return df: dataframe contains headers

        Note:
        =====
        This method calls the retrieveMergedMODEL() method of externaldb.py
        """
        # load data
        model = self.getMODEL(db='CCLE')
        # check if colStr in model
        if not colStr in model.columns:
            print( 'ERROR: {:} not in model, options={:}'.format(colStr, '[model_id, model_name, BROAD_ID, RRID, COSMIC_ID]') )
            sys.exit(1)
        # return
        found_in_model = list(set(cellList) & set(model[colStr]))
        if len(found_in_model) != len(cellList):
            not_found_in_model = list(set(cellList) - set(model[colStr]))
            print( 'WARNING: Not Found={:},\n    {:}'.format(len(not_found_in_model), not_found_in_model) )
            df = model.loc[model[colStr].isin(cellList)]
        else:
            df = model.loc[model[colStr].isin(cellList)]
        return df

    def getMUT(self, cellList=None, to_dict=False):
        """
        Retrieve cell:mutation dictionary for the given cells

        :param cellList: list containing a list of cells (i.e., DepMap_ID)
        :param to_dict: boolean indicating to output dictionary or dataframe
        :return cell_mutList_dict: dictionary containing cell: mutation_list pairs

        Note:
        =====
        if cellList=None, the program will return all cells available in parseMUT()
        """
        # load data
        df = self.parseMUT()
        # clean data
        use_cols = ['DepMap_ID', 'Hugo_Symbol']
        df = df[use_cols].copy()
        df = df.drop_duplicates(keep='first')
        # get cellList
        if cellList == None:
            cellList = df['DepMap_ID'].values.tolist()
        else:
            cellList = sorted(list( set(cellList) & set(df['DepMap_ID']) ))
        # subsetting to include only cellList
        use_df = df.loc[df['DepMap_ID'].isin(cellList)]
        # create dict
        if to_dict == True:
            cell_mutList_dict = { cid: gnm['Hugo_Symbol'].values.tolist() for cid, gnm in use_df.groupby('DepMap_ID') }
            outputs = cell_mutList_dict
        else: 
            outputs = use_df
        # return
        return outputs

    def getCNV(self, cellList=None, categorical=False):
        """
        Retrieve copy number for the given cells

        :param cellList: list containing a list of cells (i.e., DepMap_ID)
        :param categorical: boolean indicating whether to return log2 CNV or one-hot CNV 
        :return df: dataframe with DepMap_ID by Gene

        Note:
        =====
        1. if cellList=None, the program will return all cells available in parseCNV()
        2. if categorical=True, the program will return one-hot matrix, that is, 
              "genes with focal CNV values between and including -0.3 and 0.3 are categorized as "neutral" (0), 
               otherwise 1. (see reference)"

        Reference:
        ==========
        https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/CNV_Pipeline/
     
        """
        # load data
        raw_cnv = self.parseCNV()
        # clean data
        gene_list = [col.split(' ')[0] for col in raw_cnv.columns]
        gid_list =  [int(col.split(' ')[1].split('(')[1].split(')')[0]) for col in raw_cnv.columns]
        raw_cnv.columns = gene_list
        # select by cell ids
        if cellList == None:
            cellList = raw_cnv.index.tolist()
        else:
            cellList = sorted(list( set(cellList) & set(raw_cnv.index.tolist()) ))
        # subsetting
        df = raw_cnv.loc[cellList]

        # return
        if categorical:
            df[(df >= -0.3) & (df <= 0.3)] = 0
            df[(df < -0.3) | (df > 0.3)] = 1
            
        return df

    def getCompoundSMILE(self, use_exdb=True, drugList=None):
        """
        Retrieve drug chemical structure: SMILE

        :param drugList: list containing  a list of drug name
        :param use_exdb: boolean indicating whether use external database or not
        :return df: dataframe contains headers=[drug, smile]

        Note:
        =====
        This method calls the retrieveSMILE() method of externaldb.py
        If drugList is None, this program will use drug name available
        """
        # load data
        resp = self.parseRESP()
        drug_list = list( resp['Compound'].unique() )
        # use external database
        drug_df = exdb.retrieveSMILE(drug_list, db='PubChem') # drug, smile
        if drug_df.shape[0] < len(drug_list):
            print( '{:}/{:} have SMILE string'.format(len(drug_df), len(drug_list)) )

        # select by drug name 
        if drugList == None:
            drugList = drug_df['drug'].values.tolist()
        df = drug_df.loc[drug_df['drug'].isin(drugList)]
        # return
        if df.shape[0] < len(drugList):
            print( 'these drugs do not have SMILE string: {:}'.format(list(set(drugList)-set(df['drug']))) )
        return df #headers=[drug, smile]

    def getCompoundTARGET(self, use_exdb=True, drugList=None):
        """
        Retrieve drug target gene/protein

        :param drugList: list containing  a list of drug name
        :param use_exdb: boolean indicating whether use external database or not
        :return df: dataframe contains headers=[drug, target]

        Note:
        =====
        This method calls the retrieveTARGET() method of externaldb.py
        If drugList is None, this program will use drug name available
        """
        # load data
        df = self.parseRESP()
        use_cols = ['Compound', 'Target']
        df = df[use_cols]
        # get targets from raw data
        df.columns = ['drug', 'target']
        target_df = df.drop_duplicates(keep='first')
        # check condition
        if use_exdb:
            # retrieve targets from external database
            map_dict = {'AEW541':'NVP-AEW541', 'ZD-6474':'Vandetanib', 'PF2341066':'Crizotinib'}
            target_df.loc[:, 'drug'] = target_df.loc[:, 'drug'].replace(to_replace=map_dict) #modify drug name to get more matches
            ex_target_df = exdb.retrieveTARGET(target_df['drug'].values.tolist(), db='DGIdb')
            # return
            drug_df = pd.concat([target_df, ex_target_df], axis=0, sort='True') # merge 
            map2_dict = { value:key for key, value in map_dict.items() } #change back to original drug name
            drug_df['drug'] = drug_df['drug'].replace(to_replace=map2_dict)
            drug_df.drop_duplicates(inplace=True) # remove duplicated drug-target pair
            # select by gene name
            if drugList == None:
                drugList = drug_df['drug'].values.tolist()
            df = drug_df.loc[drug_df['drug'].isin(drugList)]

        else:
            # return
            df = target_df
            df.drop_duplicates(inplace=True) # remove duplicated drug-target pair
            # select by gene name
            if drugList == None:
                drugList = drug_df['drug'].values.tolist()
            df = drug_df.loc[drug_df['drug'].isin(drugList)]

        df = df.dropna(axis=0)
        return df # headers = ['drug', 'target']

    def getEXP(self, use='TPM', cellList=None):
        """
        Retrieve expression for the given cells

        :param cellList: list containing a list of cells (i.e., DepMap_ID)
        :return df: dataframe with DepMap_ID by Gene

        Note:
        =====
        if cellList=None, the program will return all cells available in parseEXP()
        """
        # load data
        if use == 'TPM':
            raw_exp = self.parseEXP(use='TPM')
        elif use == 'RPKM':
            raw_exp = self.parseEXP(use='RPKM')
        else:
            print( 'ERROR: use=[TPM | RPKM], got={:}'.format(use) )
            sys.exit(1)
        # clean data
        gene_list = [col.split(' ')[0] for col in raw_exp.columns]
        gid_list =  [int(col.split(' ')[1].split('(')[1].split(')')[0]) for col in raw_exp.columns]
        raw_exp.columns = gene_list
        # select by cell ids
        if cellList == None:
            cellList = raw_exp.index.tolist()
        else:
            cellList = sorted(list( set(cellList) & set(raw_exp.index.tolist()) ))
        # return
        df = raw_exp.loc[cellList]
        return df

    def getCIDMAPDICT(self, keyStr, valueStr):
        # load data
        model = self.getMODEL()
        # check columns
        for col in [keyStr, valueStr]:
            if not col in model.columns:
                print( 'ERROR: {:} not in model'.format(col) )
                sys.exit(1)
        # get mapping dict
        map_dict = dict( zip(model[keyStr].values.tolist(), model[valueStr].values.tolist()) )
        return map_dict

    def _replace_cell_name(self, use=['resp']):
        """
        Return DepMap_ID for cell name in given data

        :param use: list containing data type
        :return cid_dict: dictionary containing model_id for each given data
        """
        # cheeting set
        resp_dict = { 'Raji':'ACH-000654', 'HCT 116':'ACH-000971', 'FaDu':'ACH-000846', 'U-87 MG':'ACH-000075', 
                      'Hep 3B2.1-7':'ACH-000625', 'Hs 739.T':'ACH-000413', 'LC-1/sq-SF':'ACH-001113', 'Hs 729':'ACH-000133', 
                      'Hs 695T':'ACH-000799', 'Hs 578T':'ACH-000148', 'Hs 852.T':'ACH-000274', 'Hs 895.T':'ACH-000043', 
                      'U-118 MG':'ACH-000040', 'KP-N-SI9s':'ACH-000446', 'Hep G2':'ACH-000739', 'Calu-6':'ACH-000264', 
                      'IST-MES2':'ACH-000331', 'Hey-A8':'ACH-000542', 'Daoy':'ACH-000211', 'CAL 27':'ACH-000832', 
                      'C2BBe1':'ACH-000009', 'L3.3':'ACH-000685', 'huH-1':'ACH-000475', 'Caki-2':'ACH-000234', 
                      'Saos-2':'ACH-000410', 'NIH:OVCAR-3':'ACH-000001', 'MIA PaCa-2':'ACH-000601', 'Panc 02.03':'ACH-000042', 
                      'P31/FUJ':'ACH-000770', 'Hs 766T':'ACH-000178', 'SW 1573':'ACH-000677', 'Malme-3M':'ACH-000477', 
                      'Panc 03.27':'ACH-000139', 'Hs 746T':'ACH-000616', 'Hs 294T':'ACH-000014', 'Calu-3':'ACH-000392', 
                      'SW 780':'ACH-000384', 'SW 900':'ACH-000669', 'Reh':'ACH-000960', 'LOX IMVI':'ACH-000750', 
                      'U-2 OS':'ACH-000364', 'SCaBER':'ACH-000839', 'HuCCT1':'ACH-000976', 'OC 315':'ACH-001144', 
                      'DU 145':'ACH-000979', 'SW 1271':'ACH-000890', 'Hs 944.T':'ACH-000632', 'COLO 201':'ACH-000253', 
                      'HuT 78':'ACH-000509', 'Capan-2':'ACH-000107', 'Ishikawa (Heraklio) 02 ER-':'ACH-000961', 
                      'TYK-nu':'ACH-000430', 'SW 1353':'ACH-000418', 'OC 316':'ACH-001145', 'HEL 92.1.7':'ACH-000005', 
                      'PLC/PRF/5':'ACH-001318', 'Hs 229.T':'ACH-000131', 'DMS 114':'ACH-000530', 'SW 1088':'ACH-000437', 
                      'TE 617.T':'ACH-000051', 'BxPC-3':'ACH-000535', 'Panc 04.03':'ACH-000235', 'Hs 683':'ACH-000067', 
                      'SK-N-BE(2)':'ACH-000312', 'OC 314':'ACH-000962', 'Hs 936.T':'ACH-000801', 'COLO 205':'ACH-001039', 
                      'AN3 CA':'ACH-000940', 'HCC2935':'ACH-000150', 'Detroit 562':'ACH-000207', 'Hs 840.T':'ACH-000284', 
                      'AsPC-1':'ACH-000222', 'MOR/CPR':'ACH-000851', 'Hs 939.T':'ACH-000814',
                      'NCIH1184':'ACH-000523', '22Rv1':'ACH-000956', 'Panc 10.05':'ACH-000060',
                      'SU.86.86':'ACH-000114', 'HPAF-II':'ACH-000094', 'HT-144':'ACH-000322',
                      'SK-MM-2':'ACH-000363', 'RPMI-8402':'ACH-000636', 'TE-9':'ACH-000694',
                      'MKN-45':'ACH-000356', 'SW 1990':'ACH-000155', 'COR-L23':'ACH-000662',
                      'SK-N-SH':'ACH-000149', '8-MG-BA':'ACH-000137', 'UM-UC-3':'ACH-000522',
                      'KYSE-450':'ACH-000865', 'MPP 89':'ACH-000319', 'JHH-2':'ACH-000577',
                      'SK-N-DZ':'ACH-000366', 'KP-3':'ACH-000108', 'UACC-257':'ACH-000579',
                      'AMO-1':'ACH-000838', 'SCC-9':'ACH-000181', 'MEC-1':'ACH-000405',
                      '786-O':'ACH-000649', 'OVCAR-4':'ACH-000617', 'GMS-10':'ACH-000102', 
                      'MES-SA':'ACH-000449', 'HT-1197':'ACH-000547', 'MFE-296':'ACH-000879',
                      'OVCAR-8':'ACH-000696', 'GB-1':'ACH-000738', 'SUP-M2':'ACH-000226',
                      'SK-MEL-2':'ACH-001190', 'MB 157':'ACH-000621', 'MeWo':'ACH-000987',
                      'COLO-205':'ACH-001039', 'CCK-81':'ACH-000963', 'MDA-MB-436':'ACH-000573',
                      'B-CPAP':'ACH-000456', 'GLC-82':'ACH-001071', 'HT-1080':'ACH-000054',
                      'OV-90':'ACH-000291', 'KYM-1':'ACH-000607', 'A-204':'ACH-000201',
                      'BGC-823':'ACH-001017', 'A-375':'ACH-000219','Ishikawa (Heraklio) 02 ER':'ACH-000961',
                      'SW-990':'ACH-000669', 'DU-145':'ACH-000979', 'SW-1353':'ACH-000418', 'CAL-27':'ACH-000832'}
        # program start
        cid_dict = {'resp':{}}
        # load data
        resp = self.parseRESP()
        model = self.getMODEL() # DepMap_ID, stripped_cell_line_name
        ## cell name in resp is a striped string 
        ## use cname to corresponding cid
        not_found_cname_list = []
        for cname_str in set(resp['Primary Cell Line Name']):
            # use rules to modify string, in order to get more matches automatically
            if '-' in cname_str:
                # try to match
                if cname_str in model['stripped_cell_line_name'].values.tolist():
                    cid_dict['resp'].update( {cname_str: model.loc[model['stripped_cell_line_name']==cname_str]['DepMap_ID'].values[0]} )
                else:
                    # use modified str to match
                    mcname_str = ''.join( cname_str.split('-') )#[0] + cname_str.split('-')[1]
                    # try to match AGAIN
                    if mcname_str in model['stripped_cell_line_name'].values.tolist():
                        cid_dict['resp'].update( {cname_str: model.loc[model['stripped_cell_line_name']==mcname_str]['DepMap_ID'].values[0]} )
                    else:
                        # use cheet dict
                        if cname_str in resp_dict.keys():
                            cid_dict['resp'].update( {cname_str:resp_dict[cname_str]} )
                        else:
                            not_found_cname_list.append(cname_str)

            elif ' ' in cname_str:
                # try to match
                if cname_str in model['stripped_cell_line_name'].values.tolist():
                    cid_dict['resp'].update( {cname_str: model.loc[model['stripped_cell_line_name']==cname_str]['DepMap_ID'].values[0]} )
                else:
                    # use modified str to match
                    mcname_str = '-'.join( cname_str.split(' ') )#[0] + '-' + cname_str.split(' ')[1]
                    # try to match AGAIN
                    if mcname_str in model['stripped_cell_line_name'].values.tolist():
                        cid_dict['resp'].update( {cname_str: model.loc[model['stripped_cell_line_name']==mcname_str]['DepMap_ID'].values[0]} )
                    else:
                        # use cheet dict
                        if cname_str in resp_dict.keys():
                            cid_dict['resp'].update( {cname_str:resp_dict[cname_str]} )
                        else:
                            not_found_cname_list.append(cname_str)

            else:
                # try to match
                if cname_str in model['stripped_cell_line_name'].values.tolist():
                    cid_dict['resp'].update( {cname_str: model.loc[model['stripped_cell_line_name']==cname_str]['DepMap_ID'].values[0]} )
                else:
                    # use cheet dict
                    if cname_str in resp_dict.keys():
                        cid_dict['resp'].update( {cname_str:resp_dict[cname_str]} )
                    else:
                        not_found_cname_list.append(cname_str)
        # send warning message
        if len(not_found_cname_list) > 0:
            print( 'WARNING: cells in resp not_found DepMap_ID={:}'.format(not_found_cname_list) )
        return cid_dict

    def getContRESPDICT(self, use='IC50', drugList=None, cellList=None):
        """
        :param use: options = [IC50 | LN_IC50 | ActArea]
        :param drugList: list containing  a list of drug name
        :param cellList: list containing a list of cells (e.g., model_id, model_name)
        :return respMat_dict: dictionary containing respMat of drug by mode_id

        """
        # load data
        resp = self.getRESP() #self._replace_cell_name(use=['resp'])
       
        # select data
        use_dict = { 'IC50': 'IC50', 'EC50':'EC50', 'ActArea': 'ActArea', 'LN_IC50': 'LN_IC50'}
        use_cols =  ['Primary Cell Line Name', 'Compound'] + [ use_dict[use] ]
        
        # indexing
        df1 = resp.set_index('Compound').sort_index()
        df2 = resp.set_index(['Compound', 'Doses (uM)']).sort_index()
        # collecting wanted drugs  by checking duplication
        # Branching based on dose range
        BASE_list = [] # If a drug has no duplicated cells
        EXP1_list = [] # If a drug-cell pair has duplicated dose range,
        EXP2_list = [] # If a drug-cell pair has different dose range
        for idx in set(df2.index):
            drug = idx[0]
            # collect drugs and different dose ranges
            drug_df = resp.loc[resp['Compound']==drug][['Primary Cell Line Name', 'Compound', 'Doses (uM)', use_dict[use]]]
            has_dup = drug_df['Primary Cell Line Name'].duplicated()
            if len(drug_df[has_dup]) == 0:
                BASE_list.append(drug)
            else:
                drug_dose_df = df2.loc[idx]
                has_dup_dose = drug_dose_df['Primary Cell Line Name'].duplicated()
                if len(has_dup_dose) == 0: 
                    #print( 'compound={:} {:} has duplicated cells, but different dose, appending to exp2 for separate set'.format(drug, idx) )
                    EXP2_list.append(drug)
                else:
                    #print( 'compound={:} {:} has duplicated cells and dose, appending to exp1 for averaging them'.format(drug, idx) )
                    EXP1_list.append(drug)
        # retrieve data and check duplication by averaging them
        print( 'BASE_list={:}, EXP1_list={:}, EXP2_list={:}'.format(len(BASE_list), len(EXP1_list), len(EXP2_list)) )
        exp_list = [EXP1_list, EXP2_list]
        exp_df_dict = {}
        for i in range(len(exp_list)):
            exp_int = i+1
            if len(exp_list[i]) != 0:  
                # load data
                wanted_list = BASE_list + exp_list[i]
                use_cols += ['Doses (uM)']
                exp_df = df1.loc[wanted_list].reset_index()[use_cols]
                # take average of duplicated data
                exp_df = exp_df.groupby(['Compound', 'Primary Cell Line Name', 'Doses (uM)']).mean().reset_index()
                exp_df = exp_df[['Compound', 'Primary Cell Line Name', use]]
                exp_df_dict.update( {'exp'+str(exp_int):exp_df} )
        # converting drug-cell matrix
        respMat_dict = {}
        for key, df in exp_df_dict.items():
            # select by drug name
            if drugList != None:
                #drugList = df.index.tolist()
                df = df.loc[df['Compound'].isin(drugList)]
            # select by cell name
            if cellList != None:
                #cellList = df.columns.tolist()
                df = df.loc[df['Primary Cell Line Name'].isin(cellList)]

            # long to wide
            respMat = df.pivot(index='Compound', columns='Primary Cell Line Name', values=use) # this could create hole in the matrix
            respMat.fillna(np.nan, inplace=True) # make holes with np.nan
            respMat.columns.name = 'DepMap_ID'
            # save to dict
            respMat_dict.update( {key:respMat} )
        return respMat_dict # key:value = exp1:respMat (i.e., drug-cell response matrix)

    def getDeltaResponse(self, respMat, by='cell'):
        """
        Return delta response of pairs

        :param respMat: matrix with Compound by DepMap_ID
        :param by: string representing key of return delta_dict, options=['cell', 'drug'], default=cell
        :return delta_dict: dictionary containing delta response of pairs

        Note:
        =====
        if by='cell', the program will calculate delta response of two drugs for each cell
        if by='drug', the program will calculate delta response of two cells for each drug
        """
        # create result dict
        delta_dict = {}
        # condition
        if by == 'cell':
            respMat = respMat.T
        elif by == 'drug':
            respMat = respMat
        else:
            print( 'ERROR: by={:}, options=[cell|drug]'.format(by) )

        # all-pair list
        pair_list = [ pair for pair in itrs.combinations(respMat.columns, 2) ]

        # looping rows and calculate delta response
        for idx in respMat.index: # if by=cell, idx will be cell id, if by=drug, idx will be drug name
            # load data
            idx_df = respMat.loc[idx] # Series
            idx_df = idx_df.to_frame(name='resp') # dataframe
            idx_df.dropna(axis=0, inplace=True)
            # calculate pairwise difference among possible column pair (i.e., cellpair, drugpair)
            delta_difference_df = pd.DataFrame(np.abs( idx_df['resp'].values - idx_df['resp'].values[:, None] ),
                                               columns = idx_df.index.tolist(), index = idx_df.index.tolist())
            # convert to long-form df
            delta_df = ut.sym2half(delta_difference_df, keepSelfPairs=False) # columns = [idxpair, similarity]
            # change column name
            if by == 'cell':
                delta_df.columns = ['drugpair', 'delta response']
            else:
                delta_df.columns = ['cellpair', 'delta response']
            # update to result
            delta_dict.update( {idx: delta_df} )
        # return
        #print(delta_dict)
        return delta_dict


    def respMat2respDf(self, respMat):
        """
        Return delta response of cell-pairs again certain drug

        :param respMat: matrix with Compound by DepMap_ID
        :return delta_dict: dictionary containing delta response of cell-pairs for each compund in respMat

        Note:
        =====
        1. delta response of cell-pairs is defined by response value of one cell - response value of the other cell
                for example: ACH-000001 has response value of 0.245932 against Compound 17-AAG, and
                             ACH-000005 has response value of 1.616860 against Compound 17-AAG, therefore
                             delta response of ACH-000001-ACH-000005 pair is 0.245932 - 1.616860 = -1.370928
        2. If any of cell-pair has one na, the cell pair will not be saved into result
        """
        # create result dict
        delta_dict = {}
        # create cell-pair list
        cellpair_list = [ pair for pair in itrs.combinations(respMat.columns, 2) ]
        # loop through drug and calculate delta response
        for drug in respMat.index:
            # load data
            drug_df = respMat.loc[drug] # Series
            drug_df = drug_df.to_frame(name='resp') # dataframe
            drug_df.dropna(axis=0, inplace=True)
            # calculate pairwise difference among possible column pair (i.e., cellpair)
            delta_difference_df = pd.DataFrame(np.abs( drug_df['resp'].values - drug_df['resp'].values[:, None] ),
                                               columns = drug_df.index.tolist(), index = drug_df.index.tolist())
            # convert to long-form df
            delta_df = ut.sym2half(delta_difference_df, keepSelfPairs=False) # columns = [cellpair, similarity]
            delta_df.columns = ['cellpair', 'delta response']
            # update to result
            delta_dict.update( {drug: delta_df} )
        # return
        return delta_dict

    def queryRESP(self, drugList, cellList, use='IC50'):
        resp = self.getRESP()
        df = resp[ resp['Compound'].isin(drugList) & resp['Primary Cell Line Name'].isin(cellList) ]
        print(df)
# execute from script for debugging
if __name__ == "__main__":
    ccle = UseCCLE() # initiate an instance
    # set choices
    test = True
    save = False
 
    if test:
        # TESTING
        #ccle.queryRESP(['17-AAG'], ['ACH-000012'])
        #df = pd.read_csv('/repo4/ytang4/PHD/db/CCLE/processed/ccle.common.drug.txt', header=0, sep="\t")
        #drug_list = df['drug_name'].values.tolist() 
        #target = ccle.getCompoundTARGET(use_exdb=True, drugList=None)
        #target.to_csv('./CCLE.TARGET.GeneName.Mat.txt',header=True, index=False, sep="\t")
        #print('missing={:}'.format(set(drug_list)-set(target['drug'])))
        #resp = ccle.getRESP()
        #resp['MaxConc'] = np.log(8)
        #use_resp = resp[['Compound', 'MaxConc']] 
        #use_resp = use_resp.loc[use_resp['Compound'].isin(drug_list)]
        #use_resp = use_resp.drop_duplicates(keep='first')
        #use_resp.columns = ['drug', 'max_concentration']
        #use_resp = use_resp.set_index('drug')
        #use_resp = use_resp.loc[drug_list]
        #print(use_resp)
        cell_list = ccle.getCELL(use=['exp', 'cnv', 'mut', 'resp'])
        cnv_c = ccle.getCNV(cellList=cell_list, categorical=True)
        df_list = []
        for idx in cnv_c.index:
            data = cnv_c.loc[[idx]].T
            ones = data[data[idx]==1]
            genes = ones.index.tolist()
            df = pd.DataFrame({'cell':[idx]*len(genes), 'gene':genes})
            df_list.append(df)
        cnv = pd.concat(df_list, axis=0)
        #cnv['gene'] = cnv['gene'].replace(gid_gnm_dict)
        print(cnv.head())
        #use_resp.to_csv('/repo4/ytang4/PHD/db/CCLE/processed/CCLE.DRUG.COMMON.LN_MaxConc.Mat.txt', header=True, index=True, sep="\t")
        DB_PATH = '/data/DR/db/CCLE'
        cnv.to_csv(DB_PATH+'/processed/CCLE.CNV.GeneName.Mat.txt',header=True, index=False, sep="\t")

        #resp_dict1 = ccle.getContRESPDICT(use='IC50', drugList=None, cellList=None)
        #resp_dict2 = ccle.getContRESPDICT(use='EC50', drugList=None, cellList=None)
        #resp_dict3 = ccle.getContRESPDICT(use='LN_IC50', drugList=None, cellList=None)
        #resp_dict4 = ccle.getContRESPDICT(use='ActArea', drugList=None, cellList=None)
        #resp_dict1['exp1'].to_csv(DB_PATH+'/CCLE/CCLE.RAW.RESP.IC50.exp1.Mat.txt',header=True, index=True, sep="\t")
        #resp_dict2['exp1'].to_csv(DB_PATH+'/CCLE/CCLE.RAW.RESP.EC50.exp1.Mat.txt',header=True, index=True, sep="\t")
        #resp_dict3['exp1'].to_csv(DB_PATH+'/CCLE/CCLE.RAW.RESP.LN_IC50.exp1.Mat.txt',header=True, index=True, sep="\t")
        #resp_dict4['exp1'].to_csv(DB_PATH+'/CCLE/CCLE.RAW.RESP.ActArea.exp1.Mat.txt',header=True, index=True, sep="\t")


    if save:
        # path
        DB_PATH = '/data/DR/db/CCLE' 
        # get cell, drug
        cell_list = ccle.getCELL(use=['exp', 'cnv', 'mut', 'resp'])
        drug_list = ccle.getDRUG(use=['target', 'smile', 'resp'])
        # generate processed data
        smile = ccle.getCompoundSMILE(use_exdb=True, drugList=drug_list)
        target = ccle.getCompoundTARGET(use_exdb=True, drugList=drug_list)
        exp = ccle.getEXP(use='TPM', cellList=cell_list)
        cnv_c = ccle.getCNV(cellList=cell_list, categorical=True)
        cnv = ccle.getCNV(cellList=cell_list, categorical=False)
        mut = ccle.getMUT(cellList=cell_list, to_dict=False)
        resp_dict1 = ccle.getContRESPDICT(use='IC50', drugList=drug_list, cellList=cell_list)
        #resp_dict2 = ccle.getContRESPDICT(use='EC50', drugList=drug_list, cellList=cell_list)
        #resp_dict3 = ccle.getContRESPDICT(use='LN_IC50', drugList=drug_list, cellList=cell_list)
        #resp_dict4 = ccle.getContRESPDICT(use='ActArea', drugList=drug_list, cellList=cell_list)
        model = ccle.getMODEL()
        # save files to DP_PATH/processed/ folder
        smile.to_csv(DB_PATH+'/processed/CCLE.TARGET.SMILE.Mat.txt',header=True, index=True, sep="\t")
        target.to_csv(DB_PATH+'/processed/CCLE.TARGET.GeneName.Mat.txt',header=True, index=False, sep="\t")
        exp.to_csv(DB_PATH+'/processed/CCLE.EXP.TPM.Mat.txt',header=True, index=True, sep="\t")
        cnv_c.to_csv(DB_PATH+'/processed/CCLE.CNV.OneHot.Mat.txt',header=True, index=True, sep="\t")
        cnv.to_csv(DB_PATH+'/processed/CCLE.CNV.Log2.Mat.txt',header=True, index=True, sep="\t")
        mut.to_csv(DB_PATH+'/processed/CCLE.MUT.GeneName.Mat.txt',header=True, index=False, sep="\t")
        resp_dict1['exp1'].to_csv(DB_PATH+'/processed/CCLE.RESP.IC50.exp1.Mat.txt',header=True, index=True, sep="\t")
        #resp_dict2['exp1'].to_csv(DB_PATH+'/processed/CCLE.RESP.EC50.exp1.Mat.txt',header=True, index=True, sep="\t")
        #resp_dict3['exp1'].to_csv(DB_PATH+'/processed/CCLE.RESP.LN_IC50.exp1.Mat.txt',header=True, index=True, sep="\t")
        #resp_dict4['exp1'].to_csv(DB_PATH+'/processed/CCLE.RESP.ActArea.exp1.Mat.txt',header=True, index=True, sep="\t")
        model.to_csv(DB_PATH+'/processed/CCLE.MODEL.Annotation.Mat.txt',header=True, index=False, sep="\t")
        # display print message
        print('save processed files to {:}'.format(DB_PATH+'/processed/'))
        print('    #common drugs (combo removed) ={:}, common cells={:}'.format(len(drug_list), len(cell_list)))
        print('    resp (drug,cell)={:}, target (drug,target)={:}, smile (drug, smile)={:}'.format(resp_dict1['exp1'].shape, target.shape, smile.shape))
        print('    exp (cell,gene)={:}, cnv (cell,gene)={:}, mut(cell,gene)={:}'.format(exp.shape, cnv.shape, mut.shape))
