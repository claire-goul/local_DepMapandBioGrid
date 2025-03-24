## WRITTEN BY CLAIRE GOUL AUG 2022
## MIT LICENSE
## Revised to fix command line argument handling and empty edgelist edge case
import numpy as np
import pandas as pd
import argparse
##get_correlations_edgelist:
#INPUTS 
#genes: (excel file) with a column titled 'Gene' with list of genes of interest
#links_filtered (excel file) from generate_corrected_coessentiality downloaded from achilles website https://depmap.org/portal/download/all/
#threshold (float in range of 0-1), e.g. 0.2: number that correlation score has to be greater than 
#corrpos (Boolean): if True, get only positive correlation genes; if false, get only negative corr genes
#num (int): number of correlated genes you want
#OUTPUT
#corr: correlation matrix with two columns, 'Gene' and 'Gene1' and their correlation scores ('corrscore')

def get_correlations_edgelist(genes,links_filtered,threshold,corrpos,num):
        links_filtered_newfinal=pd.merge(links_filtered, genes_of_interest,on=['Gene']) #filter by  genes of interest 
        
        # Check if initial merge produced empty result
        if links_filtered_newfinal.empty:
            print("Warning: No genes from the input list found in the correlation dataset.")
            # Return empty DataFrame with expected columns
            return pd.DataFrame(columns=['Gene', 'Gene1', 'corrscore'])
            
        if corrpos:
                links_filtered_newfinal=links_filtered_newfinal[links_filtered_newfinal['corrscore']>=threshold]#threshold for degree of correlation
                # Check if filtering by threshold resulted in empty DataFrame
                if links_filtered_newfinal.empty:
                    print(f"Warning: No genes have correlation scores >= {threshold}.")
                    return pd.DataFrame(columns=['Gene', 'Gene1', 'corrscore'])
                    
                grouped= links_filtered_newfinal.groupby('Gene')
                toplargestdf=pd.DataFrame()
                for var1,subdict in grouped:
                        sub=subdict[subdict['corrscore']>0]
                        # Check if positive correlations exist for this gene
                        if sub.empty:
                            print(f"Warning: No positive correlations found for gene {var1}.")
                            continue
                        sublargest=sub.nlargest(n=num,columns='corrscore') 
                        toplargestdf=pd.concat([toplargestdf,sublargest])
        else:
                links_filtered_newfinal=links_filtered_newfinal[links_filtered_newfinal['corrscore']<=threshold]#threshold for degree of correlation
                # Check if filtering by threshold resulted in empty DataFrame
                if links_filtered_newfinal.empty:
                    print(f"Warning: No genes have correlation scores <= {threshold}.")
                    return pd.DataFrame(columns=['Gene', 'Gene1', 'corrscore'])
                    
                grouped= links_filtered_newfinal.groupby('Gene')
                toplargestdf=pd.DataFrame()
                for var1,subdict in grouped:
                        sub=subdict[subdict['corrscore']<0]
                        # Check if negative correlations exist for this gene
                        if sub.empty:
                            print(f"Warning: No negative correlations found for gene {var1}.")
                            continue
                        sublargest=sub.nlargest(n=num,columns='corrscore') 
                        toplargestdf=pd.concat([toplargestdf,sublargest])
                        
        # Check if no genes passed all filters
        if toplargestdf.empty:
            print("Warning: No genes passed all correlation filters.")
            return pd.DataFrame(columns=['Gene', 'Gene1', 'corrscore'])
            
        corr=(toplargestdf.reset_index()).drop(['index'],axis=1)
        return corr

##INPUTS
#genes: (excel file) with a column titled 'Gene' with list of genes of interest
#bg = csv file of all Biogrid interactors for human genes (downloaded from https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-4.4.220/ # BIOGRID-MV-Physical)
#filters: list of filters you want. default is empty list ([]). You can look at the bg file to see what the possible filters are; for example ['psi-mi:"MI:0407"(direct interaction)']
##OUTPUT
#edgelist_biogrid: file with two columns of  biogrid interactions, 'InteractorA' and 'InteractorB' and 'tuples' column containing a tuple of those
#also make an option for biogrid only for corr and for hits only 
def get_biogrid_edgelist(genes,bg,filters,numcitations): 
        bg_df_final=bg
        if len(filters)>0:
                bg_df_final=bg_df_final[bg_df_final['Interaction Types'].isin(filters)]
        bg_df_final=bg_df_final.reset_index()
        bg_df_final=bg_df_final.drop(columns=['index'])
        A=list(bg_df_final['Aliases Interactor A'])
        B=list(bg_df_final['Aliases Interactor B'])
        Anewog=[]
        Bnewog=[]
        Aog=list(bg_df_final['Alt IDs Interactor A'])
        Bog=list(bg_df_final['Alt IDs Interactor B'])
        for ele in range(0,len(Aog)):
                if len(Aog[ele])>1:
                        Aalts=[]
                        elesnewA=Aog[ele].split('|')
                        gene=elesnewA[1].split(':')[1].split('(')[0]
                        Aalts.append(gene)
                        Anewog.append(Aalts)
                else:
                        Anewog.append([])
        for ele in range(0,len(Bog)):
                if len(Bog[ele])>1:
                        Balts=[]
                        elesnewB=Bog[ele].split('|')
                        gene=elesnewB[1].split(':')[1].split('(')[0]
                        Balts.append(gene)
                        Bnewog.append(Balts)
                else:
                        Bnewog.append([])
        Anew=[]
        Bnew=[]
        for ele in range(0,len(A)):
                if len(A[ele])>1:
                        Aaliases=[]
                        elesnewA=A[ele].split('|')
                        for i in elesnewA:
                                gene=i.split(':')[1].split('(')[0]
                                Aaliases.append(gene)#add in aliases 
                        Aaliases.append(Anewog[ele][0]) #add in the original alt ID too
                        Anew.append(Aaliases)
                else:
                        Anew.append([])
        for ele in range(0,len(B)):
                if len(B[ele])>1:
                        Baliases=[]
                        elesnewB=B[ele].split('|')
                        for i in elesnewB:
                                gene=i.split(':')[1].split('(')[0]
                                Baliases.append(gene)
                        Baliases.append(Bnewog[ele][0])
                        Bnew.append(Baliases)
                else:
                        Bnew.append([])
        #then filter by genes of interest
        intAfin=[]
        intBfin=[]
        for gene in list(genes['Gene']):
                for a in range(0,len(Anew)):
                        if gene in Anew[a]:# this bg interaction is for a gene in genes of interest
                                intAfin.append(gene)
                                intBfin.append(Bnewog[a][0])#add corresponding B interactor (use the alias gene)
        for gene in list(genes['Gene']):
               for b in range(0,len(Bnew)):
                       if gene in Bnew[b]:# this bg interaction is for a gene in genes of interest
                               intAfin.append(gene)
                               intBfin.append(Anewog[b][0])#add corresponding A interactor (use the alias gene)
        
        # Handle the case where no interactions are found
        if len(intAfin) == 0 or len(intBfin) == 0:
            print("Warning: No BioGRID interactions found matching your criteria.")
            # Return an empty DataFrame with the expected columns
            empty_df = pd.DataFrame(columns=['Gene', 'Gene1', 'tuples', 'bg'])
            return empty_df
            
        edgelist_biogrid=pd.DataFrame()
        edgelist_biogrid['Final IDs Interactor A']=pd.DataFrame(intAfin)
        edgelist_biogrid['Final IDs Interactor B']=pd.DataFrame(intBfin)
        edgelist_biogrid['tuples']=list(zip(edgelist_biogrid['Final IDs Interactor A'],edgelist_biogrid['Final IDs Interactor B']))
        
        # Check if filtering by numcitations would result in an empty DataFrame
        filtered_edgelist = edgelist_biogrid.groupby('tuples').filter(lambda x : len(x)>=numcitations)
        if filtered_edgelist.empty:
            print(f"Warning: No interactions have {numcitations} or more citations. Returning empty DataFrame.")
            empty_df = pd.DataFrame(columns=['Gene', 'Gene1', 'tuples', 'bg'])
            return empty_df
            
        edgelist_biogrid = filtered_edgelist
        edgelist_biogrid=edgelist_biogrid.reset_index()
        edgelist_biogrid=edgelist_biogrid.rename(columns={'Final IDs Interactor A':'Gene','Final IDs Interactor B':'Gene1'})
        genestuplesbiogrid=list(zip(edgelist_biogrid.Gene, edgelist_biogrid.Gene1))
        edgelist_biogrid=edgelist_biogrid.drop_duplicates(subset='tuples',keep='first')
        edgelist_biogrid=edgelist_biogrid.reset_index()
        edgelist_biogrid_final = edgelist_biogrid.drop(edgelist_biogrid[edgelist_biogrid.Gene==edgelist_biogrid.Gene1].index)
        
        # Check if removing self-interactions results in an empty DataFrame
        if edgelist_biogrid_final.empty:
            print("Warning: After removing self-interactions, no BioGRID interactions remain.")
            empty_df = pd.DataFrame(columns=['Gene', 'Gene1', 'tuples', 'bg'])
            return empty_df
            
        edgelist_biogrid_final=edgelist_biogrid_final.drop(columns=['index'])
        edgelist_biogrid_final=edgelist_biogrid_final.reset_index()
        edgelist_biogrid_final['bg']='yes'
        return edgelist_biogrid_final

def merge3(list1, list2,list3):
    merged_list = [(list1[i], list2[i],list3[i]) for i in range(0, len(list1))]
    return merged_list

##LOAD IN BIOGRID DN DEPMAP DATA
links_filtered=pd.read_excel('links_achilles.xlsx')

bg=pd.read_excel('Biogrid_MV-Physical_4.4.243_Human.xlsx')
genes_of_interest=pd.read_excel('genestest.xlsx')## can read in any excel file to filter the correlation matrix by

# Set up command-line argument parsing
parser = argparse.ArgumentParser(description='Generate network from gene correlation and BioGRID interaction data')
parser.add_argument('--threshold', type=float, default=0.2, 
                    help='Correlation threshold (float in range 0-1). Correlations must be greater than this value')
parser.add_argument('--corrpos', type=str, default='True', 
                    help='If True, get only positive correlation genes; if False, get only negative correlation genes')
parser.add_argument('--num', type=int, default=3, 
                    help='Number of correlated genes to include for each gene of interest')
# Modified to properly handle list arguments
parser.add_argument('--filters', nargs='+', default=[], 
                    help='List of filters for BioGRID interactions')
parser.add_argument('--numcitations', type=int, default=2, 
                    help='Minimum number of citations required for BioGRID interactions')

# Parse the arguments
args = parser.parse_args()

# Convert corrpos string to boolean
corrpos_bool = True if args.corrpos.lower() == 'true' else False

# Print arguments for debugging
print("Arguments:")
print(f"  threshold: {args.threshold}")
print(f"  corrpos: {corrpos_bool}")
print(f"  num: {args.num}")
print(f"  filters: {args.filters}")
print(f"  numcitations: {args.numcitations}")

#GET BIOGRID INTERACTIONS / COESSENTIAL GENES FOR GENES IN GENE LIST
corr=get_correlations_edgelist(genes_of_interest, links_filtered, threshold=args.threshold, corrpos=corrpos_bool, num=args.num)
# Pass filters directly - it's already a list
edgelist_biogrid=get_biogrid_edgelist(genes_of_interest, bg, filters=args.filters, numcitations=args.numcitations)

# Ensure empty DataFrames have the necessary columns for merge operations
if corr.empty and 'Gene' not in corr.columns:
    corr = pd.DataFrame(columns=['Gene', 'Gene1', 'corrscore'])

if edgelist_biogrid.empty and 'Gene' not in edgelist_biogrid.columns:
    edgelist_biogrid = pd.DataFrame(columns=['Gene', 'Gene1', 'tuples', 'bg'])

#OPTIONS
#CORR: GET CORR MATRIX FOR ALL GENES IN GENE LIST
# Handle empty correlation matrix when saving
if corr.empty:
    print("Warning: Empty correlation matrix. Creating a minimal file with column headers only.")
    pd.DataFrame(columns=['Gene', 'Gene1', 'corrscore']).to_excel('genes_corr.xlsx')
else:
    corr.to_excel('genes_corr.xlsx')

#EDGELIST BIOGRID: GET BIOGRID INTERACTORS FOR ALL GENES IN GENE LIST
# Handle empty edgelist when saving
if edgelist_biogrid.empty:
    print("Warning: Empty BioGRID edgelist. Creating a minimal file with column headers only.")
    pd.DataFrame(columns=['Gene', 'Gene1', 'tuples', 'bg']).to_excel('genes_bg.xlsx')
else:
    edgelist_biogrid.to_excel('genes_bg.xlsx')

#DEFAULT IS TO COMBINE:
#COMBINE BIOGRID AND CORR INTO ONE NETWORK: OVERLAY BIOGRID INTERACTIONS (FOR GENES IN genes_of_interest ONLY) ONTO COESSENTIALITY 
# Handle different combinations of empty DataFrames
if corr.empty and edgelist_biogrid.empty:
    print("Warning: Both correlation matrix and BioGRID edgelist are empty. Creating a minimal merged file.")
    pd.DataFrame(columns=['Gene', 'Gene1', 'corrscore', 'bg']).to_excel('genes_corr_bg_merge.xlsx')
elif corr.empty:
    print("Warning: Correlation matrix is empty, exporting BioGRID data only.")
    # Add empty corrscore column to BioGRID data
    edgelist_with_corr = edgelist_biogrid.copy()
    edgelist_with_corr['corrscore'] = np.nan
    edgelist_with_corr[['Gene', 'Gene1', 'corrscore', 'bg']].to_excel('genes_corr_bg_merge.xlsx')
elif edgelist_biogrid.empty:
    print("Warning: BioGRID edgelist is empty, exporting correlation data only.")
    # Add empty bg column to correlation data
    corr_with_bg = corr.copy()
    corr_with_bg['bg'] = np.nan
    corr_with_bg[['Gene', 'Gene1', 'corrscore', 'bg']].to_excel('genes_corr_bg_merge.xlsx')
else:
    # Both DataFrames have data, proceed with merge
    try:
        corrwithbgforcorr = pd.merge(corr, edgelist_biogrid, how='left', left_on=['Gene','Gene1'], right_on=['Gene','Gene1'])
        corrwithbgforcorr[['Gene', 'Gene1', 'corrscore', 'bg']].to_excel('genes_corr_bg_merge.xlsx')
    except KeyError as e:
        print(f"Warning: Merge operation failed with error: {str(e)}")
        print("Creating separate output files instead.")
        # If merge fails, create separate output with proper structure
        combined_df = pd.DataFrame(columns=['Gene', 'Gene1', 'corrscore', 'bg'])
        combined_df.to_excel('genes_corr_bg_merge.xlsx')
