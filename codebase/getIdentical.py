import pandas as pd
import numpy as np

def getIdenticalRule(rule,radius):
    df = pd.read_csv('rule_' + str(radius) + '.csv',\
            index_col=0).fillna(0)

    for col in df.columns.values:
        if df[col].equals(df[rule]) or df[col].equals(-df[rule]):
            print col

def getIdenticalMolecularSignature(met,radius):
    df = pd.read_csv('molecular_signature_' + str(radius) + '.csv',\
        index_col=0).fillna(0)

    for col in df.columns.values:
        if df[col].equals(df[met]):
            print col   

def duplicate_columns(df, return_dataframe = False, verbose = False):
    '''
        a function to detect and possibly remove duplicated columns for a pandas dataframe
    '''
    from pandas.core.common import array_equivalent
    # group columns by dtypes, only the columns of the same dtypes can be duplicate of each other
    groups = df.columns.to_series().groupby(df.dtypes).groups
    duplicated_columns = []
 
    for dtype, col_names in groups.items():
        column_values = df[col_names]
        num_columns = len(col_names)
 
        # find duplicated columns by checking pairs of columns, store first column name if duplicate exist 
        for i in range(num_columns):
            column_i = column_values.iloc[:,i].values
            for j in range(i + 1, num_columns):
                column_j = column_values.iloc[:,j].values
                if array_equivalent(column_i, column_j):
                    if verbose: 
                        print("column {} is a duplicate of column {}".format(col_names[i], col_names[j]))
                    duplicated_columns.append(col_names[i])
                    break
    if not return_dataframe:
        # return the column names of those duplicated exists
        return duplicated_columns
    else:
        # return a dataframe with duplicated columns dropped 
        return df.drop(labels = duplicated_columns, axis = 1)

if __name__ == '__main__':
    radius = 1
    rule = 'R09078'
    # rule = 'R02963'
    # rule = 'R01173'
    getIdenticalRule(rule,radius)
