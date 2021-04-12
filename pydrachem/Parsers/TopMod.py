import pandas as pd

def pmout_parser(filename):
    """
    Reads and processes the output from TopMod's pop90 program output and returns a pandas dataframe
    containing the basins, electronic populations, and multipolar breakdowns for each basin.
    Arguments:
        filename : str
            path to the file containing output from pop90
    Returns:
        cleaned_df : Pandas DataFrame object
            contains a cleaned and sorted pandas dataframe with the data listed above.
    """
    pmout = open(filename,"r")

    line = pmout.readline()
    while "basin " not in line: ## Get to start of results:
        line = pmout.readline()
        while line == "\n":
            line = pmout.readline()

    while line == "\n":
        line = pmout.readline()
    df = pd.DataFrame(columns=line.replace("std. dev.","stdev").split())
    line = pmout.readline()
    while line == "\n":
        line = pmout.readline()
    while "sum of populations" not in line:
        row_data = line.split()
        row_data[1] = str(row_data[0]+" "+row_data[1])
        row_data = row_data[1:]
        column_names=list(df)
        d = {k:[v] for k,v in zip(column_names,row_data)}
        basins = pd.DataFrame(d)
        df = df.append(basins,ignore_index=True)
        line = pmout.readline()
        while line == "\n":
            line = pmout.readline()

    sum_of_populations = float(line.split()[-1])
    print(sum_of_populations)
    line = pmout.readline()
    while line == "\n":
        line = pmout.readline() ## Gets to next section of results

    df2=pd.DataFrame(columns=line.split()[:-1])
    line = pmout.readline()
    while line == "\n":
        line = pmout.readline()

    while "First Moments" not in line:
        row_data = line.split()
        row_data[1] = str(row_data[0]+" "+row_data[1])
        row_data = row_data[1:]
        column_names=list(df2)
        d = {k:[v] for k,v in zip(column_names,row_data)}
        basins = pd.DataFrame(d)
        df2 = df2.append(basins,ignore_index=True)
        line = pmout.readline()
        while line == "\n":
            line = pmout.readline()

    df3 = pd.merge(df,df2,how="outer",on="basin")
    df = df3
    line = pmout.readline()
    while line == "\n":
        line = pmout.readline()
    column_names = line.split()
    column_names[-1] = "Norm.Dipole"
    df2 = pd.DataFrame(columns=column_names)
    line = pmout.readline()
    while "Second Moments" not in line:
        row_data = line.split()
        row_data[1] = str(row_data[0]+" "+row_data[1])
        row_data = row_data[1:]
        column_names=list(df2)
        d = {k:[v] for k,v in zip(column_names,row_data)}
        basins = pd.DataFrame(d)
        df2 = df2.append(basins,ignore_index=True)
        line = pmout.readline()
        while line == "\n":
            line = pmout.readline()

    df3 = pd.merge(df,df2,how="outer",on="basin")
    df = df3

    line = pmout.readline()
    while line == "\n":
        line = pmout.readline()
    column_names = line.split()
    column_names[-1] = "Norm.Quadrupole"
    df2 = pd.DataFrame(columns=column_names)
    line = pmout.readline()
    while line == "\n":
        line = pmout.readline()
    while "Dipolar Contributions" not in line:
        row_data = line.split()
        if row_data != []:
            row_data[1] = str(row_data[0]+" "+row_data[1])
            row_data = row_data[1:]
            column_names=list(df2)
            d = {k:[v] for k,v in zip(column_names,row_data)}
            basins = pd.DataFrame(d)
            df2 = df2.append(basins,ignore_index=True)
            line = pmout.readline()
        while line == "\n" or line.split()==[]:
            line = pmout.readline()

    df3 = pd.merge(df,df2,how="outer",on="basin")
    df = df3

    pmout.close() ## End of the job

    cleaned_df=df.drop(labels=["vol.","pab","paa","paa.","pb","pbb_x","pbb_y","sigma2_x","sigma2_y","stdev","pa","paa"],axis=1)
    return cleaned_df

def pmout_comparison(reference, subject):
    if type(reference) == str:
        reference = pmout_parser(reference)
    if type(subject) == str:
        subject = pmout_parser(subject)
    test = reference.drop("basin",axis=1).astype(float) - subject.drop("basin",axis=1).astype(float)
    test.insert (0, "basin", list(subject["basin"]))
    return test

def multiframes_to_excel(dfs,sheets,filename):
    if type(dfs) != list:
        print("Type Error:  'dfs' must be a list of pandas dataframes")
        return
    if type(sheets) != list:
        print("Type Error:  'sheets' must be a list of strings")
        return
    if len(dfs) != len(sheets):
        print("List error:  'dfs' ("+str(len(dfs))+") and 'sheets' ("+str(len(sheets))+")must be of equal lengths.")
        return
    if type(filename) != str:
        print("Type Error:  'filename' must be a string for an output filename.")
        return
    excel = pd.ExcelWriter(filename, engine='xlsxwriter') # pylint: disable=abstract-class-instantiated
    for i in range(len(dfs)):
        dfs[i].to_excel(excel,sheet_name=sheets[i],index=False)
    excel.save()


    