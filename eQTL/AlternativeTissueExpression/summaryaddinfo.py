#!/usr/bin/env python3

import scipy.stats as stats

def t_test_greater_than(mean_A, std_A, n_A, mean_B, std_B, n_B):

        if "NA" in [mean_A, std_A, n_A, mean_B, std_B, n_B]:
                return "NA"

        mean_A, std_A, n_A, mean_B, std_B, n_B = map(float,[mean_A, std_A, n_A, mean_B, std_B, n_B])

        if min(n_A,n_B) < 2:
                return "NA"
        # Calculate the t-statistic
        t_statistic = (mean_A - mean_B) / ((std_A**2 / n_A) + (std_B**2 / n_B))**0.5

        # Calculate degrees of freedom using Welch-Satterthwaite equation
        df = ((std_A**2 / n_A) + (std_B**2 / n_B))**2 / (((std_A**2 / n_A)**2 / (n_A - 1)) + ((std_B**2 / n_B)**2 / (n_B - 1)))

        # Calculate the p-value
        p_value = 1 - stats.t.cdf(t_statistic, df)

        return p_value

typefile = "/panfs/qcb-panasas/wangfeim/Build/fourthrun/PangenomeAlleles_typefix.tsv"

inputfile  = "expr_tissuetable.txt"

outputfile = "expr_tissuesummary.txt"

typename = dict()
with open(typefile , mode = 'r') as f:
        for line in f:
                line = line.strip().split()

                typename["_".join(line[0].split("_")[:-1])] = [line[1], line[4].split(":")[0]]


outtext = []
with open(inputfile , mode = 'r') as f:
        for line in f:
                line = line.strip().split()

                Atissue,A1_mean, A1_std, A2_mean, A2_std,A_df, Btissue,B1_mean,B1_std, B2_mean,B2_std, B_df = line[2:]


                p1 = t_test_greater_than(A1_mean,A1_std,A_df,B1_mean,B1_std,B_df)
                p2 = t_test_greater_than(B2_mean,B2_std,B_df,A2_mean,A2_std,A_df)

                outtext.append("\t".join([line[0],typename["_".join(line[1].split("_")[:-1])][0],typename["_".join(line[1].split("_")[:-1])][1]]  + line[1:8]+ [str(p1)] +line[8:] + [str(p2)])+"\n")


with open(outputfile, mode = 'w') as f:
        f.write("".join(outtext))
