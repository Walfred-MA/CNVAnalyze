import os

output = "./summary.txt"
folder = "results/"
allfiles = [folder+x for x in os.listdir(folder) if x.endswith("txt")]

outtext = []
for thefile in allfiles:
        with open(thefile , mode = 'r') as f:
                text = f.read().splitlines()

        if len(text) < 4:
                continue


        gene = thefile.split("/")[-1].split("%")[0]

        name = "_".join(thefile.split("/")[-1].split("%")[1].split("_")[:-1])

        pvalue = text[-1].split(" ")[-1]

        #print(text[0].strip())
        outtext.append(" ".join([gene, name, text[0].strip() , text[1].strip(), pvalue]))

with open(output,mode = 'w') as f:

        f.write("\n".join(outtext)+"\n")
