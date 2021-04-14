

secreted = set([])


with open("secreted.txt", "r") as fh:
    for line in fh:
       line = line.rstrip()
       line = line + "-T1"
       secreted.add(line)


f_out =  open("Orthogroups_secreted.txt", "w")  
with open("Orthogroups.txt", "r") as fh:
    for line in fh:
       elements = line.split()
       for gene in elements:
          if gene.startswith("GPALN"):
             if not gene in secreted:
                line = line.replace(gene, "")
       f_out.write(line)
       
f_out.close() 
