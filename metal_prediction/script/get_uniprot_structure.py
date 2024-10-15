import statistics
#wget http://www.uniprot.org/uniprot/P27487.txt
inputfile = "p0DTC1.txt"

r_list=[]
c_list=[]

for l in open(inputfile,"r").readlines():
    s = l.split(';')
    if s[0][0:2] == "DR" and s[0][5:] =="PDB":
        try:
            resolution = float(s[3][1:5])
        except:
            pass
        c_aa = s[4].split('=')[-1].split('.')[0].split('-')
        complete =int(c_aa[1])-int(c_aa[0])
        r_list.append(resolution)
        c_list.append(complete)

r_m=statistics.median(r_list)
c_m=statistics.median(c_list)
pdb_dict={}
for l in open(inputfile,"r").readlines():
    s = l.split(';')
    if s[0][0:2] == "DR" and s[0][5:] =="PDB":
        pdbid = s[1][1:]
        print(pdbid)
        try:
            resolution = float(s[3][1:5])
        except:
            pass
        c_aa = s[4].split('=')[-1].split('.')[0].split('-')
        complete =int(c_aa[1])-int(c_aa[0])
        score = (1/(resolution/r_m))+(complete/c_m)
        pdb_dict[pdbid]=score
if pdb_dict:
    d_order=sorted(pdb_dict.items(),key=lambda x:x[1],reverse=False)
    best_pdbid = d_order[-1][0]
    print(best_pdbid)
    #wget https://files.rcsb.org/download/best_pdbid.pdb
#else:
    #https://alphafold.ebi.ac.uk/files/AF-P15121-F1-model_v4.pdb

