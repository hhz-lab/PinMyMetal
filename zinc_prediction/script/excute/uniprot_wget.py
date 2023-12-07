import statistics
import wget
import os

def get_uniprot_file(uniprot):
    url = "http://www.uniprot.org/uniprot/"+str(uniprot)+".txt"
    inputfile = wget.download(url)

    
    r_list=[]
    c_list=[]
    try:
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
                try:
                    resolution = float(s[3][1:5])
                except:
                    pass
                c_aa = s[4].split('=')[-1].split('.')[0].split('-')
                complete =int(c_aa[1])-int(c_aa[0])
                score = (1/(resolution/r_m))+(complete/c_m)
                pdb_dict[pdbid]=score

        os.remove(inputfile)
        if pdb_dict:
            d_order=sorted(pdb_dict.items(),key=lambda x:x[1],reverse=False)
            best_pdbid = d_order[-1][0]
            return "https://files.rcsb.org/download/"+best_pdbid+".pdb",uniprot+'-'+best_pdbid
        else:
            return "https://alphafold.ebi.ac.uk/files/AF-"+uniprot+"-F1-model_v4.pdb",uniprot
    except:
        return "https://alphafold.ebi.ac.uk/files/AF-"+uniprot+"-F1-model_v4.pdb",uniprot
