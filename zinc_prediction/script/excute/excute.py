import argparse
import subprocess
import requests
import uniprot_wget
import os

parser = argparse.ArgumentParser()
parser.add_argument("-u", "--uniprot",help='UNIPROT ID')
parser.add_argument("-p", "--pdb",help='PDB ID')
parser.add_argument("-f", "--file",help = 'PDB file name')
args = parser.parse_args()
if(args.pdb):
    url = "https://files.rcsb.org/download/" + args.pdb + ".pdb"
    res = requests.get(url)
    save_path = '/zinc_prediction/script/excute/input_data/%s' % (args.pdb + '.ent')
    with open(save_path, 'wb') as f:
        f.write(res.content)
    with open('test.log', 'a') as log_file:
        subprocess.run(['./../script.sh', args.pdb], stdout=log_file, stderr=subprocess.STDOUT)

if(args.file):
    name = args.file.split('.')[0]
    os.rename('input_data/'+args.file,'input_data/'+name+'.ent')
    with open('test.log', 'a') as log_file:
        subprocess.run(['./../script.sh', name], stdout=log_file, stderr=subprocess.STDOUT)

if(args.uniprot):
    url, re_name = uniprot_wget.get_uniprot_file(args.uniprot)
    res = requests.get(url)
    save_path = '/zinc_prediction/script/excute/input_data/%s'%(re_name+".ent")
    with open(save_path, 'wb') as f:
        f.write(res.content)
    with open('test.log', 'a') as log_file:
        subprocess.run(['./../script.sh', re_name], stdout=log_file, stderr=subprocess.STDOUT)

