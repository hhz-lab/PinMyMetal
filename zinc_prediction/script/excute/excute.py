import argparse
import os
import requests
import uniprot_wget

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
    os.system('./../test_script.sh '+args.pdb+' &>test.log')

if(args.file):
    name = args.file.split('.')[0]
    os.rename('input_data/'+args.file,'input_data/'+name+'.ent')
    os.system('./../test_script.sh '+name+' &>test.log')

if(args.uniprot):
    url,re_name = uniprot_wget.get_uniprot_file(uniprot)
    res = requests.get(url)
    save_path = '/zinc_prediction/script/excute/input_data/%s'%(re_name+".ent")
    with open(save_path, 'wb') as f:
        f.write(res.content)
    os.system('./../test_script.sh '+re_name+' &>test.log')