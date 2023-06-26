import gc
from transformers import EsmTokenizer, EsmModel
import torch
import pandas as pd
import requests as r
from Bio import SeqIO
from io import StringIO
import os
from numba import cuda
import csv


gc.collect()
torch.cuda.empty_cache()
f = open("esm2_output.csv", "a")
writer = csv.writer(f)

#Load Abbasov paper data
df = pd.read_excel("../lysine_data.xlsx", sheet_name=2)

def get_full_seq(cID):
    baseUrl="http://www.uniprot.org/uniprot/"
    currentUrl=baseUrl+cID+".fasta"
    response = r.post(currentUrl)
    cData=''.join(response.text)

    Seq=StringIO(cData)
    pSeq=list(SeqIO.parse(Seq,'fasta'))
    if not pSeq:
        return("")
    else:
        return(str(pSeq[0].seq))

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
#device = torch.device("cpu")
model_name = "facebook/esm2_t33_650M_UR50D"
tokenizer = EsmTokenizer.from_pretrained(model_name)
model = EsmModel.from_pretrained(model_name).to(device)

model.eval().to(device)

ids = []
lys_esm2 = []
lys_pos = []
for i in range(0, df.shape[0]):
    #print(i)
    protID = df.loc[i, 'uniprot']
    seq = get_full_seq(protID)
    idx_lys = df.loc[i, 'probe P1-labeled site']
    if seq and idx_lys < len(seq):
        with torch.no_grad():
            inputs = tokenizer(seq, return_tensors="pt").to(device)
            outputs = model(**inputs)
        last_hidden_states = outputs.last_hidden_state
        residue_rep = last_hidden_states.cpu().detach().numpy()
        lys_rep = residue_rep[0][idx_lys]
        #ids.append(protID)
        #lys_esm2.append(lys_rep)
        #lys_pos.append(idx_lys)
        data = [protID, idx_lys]
        for j in lys_rep:
            data.append(j)
        writer.writerow(data)
        del inputs, outputs, last_hidden_states
        torch.cuda.empty_cache()
        gc.collect()

#data = {'PDB code': ids,
        #'Lysine posiion': lys_pos,
        #'ESM2 Representation of lysine': lys_esm2}

#df = pd.DataFrame(data)
#df.to_csv('esm2_info2.csv', sep='\t', index=False, encoding='utf-8')
f.close()

#INITIAL TRIAL
#inputs = tokenizer("K.YGIEPTMVVQGVKMLYVPVMPGHAK.R", return_tensors="pt").to(device)
#outputs = model(**inputs)

#last_hidden_states = outputs.last_hidden_state

#31 is length of inputs array
#token_len = inputs['input_ids'].numpy().shape[1]

#print(last_hidden_states)
#print(inputs)


