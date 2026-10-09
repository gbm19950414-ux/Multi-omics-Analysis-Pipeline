"""Synthetic RNA count merge and lipid PQN, using the original R programs."""
from pathlib import Path
import csv, argparse, subprocess
p=argparse.ArgumentParser();p.add_argument('--rscript',default='Rscript');args=p.parse_args()
ROOT=Path(__file__).resolve().parents[1];out=ROOT/'demo-output';out.mkdir(exist_ok=True)
processed=out/'processed/rna'
for batch,rows in [('batch1',[['gene-a',10,20],['gene-b',30,40]]),('batch2',[['gene-a',50,60],['gene-b',70,80]])]:
 f=processed/batch/'counts/featureCounts_matrix.tsv';f.parent.mkdir(parents=True,exist_ok=True)
 with f.open('w',newline='') as s:
  w=csv.writer(s,delimiter='\t');w.writerow(['GeneID','WT_1','KO_1']);w.writerows(rows)
subprocess.run([args.rscript,str(ROOT/'scripts/rna/07_merge_counts.R'),str(processed),str(out/'merged_counts.tsv')],check=True)
with (out/'merged_counts.tsv').open() as s:rows=list(csv.reader(s,delimiter='\t'))
assert len(rows)==3 and len(rows[0])==5,rows
assert {tuple(row[1:]) for row in rows[1:]}=={('10','20','50','60'),('30','40','70','80')},rows
tables=out/'lipid-workspace/results/lipid/tables';tables.mkdir(parents=True,exist_ok=True)
data=[]
for batch in ['batch1','batch2']:
 for sample,factor,qc in [('QC-1',1,True),('QC-2',1,True),('sample-x',2,False)]:
  for lipid,base in [('lipid-a',10),('lipid-b',20),('lipid-c',30)]:data.append([batch,sample,lipid,base*factor,str(qc).upper(),'TRUE'])
for name in ['lipid_long_dedup_with_qcflag.tsv','lipid_long_for_batch_correction.tsv']:
 with (tables/name).open('w',newline='') as s:
  w=csv.writer(s,delimiter='\t');w.writerow(['batch','sample_id','lipidName','intensity','is_qc','keep_for_correction']);w.writerows(data)
subprocess.run([args.rscript,str(ROOT/'scripts/lipid/04_within_batch_pqn.R')],cwd=out/'lipid-workspace',check=True)
with (tables/'lipid_long_dedup_pqn.tsv').open() as s:rows=list(csv.DictReader(s,delimiter='\t'))
assert len(rows)==18
for row in rows:
 expected={'lipid-a':10,'lipid-b':20,'lipid-c':30}[row['lipidName']]
 assert abs(float(row['intensity_pqn'])-expected)<1e-8,row
print('Synthetic count merge and known-factor PQN checks passed. No real data or scientific result was evaluated.')
