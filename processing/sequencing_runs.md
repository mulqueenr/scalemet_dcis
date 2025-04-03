Overview:

located in /data/rmulqueen/projects/scalebio_dcis/seq
Preliminary 1:
(Plates 1 and 2 of alpha kit)
MCF10A, MCF7, 231
HBCA BCMH13CA16R
HBCA BCMH13CA83L
https://mdandersonorg-my.sharepoint.com/personal/rmulqueen_mdanderson_org/_layouts/OneNote.aspx?id=%2Fpersonal%2Frmulqueen_mdanderson_org%2FDocuments%2FmetACT&wd=target%28scalebio%20sciMETv2.one%7CD3F046A2-B151-0443-938E-82A415D420EB%2F240109%20Scale%20Met%20Alpha%20Test%20Kit%7C2C45CEED-7824-7C4B-8CFC-697EE8D6A947%2F%29
- 240111_VH00219_563_AAFFFCCM5 #initial test on cell lines + HBCA
- PM2517 #240202_A01819_0450_AHNTKYDMXY high coverage run on cell lines + HBCA

Preliminary 2:
(Plate 3 of alpha kit)
BCMDCIS41T
BCMDCIS66T
BCMDCIS88T
https://mdandersonorg-my.sharepoint.com/personal/rmulqueen_mdanderson_org/_layouts/OneNote.aspx?id=%2Fpersonal%2Frmulqueen_mdanderson_org%2FDocuments%2FmetACT&wd=target%28scalebio%20sciMETv2.one%7CD3F046A2-B151-0443-938E-82A415D420EB%2F240422%20Scale%20Met%20Alpha%20Test%20Kit%202%7CC8153809-4616-E444-A791-2800DD23B717%2F%29

- 240523_VH00219_594_AAFLYGNM5

Preliminary 3:
- 241004_A01819_0637_BHY5MJDMXY #second preliminary run (HBCA+DCIS)

Batch 1 4 initial plates:
Collating from other servers:

```bash
#from sysbio sftp to geo and get
cd /data/rmulqueen/projects/scalebio_dcis/seq
sftp mulqueen@qcprpgeo.mdanderson.edu
#get -R /volumes/seq/flowcells/MDA/nextseq2000/2024/240111_VH00219_563_AAFFFCCM5 #initial test on cell lines + HBCA
#get -R /volumes/seq/flowcells/MDA/NovaSeq/2024/PM2517 #240202_A01819_0450_AHNTKYDMXY high coverage run on cell lines + HBCA
#get -R /volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/241004_A01819_0637_BHY5MJDMXY #novaseq of additional samples (second prelim run)
```

```bash
#from seadragon to sysbio
bsub -Is -W 4:00 -q transfer -n 10 -M 10 -R rusage[mem=100] /bin/bash #small interactive node for bcl-convert 
cd /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240526_RMMM_scalebio_dcis
sftp mulqueen@qcprpsysbio.mdanderson.edu
#put -R 240523_VH00219_594_AAFLYGNM5 #high coverage first prelim run HBCA + DCIS 
```

Just currently missing batch 1 nextseq run
```bash
cd /data/rmulqueen/projects/scalebio_dcis/seq
sftp mulqueen@qcprpneo.mdanderson.edu

```
