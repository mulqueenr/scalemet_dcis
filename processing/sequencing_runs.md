# Sequenced runs
All bcl data is located in /data/rmulqueen/projects/scalebio_dcis/seq

## Preliminary 1:
(Plates 1 and 2 of alpha kit)
MCF10A, MCF7, 231
HBCA BCMH13CA16R
HBCA BCMH13CA83L
https://mdandersonorg-my.sharepoint.com/personal/rmulqueen_mdanderson_org/_layouts/OneNote.aspx?id=%2Fpersonal%2Frmulqueen_mdanderson_org%2FDocuments%2FmetACT&wd=target%28scalebio%20sciMETv2.one%7CD3F046A2-B151-0443-938E-82A415D420EB%2F240109%20Scale%20Met%20Alpha%20Test%20Kit%7C2C45CEED-7824-7C4B-8CFC-697EE8D6A947%2F%29
- 240111_VH00219_563_AAFFFCCM5 #initial test on cell lines + HBCA
- PM2517 #240202_A01819_0450_AHNTKYDMXY high coverage run on cell lines + HBCA

## Preliminary 2:
(Plate 3 of alpha kit)
BCMDCIS41T
BCMDCIS66T
BCMDCIS88T
https://mdandersonorg-my.sharepoint.com/personal/rmulqueen_mdanderson_org/_layouts/OneNote.aspx?id=%2Fpersonal%2Frmulqueen_mdanderson_org%2FDocuments%2FmetACT&wd=target%28scalebio%20sciMETv2.one%7CD3F046A2-B151-0443-938E-82A415D420EB%2F240422%20Scale%20Met%20Alpha%20Test%20Kit%202%7CC8153809-4616-E444-A791-2800DD23B717%2F%29
- 240523_VH00219_594_AAFLYGNM5

## Preliminary 3:
(Small kit 1, production kit)
- HBCA-19R 
- HBCA-17R 
- DCIS-92T 
- DCIS-66T 
- DCIS-79T 
- IDC-79T 
[BCMDCIS79_IDC_T 8/22/2023 1M 86%]
https://mdandersonorg-my.sharepoint.com/personal/rmulqueen_mdanderson_org/_layouts/OneNote.aspx?id=%2Fpersonal%2Frmulqueen_mdanderson_org%2FDocuments%2FmetACT&wd=target%28scalebio%20sciMETv2.one%7CD3F046A2-B151-0443-938E-82A415D420EB%2F240910%20ScaleBio%20DCIS%20Samples%7C30534461-040E-C54F-BB40-7D53F8115495%2F%29
- 241004_A01819_0637_BHY5MJDMXY 

## Batch 1 
(Small kit 2, production kit)
BCMDCIS41T
BCMHBCA03R
BCMDCIS07T
BCMHBCA10L_3h
BCMHBCA15R_3h
BCMDCIS22T
BCMHBCA26L_24hTis_4h
BCMDCIS28T
BCMDCIS32T
ECIS36T
BCMHBCA38L_3h
BCMDCIS52T
BCMDCIS66T
BCMDCIS74T
BCMDCIS80T
BCMDCIS94T_24hTis
BCMDCIS99T

### 4 initial plates:
https://mdandersonorg-my.sharepoint.com/personal/rmulqueen_mdanderson_org/_layouts/OneNote.aspx?id=%2Fpersonal%2Frmulqueen_mdanderson_org%2FDocuments%2FmetACT&wd=target%28scalebio_dcis_processing.one%7CA0170AC0-CBA6-B848-B58D-A9AF996E46F1%2F250311%20Batch%201%20Plates%201%2C%206%20Processing%20Notes%20Scalebio%20Kit%7C3BFEC095-CE37-DC4C-A08C-E6F1F4B4067B%2F%29

https://mdandersonorg-my.sharepoint.com/personal/rmulqueen_mdanderson_org/_layouts/OneNote.aspx?id=%2Fpersonal%2Frmulqueen_mdanderson_org%2FDocuments%2FmetACT&wd=target%28scalebio_dcis_processing.one%7CA0170AC0-CBA6-B848-B58D-A9AF996E46F1%2F250324%20Batch%201%20Plates%202%2C7%20Processing%20Notes%20Homebrew%7C59F268C2-A396-ED4B-B87E-AB0927640532%2F%29


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
