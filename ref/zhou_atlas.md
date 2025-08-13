
#Downloading from https://www.biorxiv.org/content/10.1101/2025.03.23.644697v1.full
#https://huggingface.co/datasets/zhoujt1994/HumanCellEpigenomeAtlas_sc_allc/

```bash
 pip install -U "huggingface_hub[cli]"

huggingface-cli download zhoujt1994/HumanCellEpigenomeAtlas_sc_allc \
--repo-type dataset \
--include "B_IOBHL/*tsv.gz" \
--local-dir /data/rmulqueen/projects/scalebio_dcis/ref/zhou_atlas/B_IOBHL \
--max-workers 100 \
--force-download 

huggingface-cli download zhoujt1994/HumanCellEpigenomeAtlas_sc_allc \
--repo-type dataset \
--include "B_JF1NV/*tsv.gz" \
--local-dir /data/rmulqueen/projects/scalebio_dcis/ref/zhou_atlas/B_JF1NV \
--max-workers 100 \
--force-download 


```