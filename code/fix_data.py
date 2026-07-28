import urllib.request
import os

studies = {
    'Human': {
        'GSE249812': 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE249nnn/GSE249812/suppl/GSE249812_genes.fpkm.anno.txt.gz',
        'GSE103182': 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103182/suppl/GSE103182_FPKM_Matrix.txt.gz',
        'GSE308783': 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE308nnn/GSE308783/suppl/GSE308783_Processed_data_files.txt.gz'
    },
    'Mouse': {
        'GSE267256': 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE267nnn/GSE267256/suppl/GSE267256_MusMusculus_Counts.txt.gz',
        'GSE285626': 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE285nnn/GSE285626/suppl/GSE285626_FPKM_allsamples.txt.gz',
        'GSE304427': 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE304nnn/GSE304427/suppl/GSE304427_fpkm.csv.gz'
    }
}

os.makedirs("data/raw", exist_ok=True)
for sp, d in studies.items():
    for gse, url in d.items():
        filename = url.split('/')[-1]
        print(f"Downloading {gse}...")
        try:
            urllib.request.urlretrieve(url, f"data/raw/{gse}_{filename}")
            print(f"Downloaded {gse}")
        except Exception as e:
            print(f"Failed {gse}: {e}")
import urllib.request

# Ensure we have a third human dataset in the automated scripts
try:
    print("Downloading GSE154294...")
    urllib.request.urlretrieve('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154294/suppl/GSE154294_combined.count.tab.gz', 'data/raw/GSE154294.tab.gz')
except Exception as e:
    print(f"Failed GSE154294: {e}")
