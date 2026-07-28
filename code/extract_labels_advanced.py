import GEOparse
import pandas as pd
import glob
import os

studies = ['GSE103182', 'GSE154294', 'GSE249812', 'GSE267256', 'GSE285626', 'GSE304427']
os.makedirs('data/metadata', exist_ok=True)

for gse_id in studies:
    print(f"Fetching metadata for {gse_id}...")
    try:
        gse = GEOparse.get_GEO(geo=gse_id, destdir="data/raw")
        labels = {}
        for gsm_name, gsm in gse.gsms.items():
            chars = gsm.metadata.get('characteristics_ch1', [])
            title = gsm.metadata.get('title', [''])[0].lower()
            source = gsm.metadata.get('source_name_ch1', [''])[0].lower()

            is_mi = False
            is_control = False

            all_text = " ".join(chars).lower() + " " + title + " " + source

            if 'mi' in all_text or 'myocardial infarction' in all_text or 'infarct' in all_text or 'disease' in all_text or 'tac' in all_text or 'ischemia' in all_text or 'ami' in all_text or 'remodeling' in all_text:
                is_mi = True
            if 'control' in all_text or 'sham' in all_text or 'healthy' in all_text or 'normal' in all_text or 'wt' in all_text or 'wild type' in all_text or 'baseline' in all_text:
                is_control = True

            if is_mi and not is_control:
                labels[gsm_name] = 'MI'
            elif is_control and not is_mi:
                labels[gsm_name] = 'Control'
            elif is_control and is_mi:
                if 'sham' in all_text or 'control' in all_text or 'healthy' in all_text:
                    labels[gsm_name] = 'Control'
                else:
                    labels[gsm_name] = 'MI'
            else:
                labels[gsm_name] = 'Unknown'

        df = pd.DataFrame(list(labels.items()), columns=['Sample', 'Condition'])
        df.to_csv(f'data/metadata/{gse_id}_labels.csv', index=False)
        print(df['Condition'].value_counts())
    except Exception as e:
        print(f"Failed {gse_id}: {e}")
