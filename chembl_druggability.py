import requests
import json
import pandas as pd
import time
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry

# Target UniProt IDs and corresponding gene names
target_uniprot_ids = ["Q5S007", "P37840", "P04062", "P38137", "Q92769", "P43844", "Q13221", "O00257", "Q8N9N0", "Q9HBF4"] 
gene_names = ["LRRK2", "SNCA", "GBA", "PARK7", "HDAC4", "EPHA4", "IFT88", "TMEM67", "BBS4", "CEP290"]

# Create a session with retry logic
def create_session():
    session = requests.Session()
    retry = Retry(
        total=3,
        backoff_factor=1,
        status_forcelist=[429, 500, 502, 503, 504],
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    return session

session = create_session()

print(f"Querying ChEMBL API for {len(target_uniprot_ids)} targets...")
print("Using extended timeouts and retry logic...\n")

results_list = []

for idx, (uniprot_id, gene_name) in enumerate(zip(target_uniprot_ids, gene_names), 1):
    print(f"[{idx}/{len(target_uniprot_ids)}] Processing {gene_name} ({uniprot_id})...")
    
    try:
        # Step 1: Get ChEMBL target ID from UniProt accession
        target_url = f"https://www.ebi.ac.uk/chembl/api/data/target.json?target_components__accession={uniprot_id}"
        
        response = session.get(target_url, timeout=30)
        response.raise_for_status()
        target_data = response.json()
        
        # Check if target was found
        if target_data.get('targets') and len(target_data['targets']) > 0:
            chembl_target_id = target_data['targets'][0]['target_chembl_id']
            target_type = target_data['targets'][0].get('target_type', 'N/A')
            print(f"  ✓ ChEMBL ID: {chembl_target_id} (Type: {target_type})")
            
            # Step 2: Get activity summary instead of all activities
            # This is much faster - gets aggregated data
            time.sleep(1)  # Be nice to the API
            
            # Try to get activities with a smaller limit first
            activity_url = f"https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id={chembl_target_id}&limit=100"
            
            activity_response = session.get(activity_url, timeout=30)
            activity_response.raise_for_status()
            activity_data = activity_response.json()
            
            activities = activity_data.get('activities', [])
            total_activities = activity_data.get('page_meta', {}).get('total_count', len(activities))
            
            # Count potent compounds from the sample
            potent_count = 0
            for activity in activities:
                std_type = activity.get('standard_type')
                std_value = activity.get('standard_value')
                std_units = activity.get('standard_units')
                
                if std_type in ['IC50', 'Ki', 'Kd'] and std_value and std_units == 'nM':
                    try:
                        value = float(std_value)
                        if value < 1000:  # Less than 1 µM
                            potent_count += 1
                    except (ValueError, TypeError):
                        continue
            
            # Estimate total potent compounds if we only sampled
            if total_activities > 100:
                estimated_potent = int((potent_count / len(activities)) * total_activities) if activities else 0
                potent_display = f"~{estimated_potent}"
                print(f"  ✓ Estimated {estimated_potent} potent compounds from {total_activities} total activities (sampled 100)")
            else:
                potent_display = potent_count
                print(f"  ✓ Found {potent_count} potent compounds from {total_activities} activities")
            
            results_list.append({
                "Gene": gene_name,
                "UniProt_ID": uniprot_id,
                "ChEMBL_ID": chembl_target_id,
                "Target_Type": target_type,
                "Potent_Compounds": potent_display,
                "Total_Activities": total_activities
            })
            
        else:
            print(f"  ✗ No ChEMBL target found")
            results_list.append({
                "Gene": gene_name,
                "UniProt_ID": uniprot_id,
                "ChEMBL_ID": "Not Found",
                "Target_Type": "N/A",
                "Potent_Compounds": 0,
                "Total_Activities": 0
            })
            
    except requests.exceptions.Timeout:
        print(f"  ✗ Timeout (>30s)")
        results_list.append({
            "Gene": gene_name,
            "UniProt_ID": uniprot_id,
            "ChEMBL_ID": "Timeout",
            "Target_Type": "N/A",
            "Potent_Compounds": 0,
            "Total_Activities": 0
        })
    except requests.exceptions.RequestException as e:
        print(f"  ✗ Error: {str(e)[:100]}")
        results_list.append({
            "Gene": gene_name,
            "UniProt_ID": uniprot_id,
            "ChEMBL_ID": "Error",
            "Target_Type": "N/A",
            "Potent_Compounds": 0,
            "Total_Activities": 0
        })
    
    print()  # Blank line for readability

# Create DataFrame and display results
df = pd.DataFrame(results_list)

print("\n" + "="*90)
print("ChEMBL DRUGGABILITY ASSESSMENT")
print("="*90)
print(f"\nPotent compounds: IC50/Ki/Kd < 1 µM (1000 nM)")
print(f"Note: '~' indicates estimate based on sample of 100 activities\n")
print(df.to_string(index=False))

# Summary statistics
print("\n" + "="*90)
print("SUMMARY")
print("="*90)

# Convert potent compounds to numeric for summary (remove ~ prefix)
df['Potent_Numeric'] = df['Potent_Compounds'].apply(
    lambda x: int(str(x).replace('~', '')) if x != 0 else 0
)

druggable = df[df['Potent_Numeric'] > 0]
print(f"Targets found in ChEMBL: {len(df[df['ChEMBL_ID'].str.startswith('CHEMBL', na=False)])}/{len(df)}")
print(f"Targets with potent compounds: {len(druggable)}/{len(df)}")
print(f"Total potent compounds found: {df['Potent_Numeric'].sum()}")
if len(df) > 0:
    print(f"Average potent compounds per target: {df['Potent_Numeric'].mean():.1f}")

# Show most druggable targets
if len(druggable) > 0:
    print("\nMost druggable targets (by compound count):")
    top_targets = druggable.nlargest(5, 'Potent_Numeric')[['Gene', 'ChEMBL_ID', 'Potent_Compounds']]
    print(top_targets.to_string(index=False))
