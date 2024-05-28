import os
import pandas as pd
import numpy as np
import multiprocessing as mp
from rdkit import Chem
from rdkit.Chem import Descriptors
from tqdm import tqdm
import warnings
from propy import PyPro
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import time
import sqlite3

# This is surely not causing a headache later on...
warnings.filterwarnings("ignore")

# The following statement was the work of 4h of debugging....
file_lock = mp.Lock()

# Function to extract all features and return a DataFrame
def extract_all_features(peptide):
    try:
        pro = PyPro.GetProDes(peptide)
        
        # Amino acid composition
        aa_comp = pro.GetAAComp()
        aa_comp_names = [f"AAComp_{aa}" for aa in aa_comp.keys()]
        aa_comp_values = list(aa_comp.values())
        
        # Dipeptide composition
        dp_comp = pro.GetDPComp()
        dp_comp_names = [f"DPComp_{dp}" for dp in dp_comp.keys()]
        dp_comp_values = list(dp_comp.values())
        
        # Tripeptide composition
        tp_comp = pro.GetTPComp()
        tp_comp_names = [f"TPComp_{tp}" for tp in tp_comp.keys()]
        tp_comp_values = list(tp_comp.values())
        
        # Moreau-Broto autocorrelation descriptors
        moreau_broto = pro.GetMoreauBrotoAuto()
        moreau_broto_names = [f"MoreauBroto_{i}" for i in range(len(moreau_broto))]
        moreau_broto_values = list(moreau_broto.values())
        
        # Moran autocorrelation descriptors
        moran = pro.GetMoranAuto()
        moran_names = [f"Moran_{i}" for i in range(len(moran))]
        moran_values = list(moran.values())
        
        # Geary autocorrelation descriptors
        geary = pro.GetGearyAuto()
        geary_names = [f"Geary_{i}" for i in range(len(geary))]
        geary_values = list(geary.values())
        
        # Quasi-sequence order descriptors
        qso = pro.GetQSO()
        qso_names = [f"QSO_{i}" for i in range(len(qso))]
        qso_values = list(qso.values())
        
        # Calculate additional physicochemical properties using Bio.SeqUtils.ProtParam
        analysed_seq = ProteinAnalysis(peptide)
        physchem_names = ['Molecular_Weight', 'Isoelectric_Point', 'Instability_Index', 'GRAVY']
        physchem_values = [
            analysed_seq.molecular_weight(),  # Molecular weight
            analysed_seq.isoelectric_point(),  # Isoelectric point (pI)
            analysed_seq.instability_index(),  # Instability index
            analysed_seq.gravy(),  # Hydrophobicity (GRAVY)
        ]
        
        # Combine all feature names and values
        feature_names = (aa_comp_names + dp_comp_names + tp_comp_names + moreau_broto_names +
                         moran_names + geary_names + qso_names + physchem_names)
        feature_values = (aa_comp_values + dp_comp_values + tp_comp_values + moreau_broto_values +
                          moran_values + geary_values + qso_values + physchem_values)
        
        # Create a DataFrame with the feature names and values
        features_df = pd.DataFrame([feature_values], columns=feature_names)
        
        return features_df
    
    except ZeroDivisionError:
        # Handle division by zero gracefully
        print(f"Error processing sequence {peptide}: Division by zero")
        return pd.DataFrame()
    except AttributeError:
        # Handle attribute error gracefully
        print(f"Error processing sequence {peptide}: Empty sequence")
        return pd.DataFrame()
    except Exception as e:
        # Handle other exceptions gracefully
        print(f"Error processing sequence {peptide}: {e}")
        return pd.DataFrame()

def calculate_descriptors(row):
    mol = Chem.MolFromSequence(row['seq'])
    descriptor_names = [desc_name for desc_name, _ in Descriptors._descList]
    descriptors = {}
    for desc_name in descriptor_names:
        try:
            descriptor_value = getattr(Descriptors, desc_name)(mol)
            descriptors[desc_name] = descriptor_value
        except Exception as e:
            descriptors[desc_name] = None
            print(f"Error calculating descriptor {desc_name}: {e}")

    return pd.concat([row.to_frame().T.reset_index(drop=True), pd.DataFrame([descriptors]).reset_index(drop=True)], axis=1)

def save_to_hdf(df, file_path):
    df = df.loc[:, ~df.columns.duplicated()]
    with file_lock:
        df.to_hdf(file_path, key='data', mode='a', format='table', append=True, complib='blosc', complevel=9)

def calc_and_save(packet, file_path, save_interval):
    process_id = mp.current_process().pid
    q.put(f"ID: {process_id:<5}, Packetsize: {len(packet['seq']):<4}, Calculation started")
    result = pd.concat([calculate_descriptors(packet.iloc[0]), extract_all_features(packet.iloc[0]['seq'])], axis=1)
    for i in range(1, len(packet['seq'])):
        q.put(f"ID: {process_id:<5}, progress: {i / len(packet['seq']) *  100:2.2f} % , calculating descriptors for seq: {packet.iloc[i]['seq']}")
        new_result = pd.concat([calculate_descriptors(packet.iloc[i]), extract_all_features(packet.iloc[i]['seq'])], axis=1)
        result = pd.concat([result, new_result], ignore_index=True)
        
        # Save at intervals
        if i % save_interval == 0 and i != 0:
            q.put(f"ID: {process_id:<5}, saving ...")
            save_to_hdf(result, file_path)
            result = result[0:0]
            q.put(f"ID: {process_id:<5}, continuing")
    
    # Save any remaining results
    if not result.empty:
        q.put(f"ID: {process_id:<5}, saving ...")
        save_to_hdf(result, file_path)
        q.put(f"ID: {process_id:<5}, finished")

def queue_printer(queue, stop_event, filename, total_messages):
    with open(filename, 'w') as f, tqdm(total=total_messages, desc="Processing") as pbar:
        while not stop_event.is_set() or not queue.empty():
            try:
                item = queue.get(timeout=0.1)
                f.write(f"{item}\n")
                f.flush()
                if "saving" not in item and "continuing" not in item and not "started" in item:
                    pbar.update(1)
            except:
                continue

def init(queue):
    global q
    q = queue

if __name__ == "__main__":
    # Get number of cores
    num_cores = os.cpu_count()
    print(f"CPU cores: {num_cores}")

    # DB import
    con = sqlite3.connect('unified.db')
    df = pd.read_sql_query("SELECT * FROM prod", con)
    con.close()
    print(f"Found {len(df['seq'])} sequences.") 

    # For interprocess communication
    msg_queue = mp.Queue()

    # Calculate the total number of messages (excluding "saving" and "continuing")
    total_messages = len(df['seq'])

    # Start the queue printer process
    stop_event = mp.Event()
    log_filename = "process_log.txt"
    printer_process = mp.Process(target=queue_printer, args=(msg_queue, stop_event, log_filename, total_messages))
    printer_process.start()

    # Perform calculation
    print("Starting Calculation")
    output_file = "output.h5"
    
    # Initialize the HDF5 file
    with pd.HDFStore(output_file, mode='w', complib='blosc', complevel=9) as store:
        pass

    df_split = np.array_split(df, num_cores)
    save_intervals = np.random.randint(1, 10, size=num_cores)  # Variable save intervals for each process
    pool = mp.Pool(num_cores, initializer=init, initargs=(msg_queue,))
    pool.starmap(calc_and_save, [(df_split[i], output_file, save_intervals[i]) for i in range(num_cores)])
    pool.close()
    pool.join()

    # Signal the process to stop and wait for it to finish
    stop_event.set()
    printer_process.join()

    print("Calculation complete.")
