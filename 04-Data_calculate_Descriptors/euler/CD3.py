import os
import sqlite3
import pandas as pd
import numpy as np
import multiprocessing as mp
from rdkit import Chem
from rdkit.Chem import Descriptors
from tqdm import tqdm
import warnings

# This is surely not causing a headache later on...
warnings.filterwarnings("ignore")

# The following statement was the work of 4h of debugging....
db_lock = mp.Lock()

def get_sql_con():
    con = sqlite3.connect('unified.db')
    return con

def reset_sql_tables(row, table_name="prod_desc"):
    print("Setting up new sql table")
    con = get_sql_con()
    cur = con.cursor()
    df = calculate_descriptors(row)
    # Create table if it does not exist
    columns = ', '.join([f"{col} REAL" for col in df.columns if col not in ['seq', 'seq_length', 'id', 'valid', 'name', 'source', 'description', 'OX']])
    create_table_query = f"""
    CREATE TABLE IF NOT EXISTS {table_name} (
        id TEXT,
        seq TEXT PRIMARY KEY,
        seq_length INTEGER,
        valid TEXT,
        name TEXT,
        source TEXT,
        description TEXT,
        OX TEXT,
        {columns}
    )
    """
    df.to_sql(table_name, con, if_exists="replace")
    cur.execute(create_table_query)
    
    con.commit()
    con.close()

def save_to_sql(df, table_name="prod_desc"):
    with db_lock:
        con = get_sql_con()
        df.to_sql(table_name, con, if_exists='append', index=False)
        con.commit()
        con.close()

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

def calc_and_save(packet, save_interval=10):
    process_id = mp.current_process().pid
    q.put(f"ID: {process_id:<5}, Packetsize: {len(packet['seq']):<4}, Calculation started")
    result = calculate_descriptors(packet.iloc[0])
    for i in range(1, len(packet['seq'])):
        q.put(f"ID: {process_id:<5}, progress: {i / len(packet['seq']) *  100:2.2f} % , calculating descriptors for seq: {packet.iloc[i]['seq']}")
        result = pd.concat([result, calculate_descriptors(packet.iloc[i])], axis=0)
        # Save interval -> memory management
        if i % save_interval == save_interval - 1:
            q.put(f"ID: {process_id:<5}, saving ...")
            save_to_sql(result)
            result = result[0:0]
            q.put(f"ID: {process_id:<5}, continuing")

def queue_printer(queue, stop_event, filename, total_messages):
    with open(filename, 'w') as f, tqdm(total=total_messages, desc="Processing") as pbar:
        while not stop_event.is_set() or not queue.empty():
            try:
                item = queue.get(timeout=0.1)
                f.write(f"{item}\n")
                f.flush()
                if "saving" not in item and "continuing" not in item:
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
    con = get_sql_con()
    df = pd.read_sql_query("SELECT * FROM prod", con)
    con.close()
    print(f"Found {len(df['seq'])} sequences.") 

    # For interprocess communication
    msg_queue = mp.Queue()

    # Setup output sql table
    reset_sql_tables(row=df.iloc[0])

    # Calculate the total number of messages (excluding "saving" and "continuing")
    total_messages = len(df['seq'])

    # Start the queue printer process
    stop_event = mp.Event()
    log_filename = "process_log.txt"
    printer_process = mp.Process(target=queue_printer, args=(msg_queue, stop_event, log_filename, total_messages))
    printer_process.start()

    # Perform calculation
    print("Starting Calculation")
    df_split = np.array_split(df, num_cores)
    pool = mp.Pool(num_cores, initializer=init, initargs=(msg_queue,))
    pool.map(calc_and_save, df_split)
    pool.close()
    pool.join()

    # Signal the process to stop and wait for it to finish
    stop_event.set()
    printer_process.join()

    print("Calculation complete.")
