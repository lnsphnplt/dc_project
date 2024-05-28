import h5py
import numpy as np
import pandas as pd
from propy import PyPro
from rdkit import Chem
from rdkit.Chem import Descriptors
import os
import multiprocessing as mp
import sqlite3
import warnings

# The following statement was the work of 4h of debugging....
file_lock = mp.Lock()

# This is surely not causing a headache later on...
warnings.filterwarnings("ignore")

def calculate_rdkit_descriptors(seq: str):
    """ Calculates all descriptors from rdkit and returns a df with one row """
    mol = Chem.MolFromSequence(seq)
    return pd.DataFrame([Descriptors.CalcMolDescriptors(mol, missingVal=None, silent=True)])

def calculate_propy_descriptors(seq: str):
    """ Calculates relevant descriptors from propy3 and returns a df with one row """
    total_descriptors = 10000  # Ensure this matches the total number of descriptors being calculated
    df = pd.DataFrame(np.nan, index=[0], columns=np.arange(total_descriptors))

    def safe_update(method, df, start_index, stop_index):
        """ Helper function to this without wrapping each statement in a try catch clause """
        try:
            result = method()
            result_values = list(result.values())
            df.iloc[0, start_index:stop_index] = result_values[:stop_index - start_index]
        except Exception as e:
            print(f"Error: {e}", end="")

    pro = PyPro.GetProDes(seq)

    # Get composition (20)
    safe_update(pro.GetAAComp, df, 0, 20)

    # Get dipeptide composition (400)
    safe_update(pro.GetDPComp, df, 20, 420)

    # Get tripeptide sequence (8000)
    safe_update(pro.GetTPComp, df, 420, 8420)

    # Get Composition Transition Distribution descriptors (147)
    safe_update(pro.GetCTD, df, 8420, 8567)

    # Get Geary autocorrelation descriptors (240)
    safe_update(pro.GetGearyAuto, df, 8567, 8807)

    # Get Moran autocorrelation descriptors (240)
    safe_update(pro.GetMoranAuto, df, 8807, 9047)

    # Get Normalized Moreau-Broto autocorrelation descriptors (240)
    safe_update(pro.GetMoreauBrotoAuto, df, 9047, 9287)

    # Get Quasi sequence order descriptors (default is 50)
    safe_update(pro.GetQSO, df, 9287, 9337)

    # Get Sequence order coupling numbers (default is 45)
    safe_update(pro.GetSOCN, df, 9337, 9382)

    # Get Type I Pseudo amino acid composition descriptors (default is 30)
    safe_update(pro.GetPAAC, df, 9382, 9412)

    # Get Amphiphilic (Type II) Pseudo amino acid composition descriptors (30)
    safe_update(pro.GetAPAAC, df, 9412, 9442)

    return df

def create_hdf5_file(file="data2.h5"):
    if not os.path.exists(file):
        with h5py.File(file, 'w') as f:
            f.create_dataset("headers", (0,), maxshape=(None,), dtype=h5py.string_dtype())
            f.create_dataset("data", (0, 10000), maxshape=(None, 10000), dtype='f')

def append_to_hdf5(df, file="data2.h5", dataset_name='df'):
    with h5py.File(file, 'a') as f:
        headers = list(df.columns.values)
        if 'headers' in f and len(f['headers']) > 0:
            existing_headers = f['headers'][:]
            all_headers = list(existing_headers)
        else:
            f.create_dataset('headers', data=headers, maxshape=(None,))
            all_headers = headers
        
        data_to_append = df.reindex(columns=all_headers, fill_value=np.nan).values
        data_shape = f['data'].shape
        f['data'].resize((data_shape[0] + data_to_append.shape[0], data_shape[1]))
        f['data'][-data_to_append.shape[0]:, :data_to_append.shape[1]] = data_to_append

def calculate_descriptors(seq):
    return pd.concat([calculate_propy_descriptors(seq), calculate_rdkit_descriptors(seq)], axis=1)

def read_sqlite():
    con = sqlite3.connect('04-Data_calculate_Descriptors/euler/unified.db')
    df = pd.read_sql_query("SELECT * FROM prod", con)
    con.close()
    return df

def worker(df, save_interval=30):
    print("Worker starting")
    tmp = pd.DataFrame()
    for i in df.index:
        tmp = pd.concat([tmp, calculate_descriptors(df.iloc[i]['seq'])], axis=0, join='outer')
        if not i == 0 and i % save_interval == 0:
            append_to_hdf5(tmp)
            tmp = tmp[0:0]
            print("saving")
    append_to_hdf5(tmp)  
    print("Worker finished")

if __name__ == "__main__":
    print("Setting up env")
    create_hdf5_file()

    num_cores = os.cpu_count()
    num_cores = 10 # automatic core detection does not always work

    # import sequences
    print("Reading data")
    df = read_sqlite()

    # split data
    print("Spliting data")
    df_split = np.array_split(df, num_cores)

    # Create processes and start them
    processes = []
    for i in range(num_cores):
        p = mp.Process(target=worker, args=(df_split[i], 3))
        processes.append(p)
        p.start()

    # Ensure all processes have finished execution
    for p in processes:
        p.join()

    print("Processing complete")
