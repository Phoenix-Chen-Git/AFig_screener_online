import argparse
from utils import *
def main():

    # parser add arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-ori','--origin_pdb_dir', type=str, required=True,help='pdbs_before_AFig')
    parser.add_argument('-g','--generated_pdb_dir', type=str, required=True,help='pdbs_after_AFig')
    parser.add_argument('-p','--project_name', type=str, required=True,help='project_name')
    parser.add_argument('-iter','--iteration', type=int, required=True,help='iteration')
    parser.add_argument('-a','--af2igpae_dir', type=str, required=True,help='af2igpae_dir')

    # args
    args = parser.parse_args()
    origin_pdb_dir = args.origin_pdb_dir
    generated_pdb_dir = args.generated_pdb_dir
    project_name = args.project_name
    iteration = args.iteration
    next_iter = iteration + 1
    raw_results_csv = project_name + '_Aligned_' + str(iteration) + '.csv'
    af2igpae_dir = args.af2igpae_dir
    screened_results_csv = project_name + '_Screened_' + str(iteration) + '.csv'
    output_dir = project_name + '_' + str(next_iter)+'candidates'
    binder_sequence_csv = project_name + '_binder_sequences_' + str(iteration) + '.csv'
    full_score_csv = project_name + '_full_score_' + str(iteration) + '.csv'

    # read pdbs
    origin_pdb = os.listdir(origin_pdb_dir)
    generated_pdb = os.listdir(generated_pdb_dir)
    df_origin = pd.DataFrame({
        'origin_pdb': origin_pdb,
        'key': [get_key(x) for x in origin_pdb]
    })


    df_generated = pd.DataFrame({
        'generated_pdb': generated_pdb,
        'key': [get_key_AF(x) for x in generated_pdb]
    })
    print(df_origin.head())
    print(df_generated.head())

    # merge generated pdb and original pdbs
    merged_df = pd.merge(df_origin, df_generated, on='key', how='inner')
    print(merged_df.shape)

    # calculate TMscore(determine if off target)
    futures = {}
    with ThreadPoolExecutor(max_workers=14) as executor:
        for pdb1, pdb2 in zip(
                origin_pdb_dir + merged_df["origin_pdb"],
                generated_pdb_dir + merged_df["generated_pdb"]
        ):
            future = executor.submit(get_tmscore, pdb1, pdb2)
            futures[future] = (pdb1, pdb2)  # 保存任务对应的文件名

        results = []
        for f in tqdm(as_completed(futures), total=len(futures), desc="Running USalign"):
            pdb1, pdb2 = futures[f]
            tm = f.result()
            results.append({
                "origin_pdb": pdb1,
                "generated_pdb": pdb2,
                "TMscore": tm
            })

    # save alignment result
    results = pd.DataFrame(results)
    results.to_csv(raw_results_csv, index=False)

    # screen out off-targets
    data = pd.read_csv(af2igpae_dir, delim_whitespace=True)
    cut_off = data['pae_interaction'].quantile(0.05)
    print(cut_off)
    results = results[results["TMscore"] > 0.95]
    data['generated_pdb'] = data['description'].apply(lambda x: generated_pdb_dir + x + '.pdb')
    results_data = results.merge(data, how='left', on='generated_pdb')
    results_data.to_csv(full_score_csv, index=False)


    # screen out weak binders
    results_data = results_data[results_data['pae_interaction'] <= cut_off]

    # save screen results
    results_data.to_csv(screened_results_csv, index=False)

    # save candidates for next iteration
    os.makedirs(output_dir, exist_ok=True)
    for i, row in results_data.iterrows():
        src = row['generated_pdb']
        filename = os.path.basename(src)
        dst = os.path.join(output_dir, filename)
        print(f"Copying {src} -> {dst}")
        shutil.copy(src, dst)

    # save candidates' sequences for validation
    parser = PDB.PDBParser(QUIET=True)
    ppb = PDB.PPBuilder()
    results = []
    for i, row in results_data.iterrows():
        pdb_path = row['generated_pdb']
        if not os.path.exists(pdb_path):
            print(f"File not found: {pdb_path}")
            continue
        structure = parser.get_structure("model", pdb_path)
        chainA = structure[0]['A'] if 'A' in structure[0] else None
        if chainA is None:
            print(f"No chain A found in {pdb_path}")
            continue
        seq = ''.join([str(pp.get_sequence()) for pp in ppb.build_peptides(chainA)])
        results.append({
            "pdb_name": os.path.basename(pdb_path),
            "chain": "A",
            "sequence": seq
        })
    df = pd.DataFrame(results)
    df.to_csv(binder_sequence_csv, index=False)
    print("Results saved to chainA_sequences.csv")


if __name__ == '__main__':
    main()