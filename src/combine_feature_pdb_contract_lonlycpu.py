import pandas as pd
import time

def process_files(contract_file_path, updated_file_path, chunksize=100):
    # 完整读取 contract_ligand_output.csv 文件并进行处理
    contract_df = pd.read_csv(contract_file_path)
    contract_df = contract_df.rename(columns={'Unnamed: 0': 'filename'})
    contract_df['outputn'] = contract_df['filename'].apply(lambda x: x.split('_')[0])

    # 将数据转换成更高效的格式，例如字典，用于快速匹配
    grouped_data = contract_df.groupby('outputn').apply(lambda df: df.set_index('filename').to_dict(orient='index')).to_dict()

    # 初始化一个空的 DataFrame 用于存储最终结果
    final_result_df = pd.DataFrame()

    # 以分块方式读取 updated_ligand_rmsd.csv 文件
    for chunk_index, chunk in enumerate(pd.read_csv(updated_file_path, chunksize=chunksize)):
        def match_and_merge(row):
            outputn = row['filename'].split('_')[0]
            if outputn in grouped_data and row['filename'] in grouped_data[outputn]:
                return pd.Series({**row, **grouped_data[outputn][row['filename']]})
            return row

        chunk_result_df = chunk.apply(match_and_merge, axis=1)
        print(chunk_index)
        # 首个块写入表头，之后的块追加数据
        if chunk_index == 0:
            chunk_result_df.to_csv("final_output.csv", index=False)
        else:
            chunk_result_df.to_csv("final_output.csv", mode='a', index=False, header=False)


print(time.ctime())
# 调用函数处理文件
output_file = process_files(contract_file_path = '/home/lbbe02/ZSD/pdb/contract_ligand_output2.csv', updated_file_path = '/home/lbbe02/ZSD/pdb/try1.csv')
print(time.ctime())