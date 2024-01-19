import torch
import os
import pandas as pd

def load_embeddings(file_path):
    """
    加载.pt文件中的embeddings
    :param file_path: .pt文件的路径
    :return: embeddings['mean_representations'][33]的内容
    """
    try:
        embeddings = torch.load(file_path)
        return embeddings['mean_representations'][33].numpy()
    except Exception as e:
        print(f"Error occurred while loading embeddings from {file_path}: {e}")
        return None

def gather_embeddings_to_csv(folder_path, output_csv):
    """
    遍历指定文件夹中的.pt文件，提取embeddings并将其保存到CSV文件中
    :param folder_path: 包含.pt文件的文件夹路径
    :param output_csv: 输出CSV文件的路径
    """
    embeddings_list = []
    for file_name in os.listdir(folder_path):
        if file_name.endswith('.pt'):
            file_path = os.path.join(folder_path, file_name)
            embedding = load_embeddings(file_path)
            if embedding is not None:
                # 将文件名作为第一个元素加入列表
                protein_name = file_name.replace('_.pt','')
                embeddings_list.append([protein_name] + list(embedding))

    # 将数据转换为DataFrame并保存到CSV
    df = pd.DataFrame(embeddings_list)
    # 设置第一列的列名为protein_name
    df.columns = ['protein_name'] + ['feature_' + str(i) for i in range(df.shape[1]-1)]
    df.to_csv(output_csv, index=False)
    print(f"Embeddings saved to {output_csv}")



# 文件夹路径和输出CSV文件路径
folder_path = '/mnt/e/SDT/input/protein_embedding/'
output_csv = '/mnt/e/SDT/input/protein_embedding.csv'  # 修改为您希望保存CSV文件的路径

import pandas as pd

def merge_data(actives_csv, protein_csv, output_csv):
    """
    合并actives_output.csv和protein_embedding.csv文件
    :param actives_csv: actives_output.csv文件的路径
    :param protein_csv: protein_embedding.csv文件的路径
    :param output_csv: 输出CSV文件的路径
    """
    # 读取CSV文件
    actives_df = pd.read_csv(actives_csv)
    actives_df = actives_df.rename(columns={'Unnamed: 0': 'filename'})
    protein_df = pd.read_csv(protein_csv)

    # 提取protein_name
    actives_df['protein_name'] = actives_df['filename'].apply(lambda x: x.split('_')[0])

    # 合并数据
    merged_df = actives_df.merge(protein_df, on='protein_name', how='left')

    # 删除额外的protein_name列
    # merged_df.drop(columns=['protein_name'], inplace=True)

    # 保存到新的CSV文件
    merged_df.to_csv(output_csv, index=False)
    print(f"Merged data saved to {output_csv}")

if __name__ == "__main__":
    # 文件路径
    actives_csv = '/mnt/e/SDT/input/actives_output.csv'
    protein_csv = '/mnt/e/SDT/input/protein_embedding.csv'
    output_csv = '/mnt/e/SDT/input/actives_combine.csv'  # 修改为您希望保存合并后CSV文件的路径

    merge_data(actives_csv, protein_csv, output_csv)

 