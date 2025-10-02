import pandas as pd

def align_csv_order(file_a, file_b, output_file):
    df_a = pd.read_csv(file_a)
    df_b = pd.read_csv(file_b)

    col_a = df_a.columns[0]
    col_b = df_b.columns[0]
    
    print(f"文件A的第一列名称: {col_a}")
    print(f"文件B的第一列名称: {col_b}")

    if col_a != col_b:
        print(f"警告: 两个文件的第一列名称不同 ({col_a} vs {col_b})")
        print("将使用文件B的第一列名称作为参考")

    order_mapping = {sample: idx for idx, sample in enumerate(df_b[col_b])}

    missing_in_b = set(df_a[col_a]) - set(df_b[col_b])
    if missing_in_b:
        print(f"警告: A文件中有{len(missing_in_b)}个样本在B文件中不存在")
        print("这些样本将被放置在排序后的文件末尾")

    max_index = len(df_b)
    df_a['temp_order'] = df_a[col_a].map(lambda x: order_mapping.get(x, max_index + 1))

    df_a_sorted = df_a.sort_values('temp_order')

    df_a_sorted = df_a_sorted.drop(columns=['temp_order'])

    df_a_sorted.to_csv(output_file, index=False)
    
    print(f"文件已保存至: {output_file}")
    print(f"原始A文件样本数: {len(df_a)}")
    print(f"调整后A文件样本数: {len(df_a_sorted)}")


if __name__ == "__main__":
    align_csv_order('self_attention_omics1.csv', 'tmb1.csv', 'tmb_aligned.csv')
    align_csv_order('self_attention_omics2.csv', 'wei1.csv', 'wei_aligned.csv')
    align_csv_order('self_attention_omics3.csv', 'rna1.csv', 'rna_aligned.csv')
    align_csv_order('self_attention_omics4.csv', 'dna1.csv', 'dna_aligned.csv')
    align_csv_order('fusion_features.csv', 'dna1.csv', 'fus_aligned.csv')