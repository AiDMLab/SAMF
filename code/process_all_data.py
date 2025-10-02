import os
import numpy as np
import pandas as pd
import tensorflow as tf
from keras.models import load_model
import gc
from sklearn.decomposition import PCA

def load_and_preprocess_data():
    print("加载原始数据...")
    
    tmb_df = pd.read_csv('tmb1.csv', index_col=0)
    micro_df = pd.read_csv('wei1.csv', index_col=0)
    mirna_df = pd.read_csv('rna1.csv', index_col=0)
    gene_df = pd.read_csv('dna1.csv', index_col=0)
    
    common_samples = list(set(tmb_df.index) & set(micro_df.index) & set(mirna_df.index) & set(gene_df.index))
    tmb_df = tmb_df.loc[common_samples]
    micro_df = micro_df.loc[common_samples]
    mirna_df = mirna_df.loc[common_samples]
    gene_df = gene_df.loc[common_samples]

    micro_df = micro_df.fillna(0)
    mirna_df = mirna_df.fillna(0)
    gene_df = gene_df.fillna(0)
    
    sample_ids = common_samples
    
    tmb_data = np.repeat(tmb_df.values, 10, axis=1)
    
    micro_data = micro_df.values 

    mirna_data = mirna_df.values
    pca1 = PCA(n_components=164)
    mirna_data = pca1.fit_transform(mirna_data)
    
    gene_data=gene_df.values
    pca2 = PCA(n_components=164)
    gene_data = pca2.fit_transform(gene_data)
    
    print(f"数据加载完成: TMB={tmb_data.shape}, Micro={micro_data.shape}, "
          f"miRNA={mirna_data.shape}, Gene={gene_data.shape}")
    
    return [tmb_data, micro_data, mirna_data, gene_data], sample_ids

def load_trained_models():
    print("加载训练好的模型...")
    
    models_dir = "models"
    
    autoencoder_path = os.path.join(models_dir, "multi_omic_autoencoder.h5")
    if not os.path.exists(autoencoder_path):
        raise FileNotFoundError(f"自编码器模型未找到: {autoencoder_path}")
    
    autoencoder = load_model(autoencoder_path, compile=False)
    print("自编码器模型加载完成")
    
    fusion_encoder_path = os.path.join(models_dir, "fusion_feature_extractor.h5")
    self_att_encoder_path = os.path.join(models_dir, "self_attention_extractor.h5")
    cross_att_encoder_path = os.path.join(models_dir, "cross_attention_extractor.h5")
    
    fusion_encoder = load_model(fusion_encoder_path, compile=False) if os.path.exists(fusion_encoder_path) else None
    self_att_encoder = load_model(self_att_encoder_path, compile=False) if os.path.exists(self_att_encoder_path) else None
    cross_att_encoder = load_model(cross_att_encoder_path, compile=False) if os.path.exists(cross_att_encoder_path) else None
    
    print("特征提取器模型加载完成")
    
    return autoencoder, fusion_encoder, self_att_encoder, cross_att_encoder

def process_data_with_models(models, data, sample_ids, batch_size=32):
    autoencoder, fusion_encoder, self_att_encoder, cross_att_encoder = models
    
    print("开始处理数据...")

    print("获取重建数据...")
    reconstructed_data = autoencoder.predict(data, batch_size=batch_size, verbose=1)

    print("获取融合特征...")
    if fusion_encoder:
        fusion_features = fusion_encoder.predict(data, batch_size=batch_size, verbose=1)
    else:
        fusion_features = None
        print("警告: 融合编码器未找到")

    print("获取自注意力特征...")
    if self_att_encoder:
        self_att_features = self_att_encoder.predict(data, batch_size=batch_size, verbose=1)
    else:
        self_att_features = None
        print("警告: 自注意力编码器未找到")

    print("获取跨模态注意力特征...")
    if cross_att_encoder:
        cross_att_features = cross_att_encoder.predict(data, batch_size=batch_size, verbose=1)
    else:
        cross_att_features = None
        print("警告: 跨模态注意力编码器未找到")
    
    return {
        'reconstructed': reconstructed_data,
        'fusion': fusion_features,
        'self_attention': self_att_features,
        'cross_attention': cross_att_features
    }, sample_ids

def save_features_to_csv(features_dict, sample_ids, output_dir="processed_features"):
    print("保存特征为CSV文件...")

    os.makedirs(output_dir, exist_ok=True)

    if features_dict['reconstructed']:
        reconstructed_data = features_dict['reconstructed']
        omics_names = ["TMB", "Microbiome", "miRNA", "Gene"]
        
        for i, (omics_name, data) in enumerate(zip(omics_names, reconstructed_data)):
            df = pd.DataFrame(data, index=sample_ids)
            feature_cols = [f"{omics_name}_feature_{j}" for j in range(data.shape[1])]
            df.columns = feature_cols
            df.index.name = "sample_id"
            
            filename = f"reconstructed_{omics_name.lower()}.csv"
            filepath = os.path.join(output_dir, filename)
            df.to_csv(filepath)
            print(f"已保存: {filename}")

    if features_dict['fusion'] is not None:
        fusion_data = features_dict['fusion']
        df = pd.DataFrame(fusion_data, index=sample_ids)
        feature_cols = [f"fusion_feature_{j}" for j in range(fusion_data.shape[1])]
        df.columns = feature_cols
        df.index.name = "sample_id"
        
        filename = "fusion_features.csv"
        filepath = os.path.join(output_dir, filename)
        df.to_csv(filepath)
        print(f"已保存: {filename}")

    if features_dict['self_attention'] is not None:
        self_att_data = features_dict['self_attention']

        for key, value in self_att_data.items():
            df = pd.DataFrame(value, index=sample_ids)
            feature_cols = [f"self_att_{key}_feature_{j}" for j in range(value.shape[1])]
            df.columns = feature_cols
            df.index.name = "sample_id"
            
            filename = f"self_attention_{key}.csv"
            filepath = os.path.join(output_dir, filename)
            df.to_csv(filepath)
            print(f"已保存: {filename}")

    if features_dict['cross_attention'] is not None:
        cross_att_data = features_dict['cross_attention']

        for key, value in cross_att_data.items():
            df = pd.DataFrame(value, index=sample_ids)
            feature_cols = [f"cross_att_{key}_feature_{j}" for j in range(value.shape[1])]
            df.columns = feature_cols
            df.index.name = "sample_id"

            safe_key = key.replace("-", "_to_")
            filename = f"cross_attention_{safe_key}.csv"
            filepath = os.path.join(output_dir, filename)
            df.to_csv(filepath)
            print(f"已保存: {filename}")
    
    print("所有特征已保存为CSV文件")

def main():
    print("开始处理所有数据...")
    
    try:
        data, sample_ids = load_and_preprocess_data()
        
        models = load_trained_models()
        
        features_dict, sample_ids = process_data_with_models(models, data, sample_ids)

        save_features_to_csv(features_dict, sample_ids)
        
        print("数据处理和保存完成!")
        
    except Exception as e:
        print(f"处理过程中发生错误: {str(e)}")
        import traceback
        traceback.print_exc()
    
    finally:
        gc.collect()
        print("内存清理完成")

if __name__ == "__main__":
    gpus = tf.config.experimental.list_physical_devices('GPU')
    if gpus:
        try:
            for gpu in gpus:
                tf.config.experimental.set_memory_growth(gpu, True)
        except RuntimeError as e:
            print(e)
    
    main()