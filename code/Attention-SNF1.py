import os
import random
import gc
import numpy as np
import pandas as pd
import tensorflow as tf
import seaborn as sns
from keras.models import Model
from keras.layers import (
    Input, Dense, Dropout, BatchNormalization, 
    MultiHeadAttention, concatenate, LayerNormalization,
    Lambda, Flatten
)
from keras.optimizers import adam_v2
from keras.regularizers import l1_l2
from keras.callbacks import EarlyStopping, ReduceLROnPlateau
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
import json
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
    try:
        for gpu in gpus:
            tf.config.experimental.set_memory_growth(gpu, True)
    except RuntimeError as e:
        print(e)

def reset_random_seeds(seed):
    os.environ['PYTHONHASHSEED'] = str(seed)
    tf.random.set_seed(seed)
    np.random.seed(seed)
    random.seed(seed)

def load_data():
    # omics1 = pd.read_csv("omics1_data.csv").drop("Unnamed: 0", axis=1).values
    # omics2 = pd.read_csv("omics2_data.csv").drop("Unnamed: 0", axis=1).values
    # omics3 = pd.read_csv("omics3_data.csv").drop("Unnamed: 0", axis=1).values
    # omics4 = pd.read_csv("omics4_data.csv").drop("Unnamed: 0", axis=1).values

    tmb_target_dim=10
    mirna_target_dim=164
    gene_target_dim=164

    tmb_df = pd.read_csv('tmb1.csv', index_col=0)#omics1
    micro_df = pd.read_csv('wei1.csv', index_col=0)#omics2
    mirna_df = pd.read_csv('rna1.csv', index_col=0)#omics3
    gene_df = pd.read_csv('dna1.csv', index_col=0)#omics4

    common_samples = list(set(tmb_df.index) & set(micro_df.index) & set(mirna_df.index) & set(gene_df.index))
    tmb_df = tmb_df.loc[common_samples]
    micro_df = micro_df.loc[common_samples]
    mirna_df = mirna_df.loc[common_samples]
    gene_df = gene_df.loc[common_samples]

    micro_df = micro_df.fillna(0)
    mirna_df = mirna_df.fillna(0)
    gene_df = gene_df.fillna(0)

    tmb_data = np.repeat(tmb_df.values, tmb_target_dim, axis=1)
    
    micro_data = micro_df.values 
    
    mirna_data = mirna_df.values
    pca1 = PCA(n_components=mirna_target_dim)
    mirna_data = pca1.fit_transform(mirna_data)
    
    gene_data=gene_df.values
    pca2 = PCA(n_components=gene_target_dim)
    gene_data = pca2.fit_transform(gene_data)
    
    print(f"预处理后维度: TMB={tmb_data.shape[1]}, 微生物={micro_data.shape[1]}, "
          f"miRNA={mirna_data.shape[1]}, 基因={gene_data.shape[1]}")
    
    return tmb_data, micro_data, mirna_data, gene_data

def create_adaptive_encoder(input_dim, output_dim, name):
    inputs = Input(shape=(input_dim,), name=f"{name}_input")
    
    if input_dim < 20:  
        x = Dense(64, activation="relu", 
                 kernel_regularizer=l1_l2(0.001, 0.001))(inputs)
        x = BatchNormalization()(x)
        x = Dropout(0.3)(x)
        
        x = Dense(128, activation="relu")(x)
        x = BatchNormalization()(x)
        x = Dropout(0.2)(x)
        
    elif input_dim < 500: 
        hidden_units = max(128, min(512, input_dim))
        x = Dense(hidden_units, activation="relu", 
                 kernel_regularizer=l1_l2(0.001, 0.001))(inputs)
        x = BatchNormalization()(x)
        x = Dropout(0.4)(x)
        
        x = Dense(hidden_units//2, activation="relu")(x)
        x = BatchNormalization()(x)
        x = Dropout(0.3)(x)
        
    else:
        x = Dense(1024, activation="relu", 
                 kernel_regularizer=l1_l2(0.001, 0.001))(inputs)
        x = BatchNormalization()(x)
        x = Dropout(0.5)(x)
        
        x = Dense(512, activation="relu")(x)
        x = BatchNormalization()(x)
        x = Dropout(0.4)(x)
        
        x = Dense(256, activation="relu")(x)
        x = BatchNormalization()(x)
        x = Dropout(0.3)(x)
    
    encoded = Dense(output_dim, activation="relu", name=f"{name}_encoded")(x)
    return Model(inputs, encoded, name=f"{name}_encoder")

def create_decoder(input_dim, output_dim, name):
    inputs = Input(shape=(input_dim,), name=f"{name}_latent_input")
    
    hidden_units = max(128, min(1024, output_dim // 2))
    
    x = Dense(max(64, hidden_units // 2), activation="relu")(inputs)
    x = BatchNormalization()(x)
    x = Dropout(0.2)(x)
    
    x = Dense(hidden_units, activation="relu")(x)
    x = BatchNormalization()(x)
    x = Dropout(0.3)(x)
    
    decoded = Dense(output_dim, activation="sigmoid", name=f"{name}_decoded")(x)
    return Model(inputs, decoded, name=f"{name}_decoder")


def self_attention(x, head_size=32, num_heads=4, name=None):

    x = LayerNormalization()(x)
    
    x_exp = tf.expand_dims(x, axis=1)
    
    attention = MultiHeadAttention(num_heads=num_heads, key_dim=head_size, 
                                  name=f"{name}_self_att")(x_exp, x_exp)
    
    return attention[:, 0, :], attention 

def create_fusion_attention(omics_encoded, mode='MM_SA'):

    omics1, omics2, omics3, omics4 = omics_encoded
    
    self_att_features = {}
    cross_att_features = {}
    
    if mode == 'MM_SA':
        omics1_att, att1 = self_attention(omics1, name="omics1")
        omics2_att, att2 = self_attention(omics2, name="omics2")
        omics3_att, att3 = self_attention(omics3, name="omics3")
        omics4_att, att4 = self_attention(omics4, name="omics4")
        
        self_att_features = {
            "omics1": omics1_att,
            "omics2": omics2_att,
            "omics3": omics3_att,
            "omics4": omics4_att
        }
        
        fused = concatenate([omics1_att, omics2_att, omics3_att, omics4_att])
        
    else:
        fused = concatenate([omics1, omics2, omics3, omics4])
    
    fused = Dense(512, activation='relu', name='fusion_dense1')(fused)
    fused = BatchNormalization(name='fusion_bn1')(fused)
    fused = Dropout(0.4, name='fusion_dropout1')(fused)
    
    fused = Dense(256, activation='relu', name='fusion_dense2')(fused)
    fused = BatchNormalization(name='fusion_bn2')(fused)
    fused = Dropout(0.3, name='fusion_dropout2')(fused)
    
    fused = Dense(128, activation='relu', name='fused_features')(fused)
    
    return fused, self_att_features, cross_att_features

def calculate_dynamic_loss_weights(omics_data, base_weight=0.25):
    variances = [np.var(data, axis=0).mean() for data in omics_data]
    total_variance = sum(variances)
    
    weights = [base_weight * (total_variance / v) for v in variances]

    weight_sum = sum(weights)
    normalized_weights = [w/weight_sum for w in weights]
    
    print(f"动态损失权重: {normalized_weights}")
    return normalized_weights

def build_multi_omic_autoencoder(omics_dims, latent_dim=64, mode='MM_SA_BA'):
    inputs = []
    for i, dim in enumerate(omics_dims):
        inputs.append(Input(shape=(dim,), name=f'omics{i+1}_input'))
    
    encoders = []
    encoded_features = []
    for i, (inp, dim) in enumerate(zip(inputs, omics_dims)):
        encoder = create_adaptive_encoder(dim, latent_dim, f'omics{i+1}')
        encoded = encoder(inp)
        encoders.append(encoder)
        encoded_features.append(encoded)
    
    fused_features, self_att_features, cross_att_features = create_fusion_attention(
        encoded_features, mode
    )
    
    decoders = []
    decoded_outputs = []
    for i, (dim, encoder) in enumerate(zip(omics_dims, encoders)):
        decoder = create_decoder(128, dim, f'omics{i+1}') 
        decoded = decoder(fused_features)
        decoders.append(decoder)
        decoded_outputs.append(decoded)
    
    autoencoder = Model(
        inputs=inputs,
        outputs=decoded_outputs,
        name='multi_omic_autoencoder'
    )
    
    fusion_encoder = Model(
        inputs=inputs,
        outputs=fused_features,
        name='fusion_feature_extractor'
    )
    
    self_att_encoder = Model(
        inputs=inputs,
        outputs=self_att_features,
        name='self_attention_extractor'
    )
    
    cross_att_encoder = Model(
        inputs=inputs,
        outputs=cross_att_features,
        name='cross_attention_extractor'
    )
    
    return autoencoder, fusion_encoder, self_att_encoder, cross_att_encoder

def train_autoencoder(autoencoder, train_data, val_data, epochs=150, batch_size=64, lr=0.0005):

    loss_weights = calculate_dynamic_loss_weights(train_data)
    
    callbacks = [
        EarlyStopping(monitor='val_loss', patience=20, restore_best_weights=True),
        ReduceLROnPlateau(monitor='val_loss', factor=0.5, patience=10, min_lr=1e-6)
    ]
    
    autoencoder.compile(
        optimizer=adam_v2.Adam(learning_rate=lr),
        loss=['mse', 'mse', 'mse', 'mse'],
        loss_weights=loss_weights, 
        metrics=['mae']
    )
    
    history = autoencoder.fit(
        train_data,
        train_data,  
        epochs=epochs,
        batch_size=batch_size,
        validation_data=(val_data, val_data),
        callbacks=callbacks,
        verbose=1
    )
    
    plt.figure(figsize=(12, 6))
    plt.plot(history.history['loss'], label='Training Loss')
    plt.plot(history.history['val_loss'], label='Validation Loss')
    plt.title('Autoencoder Training History')
    plt.ylabel('Loss')
    plt.xlabel('Epoch')
    plt.legend()
    plt.savefig('autoencoder_training_history.png')
    plt.close()

    with open('training_history.json', 'w') as f:
        history_dict = {}
        for key, value in history.history.items():
            if isinstance(value, list) and len(value) > 0:
                history_dict[key] = [float(v) if hasattr(v, 'dtype') else v for v in value]
            else:
                history_dict[key] = value
        
        json.dump(history_dict, f)
    
    return history

def save_attention_features(encoder, data, output_prefix, feature_names=None):
    features = encoder.predict(data)
    
    if isinstance(features, dict):
        for key, value in features.items():
            np.save(f"{output_prefix}_{key}.npy", value)
            print(f"Saved {key} features to {output_prefix}_{key}.npy")
    else:
        if feature_names:
            for i, name in enumerate(feature_names):
                np.save(f"{output_prefix}_{name}.npy", features[:, i])
        else:
            np.save(f"{output_prefix}.npy", features)
            print(f"Saved features to {output_prefix}.npy")
    
    return features

def visualize_attention(attention_matrix, omics_names, sample_idx=0):
    plt.figure(figsize=(10, 8))
    
    if len(attention_matrix.shape) == 3:
        attention_matrix = attention_matrix[0]

    sns.heatmap(
        attention_matrix[sample_idx, :, :],
        annot=True,
        fmt=".2f",
        cmap="viridis",
        xticklabels=omics_names,
        yticklabels=omics_names
    )
    
    plt.title("Attention Matrix")
    plt.xlabel("Keys")
    plt.ylabel("Queries")
    plt.tight_layout()
    plt.savefig('attention_matrix.png')
    plt.close()

def visualize_reconstruction(original, reconstructed, omics_index, n_samples=5):
    plt.figure(figsize=(15, 8))
    
    sample_indices = np.random.choice(range(len(original[0])), n_samples, replace=False)
    
    for i, idx in enumerate(sample_indices):
        plt.subplot(2, n_samples, i + 1)
        plt.bar(range(len(original[omics_index][idx])), original[omics_index][idx])
        plt.title(f"Original Omics {omics_index+1}")
        plt.ylim(0, 1)
        
        plt.subplot(2, n_samples, i + 1 + n_samples)
        plt.bar(range(len(reconstructed[omics_index][idx])), reconstructed[omics_index][idx])
        plt.title(f"Reconstructed Omics {omics_index+1}")
        plt.ylim(0, 1)
    
    plt.tight_layout()
    plt.savefig(f'reconstruction_omics{omics_index+1}.png')
    plt.close()

def main():
    seed = 42
    reset_random_seeds(seed)
    

    omics1, omics2, omics3, omics4 = load_data()

    omics1 = omics1.astype('float32') / np.max(omics1, axis=0, keepdims=True)
    omics2 = omics2.astype('float32') / np.max(omics2, axis=0, keepdims=True)
    omics3 = omics3.astype('float32') / np.max(omics3, axis=0, keepdims=True)
    omics4 = omics4.astype('float32') / np.max(omics4, axis=0, keepdims=True)
    
    train_data = []
    val_data = []
    
    train_idx, val_idx = train_test_split(
        np.arange(len(omics1)), test_size=0.1, random_state=seed
    )
    
    for omics in [omics1, omics2, omics3, omics4]:
        train_data.append(omics[train_idx])
        val_data.append(omics[val_idx])
    
    omics_dims = [omics1.shape[1], omics2.shape[1], 
                 omics3.shape[1], omics4.shape[1]]
    
    print(f"Omics dimensions: {omics_dims}")
    print(f"Training samples: {len(train_idx)}")
    print(f"Validation samples: {len(val_idx)}")
    
    autoencoder, fusion_encoder, self_att_encoder, cross_att_encoder = build_multi_omic_autoencoder(
        omics_dims=omics_dims,
        latent_dim=64, 
        mode='MM_SA'
    )
    
    autoencoder.summary()
    
    history = train_autoencoder(
        autoencoder,
        train_data,
        val_data,
        epochs=150,
        batch_size=64,
        lr=0.0005
    )
    
    fusion_features = save_attention_features(
        fusion_encoder, val_data, "features/fusion"
    )
    
    self_att_features = save_attention_features(
        self_att_encoder, val_data, "features/self_attention"
    )
    
    cross_att_features = save_attention_features(
        cross_att_encoder, val_data, "features/cross_attention"
    )
    
    sample_idx = 0
    omics_names = ["omics1", "omics2", "omics3", "omics4"]
    
    reconstructed = autoencoder.predict(val_data)
    visualize_reconstruction(val_data, reconstructed, omics_index=0, n_samples=5)
    
    autoencoder.save("models/multi_omic_autoencoder.h5")
    fusion_encoder.save("models/fusion_feature_extractor.h5")
    self_att_encoder.save("models/self_attention_extractor.h5")
    #cross_att_encoder.save("models/cross_attention_extractor.h5")
    print("Models saved successfully")
    
    del autoencoder, fusion_encoder, self_att_encoder, cross_att_encoder
    gc.collect()

if __name__ == "__main__":
    os.makedirs("features", exist_ok=True)
    os.makedirs("models", exist_ok=True)
    main()