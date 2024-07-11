from typing import List
import numpy as np
import tqdm
import torch
import category_encoders

# local import
def plant_seed(seed: int, USE_GPU:bool = True) -> None:
    """Set seed for reproducibility purposes.
    
    Applies to the `random`, `numpy` and `torch` packages.
    
    Args:
        seed: Seed state.
    """
    np.random.seed(seed)
    torch.manual_seed(seed)
    if USE_GPU:
        torch.cuda.manual_seed_all(seed)


class MultiOutputTargetEncoder:
    """Multi-output target encoder.
    
    Each input (categorical) feature will be encoded based on each (continuous) output variable.
    
    Attributes:
        encoders: List of encoders, of shape `n_genes`.
    """
    
    def __init__(self):
        self.encoders: List[category_encoders.leave_one_out.LeaveOneOutEncoder] = []
        
    @staticmethod
    def new_encoder() -> category_encoders.leave_one_out.LeaveOneOutEncoder:
        """Instantiates a new gene-specific target encoder.
        
        Returns:
            New instance of a target encoder.
        """
        return category_encoders.leave_one_out.LeaveOneOutEncoder(return_df=False)
    
    def fit(self, X: np.ndarray, y: np.ndarray) -> None:
        """Fit the encoders for each input feature and output variable.
        
        Args:
            X: Array of shape `(n, n_features)` containing categories as strings or integers.
                Typical, `n_features` is equal to 2 (cell type + compound).
            y: Array of shape `(n, n_genes)` containing the DE values for all the genes.
        """
        self.encoders = []
        for j in tqdm.tqdm(range(y.shape[1]), desc='fit LOO encoders'):
            self.encoders.append(MultiOutputTargetEncoder.new_encoder())
            self.encoders[-1].fit(X, y[:, j])
    
    def transform(self, X: np.ndarray) -> np.ndarray:
        """Encodes the categories. Assumes the encoders have already been fitted.
        
        Args:
            X: Array of shape `(n, n_features)` containing categories as strings or integers.
        
        Returns:
            Array of shape `(n, n_genes, n_features)` containing the encoding of each input
                feature with respect to each output variable.
        """
        Z = []
        for encoder in tqdm.tqdm(self.encoders, desc='transform LOO encoders'):
            y_hat = encoder.transform(X)
            Z.append(y_hat)
        Z = np.asarray(Z)
        return np.transpose(Z, (1, 0, 2))

class NN(torch.nn.Module):
    """Deep learning architecture composed of 2 modules: a sample-centric MLP
    and a gene-centric MLP.
    
    Attributes:
        n_genes: Total number of genes
        n_input_channels: Number of input channels. When using target encoding for both
            cell types and compounds, this results in one channel for the cell type and
            one for the compound.
        n_output_channels: Number of channels outputed by the sample-centric MLP.
        mlp1: Sample-centric MLP.
        mlp2: Gene-centric MLP.
        net1: Sparse MLP defined by the gene-pathway network.
        net2: Sparse MLP defined by the co-accessibility network.
        net3: Sparse MLP defined by the gene regulatory network.
    """
    
    def __init__(self, n_genes: int, n_input_channels: int):
        torch.nn.Module.__init__(self)
        
        # Total number of genes
        self.n_genes: int = n_genes
            
        # Number of input channels = 2 (cell type + compound)
        self.n_input_channels: int = n_input_channels
            
        # Number of channels outputed by the first MLP = 2 (an arbitrary number)
        self.n_output_channels: int = 2
            
        # No bias term is used in the full-connected layers, to encourage the model
        # to learn an output distribution for which the modal value is 0 (for most genes).
        bias = False
        
        # First MLP module to treat each row of the original data matrix as individual observation.
        # In other words,the MLP processes each (cell_type, compound) combination separately.
        # Shape: (batch_size, n_input_channels*n_genes) -> (batch_size, n_output_channels*n_genes)
        self.mlp1 = torch.nn.Sequential(

            torch.nn.Linear(self.n_genes * self.n_input_channels, 16, bias=bias),
            torch.nn.PReLU(num_parameters=16),
            
            torch.nn.Linear(16, 128, bias=bias),
            torch.nn.PReLU(num_parameters=128),
            
            torch.nn.Linear(128, 128, bias=bias),
            torch.nn.PReLU(num_parameters=128),
            
            torch.nn.Linear(128, 64, bias=bias),
            torch.nn.PReLU(num_parameters=64),
            
            torch.nn.Linear(64, 32, bias=bias),
            torch.nn.PReLU(num_parameters=32),

            torch.nn.Linear(32, 16, bias=bias),
            torch.nn.PReLU(num_parameters=16),
            
            torch.nn.Linear(16, self.n_genes * self.n_output_channels, bias=bias)
        )

        # Second MLP module treats each sample and gene as individual observation.
        # In other words,the MLP processes each (cell_type, compound, gene_name) combination separately.
        # Shape: (batch_size, n_channels) -> (batch_size, 1)
        h = 12
        self.mlp2 = torch.nn.Sequential(

            torch.nn.Linear(self.n_input_channels + self.n_output_channels, h, bias=bias),
            torch.nn.PReLU(num_parameters=h),
            
            torch.nn.Linear(h, h, bias=bias),
            torch.nn.PReLU(num_parameters=h),
            
            torch.nn.Linear(h, h, bias=bias),
            torch.nn.PReLU(num_parameters=h),
            
            torch.nn.Linear(h, h, bias=bias),
            torch.nn.PReLU(num_parameters=h),
            
            torch.nn.Linear(h, 1, bias=bias)
        )
        
        self.net1 = None
        self.net2 = None
        self.net3 = None
        
        # Xavier initialization
        def init_weights(m):
            if isinstance(m, torch.nn.Linear):
                torch.nn.init.xavier_uniform_(m.weight)
                if m.bias is not None:
                    m.bias.data.fill_(0.001)
        self.mlp1.apply(init_weights)
        self.mlp2.apply(init_weights)
        
    def forward(self, X: torch.Tensor) -> torch.Tensor:
        """Predict DE values from the categorical target encoding of DE values.
        
        Args:
            X: Categorical target encoding of cell types and compounds based on DE values.
                `X` is a tensor of shape `(batch_size, n_genes, n_input_channels)`,
                where `n_input_channels` is the number of target-encoded features.
                In the present case, we consider only the cell type and the compound.
        
        Returns:
            Estimated DE values.
        """
        
        # Reshape from (batch_size, n_genes, n_input_channels) to 
        # (batch_size, n_genes*n_input_channels) and process with sample-centric MLP
        Y = self.mlp1(X.reshape(len(X), -1))
        
        # Reshape from (batch_size, n_genes*n_output_channels) to
        # (batch_size*n_genes, n_output_channels)
        Y = Y.reshape(-1, self.n_output_channels)
        
        # Concatenate input and output channels into (batch_size, n_genes*n_channels)
        # and process with gene-centric MLP
        Y = torch.cat((X.reshape(-1, self.n_input_channels), Y), dim=1)
        Y = self.mlp2(Y)
        
        # Output is of shape (batch_size, n_genes)
        Y = Y.reshape(len(X), -1)

        return Y

def background_noise(
    *size: int,
    cutoff: float = 0.001,
    device: str = 'gpu',
    generator: torch.Generator = None) -> torch.Tensor:
    """Generates random DE values under the null hypothesis.

    In the absence of any biological signal, p-values can be expected to be 
    uniformly distributed. Also, positive and negative DE values are equally likely.

    Args:
        size: shape of the output tensor.
        cutoff: Significance threshold used to make sure we don't introduce huge outliers in the data.
            This cutoff does not have a real statistical meaning, and is only meant for numerical stability.
            P-values will be randomly and uniformly sampled from the ``[cutoff, 1]`` interval.
        device: Device where to do the computations (cpu or gpu).
        generator: RNG used for sampling.

    Returns:
        Random DE values, stored in a tensor of the desired shape.
    """
    sign = 2 * torch.randint(0, 2, size, device=device) - 1

    return sign * torch.log10(cutoff +  torch.rand(*size, generator=generator, device=device) * (1. - cutoff))


def train(
        X: torch.Tensor,
        Y: torch.Tensor,
        idx_train: np.ndarray,
        seed: int,
        n_iter: int = 40,  # 200
        USE_GPU: bool = True,
        **kwargs) -> NN:
    """Trains a neural network and returns it.
    
    Args:
        X: Input features. A tensor of shape `(n, n_genes, n_input_channels)`, where `n`
            is the total number of rows in the original data matrix.
        Y: Output variables. A tensor of shape `(n, n_genes)`.
        idx_train: Indices of the training data.
        seed: Seed for reproducibility purposes.
        n_iter: Number of epochs.
        kwargs: Optional arguments to pass to the constructor of the model.
    
    Returns:
        Trained model.
    """
    
    idx_train = np.copy(idx_train)
    
    # Initialize RNG
    rng = np.random.default_rng(seed=seed)
    
    plant_seed(seed)
    
    # Initialize NN
    n_input_channels = X.size()[2]
    model = NN(X.shape[1], n_input_channels, **kwargs)
    if USE_GPU:
        model.cuda()
        device='cuda'
    else:
        device=None
    generator = torch.Generator(device=device).manual_seed(seed)
    # Initialize optimizer and scheduler
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3, eps=1e-7)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimizer, mode='min', factor=0.9, patience=5, verbose=False, threshold=0.0001,
            threshold_mode='rel', cooldown=2, min_lr=1e-5, eps=1e-08)
    
    # Place data on GPU
    if USE_GPU:
        X = X.cuda()
        Y = Y.cuda()
    
    pbar = tqdm.tqdm(range(n_iter))
    for epoch in pbar:
        total_loss, baseline_total_loss = 0, 0
        rng.shuffle(idx_train)
        for it_idx in np.array_split(idx_train, 200):  # 120
            
            # Reset gradients
            optimizer.zero_grad()
            
            # Create batch
            Y_target = Y[it_idx, :]
            if USE_GPU:
                Y_target = Y_target.cuda()
            
            # Sample from the null distribution and add background noise to the targets
            Y_target = Y_target + 0.5 * background_noise(*Y_target.size(), device=X.device, generator=generator)
            
            # Forward pass
            x = X[it_idx]
            if USE_GPU:
                x = x.cuda()
            Y_pred = model.forward(x)
            
            # Compute loss (mean absolute error) in a differentiable fashion
            loss = torch.mean(torch.abs(Y_target - Y_pred))
            
            # Backward pass
            loss.backward()
            
            # Update NN parameters
            optimizer.step()
            
            # Update training MRRMSE, and also update MRRMSE for the baseline (zero) predictor
            # for comparison purposes
            mrrmse = torch.sum(torch.sqrt(torch.mean(torch.square(Y_target - Y_pred), dim=1)))
            baseline_mrrmse = torch.sum(torch.sqrt(torch.mean(torch.square(Y_target), dim=1)))
            total_loss += mrrmse.item()
            baseline_total_loss += baseline_mrrmse.item()
            
        # Compute relative error. Relative error is < 1 when the model performs better than
        # the baseline on the training data.
        rel_error = total_loss / baseline_total_loss
        
        # Update learning rate
        scheduler.step(rel_error)
        
        # Show relative error
        pbar.set_description(f'{rel_error:.3f}')
    return model
