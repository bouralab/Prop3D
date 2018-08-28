import torch
import skorch
from skorch.callbacks import LRScheduler

from molmimic.torch_model.torch_model import UNet3D, ResNetUNet
from molmimic.torch_model.torch_loader import IBISDataset, sparse_collate

p = 0.01 #Update with true values
true_distribution = torch.FloatTensor([1-p,p])

class 3DUnet(skorch.NeuralNet):
    number_of_features = 70
    number_of_classes = 1

    def __init__(self,
      max_epochs=100,
      lr=0.0001,
      optimizer=torch.optim.SGD,
      optimizer__momentum = 0.999,
      optimizer__weight_decay = 1e-4,
      optimizer__nesterov = True,
      module__dropout_depth = True,
      module__dropout_width = True,
      module__dropout_p = 0.5,
      **kwds):
        args = [ResNetUNet, torch.nn.CrossEntropyLoss]
        kwds["module__nInputFeatures"]    = self.number_of_features
        kwds["module__nClasses"]          = self.number_of_classes
        kwds["criterion__weight"]         = true_distribution
        kwds["iterator_train"]            = MMDBDataset
        kwds["iterator_train_collate_fn"] = sparse_collate
        kwds["iterator_valid"]            = MMBDDataset
        kwds["iterator_valid_collate_fn"] = sparse_collate

        scheduler = LRScheduler(policy="LambdaLR", lr_lambda=lambda epoch: math.exp((1 - epoch) * lr_decay))
        if "callbacks" in kwds and isinstance(kwds["callbacks"], list):
            callbacks.append(scheduler)
        else:
            kwds["callbacks"] = [scheduler]

        super().__init__(*args, **kwargs)

class ProteinAutoEncoder(skorch.NeuralNet):
    number_of_features = 70
    number_of_classes = 70

class BindingSitePredictor(skorch.NeuralNet):
    number_of_features = 70
    number_of_classes = 2

def Molmimic(dataset_name, cdd):
    pae = ProteinAutoEncoder()
    pae.fit()

    bsp = BindingSitePredictor()

    pae_weights = pae.module.state_dict()
    bsp_weights = bsp.module.state_dict()
    pae_weights["linear.weight"] = bsp_weights["linear.weight"]
    pae_weights["linear.bias"] = bsp_weights["linear.bias"]
    bsp.module.load_state_dict(weights)

    bsp.fit()
