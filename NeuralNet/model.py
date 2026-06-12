import torch
import torch.nn as nn
import torch.nn.functional as F

class BinaryClassifier(nn.Module):

    def __init__(self, input_size):
        super().__init__()

        self.network = nn.Sequential(
                nn.Linear(input_size, 128),
                nn.GELU(),
                nn.Linear(128, 64),
                nn.GELU(),
                nn.Linear(64, 1)
        )

    def forward(self, x):
        return self.network(x)
