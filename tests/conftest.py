import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
import torch

@pytest.fixture
def device():
    return torch.device('cpu')

@pytest.fixture
def batch_size():
    return 2
