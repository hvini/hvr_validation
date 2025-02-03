import numpy as np
import torch
from torch import nn
from torch.autograd import Variable
from torch.nn import functional as F
import math


class SimpleFCN(nn.Module):

    def __init__(self, params, layers):
        super(SimpleFCN, self).__init__()

        self.patch_sz  = params.patch_sz
        self.n_cls     = params.n_cls
        norm_fn        = nn.InstanceNorm3d if params.instance_norm else nn.BatchNorm3d
        
        conv_layers=[]
        in_chan=1
        for i,j in enumerate(layers):
          (kernel,out_chan)=j
          _pad=int((kernel-1)/2)
          
          conv_layers.append(nn.ReplicationPad3d(_pad))
          conv_layers.append(nn.Conv3d(in_chan, out_chan, kernel))
          
          if i > 0:
              conv_layers.append(norm_fn(out_chan, affine=True))
            
          conv_layers.append(nn.LeakyReLU(0.2, inplace=True))
          in_chan=out_chan
          
        # last layer is going to be linear (?)
        self.conv_layers = nn.Sequential(*conv_layers)
        
    def forward(self, x):
        assert x.size()[1:] == (1, self.patch_sz, self.patch_sz, self.patch_sz)
        
        # 1. pass through the regular convolution stages
        out = self.conv_layers(x)
        # 1.5 TODO: add cropping ?
        # 1.6 TODO: add dropout ?
        # 2. transpose, to move the classes to the last dimensions
        out=out.transpose(1,2)
        out=out.transpose(2,3)
        out=out.transpose(3,4)
        # make contiguous (transpose is a culprit)
        out=out.contiguous()
        # reformat for nnls
        out=out.view(-1, self.n_cls)
        # finally apply log-softmax along the feature dimension
        return F.log_softmax(out, dim=1)
