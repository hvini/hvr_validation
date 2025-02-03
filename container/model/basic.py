import torch
from torch import nn

# local stuff
from .util import load_model
from .util import PreprocessModule
from torch.cuda.amp import autocast

class InceptionModule(nn.Module):
    def __init__(self, in_chan, layers, out_chan, 
                 norm_fn = 'batch', instance_norm=None, 
                 leaky=0.2, groups=8, rectify=True):
        super(InceptionModule, self).__init__()
        
        if  norm_fn == 'instance' or instance_norm == True:
            _norm_fn = nn.InstanceNorm3d
        elif norm_fn == 'layer' :
            _norm_fn = nn.LayerNorm
        else:
            _norm_fn = nn.BatchNorm3d

        self.mods      = nn.ModuleList()
        
        out_module_channels = 0
        
        for i in layers:
            (kernel, m_out_chan)=i
            _pad=int((kernel-1)/2)
            conv_layers = []
            # TODO: check if _pad > 0
            if _pad>0:
                conv_layers.append(nn.ReplicationPad3d(_pad))
            conv_layers.append(nn.Conv3d(in_chan, m_out_chan, kernel))
            
            if rectify:
                if norm_fn == 'group' and (m_out_chan % groups) == 0:
                    conv_layers.append(nn.GroupNorm(groups, m_out_chan, affine=True))
                else:
                    conv_layers.append(_norm_fn(m_out_chan, affine=True))
                
                if leaky>0:
                    conv_layers.append(nn.LeakyReLU(leaky, inplace=True))
                else:
                    conv_layers.append(nn.ReLU(inplace=True))
            out_module_channels += m_out_chan

            self.mods.append(nn.Sequential(*conv_layers))

        if out_module_channels == out_chan and len(self.mods) == 1:
            self.merge = None
        else:
            merge_layers = []
            merge_layers.append(nn.Conv3d(out_module_channels, out_chan, 1))
            #merge_layers.append(_norm_fn(out_chan, affine=True))
            #merge_layers.append(nn.LeakyReLU(leaky, inplace=True))
            self.merge = nn.Sequential(*merge_layers)

    @autocast()
    def forward(self, x):
        v = []
        for i in self.mods:
            v.append( i.forward(x) ) 

        x = torch.cat(v, 1) # merge feature dimension
        
        if self.merge is not None:
            x = self.merge.forward(x)
        
        return x



class InceptionModule_v2(nn.Module):
    def __init__(self, in_chan, layers, out_chan, norm = True, leaky=0.2, groups=8):
        super(InceptionModule_v2, self).__init__()
        
        self.mods      = nn.ModuleList()
        
        out_module_channels = 0
        
        for i in layers:
            (kernel, m_out_chan)=i
            _pad=int((kernel-1)/2)
            conv_layers = []
            if norm:
                if in_chan % groups == 0:
                    conv_layers.append(nn.GroupNorm(groups, in_chan, affine=True))
                else:
                    conv_layers.append(nn.BatchNorm3d(in_chan, affine=True))
            
            if leaky>0.0:
                conv_layers.append(nn.LeakyReLU(leaky, inplace=True))
            else:
                conv_layers.append(nn.ReLU(inplace=True))

            if _pad>0:
                conv_layers.append(nn.ReplicationPad3d(_pad))
            
            conv_layers.append(nn.Conv3d(in_chan, m_out_chan, kernel))
            out_module_channels += m_out_chan

            self.mods.append(nn.Sequential(*conv_layers))

        if out_module_channels == out_chan and len(self.mods) == 1:
            self.merge = None
        else:
            merge_layers = []
            if norm:
                if in_chan % groups==0:
                    conv_layers.append(nn.GroupNorm(groups, in_chan, affine=True))
                else:
                    conv_layers.append(nn.BatchNorm3d(in_chan, affine=True))
            if leaky>0.0:
                conv_layers.append(nn.LeakyReLU(leaky, inplace=True))
            else:
                conv_layers.append(nn.ReLU(inplace=True))
            
            merge_layers.append(nn.Conv3d(out_module_channels, out_chan, 1))

            self.merge = nn.Sequential(*merge_layers)

    @autocast()
    def forward(self, x):
        v = []
        for i in self.mods:
            v.append( i.forward(x) ) 

        x = torch.cat(v, 1) # merge feature dimension
        
        if self.merge is not None:
            x = self.merge.forward(x)
        
        return x
