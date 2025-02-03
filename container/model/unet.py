import numpy as np
import torch
from torch import nn
from torch.autograd import Variable
from torch.nn import functional as F
import math

class InceptionModule(nn.Module):
    def __init__(self, in_chan, layers, out_chan, params):
        super(InceptionModule, self).__init__()
        
        norm_fn        = nn.InstanceNorm3d if params.instance_norm else nn.BatchNorm3d
        self.mods      = nn.ModuleList()
        
        out_module_channels=0
        
        for i in layers:
            (kernel,m_out_chan)=i
            conv_layers=[]
            _pad=int((kernel-1)/2)
            
            conv_layers=[]
            conv_layers.append(nn.ReplicationPad3d(_pad))
            conv_layers.append(nn.Conv3d(in_chan, m_out_chan, kernel))
            conv_layers.append(norm_fn(m_out_chan, affine=True))
            conv_layers.append(nn.LeakyReLU(params.leaky, inplace=True))
            out_module_channels+=m_out_chan
            
            self.mods.append(nn.Sequential(*conv_layers))
            
        if out_module_channels==out_chan and len(self.mods)==1:
            self.merge=None
        else:
            merge_layers=[]
            merge_layers.append(nn.Conv3d(out_module_channels, out_chan, 1))
            merge_layers.append(norm_fn(out_chan, affine=True))
            merge_layers.append(nn.LeakyReLU(params.leaky, inplace=True))
            self.merge=nn.Sequential(*merge_layers)
            

    def forward(self, x):
        
        v=[]
        for i in self.mods:
            v.append(i.forward(x))
        x=torch.cat(v,1) # merge feature dimension
        
        if self.merge is not None:
            x=self.merge.forward(x)
        
        return x

class InceptionNet(nn.Module):
    def __init__(self, params, layers):
        super(InceptionNet, self).__init__()

        self.patch_sz  = params.patch_sz
        self.n_cls     = params.n_cls
        
        norm_fn        = nn.InstanceNorm3d if params.instance_norm else nn.BatchNorm3d
        in_chan=1
        
        inception_modules=[]
        
        for i in layers:
            (inception_layers,out_chan)=i
            inception_modules.append(InceptionModule(in_chan, inception_layers, out_chan, params))
            in_chan=out_chan
            
        self.inc_modules=nn.Sequential(*inception_modules)
        
        # not trainable parameters
        self.sample_mean=0.0
        self.sample_sd=1.0
        
        # additial data augmentation
        self.intensity_variance=None
        self.noise=None
        
        
    def init_prefilter(self, sampler ):
        # set sample_mean and sample_sd from the sampler
        self.sample_mean=sampler.sample_mean
        self.sample_sd=sampler.sample_sd

    def forward(self, x):
        assert x.size()[1:] == (1, self.patch_sz, self.patch_sz, self.patch_sz)
        
        if self.training:
            if self.intensity_variance is not None:
                x=x.mul(torch.normal(torch.FloatTensor([1.0]),std=torch.FloatTensor([self.intensity_variance])))
            if self.noise is not None:
                x=x.add(x.clone().normal_(0.0,self.noise))
        
        # apply pre-normalization step
        x=x.add(-self.sample_mean).mul(1.0/self.sample_sd)
        
        # 1. pass through the regular convolution stages
        out = self.inc_modules(x)
        
        # 1.5 TODO: add cropping ?
        # 1.6 TODO: add dropout ?
        # 2. transpose, to move the classes to the last dimensions
        out = out.transpose(1,2)
        out = out.transpose(2,3)
        out = out.transpose(3,4)
        # make contiguous (transpose is a culprit)
        out = out.contiguous()
        # reformat for nnls
        out = out.view(-1, self.n_cls)
        # finally apply log-softmax along the feature dimension
        return F.log_softmax(out, dim=1)

class InceptionNetAE(nn.Module):
    def __init__(self, params, layers):
        super(InceptionNetAE, self).__init__()

        self.patch_sz  = params.patch_sz
        
        norm_fn        = nn.InstanceNorm3d if params.instance_norm else nn.BatchNorm3d
        in_chan=1
        
        inception_modules=[]
        
        for i in layers:
            (inception_layers,out_chan)=i
            inception_modules.append(InceptionModule(in_chan, inception_layers, out_chan, params))
            in_chan=out_chan
            
        self.inc_modules=nn.Sequential(*inception_modules)
        
        # not trainable parameters
        self.sample_mean=0.0
        self.sample_sd=1.0
        
        # additial data augmentation
        self.intensity_variance=None
        self.noise=None
        
        
    def init_prefilter(self, sampler ):
        # set sample_mean and sample_sd from the sampler
        self.sample_mean=sampler.sample_mean
        self.sample_sd=sampler.sample_sd

    def forward(self, x):
        assert x.size()[1:] == (1, self.patch_sz, self.patch_sz, self.patch_sz)
        
        if self.training:
            if self.intensity_variance is not None:
                x=x.mul(torch.normal(torch.FloatTensor([1.0]),std=torch.FloatTensor([self.intensity_variance])))
            if self.noise is not None:
                x=x.add(x.clone().normal_(0.0,self.noise))
        
        # apply pre-normalization step
        x=x.add(-self.sample_mean).mul(1.0/self.sample_sd)
        
        # 1. pass through the regular convolution stages
        out = self.inc_modules(x)
        
        # 1.5 TODO: add cropping ?
        # 1.6 TODO: add dropout ?
        # 2. transpose, to move the classes to the last dimensions
        out = out.transpose(1,2)
        out = out.transpose(2,3)
        out = out.transpose(3,4)
        # make contiguous (transpose is a culprit)
        out = out.contiguous()
        # reformat for nnls
        out = out.view(-1)
        # finally apply log-softmax along the feature dimension
        return out


class UNet(nn.Module):
    def __init__(self, params, layers_down, layers_up, layers_add, in_chan=1 ):
        super(UNet, self).__init__()
        # make sure we have equal number of layers
        assert len(layers_down) == len(layers_up)
        
        self.patch_sz  = params.patch_sz
        self.n_cls     = params.n_cls
        
        norm_fn        = nn.InstanceNorm3d if params.instance_norm else nn.BatchNorm3d
        downsample     = nn.AvgPool3d
        upsample       = nn.Upsample # TODO: replace with deconv?
        in_layer_chan  = in_chan
        
        self.mods_down  = nn.ModuleList()
        self.mods_up    = nn.ModuleList()
        self.mods_add   = nn.ModuleList()
        self.mods_ds    = nn.ModuleList()
        self.mods_us    = nn.ModuleList()
        #
        self.channels_down = []
        self.channels_up   = []
        self.channels_add  = []
        
        self.levels        = len(layers_down)
        
        
        # downstream
        for i in layers_down:
            self.channels_down.append(in_layer_chan)
            (inception_layers, out_layer_chan)=i
            self.mods_down.append(InceptionModule(in_layer_chan, inception_layers, out_layer_chan, params))
            in_layer_chan = out_layer_chan
            self.mods_ds.append(downsample(2))
        
        # upstream
        for j,i in enumerate(layers_up):
            (inception_layers, out_layer_chan)=i
            self.mods_up.append(InceptionModule(in_layer_chan, inception_layers, out_layer_chan, params))
            self.mods_us.append(upsample(scale_factor=2))
            # concat not downsampled channels here
            in_layer_chan = out_layer_chan + self.channels_down[-j-1]
            self.channels_up.append(in_layer_chan)

        # addon
        for i in layers_add:
            self.channels_add.append(in_layer_chan)
            (inception_layers, out_layer_chan)=i
            self.mods_add.append(InceptionModule(in_layer_chan, inception_layers, out_layer_chan, params))
            in_layer_chan = out_layer_chan

        # not trainable parameters
        self.sample_mean=0.0
        self.sample_sd=1.0
        
        
    def init_prefilter(self, sampler):
        # set sample_mean and sample_sd from the sampler
        self.sample_mean=sampler.sample_mean
        self.sample_sd=sampler.sample_sd
    
    def forward(self, x):
        
        # technically we don't need this
        #assert x.size()[1:] == (1, self.patch_sz, self.patch_sz, self.patch_sz)
        
        # apply pre-normalization step
        #x=x.add(-self.sample_mean).mul(1.0/self.sample_sd)
        layer_in=[]
        #pass throgh the U-net down
        for i,j in enumerate(self.mods_down):
            layer_in.append(x)
            # pass the filter
            x=j(x)
            # pass downsampler
            x=self.mods_ds[i](x)

        #pass throgh the U-net up
        for i,j in enumerate(self.mods_up):
            
            x=j(x)
            # pass through the module
            x=self.mods_us[i](x)
            # concatenate with input 
            x=torch.cat( [x, layer_in[-i-1] ], 1)
        
        #pass throgh addon modules
        for i,j in enumerate(self.mods_add):
            x=self.mods_add[i](x)
            
        # 1.5 TODO: add cropping ?
        # 1.6 TODO: add dropout ?
        # 2. transpose, to move the classes to the last dimensions
        x=x.transpose(1,2)
        x=x.transpose(2,3)
        x=x.transpose(3,4)
        # make contiguous (transpose is a culprit)
        x=x.contiguous()
        # reformat for nnls
        x=x.view(-1, self.n_cls)
        # finally apply log-softmax along the feature dimension
        return F.log_softmax(x, dim=1)
    
class UNetAE(nn.Module):
    def __init__(self, params, layers_down, layers_up, layers_add, in_chan=1 ):
        super(UNetAE, self).__init__()
        # make sure we have equal number of layers
        assert len(layers_down) == len(layers_up)
        
        self.patch_sz  = params.patch_sz
        
        norm_fn        = nn.InstanceNorm3d if params.instance_norm else nn.BatchNorm3d
        downsample     = nn.AvgPool3d
        upsample       = nn.Upsample # TODO: replace with deconv?
        in_layer_chan  = in_chan
        
        self.mods_down  = nn.ModuleList()
        self.mods_up    = nn.ModuleList()
        self.mods_add   = nn.ModuleList()
        self.mods_ds    = nn.ModuleList()
        self.mods_us    = nn.ModuleList()
        #
        self.channels_down = []
        self.channels_up   = []
        self.channels_add  = []
        
        self.levels        = len(layers_down)
        
        
        # downstream
        for i in layers_down:
            self.channels_down.append(in_layer_chan)
            (inception_layers, out_layer_chan)=i
            self.mods_down.append(InceptionModule(in_layer_chan, inception_layers, out_layer_chan, params))
            in_layer_chan = out_layer_chan
            self.mods_ds.append(downsample(2))
        
        # upstream
        for j,i in enumerate(layers_up):
            (inception_layers, out_layer_chan)=i
            self.mods_up.append(InceptionModule(in_layer_chan, inception_layers, out_layer_chan, params))
            self.mods_us.append(upsample(scale_factor=2))
            # concat not downsampled channels here
            in_layer_chan = out_layer_chan + self.channels_down[-j-1]
            self.channels_up.append(in_layer_chan)

        # addon
        for i in layers_add:
            self.channels_add.append(in_layer_chan)
            (inception_layers, out_layer_chan)=i
            self.mods_add.append(InceptionModule(in_layer_chan, inception_layers, out_layer_chan, params))
            in_layer_chan = out_layer_chan

        # not trainable parameters
        self.sample_mean=0.0
        self.sample_sd=1.0
        
        
    def init_prefilter(self, sampler):
        # set sample_mean and sample_sd from the sampler
        self.sample_mean=sampler.sample_mean
        self.sample_sd=sampler.sample_sd
    
    def forward(self, x):
        
        # technically we don't need this
        #assert x.size()[1:] == (1, self.patch_sz, self.patch_sz, self.patch_sz)
        
        # apply pre-normalization step
        #x=x.add(-self.sample_mean).mul(1.0/self.sample_sd)
        layer_in=[]
        #pass throgh the U-net down
        for i,j in enumerate(self.mods_down):
            layer_in.append(x)
            # pass the filter
            x=j(x)
            # pass downsampler
            x=self.mods_ds[i](x)

        #pass throgh the U-net up
        for i,j in enumerate(self.mods_up):
            
            x=j(x)
            # pass through the module
            x=self.mods_us[i](x)
            # concatenate with input 
            x=torch.cat( [x, layer_in[-i-1] ], 1)
        
        #pass throgh addon modules
        for i,j in enumerate(self.mods_add):
            x=self.mods_add[i](x)
            
        # 1.5 TODO: add cropping ?
        # 1.6 TODO: add dropout ?
        # 2. transpose, to move the classes to the last dimensions
        x=x.transpose(1,2)
        x=x.transpose(2,3)
        x=x.transpose(3,4)
        # make contiguous (transpose is a culprit)
        x=x.contiguous()
        x=x.view(-1)
        return x



class CriticAE(nn.Module):
    def __init__(self, params, layers):
        super(CriticAE, self).__init__()

       
        norm_fn        = nn.InstanceNorm3d if params.instance_norm else nn.BatchNorm3d
        in_chan=1
        
        inception_modules=[]
        
        for i in layers:
            (inception_layers,out_chan)=i
            inception_modules.append(InceptionModule(in_chan, inception_layers, out_chan, params))
            in_chan=out_chan
            
        self.inc_modules=nn.Sequential(*inception_modules)

        self.project=nn.Sequential(
            nn.Conv3d(out_chan, 1, 1),
            nn.Sigmoid()
        )
        
    def forward(self, x):
        
        # 1. pass through the regular convolution stages
        out = self.inc_modules(x)
        out = self.project(out)

        return out.view(x.size(0),-1).mean(1).view(x.size(0))


