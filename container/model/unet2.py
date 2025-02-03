import torch
from torch import nn

# local stuff
from .util import load_model
from .util import PreprocessModule
from .basic import InceptionModule

from torch.cuda.amp import autocast

class UNet(nn.Module):
    """
    Unet implimentation without transposition
    for GANN
    """
    def __init__(self, layers_down, layers_up, layers_add, in_chan=1, patch_sz=32, n_cls=2,
                 norm_fn='batch',
                 instance_norm=False, leaky=0.2, dropout=False, upsample_mode='nearest'):
        super(UNet, self).__init__()
        # make sure we have equal number of layers
        assert len(layers_down) == len(layers_up)
        
        self.patch_sz  = patch_sz
        self.n_cls     = n_cls
        
        downsample     = nn.AvgPool3d
        in_layer_chan  = in_chan
        
        self.mods_down  = nn.ModuleList()
        self.mods_up    = nn.ModuleList()
        self.upsample_mode = upsample_mode
        self.mods_ds    = nn.ModuleList()
        #self.mods_us    = list() # nn.ModuleList()
        #
        self.channels_down = []
        self.channels_up   = []
        self.channels_add  = []
        
        self.levels        = len(layers_down)
        
        
        # downstream
        for i in layers_down:
            self.channels_down.append(in_layer_chan)
            (inception_layers, out_layer_chan) = i
            self.mods_down.append(InceptionModule(in_layer_chan, inception_layers, out_layer_chan,
                                                  leaky=leaky,norm_fn=norm_fn,
                                                  instance_norm=instance_norm))
            in_layer_chan = out_layer_chan
            self.mods_ds.append(downsample(2))
        
        # upstream
        for j,i in enumerate(layers_up):
            (inception_layers, out_layer_chan) = i
            self.mods_up.append(InceptionModule(in_layer_chan, inception_layers, out_layer_chan,
                                                leaky=leaky,norm_fn=norm_fn,
                                                instance_norm=instance_norm))
            #self.mods_us.append(upsample(2)) 
            #self.mods_us.append(2)
            # concat not downsampled channels here
            in_layer_chan = out_layer_chan + self.channels_down[-j-1]
            self.channels_up.append(in_layer_chan)

        # addon
        _add = []

        if dropout:
            _add.append(nn.Dropout3d(0.5))

        for i in layers_add:
            self.channels_add.append(in_layer_chan)
            (inception_layers, out_layer_chan) = i
            _add.append(InceptionModule(in_layer_chan, inception_layers, out_layer_chan,
                                        leaky=leaky,norm_fn=norm_fn,
                                        instance_norm=instance_norm))
            in_layer_chan = out_layer_chan

        self.mods_add = nn.Sequential(*_add)
    
    @autocast()
    def forward(self, x):
        layer_in = []
        #pass throgh the U-net down
        for i, j in enumerate(self.mods_down):
            layer_in.append(x)
            # pass the filter
            x = j(x)
            # pass downsampler
            x = self.mods_ds[i](x)

        # pass through the U-net up,
        # concatenating skip connections
        for i, j in enumerate(self.mods_up):
            x = j(x)
            # pass through the module
            #x = self.mods_us[i](x)
            x = nn.functional.interpolate(x, scale_factor=2, mode=self.upsample_mode)
            # concatenate with input 
            x = torch.cat([x, layer_in[-i-1]], 1)
        
        #pass throgh addon stages
        x = self.mods_add(x)

        return x


class LNet(nn.Module):
    """
    Half of Unet for classifications
    """

    def __init__(self, layers_down, layers_add, in_chan=1, patch_sz=32, n_cls=2,
                 instance_norm=False, leaky=0.2, dropout=False):
        super(LNet, self).__init__()
        # make sure we have equal number of layers

        self.patch_sz = patch_sz
        self.n_cls = n_cls

        norm_fn = nn.InstanceNorm3d if instance_norm else nn.BatchNorm3d
        downsample = nn.AvgPool3d
        in_layer_chan = in_chan

        self.mods_down = nn.ModuleList()
        self.mods_ds = nn.ModuleList()
        #
        self.channels_down = []
        self.channels_add = []

        self.levels = len(layers_down)

        # downstream
        for i in layers_down:
            self.channels_down.append(in_layer_chan)
            (inception_layers, out_layer_chan) = i
            self.mods_down.append(InceptionModule(in_layer_chan, inception_layers, out_layer_chan,
                                                  leaky=leaky,
                                                  instance_norm=instance_norm))
            in_layer_chan = out_layer_chan
            self.mods_ds.append(downsample(2))

        # addon
        _add = []

        if dropout:
            _add.append(nn.Dropout3d(0.5))

        for i in layers_add:
            self.channels_add.append(in_layer_chan)
            (inception_layers, out_layer_chan) = i
            _add.append(InceptionModule(in_layer_chan, inception_layers, out_layer_chan,
                                        leaky=leaky,
                                        instance_norm=instance_norm))
            in_layer_chan = out_layer_chan

        self.mods_add = nn.Sequential(*_add)
        # run log soft max across all features
        #self.project = nn.LogSoftmax(1)

    def forward(self, x):
        layer_in = []
        # pass throgh the Left U-net down
        for i, j in enumerate(self.mods_down):
            layer_in.append(x)
            # pass the filter
            x = j(x)
            # pass downsampler
            x = self.mods_ds[i](x)

        # pass throgh addon stages
        x = self.mods_add(x)
        #x = self.project(x)

        # now we have to average across all voxels in the patch
        # TODO: figure out if it's better to average first
        #if avg:
        #    x = x.view(x.size(0),x.size(1), -1).mean(2).view(x.size(0), x.size(1))

        return x


class CriticSEG(nn.Module):
    """
    Segmentation Critic
    input: N+1 channels (scan+1-hot encoded output of segmentation)
    output: yes/no label
    """
    def __init__(self, layers, in_chan=3, instance_norm=False, leaky=0.2, dropout=False):
        super(CriticSEG, self).__init__()
       
        downsample     = nn.AvgPool3d
        
        inception_modules = []
        for i in layers:
            inception_layers, out_chan = i
            inception_modules.append(InceptionModule(in_chan, inception_layers, out_chan,
                                                     leaky=leaky, instance_norm=instance_norm))
            # figure out if downsampling improves things
            inception_modules.append(downsample(2))
            in_chan = out_chan

        if dropout:
            inception_modules.append(nn.Dropout3d(0.5))
        
        self.inc_modules = nn.Sequential(*inception_modules)

        self.project = nn.Sequential(
            nn.Conv3d(out_chan, 1, 1),
            nn.Sigmoid())

    def forward(self, x):
        # 1. pass through the regular convolution stages
        out = self.inc_modules(x)
        out = self.project(out)

        return out.view(x.size(0), -1).mean(1).view(x.size(0))


def init_model(desc,parallel=False):
    """
    Initialize a single model based on description
    :param desc: a dict
    :param n_cls: number of classes
    :return: a model
    """
    _version = desc.get('version', None)
    _type = desc.get('type', 'unet')

    _leaky = desc.get('leaky', 0.2)
    _instance_norm = desc.get('instance_norm', False)
    _dropout = desc.get('dropout', False)
    _upsample_mode = desc.get('upsample_mode','nearest')

    _input_n = desc.get('input')
    _output_n = desc.get('output')

    print(repr(desc))

    # legacy cases
    if  _version== '20180621' and _type == 'unet':
        # four stage U-net (4 downscaling -> 4 upscaling )
        model = UNet( [[[[5, 4], [3, 4], [1, 4]], 12],  # downstream
                       [[[5, 8], [3, 8], [1, 8]], 24],
                       [[[5, 8], [3, 8], [1, 8]], 24],
                       [[[1, 32]], 32],
                      ],
                     [[[[5, 16], [3, 16], [1, 16]], 48],  # upstream
                      [[[3, 16], [1, 16]], 64],
                      [[[3, 16], [1, 16]], 64],
                      [[[1, 32]], 64],
                     ],
                     [[[[3, 4], [1, 4]], 48],  # addon
                      [[[1, _output_n]], _output_n]],
                      leaky=_leaky, instance_norm=_instance_norm, dropout=_dropout, in_chan=_input_n)

    elif _version == '20180615' and _type == 'unet':
        model = UNet([
                    [[[5, 4], [3, 4], [1, 4]], 12],  # downstream
                    [[[5, 8], [3, 8], [1, 8]], 24],
                    [[[5, 8], [3, 8], [1, 8]], 24],
                    [[[1, 32]], 32],
                    ],
                     [[[[5, 16], [3, 16], [1, 16]], 48],  # upstream
                      [[[3, 16], [1, 16]], 32],
                      [[[3, 16], [1, 16]], 32],
                      [[[1, 32]], 32],
                      ],
                     [[[[3, 4], [1, 4]], 48],  # addon
                      [[[1, _output_n]], _output_n]],
                    leaky=_leaky, instance_norm=_instance_norm, dropout=_dropout, in_chan=_input_n)

    elif _version == '20180621' and _type == 'critic':
        model = CriticSEG([[[[5, 4], [3, 4], [1, 4]], 12],
                                  [[[5, 8], [3, 8], [1, 8]], 24],
                                  [[[1, 32]], 32]], in_chan=1 + n_cls,
                          leaky=_leaky, instance_norm=_instance_norm, dropout=_dropout)

    elif _version == '20180615' and _type == 'critic':
        model = CriticSEG([[[[5, 4], [3, 4], [1, 4]], 12],
                                          [[[5, 8], [3, 8], [1, 8]], 24],
                                          [[[1, 32]], 32]], in_chan=1 + n_cls,
                          leaky=_leaky, instance_norm=_instance_norm, dropout=_dropout)

    elif _version == '20180706' and _type == 'unet':
        model = UNet([
            [[[5, 4], [3, 4], [1, 4]], 12],  # downstream
            [[[5, 8], [3, 8], [1, 8]], 24],
            [[[5, 8], [3, 8], [1, 8]], 24],
            [[[1, 32]], 32],
                  ],
                     [[[[5, 16], [3, 16], [1, 16]], 48],  # upstream
                      [[[3, 16], [1, 16]], 32],
                      [[[3, 16], [1, 16]], 32],
                      [[[1, 32]], 32],
                      ],
                     [[[[3, 4], [1, 4]], 48],  # addon
                      [[[1, _output_n]], _output_n]],
            leaky=_leaky, instance_norm=_instance_norm, dropout=_dropout, in_chan=_input_n)

    elif _version == '20180706' and _type == 'critic':
        # TODO : pass this outside somehow (?)
        feed_only_output = True
        # only analyze segmentations
        model = CriticSEG( [[[[5, 4], [3, 4], [1, 4]], 12],
                            [[[5, 8], [3, 8], [1, 8]], 24],
                            [[[1, 32]], 32]], in_chan=_input_n,
                           leaky=_leaky, instance_norm=_instance_norm, dropout=_dropout)
    # generic versions
    elif _type == 'unet':
        model = UNet(desc.get('downstream', [
                        [[[5, 4], [3, 4], [1, 4]], 12],  # downstream
                        [[[5, 8], [3, 8], [1, 8]], 24],
                        [[[5, 8], [3, 8], [1, 8]], 24],
                        [[[1, 32]], 32],
                      ]),
                     desc.get('upstream',[
                         [[[5, 16], [3, 16], [1, 16]], 48],  # upstream
                         [[[3, 16], [1, 16]], 32],
                         [[[3, 16], [1, 16]], 32],
                         [[[1, 32]], 32],
                      ]),
                     desc.get('addon',[[[[3, 4], [1, 4]], 48],  # addon
                                       [[[1, _output_n]], _output_n]]),
                     leaky=_leaky, instance_norm=_instance_norm, dropout=_dropout, in_chan=_input_n,upsample_mode=_upsample_mode)
    # TODO: add inception net type
    elif _type == 'critic':
        model = CriticSEG(desc.get('downstream', [[[[5, 4], [3, 4], [1, 4]], 12],
                            [[[5, 8], [3, 8], [1, 8]], 24],
                            [[[1, 32]], 32]]), in_chan=_output_n,
                          leaky=_leaky, instance_norm=_instance_norm, dropout=_dropout)
    elif _type == 'cls' or _type == 'lnet':
        # TODO: add l-net style classifier
        model = LNet(desc.get('downstream', [
                        [[[5, 4], [3, 4], [1, 4]], 12],  # downstream
                        [[[5, 8], [3, 8], [1, 8]], 24],
                        [[[5, 8], [3, 8], [1, 8]], 24],
                        [[[1, 32]], 32],
                      ]),
                     desc.get('addon', [[[[3, 4], [1, 4]], 48],  # addon
                                        [[[1, _output_n]], _output_n]]),
                     leaky=_leaky, instance_norm=_instance_norm, dropout=_dropout, in_chan=_input_n)
    else:
        # this should never happen
        model = None

    model.cuda()

    if 'load' in desc:
        load_model(model, desc['load'])
    
    if parallel:
      model=nn.DataParallel(model)

    return model


def init_models(desc, n_cls=2):
    """
    Initialize training and critic model (for GANN)
    TODO: add task-specific models
    :param desc: dictionary with entry "model" and "critic_model"
    :param n_cls: number of classes
    :return: tuple of two models
    """
    model = init_model(desc['model'], n_cls)
    if 'model_critic' in desc:
        model_critic = init_model(desc['model_critic'], n_cls)
    else:
        model_critic = None
    return model, model_critic
