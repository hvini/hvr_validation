import torch
from torch import nn

# local stuff
from .util import load_model,import_parameters
from .util import PreprocessModule
#from .basic import InceptionModule,InceptionModule_v2

from .basic2 import InceptionModule_,ResidualModule
from torch.cuda.amp import autocast

class AttnNet(nn.Module):
    """
    Attention network for classifications/survival/regression
    """
    def __init__(self, 
                 layers_att, layers_out, 
                 in_chan=1, view_sz=32, 
                 norm_fn='batch', leaky=0.2, 
                 dropout=False,
                 in_context_chan=0,
                 use_softmax=True,
                 use_resnet=False,
                 group_norm_channels=8,
                 in_addon_chan=0,
                 regression_mode=False ):
        super(AttnNet, self).__init__()

        if isinstance(view_sz,list):
            self.view_sz = view_sz
        else:
            self.view_sz = [view_sz,view_sz,view_sz]
        self.context_chan = in_context_chan
        # add context channels
        in_layer_chan = in_chan + in_context_chan
        print("ATTN in_context_chan:",in_context_chan)

        self.mods_att = nn.ModuleList()
        #
        self.channels_att = []
        self.channels_out = []
        self.use_softmax = use_softmax

        self.levels = len( layers_att )

        # generate attention map
        for i in layers_att:
            self.channels_att.append(in_layer_chan)
            (_layers, out_layer_chan) = i

            if use_resnet:
                self.mods_att.append(ResidualModule(in_layer_chan, _layers, out_layer_chan,
                                    leaky=leaky, norm_fn='group',
                                    groups = group_norm_channels ))
            else:
                self.mods_att.append(InceptionModule_(in_layer_chan, _layers, out_layer_chan,
                                    leaky = leaky, norm_fn = 'group',
                                    groups = group_norm_channels ))

            # self.mods_att.append( InceptionModule_(in_layer_chan, inception_layers, out_layer_chan,
            #                 leaky=leaky, norm_fn=norm_fn))
            in_layer_chan = out_layer_chan


        self.att_channels_in = in_layer_chan
        # this layer will generate attention map (weight map)        
        self.attention_map = nn.Conv3d(in_layer_chan, 1, 1) # HACK

        # here we will have signal after passing through attention detection
        # the context will be stripped from layers before passing through the final classifier
        in_layer_chan = in_chan + in_addon_chan
        # out
        _out = []
        # TODO: figure out how to deal with addon channels and dropout
        if dropout:
            _out.append(nn.Dropout3d(0.5))

        for (ii,i) in enumerate(layers_out):
            # don't enable rectify 
            rectify=(not regression_mode) or ii<(len(layers_out)-1)
            self.channels_out.append(in_layer_chan)
            (_layers, out_layer_chan) = i
            if use_resnet:
                _out.append(ResidualModule(in_layer_chan, _layers, out_layer_chan,
                                    leaky=leaky, norm_fn='group',
                                    groups = group_norm_channels,
                                    rectify=rectify ))
            else:
                _out.append(InceptionModule_(in_layer_chan, _layers, out_layer_chan,
                                        leaky=leaky, norm_fn=norm_fn,
                                        rectify=rectify ))
            in_layer_chan = out_layer_chan

        self.mods_out = nn.Sequential(*_out)

    @autocast()
    def forward(self, x, addon_scalar=None):
        # pass throgh attention layers
        x_ = x
        for _, j in enumerate(self.mods_att):
            # pass the filter
            x_ = j(x_)
        # xd = x.size()
        # generate attention weights
        a = self.attention_map(x_).reshape(-1, self.att_channels_in*self.view_sz[0]*self.view_sz[1]*self.view_sz[2])

        # TODO: think about replacing this with nn.functional.scaled_dot_product_attention
        
        if self.use_softmax:
            # normalize , so that they sum to 1.0
            a = nn.functional.softmax(a, dim=1)
        else:
            a = torch.sigmoid(a)
        
        # this is again 3D attention map
        a = a.reshape(-1, 1 ,self.view_sz[0], self.view_sz[1], self.view_sz[2])

        # removing context channel(s) if present
        if self.context_chan>0:
            x=x[:,0:-self.context_chan,:,:,:]

        # apply attention
        x = x*a
        # now pool
        x = nn.functional.avg_pool3d(x, tuple(self.view_sz ) )

        #TODO: figure out what to do with dropout here
        if addon_scalar is not None:
            x = torch.cat( [x, addon_scalar.unsqueeze(-1).unsqueeze(-1).unsqueeze(-1) ], dim=1 )

        # and now finally apply output layers
        for _, j in enumerate( self.mods_out ):
            x = j(x)
        
        # here x supposed to be: (batch,n_cls,1,1,1)
        return x,a


def init_model(desc, view_sz=32, cpu=False, parallel=False):
    """
    Initialize a single model based on description
    :param desc: a dict
    :return: a model
    """
    _version = desc.get('version', None)
    _type    = desc.get('type', 'attn')

    _leaky   = desc.get('leaky', 0.1)
    _dropout = desc.get('dropout', None)

    _input_n    = desc.get('input')
    _output_n   = desc.get('output')

    _context_chan = desc.get('context_chan', 0)
    _softmax      = desc.get('softmax', True)
    _use_resnet = desc.get('use_resnet',False)

    _in_addon_chan= desc.get('addon_chan', 0)
    _regression_mode=desc.get('regression_mode', False)

    model = AttnNet( desc['attn'],
                     desc['out'],
                     leaky = _leaky,
                     dropout = _dropout,
                     in_chan = _input_n,
                     in_context_chan = _context_chan,
                     in_addon_chan = _in_addon_chan,
                     view_sz = view_sz,
                     use_resnet = _use_resnet,
                     use_softmax = _softmax,
                     regression_mode=_regression_mode )
    if not cpu:
        model.cuda()
    if 'load' in desc:
        load_model(model, desc['load'], map_location=next(model.parameters()).device)
    if 'import' in desc:
        import_parameters(model, desc['import'], map_location=next(model.parameters()).device)
    if parallel:
      #model = nn.parallel.DistributedDataParallel(model)
      pass # HACK - move parallel out of model initialization

    print("ATTN:",model)
    return model
