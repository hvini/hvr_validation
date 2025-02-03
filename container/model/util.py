import os
import re
import pickle
import random
import torch
import datetime
import inspect
import math
import copy

import torch
from torch import nn
# import torchvision

import csv
import json
from torch import optim
from torch.nn.parallel import DistributedDataParallel as DDP

# scheduler
from torch.optim       import lr_scheduler

from torch.cuda.amp import autocast

class PreprocessModule(nn.Module):
    """
    Simple preprocessing module
    will subtract mean and divide by sd
    """
    def __init__(self, sample_mean, sample_sd):
        """
        Arguments:
            sample_mean - torch tensor broadcastable to the number of input channels
            sample_sd - torch tensor broadcastable to the number of input channels
        """
        super(PreprocessModule, self).__init__()
        sample_mean=torch.nn.Parameter(sample_mean, requires_grad=False)
        sample_sd=torch.nn.Parameter(sample_sd, requires_grad=False)
        self.register_parameter("sample_mean",sample_mean)
        self.register_parameter("sample_sd",sample_sd)

    def forward(self, x):
        x=x.sub(self.sample_mean)
        x=x.div(self.sample_sd)
        return x


# based on code from facebookresearch's FaderNetwork: https://github.com/facebookresearch/FaderNetworks
def get_optimizer(model, s):
    """
    Parse optimizer parameters.
    Input should be of the form:
        - "sgd,lr=0.01"
        - "adagrad,lr=0.1,lr_decay=0.05"
    """
    if "," in s:
        method = s[:s.find(',')]
        optim_params = {}
        for x in s[s.find(',') + 1:].split(','):
            split = x.split('=')
            assert len(split) == 2
            assert re.match( r"^[+-]?(\d+(\.\d*)?|\.\d+)$", split[1]) is not None
            optim_params[split[0]] = float(split[1])
    else:
        method = s
        optim_params = {}

    if method == 'adadelta':
        optim_fn = optim.Adadelta
    elif method == 'adagrad':
        optim_fn = optim.Adagrad
    elif method == 'adam':
        optim_fn = optim.Adam
        optim_params['betas'] = (optim_params.get('beta1', 0.5), optim_params.get('beta2', 0.999))
        optim_params.pop('beta1', None)
        optim_params.pop('beta2', None)
    elif method == 'adamw':
        optim_fn = optim.AdamW
        optim_params['betas'] = (optim_params.get('beta1', 0.5), optim_params.get('beta2', 0.999))
        optim_params.pop('beta1', None)
        optim_params.pop('beta2', None)
    elif method == 'adamax':
        optim_fn = optim.Adamax
    elif method == 'asgd':
        optim_fn = optim.ASGD
    elif method == 'rmsprop':
        optim_fn = optim.RMSprop
    elif method == 'rprop':
        optim_fn = optim.Rprop
    elif method == 'sgd':
        optim_fn = optim.SGD
        assert 'lr' in optim_params
    else:
        raise Exception('Unknown optimization method: "%s"' % method)

    # check that we give good parameters to the optimizer
    # TODO: figure out if we need to use kwonlyargs also
    expected_args = inspect.getfullargspec(optim_fn.__init__)[0]
    assert expected_args[:2] == ['self', 'params']
    if not all(k in expected_args[2:] for k in optim_params.keys()):
        raise Exception('Unexpected parameters: expected "%s", got "%s"' % (
            str(expected_args[2:]), str(optim_params.keys())))

    return optim_fn(model.parameters(), **optim_params)


# depricated in favor of torch.nn.utils.clip_grad_norm_
# def clip_grad_norm(parameters, max_norm, norm_type=2):
#     """Clips gradient norm of an iterable of parameters.
#     The norm is computed over all gradients together, as if they were
#     concatenated into a single vector. Gradients are modified in-place.
#     Arguments:
#         parameters (Iterable[Variable]): an iterable of Variables that will have
#             gradients normalized
#         max_norm (float or int): max norm of the gradients
#         norm_type (float or int): type of the used p-norm. Can be ``'inf'`` for infinity norm.
#     Returns: grad norm before clipping (for logging mostly)
#     """
#     parameters = list(parameters)
#     max_norm = float(max_norm)
#     norm_type = float(norm_type)
#     if norm_type == float('inf'):
#         total_norm = max(p.grad.abs().max() for p in parameters)
#     else:
#         total_norm = 0
#         for p in parameters:
#             if p.grad is not None:
#                 param_norm = p.grad.norm(norm_type)
#                 total_norm += param_norm ** norm_type

#         total_norm = total_norm ** (1. / norm_type)
#     clip_coef = max_norm / (total_norm + 1e-6)
#     if clip_coef >= 1:
#         return total_norm
#     for p in parameters:
#         if p.grad is not None:
#             p.grad.mul_(clip_coef)

#     return total_norm


def get_model_grad_norm(model,norm_type=2):
    parameters = model.parameters()
    parameters = list(parameters)

    if norm_type == float('inf'):
        total_norm = max(p.grad.abs().max() for p in parameters)
    else:
        total_norm = 0.0
        for p in parameters:
            if p.grad is None:
                # should not happen ?
                return None
            param_norm = p.grad.norm(norm_type)
            total_norm += param_norm ** norm_type

    return float(total_norm)

def model_param_norm(model,norm_type=2):
    # based on https://discuss.pytorch.org/t/how-does-one-implement-weight-regularization-l1-or-l2-manually-without-optimum/7951
    parameters = model.parameters()

    total_norm = None
    for p in model.parameters():
        param_norm = p.norm(norm_type)
        if total_norm is not None:
            total_norm += param_norm ** norm_type
        else:
            total_norm = param_norm ** norm_type

    return total_norm


def get_grad_norms(model, norm_type=2):
    """
    Get grad norms
    :param model:  torch.nn.ModuleDict
    :param norm_type: norm type, default L2
    :return: dict with norms
    """
    r={}
    for m in model:
        r.update({m:get_model_grad_norm(model[m])})

    return r

def print_accuracies(values):
    """
    Pretty plot of accuracies.
    """
    assert all(len(x) == 2 for x in values)
    for name, value in values:
        logger.info('{:<20}: {:>6}'.format(name, '%.3f%%' % (100 * value)))
    logger.info('')



def calculate_genkappa_inter(target,output):
    """
    calculate generalized kappa overlap
    target - integer tensor
    output - one-hot tensor
    """
    _,_output = output.max(1) # last dimension contains the outputs

    total_volume     = float(torch.sum(_output.gt(0).long())+torch.sum(target.gt(0).long()))
    intersect_volume = float(torch.sum( torch.mul( torch.mul( _output.gt(0).long() , target.gt(0).long() ),torch.eq(_output,target).long()) ))

    if total_volume>0.0:
        return 2.0*intersect_volume/total_volume
    else:
        return 0.0

def calculate_genkappa_inter_3d(target,output):
    """
    calculate generalized kappa overlap
    target - integer tensor: Bx  ZxYxX
    output - one-hot tensor: BxCxZxYxX
    """
    _output = output.max(1)[1].view(-1) #
    _target = target.view(-1)

    total_volume     = float(torch.sum(_output.gt(0).long())+torch.sum(_target.gt(0).long()))
    intersect_volume = float(torch.sum( torch.mul( torch.mul( _output.gt(0).long() , _target.gt(0).long() ),torch.eq(_output,_target).long()) ))

    if total_volume>0.0:
        return 2.0*intersect_volume/total_volume
    else:
        return 0.0

def calculate_accuracy_inter(target, output):
    """
    calculate accuracy
    target - integer tensor
    output - one-hot tensor
    """
    _,_output = output.max(1) # last dimension contains the outputs

    return float( torch.sum( torch.eq(_output, target).long()) )/float(target.nelement())

def calculate_genkappa(target, output):
    """
    calculate generalized kappa overlap
    target - integer tensor
    output - integer tensor
    """
    _target=target.view(-1).long()
    _output=output.view(-1).long()
    _t = _target.gt(0).long()
    _o = _output.gt(0).long()
    total_volume     = float(torch.sum(_o)+torch.sum(_t))
    intersect_volume = float(torch.sum( torch.mul( torch.mul( _o , _t ),torch.eq(_output.long(),_target).long()) ))

    if total_volume>0.0:
        return 2.0*intersect_volume/total_volume
    else:
        return 0.0

def calculate_kappa_by_class(target, output, labels=None, n_cls=None):
    """
    calculate kappa overlap per class, for values above 0

    target - integer tensor
    output - integer tensor
    labels - (optional) names of classes, first one corresponds to background (not used)
    n_cls  -
    returns: dict of kappa_{labels[i]}:kappa

    """
    _target=target.view(-1).long()
    _output=output.view(-1).long()
    if n_cls is None:
        n_cls = _target.max()+1

    res={}
    for c in range(1,n_cls): # not checking kappa for background class
        _t=_target.eq(c).long()
        _o=_output.eq(c).long()
        total_volume     = float(torch.sum(_o)+torch.sum(_t))
        intersect_volume = float(torch.sum( torch.mul( torch.mul( _o , _t ),torch.eq(_o,_t).long()) ))

        if total_volume>0.0:
            _k = 2.0*intersect_volume/total_volume
        else:
            _k = 0.0

        if labels is None:
            res.update({'kappa_{}'.format(c):_k})
        else:
            res.update({'kappa_{}'.format(labels[c]):_k})

    return res


def calculate_hds(seg, dist):
    #
    # convert back to one-hot
    n_cls=seg.shape[0]
    one_hot = torch.zeros_like(dist)
    one_hot.scatter_(1, seg, 1)
    one_hot *= dist
    return one_hot.reshape(n_cls,-1).max(dim=1)[0]


def calculate_accuracy(target, output):
    """
    calculate accuracy (proportion of matching elements)
    target - integer tensor
    output - integer tensor
    """
    _target=target.view(-1).long()
    _output=output.view(-1).long()

    return float( torch.sum( torch.eq(_output, _target).long()) )/float(_target.nelement())


def import_parameters(model,to_load,map_location=None):
    """
    load parameters from previously trained model,
    skip initializing missing or non-matching ones
    main purpose: to use pre-trained models
    """
    loaded = torch.load(to_load, map_location=map_location)
    #assert() TODO: make sure it's a state dict
    if isinstance(model,DDP):
        _model=model.module
    else:
        _model=model

    # debug
    cnt=0
    for k in _model.state_dict().keys():
        if k in loaded.keys() and \
            _model.state_dict()[k].size() == loaded[k].size():
            _model.state_dict()[k].copy_(loaded[k])
            cnt+=1
        else:
            #### DEBUG
            print("model import:ignored key:",k)
    print(f"Imported {cnt} of {len(_model.state_dict().keys())} parameters")
    return model

def load_model(model, to_load, map_location=None):
    """
    load a previously trained model.
    """
    to_load = torch.load(to_load, map_location=map_location)

    if isinstance(to_load, nn.Module):
        #old style , saving the whole model
        # check parameters sizes
        model_params   = set(model.state_dict().keys())
        to_load_params = set(to_load.state_dict().keys())

        assert model_params == to_load_params, (model_params - to_load_params, to_load_params - model_params)
        ### TEMPORARY HACK

        # copy saved parameters
        for k in model.state_dict().keys():
        # HACK
        #if k in to_load.state_dict():
            if model.state_dict()[k].size() != to_load.state_dict()[k].size():
                raise Exception("Expected tensor {} of size {}, but got {}".format(
                    k, model.state_dict()[k].size(),
                    to_load.state_dict()[k].size()
                ))
            model.state_dict()[k].copy_(to_load.state_dict()[k])
    else:
        # it's a state_dict
        print("Loading model state dict...")
        if isinstance(model,DDP):
            model.module.load_state_dict(to_load)
        else:
            model.load_state_dict(to_load)


def save_model(model, name, params, scaler=None):
    """
    Save the model.
    """
    if not os.path.exists(params.output):
        os.makedirs(params.output)

    path = os.path.join(params.output, '{}.pth'.format( name) )
    scaler_path = os.path.join(params.output, '{}_scaler.pth'.format( name) )
    print('Saving the model parameters to {} ...' .format( path))
    if isinstance(model,DDP):
        torch.save(model.module.state_dict(), path)
    else:
        torch.save(model.state_dict(), path)

    if scaler is not None:
        torch.save(scaler.state_dict(),scaler_path)

def save_testing_results(res, name, params):
    if not os.path.exists(params.output):
        os.makedirs(params.output)

    path = os.path.join(params.output, '{}.csv'.format( name) )

    with open(path, 'w', newline='') as csvfile:
        wrt = csv.writer(csvfile, quoting=csv.QUOTE_MINIMAL)
        wrt.writerow(['sample','gkappa'])
        for i in res:
            wrt.writerow(i)

def save_testing_results_json(res, name, params):
    if not os.path.exists(params.output):
        os.makedirs(params.output)

    path = os.path.join(params.output, '{}.json'.format(name))
    print(f"Saving results to {path}")
    with open(path, 'w') as f:
        json.dump(res,f)


def segment_with_patches(dataset, model, use_cuda=True, crop=0, bck=0, old=True):
    """
    Apply model to dataset of arbitrary size
    """
    dsize  = dataset.size()

    output_dataset = torch.LongTensor( dataset.size() )
    output_dataset.fill_( bck )

    patch_sz = model.patch_sz
    patch_sz_= patch_sz - crop*2

    out_roi=[dsize[2]-crop*2,dsize[3]-crop*2,dsize[4]-crop*2]

    # TODO:
    with torch.no_grad():
        for k in range(math.ceil( out_roi[0]/patch_sz_ )):
            for l in range(math.ceil( out_roi[1]/patch_sz_ )):
                for m in range(math.ceil( out_roi[2]/patch_sz_ )):

                    c  =  [ k*patch_sz_ + crop, l*patch_sz_ + crop, m*patch_sz_ + crop]

                    for i in range(3):
                        c[i]=min(c[i], dsize[i+2] - patch_sz + crop)

                    # extract a patch
                    in_data = dataset.narrow(2, c[0]-crop, patch_sz).narrow(3,c[1]-crop, patch_sz).narrow(4,c[2]-crop, patch_sz)

                    if use_cuda:
                        in_data=in_data.cuda()

                    output = model.forward(in_data)

                    _output = output.max(1)[1].cpu()

                    _output_data=_output.view(1,1, patch_sz, patch_sz, patch_sz)

                    # paste data into the ouput buffer
                    output_dataset.narrow(2, c[0], patch_sz_).narrow(3, c[1],patch_sz_ ).narrow(4, c[2], patch_sz_).copy_(
                        _output_data.narrow(2, crop, patch_sz_).narrow(3, crop, patch_sz_).narrow(4, crop, patch_sz_)
                    )
    #
    return output_dataset


def segment_with_patches2(dataset, model, use_cuda=True, crop=0, bck=0, patch_sz=None):
    """
    Apply model to dataset of arbitrary size
    """
    dsize = dataset.size()
    output_size = list(dsize)
    output_size[1] = 1 # strip additional dimensions

    output_dataset = torch.LongTensor(*output_size)
    output_dataset.fill_(bck)

    if patch_sz is None: # backward compatibility
        patch_sz = model.patch_sz

    patch_sz_ = patch_sz - crop*2

    out_roi = [ dsize[2]-crop*2, dsize[3]-crop*2, dsize[4]-crop*2 ]

    # TODO:
    with torch.no_grad():
        for k in range(math.ceil( out_roi[0]/patch_sz_ )):
            for l in range(math.ceil( out_roi[1]/patch_sz_ )):
                for m in range(math.ceil( out_roi[2]/patch_sz_ )):

                    c = [k*patch_sz_ + crop, l*patch_sz_ + crop, m*patch_sz_ + crop]

                    for i in range(3):
                        c[i] = min(c[i], dsize[i+2] - patch_sz + crop-1)

                    # extract a patch
                    in_data = dataset[:, :, c[0]-crop: c[0]-crop+patch_sz, c[1]-crop: c[1]-crop+patch_sz, c[2]-crop: c[2]-crop+patch_sz]

                    if use_cuda:
                        in_data = in_data.cuda()

                    # TODO: add VAE output
                    output = model.forward(in_data)

                    _output = output.max(1)[1].cpu().unsqueeze_(1)

                    # paste data into the output buffer
                    output_dataset[:, :, c[0]: c[0]+patch_sz_, c[1]: c[1]+patch_sz_, c[2]: c[2]+patch_sz_] = \
                           _output[:, :, crop: crop+patch_sz_, crop: crop+patch_sz_, crop: crop+patch_sz_]

    return output_dataset


def segment_with_patches_fuzzy(
        dataset, model,
        use_cuda = True,
        crop     = 0,
        patch_sz = None,
        stride   = None,
        use_max  = False,
        n_cls    = None,
        bck      = 0,
        prefilter=None,
        out_prefiltered = False,
        use_group = False,
        use_seg = False,
        use_vae = False, # TODO: implement this
        group_downsample = None,
        group_fuzzy = False,
        groups = 1 ):
    """
    Apply model to dataset of arbitrary size
    Result is fuzzy, can be convert to discrete by .max(1)[1]
    """
    dsize = dataset.size()
    output_size = list(dsize)
    output_size[1] = 1 # strip additional dimensions

    output_dataset = torch.LongTensor( *output_size )
    output_dataset.fill_(bck)

    if use_group:
        if group_downsample is None:
            if isinstance(model, DDP):
                group_downsample = model.module.group_downsample
            else:
                group_downsample = model.group_downsample
        if group_fuzzy:
            output_group_size = [output_size[0], groups, output_size[2]//group_downsample, output_size[3]//group_downsample, output_size[4]//group_downsample]
            output_group_dataset = torch.Tensor( *output_group_size )
        else:
            output_group_size = [output_size[0], output_size[1], output_size[2]//group_downsample, output_size[3]//group_downsample, output_size[4]//group_downsample]
            output_group_dataset = torch.LongTensor( *output_group_size )
            # TODO: use different default group?
            output_group_dataset.fill_(bck)
    else:
        output_group_dataset = None

    if out_prefiltered : # assume that prefilter segments into n_cls too
        output_dataset_pre = torch.LongTensor( *output_size )
        output_dataset_pre.fill_(bck) # TODO: group bck?

    if patch_sz is None: # backward compatibility
        if isinstance(model, DDP):
            patch_sz = model.module.patch_sz
        else:
            patch_sz = model.patch_sz

    patch_sz_ = patch_sz - crop*2

    out_roi = [ dsize[2]-crop*2, dsize[3]-crop*2, dsize[4]-crop*2 ]

    # TODO:
    with torch.no_grad():
        for k in range(math.ceil( out_roi[0]/patch_sz_ )):
            for l in range(math.ceil( out_roi[1]/patch_sz_ )):
                for m in range(math.ceil( out_roi[2]/patch_sz_ )):

                    c = [k*patch_sz_ + crop, l*patch_sz_ + crop, m*patch_sz_ + crop]

                    for i in range(3):
                        c[i] = max(min(c[i], dsize[i+2] - patch_sz + crop-1),0)

                    # extract a patch
                    in_data = dataset[:, :, c[0]-crop: c[0]-crop+patch_sz, c[1]-crop: c[1]-crop+patch_sz, c[2]-crop: c[2]-crop+patch_sz]

                    if use_cuda:
                        in_data = in_data.cuda()

                    if prefilter is not None: in_data = prefilter.forward(in_data)['seg']

                    if out_prefiltered: # HACK, assumes prefilter outpus last n_cls of segmenttions
                        _output_pre = in_data[:,-n_cls:,:,:,:].max(1)[1].cpu().unsqueeze_(1)
                        # paste data into the output buffer
                        output_dataset_pre[:, :, c[0]: c[0]+patch_sz_, c[1]: c[1]+patch_sz_, c[2]: c[2]+patch_sz_] = \
                            _output_pre[:, :, crop: crop+patch_sz_, crop: crop+patch_sz_, crop: crop+patch_sz_]

                    #print(f"{c=} {k=} {l=} {m=}")
                    output = model.forward( in_data ) #TODO: add skip_vae

                    if use_seg:
                        _output = output['seg'].max(1)[1].cpu().unsqueeze_(1)

                        # paste data into the output buffer
                        output_dataset[:, :, c[0]: c[0]+patch_sz_, c[1]: c[1]+patch_sz_, c[2]: c[2]+patch_sz_] = \
                               _output[:, :, crop: crop+patch_sz_, crop: crop+patch_sz_, crop: crop+patch_sz_]

                    # assemble group prediction
                    if use_group and group_downsample is not None:
                        if group_fuzzy:
                            output_group_dataset[:,:, c[0]//group_downsample : (c[0]+patch_sz)//group_downsample,
                                                      c[1]//group_downsample : (c[1]+patch_sz)//group_downsample,
                                                      c[2]//group_downsample : (c[2]+patch_sz)//group_downsample] = \
                                output['group'].cpu()
                        else:
                            output_group_dataset[:,:, c[0]//group_downsample : (c[0]+patch_sz)//group_downsample,
                                                      c[1]//group_downsample : (c[1]+patch_sz)//group_downsample,
                                                      c[2]//group_downsample : (c[2]+patch_sz)//group_downsample] = \
                                output['group'].max(1)[1].cpu().unsqueeze_(1)

                    # TODO: store VAE output
    if out_prefiltered:
        return output_dataset, output_dataset_pre, output_group_dataset
    else:
        return output_dataset, output_group_dataset

def segment_with_patches_overlap(
        dataset, model,
        use_cuda=True,
        crop=0,
        patch_sz = None,
        stride = None,
        bck = 0,
        out_fuzzy=False,
        out_vae  =False,
        out_latent_vectors=False,
        loc=None,prec=None):
    """
    Apply model to dataset of arbitrary size
    Arguments:
        dataset - torch.Tensor of input data, 5D
        model - torch model
    Keyword arguments:
        use_cuda - use CUDA for inference
        crop - crop patches by this many voxels in all spatial dimensions for output
        patch_sz - size of patch to process with a model
        stride - step between patches, patches can overlap
        bck - background value, to be used for areas where model was not applied
        out_fuzzy - output fuzzy results instead of just discrete
        out_latent_vectors - output internal latent vectors
        loc  - latent vectors location
        prec - latent vectors precision matrix (inv covariance)
    Returns:
        3D segmentation , if out_fuzzy is False
        tuple: 3D segmentation, 4D fuzzy output if out_fuzzy is True
        tuple: 3D segmentation, 3D -log likelyhood if loc and prec are given
    """
    dsize = dataset.size()
    output_size = list( dsize )
    output_size_fuzzy = list( dsize )
    output_size_vae   = list( dsize )
    output_size_likelihood = list( dsize )

    output_size_likelihood[1]=1
    output_size[1] = 1

    output_fuzzy  = None
    output_vae    = None
    output_weight = torch.zeros( output_size )
    output_likelihood = None

    if patch_sz is None: # backward compatibility
        patch_sz = model.patch_sz

    patch_sz_ = patch_sz - crop*2

    if stride is None:
        stride = patch_sz_

    out_roi = [ dsize[2]-crop*2, dsize[3]-crop*2, dsize[4]-crop*2 ]
    ones = torch.ones([1,1,patch_sz_,patch_sz_,patch_sz_])
    latent_vectors = []

    if loc is not None and prec is not None:
        # convert to expected format
        _loc  = loc.reshape(1, loc.shape[0]).to(torch.double)
        _prec = prec.reshape(1,prec.shape[0],prec.shape[1]).to(torch.double)
        # https://stats.stackexchange.com/questions/97408/relation-of-mahalanobis-distance-to-log-likelihood
        # https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Likelihood_function
        # import math
        # here we replaced the determinant of covariance matrix with
        # inverse of the determinat of the precision matrix (since it't inverse of covariance matrix)
        _c = -0.5*math.log(float(torch.linalg.det(prec))) + prec.shape[0]*math.log(2*math.pi)
        print('_c',_c)


    # TODO:
    with torch.no_grad():
        for k in range(math.ceil( (out_roi[0]-crop)/stride )):
            for l in range(math.ceil( (out_roi[1]-crop)/stride )):
                for m in range(math.ceil( (out_roi[2]-crop)/stride )):

                    c = [k*stride + crop, l*stride + crop, m*stride + crop]

                    for i in range(3):
                        c[i] = max(min(c[i], dsize[i+2] - patch_sz + crop - 1),crop)
                    # extract a patch
                    in_data = dataset[:, :, c[0]-crop: c[0]-crop+patch_sz, c[1]-crop: c[1]-crop+patch_sz, c[2]-crop: c[2]-crop+patch_sz]

                    if use_cuda:
                        in_data = in_data.cuda()

                    out_ = model.forward(in_data)
                    out = out_['seg']
                    if 'vae' in out_:
                        vae =  out_['vae']
                    else:
                        vae = None

                    if out_latent_vectors and 'latent' in out_:
                        latent_vectors += [out_['latent'].cpu()]

                    patch_output = torch.log_softmax(out,1)

                    if use_cuda:
                        patch_output = patch_output.cpu()

                    if output_fuzzy is None:
                        # need to allocate depending on the number of output classses
                        output_size_fuzzy[1] = patch_output.size()[1]
                        output_fuzzy = torch.zeros( *output_size_fuzzy )

                    if out_vae and vae is not None:
                        if output_vae is None:
                            # need to allocate depending on the number of output classses
                            output_size_vae[1] = vae.size()[1]
                            output_vae = torch.zeros( *output_size_vae )

                        if use_cuda:
                            vae=vae.cpu()

                        output_vae[:, :, c[0]: c[0]+patch_sz_, c[1]: c[1]+patch_sz_, c[2]: c[2]+patch_sz_] += \
                            vae[:, :, crop: crop+patch_sz_, crop: crop+patch_sz_, crop: crop+patch_sz_]

                    if loc is not None and prec is not None:
                        if output_likelihood is None:
                            output_likelihood= torch.zeros( *output_size_likelihood )

                        if 'latent' in out_:

                            lat=out_['latent']
                            # calculate log-likelyhood of latent vector being from N-d gaussian distribution of latent vectors
                            #print("lat:",lat.shape)
                            _lat=lat.reshape(lat.shape[0],lat.shape[1],-1).transpose(1,2).reshape(-1,lat.shape[1]).to(torch.double)

                            #print("_lat:",_lat.shape,'loc:',_loc.shape)
                            _diff=(_lat-_loc).reshape(-1,1,lat.shape[1])
                            #print("_diff:",_diff.shape)
                            # now we have dist calculated by n_batch*n_vox:
                            mahalanobis_dist2=torch.matmul(torch.matmul(_diff,_prec),_diff.transpose(1,2))
                            #print("mahalanobis_dist2:",mahalanobis_dist2.shape)
                            log_likelihood = 0.5*mahalanobis_dist2 + _c
                            # TODO: fix this?
                            #log_likelihood = torch.sqrt(mahalanobis_dist2)
                            # reshape back into a patch
                            log_likelihood = log_likelihood.reshape(lat.shape[0],1,lat.shape[2],lat.shape[3],lat.shape[4]).to(torch.float)
                            # choose resample function
                            log_likelihood = nn.functional.interpolate(
                                                        log_likelihood,
                                                        size=patch_output.shape[2:5],  mode='trilinear', align_corners=False) # TODO
                            if use_cuda:
                                log_likelihood=log_likelihood.cpu()

                            output_likelihood[:, :, c[0]: c[0]+patch_sz_, c[1]: c[1]+patch_sz_, c[2]: c[2]+patch_sz_] += \
                                log_likelihood[:, :, crop: crop+patch_sz_, crop: crop+patch_sz_, crop: crop+patch_sz_]

                    # accumulate data
                    output_fuzzy[:, :, c[0]: c[0]+patch_sz_, c[1]: c[1]+patch_sz_, c[2]: c[2]+patch_sz_] += \
                         patch_output[:, :, crop: crop+patch_sz_, crop: crop+patch_sz_, crop: crop+patch_sz_]

                    output_weight[:, :, c[0]: c[0]+patch_sz_, c[1]: c[1]+patch_sz_, c[2]: c[2]+patch_sz_] += \
                         ones

        # aggregate weights
        #output_weight_save = output_weight.clone().detach()

        invalid = output_weight<1.0

        output_weight.masked_fill_(invalid, 1.0 )

        output_fuzzy /= output_weight


        output_fuzzy = nn.functional.softmax(output_fuzzy, 1)

        # set BG to 1 where mask was invalid
        output_fuzzy[:,0:1,:,:,:].masked_fill_(invalid, 1.0 )
        output_fuzzy[:,1:,:,:,:].masked_fill_(invalid, 0.0 )

        output_dataset = output_fuzzy.max(1)[1].cpu().unsqueeze_(1)
        output_dataset.masked_fill_(invalid, bck )

        if out_vae and output_vae is not None:
            output_vae /= output_weight

        if output_likelihood is not None:
            output_likelihood /= output_weight
        if out_fuzzy or out_vae or out_latent_vectors:
            return output_dataset, output_fuzzy, output_vae, latent_vectors
        elif loc is not None and prec is not None:
            return output_dataset, torch.exp(-1.0*output_likelihood)
        else:
            return output_dataset #, output_weight_save


def ae_with_patches(dataset, model, use_cuda=True, crop=0, bck=0.0, old=True, patch_sz=None):
    """
    Apply auto-encoder model to dataset of arbitrary size
    """

    dsize    = dataset.size()

    output_dataset = torch.FloatTensor( dataset.size() )
    output_dataset.fill_( bck )

    if patch_sz is None:# backward compatibility
        patch_sz = model.patch_sz

    patch_sz_= patch_sz - crop*2

    out_roi=[dsize[2]-crop*2,dsize[3]-crop*2,dsize[4]-crop*2]

    # TODO:
    with torch.no_grad():
        for k in range(math.ceil( out_roi[0]/patch_sz_ )):
            for l in range(math.ceil( out_roi[1]/patch_sz_ )):
                for m in range(math.ceil( out_roi[2]/patch_sz_ )):

                    c  =  [ k*patch_sz_ + crop, l*patch_sz_ + crop, m*patch_sz_ + crop]

                    for i in range(3):
                        c[i]=min(c[i], dsize[i+2] - patch_sz + crop)

                    # extract a patch
                    in_data = dataset.narrow(2,c[0]-crop,patch_sz).narrow(3,c[1]-crop,patch_sz).narrow(4,c[2]-crop,patch_sz)

                    if use_cuda:
                        in_data=in_data.cuda()

                    output = model.forward(in_data)

                    _output = output.cpu()

                    _output_data=_output.view(1,1, patch_sz, patch_sz, patch_sz)

                    # paste data into the ouput buffer
                    output_dataset.narrow(2, c[0], patch_sz_).narrow(3, c[1],patch_sz_ ).narrow(4, c[2], patch_sz_).copy_(
                        _output_data.narrow(2, crop, patch_sz_).narrow(3, crop, patch_sz_).narrow(4, crop, patch_sz_)
                        )
    #
    return output_dataset



def dump_minc_file(dataset,path):

        from minc2_simple import minc2_file

        image_dims=dataset.shape[2:5]

        dims=[
            { 'id':minc2_file.MINC2_DIM_X,  'length':image_dims[2],'start':0, 'step':1.0},
            { 'id':minc2_file.MINC2_DIM_Y,  'length':image_dims[1],'start':0, 'step':1.0},
            { 'id':minc2_file.MINC2_DIM_Z,  'length':image_dims[0],'start':0, 'step':1.0}
        ]

        o=minc2_file()
        if dataset.dtype.name=='int8' or dataset.dtype.name=='int16' or dataset.dtype.name=='int32':
            o.define(dims, minc2_file.MINC2_BYTE, minc2_file.MINC2_BYTE)
        else:
            o.define(dims, minc2_file.MINC2_FLOAT, minc2_file.MINC2_FLOAT)

        o.create(path)
        o.setup_standard_order()
        if dataset.dtype.name=='int8' or dataset.dtype.name=='int16' or dataset.dtype.name=='int32':
            o.save_complete_volume(dataset.astype('int8'))
        else:
            o.save_complete_volume(dataset)
        o.close()


# def qc_minibatch(x, y, n_cls):
    # s = x.size()

    # return torchvision.utils.make_grid(
           # [y.narrow(0,0,1).narrow(2, int(s[2]*0.25 ), 1).squeeze().unsqueeze(0).float()/(n_cls-1)*2-1.0,
            # x.narrow(0,0,1).narrow(2, int(s[2]*0.25 ), 1).squeeze().unsqueeze(0),
            # y.narrow(0,0,1).narrow(2, int(s[2]*0.50 ), 1).squeeze().unsqueeze(0).float()/(n_cls-1)*2-1.0,
            # x.narrow(0,0,1).narrow(2, int(s[2]*0.50 ), 1).squeeze().unsqueeze(0),
            # y.narrow(0,0,1).narrow(2, int(s[2]*0.75 ), 1).squeeze().unsqueeze(0).float()/(n_cls-1)*2-1.0,
            # x.narrow(0,0,1).narrow(2, int(s[2]*0.75 ), 1).squeeze().unsqueeze(0)
            # ], normalize=True
        # )


# def qc_minibatch_ovl(x, y, lut, alpha=0.5, flip_y=True, axis=1, samples=10, nrow=5, padding=2, norm=True):
    # """
    # Generate QC image for multi-class segmentation task
    # :param x: grayscale volume
    # :param y: segmentation volume (discrete)
    # :param lut: lookup table
    # :param alpha: alpha value for gray-scale for blending
    # :param flip_y: flip y axis (to account fo on-screen origin on top left)
    # :param axis: axis along which to extract slices (0-z, 1-y, 2-x)
    # :param samples: number of slices to extract
    # :param nrow: number of samples per row in a grid
    # :param padding: add padding to the image
    # :param norm: normalize between RGB slices
    # :return:
    # """
    # s = x.size()
    # x_ = x
    # if norm:
        # _min = x.min()
        # _max = x.max()
        # x_ = (x_-_min)/(_max-_min)

    # d_ = s[axis]/(samples+1)

    # y_ = [y.narrow(axis, i, 1).squeeze() for i in (torch.arange(samples) * d_ + d_).long()]
    # x_ = [x_.narrow(axis, i, 1).squeeze().unsqueeze(0) for i in (torch.arange(samples)*d_+d_).long()]

    # y_ = [seg_to_rgb(i, lut) for i in y_]
    # alpha_ = torch.full_like(x_[0], alpha)
    # x_ = [torch.cat((i, i, i, alpha_), 0) for i in x_]

    # b_ = [ alpha_blend(i, j).narrow(0, 0, 3) for (i,j) in zip(x_,y_)]

    # if flip_y:
         # b_ = [i.index_select(1, torch.arange(i.size(1)-1, -1, -1, dtype=torch.long)) for i in b_]

    # return torchvision.utils.make_grid(
        # b_, normalize=False, nrow=nrow, padding=padding
    # )

# def qc_minibatch_gs(x, flip_y=True, axis=1, samples=10, nrow=5, padding=2, norm=True):
    # """
    # Generate QC image for grayscale image
    # :param x: grayscale volume
    # :param flip_y: flip y axis (to account fo on-screen origin on top left)
    # :param axis: axis along which to extract slices (0-z, 1-y, 2-x)
    # :param samples: number of slices to extract
    # :param nrow: number of samples per row in a grid
    # :param padding: add padding to the image
    # :param norm: normalize between RGB slices
    # :return:
    # """
    # s = x.size()
    # x_ = x
    # if norm:
        # _min = x.min()
        # _max = x.max()
        # x_ = (x_-_min)/(_max-_min)

    # d_ = s[axis]/(samples+1)

    # x_ = [x_.narrow(axis, i, 1).squeeze().unsqueeze(0) for i in (torch.arange(samples)*d_+d_).long()]

    # b_ = x_

    # if flip_y:
         # b_ = [i.index_select(1, torch.arange(i.size(1)-1, -1, -1, dtype=torch.long)) for i in b_]

    # return torchvision.utils.make_grid(
        # b_, normalize=False, nrow=nrow, padding=padding
    # )


# def qc_ae(x, y):
    # s=x.size()
    # return torchvision.utils.make_grid(
        # [ # Z
          # x.narrow(2, s[2]//4,   1).squeeze().unsqueeze(0).cpu(),
          # y.narrow(2, s[2]//4,   1).squeeze().unsqueeze(0).cpu(),

          # x.narrow(2, s[2]//2,   1).squeeze().unsqueeze(0).cpu(),
          # y.narrow(2, s[2]//2,   1).squeeze().unsqueeze(0).cpu(),

          # x.narrow(2, 3*s[2]//4, 1).squeeze().unsqueeze(0).cpu(),
          # y.narrow(2, 3*s[2]//4, 1).squeeze().unsqueeze(0).cpu(),
          # # Y
          # x.narrow(1, s[1]//4,   1).squeeze().unsqueeze(0).cpu(),
          # y.narrow(1, s[1]//4,   1).squeeze().unsqueeze(0).cpu(),

          # x.narrow(1, s[1]//2,   1).squeeze().unsqueeze(0).cpu(),
          # y.narrow(1, s[1]//2,   1).squeeze().unsqueeze(0).cpu(),

          # x.narrow(1, 3*s[1]//4, 1).squeeze().unsqueeze(0).cpu(),
          # y.narrow(1, 3*s[1]//4, 1).squeeze().unsqueeze(0).cpu(),
          # # X
          # x.narrow(0, s[0]//4,   1).squeeze().unsqueeze(0).cpu(),
          # y.narrow(0, s[0]//4,   1).squeeze().unsqueeze(0).cpu(),

          # x.narrow(0, s[0]//2,   1).squeeze().unsqueeze(0).cpu(),
          # y.narrow(0, s[0]//2,   1).squeeze().unsqueeze(0).cpu(),

          # x.narrow(0, 3*s[0]//4, 1).squeeze().unsqueeze(0).cpu(),
          # y.narrow(0, 3*s[0]//4, 1).squeeze().unsqueeze(0).cpu()
        # ], normalize=True ,nrow=6
        # )

# def qc_ae_axial(x, y):
    # s=x.size()
    # return torchvision.utils.make_grid(
        # [ # Z
          # x.narrow(2, s[2]//4,   1).squeeze().unsqueeze(0).cpu(),
          # y.narrow(2, s[2]//4,   1).squeeze().unsqueeze(0).cpu(),

          # x.narrow(2, s[2]//2,   1).squeeze().unsqueeze(0).cpu(),
          # y.narrow(2, s[2]//2,   1).squeeze().unsqueeze(0).cpu(),

          # x.narrow(2, 3*s[2]//4, 1).squeeze().unsqueeze(0).cpu(),
          # y.narrow(2, 3*s[2]//4, 1).squeeze().unsqueeze(0).cpu()
        # ], normalize=True ,nrow=6
        # )

def one_hot_3d(input_seg, n_cls, val=1.0):
    """
    Converts categorical input into one-hot
    for 3D batched input: BxZxYxX
    outputs: BxCxZxYxX
    """
    sz = input_seg.size()

    # create in cuda?
    one_hot_out = torch.FloatTensor(sz[0], n_cls, sz[1], sz[2], sz[3]).cuda().zero_()

    # scatter 1.0
    one_hot_out.scatter_(1, input_seg.unsqueeze(1), val)
    return one_hot_out

def seg_to_rgb(seg, lut):
    """
    convert integer segmentation 3D tensor to RGB 3D tensor using LUT
    """
    sz = seg.size()
    lz = lut.size(0)

    return lut.index_select(1, seg.contiguous().view(-1)).view(lz, *sz)


def alpha_blend(bck, drawing):
    """
    Blend two images:
        bck - RGB tensor
        drawing - RGBA tensor
        weight - floating weight of drawing
    Output:
        RGB tensor
    """
    _eps = 1e-4

    alpha_bck = bck.narrow(0, 3, 1)
    alpha_drawing = drawing.narrow(0, 3, 1)

    alpha = alpha_bck + alpha_drawing - alpha_drawing*alpha_bck
    alpha_mask = alpha.le(_eps)
    out = torch.empty_like(bck)

    out.narrow(0, 0, 3).copy_((bck.narrow(0, 0, 3)*alpha_bck + drawing.narrow(0, 0, 3)*alpha_drawing - drawing.narrow(0, 0, 3)*alpha_drawing*alpha_bck)/alpha )
    out.narrow(0, 0, 3)[torch.cat([alpha_mask, alpha_mask, alpha_mask], 0)].fill_(0.0)
    out.narrow(0, 3, 1).copy_(alpha)

    return out


def gen_discrete_lut(n_cls,cmap="plasma"):
    """
    Generate discrete lut ,
    :param cmap: color map name, see https://matplotlib.org/examples/color/colormaps_reference.html
    :param n_cls: number of colours (including background)
    :return: tensor n_cls * 4
    """
    # for colour luts
    import matplotlib.pyplot as plt
    import matplotlib.cm  as cmx
    import matplotlib.colors as colors
    import numpy as np

    cmo = copy.copy(plt.get_cmap(cmap))

    # create colour map with 0 set to transparent
    cmo.set_bad('k', alpha=0.0)
    cmo.set_under('k', alpha=0.0)

    oNorm = colors.Normalize(vmin=1, vmax=n_cls)

    oscalarMap = cmx.ScalarMappable(norm=oNorm, cmap=cmo)

    lut = torch.FloatTensor(oscalarMap.to_rgba(np.linspace(0, n_cls, n_cls+1)).astype('float32')).t()

    return lut


def init_optimizer(desc, model):
    """
    Initialize optimizers and scheduler
    :param desc: dict
    :param model: model to optimize
    :return: (optimizer, scheduler, clip_grad_norm )
    """
    optimizer       = get_optimizer(model, desc['optimizer'])
    scheduler       = None
    clip_grad_norm  = 0
    loss_norm_decay = 1000
    init_lr         = None

    if 'scheduler' in desc:
        if desc['scheduler'] == 'step':
            scheduler = lr_scheduler.StepLR(optimizer, step_size=desc['scheduler_step'], gamma=desc.get('gamma',0.1))
        elif desc['scheduler'] == 'plateau':
            scheduler = lr_scheduler.ReduceLROnPlateau(optimizer, mode=desc.get('mode', 'max'))
    if 'clip_grad_norm' in desc:
        clip_grad_norm = desc['clip_grad_norm']
    if 'loss_norm_decay' in desc:
      loss_norm_decay = desc['loss_norm_decay']
    # hack to get inital learning rate
    for param_group  in optimizer.param_groups:
      init_lr = param_group.get('lr',init_lr)

    return optimizer, scheduler, clip_grad_norm, loss_norm_decay, init_lr


def sklearn_auroc(y_pred,y_true):
    """
    calculate AUROC using scikit learn
    works only for binary (two class labels)
    """
    try:
        from sklearn.metrics import roc_auc_score

        y_true_ = y_true.cpu().detach()
        y_pred_ = y_pred[:,1].cpu().detach()
        try:
            return roc_auc_score(y_true_.numpy(), y_pred_.numpy() )
        except ValueError:
            # probably only one class is present
            return 0.0
    except ModuleNotFoundError:
        # no sklearn :(
        return 0.0

# kate: space-indent on; indent-width 4; indent-mode python;replace-tabs on;word-wrap-column 80
