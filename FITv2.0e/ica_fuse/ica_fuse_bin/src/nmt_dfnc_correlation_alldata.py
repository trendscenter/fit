'''
Build a attention-based neural machine translation model
'''
import os
os.environ["THEANO_FLAGS"] = "device=cpu,floatX=float32,gpuarray.preallocate=1"
#os.environ["THEANO_FLAGS"] = "device=cpu,floatX=float32,gpuarray.preallocate=1,blas.ldflags=-L'C:\\Users\\srrac\\Anaconda2\\Library\\bin' -lmkl_core -lmkl_intel_thread -lmkl_rt"
#os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS']='1'
#os.environ['FOR_DISABLE_CONSOLE_CTRL_HANDLER'] = 'T'
import theano
import theano.tensor as tensor
#theano.sandbox.cuda.use("gpu0")
#from theano.sandbox.rng_mrg import MRG_RandomStreams as RandomStreams
from theano.tensor.shared_randomstreams import RandomStreams
import numpy
import csv
import sys

import warnings

from collections import OrderedDict

trng = RandomStreams(1234)
use_noise = theano.shared(0)
prob = 0.5

# push parameters to Theano shared variables
def zipp(params, tparams):
    for kk, vv in params.iteritems():
        tparams[kk].set_value(vv)

# pull parameters from Theano shared variables
def unzip(zipped):
    new_params = OrderedDict()
    for kk, vv in zipped.iteritems():
        new_params[kk] = vv.get_value()
    return new_params

# get the list of parameters: Note that tparams must be OrderedDict
def itemlist(tparams):
    return [vv for kk, vv in tparams.iteritems()]

# dropout
def dropout_layer(state_before):
    proj = tensor.switch(use_noise, 
            state_before * trng.binomial(state_before.shape, p=prob, n=1, dtype=state_before.dtype),
            state_before * prob)
    return proj

# make prefix-appended name
def _p(pp, name):
    return '%s_%s'%(pp, name)

# initialize Theano shared variables according to the initial parameters
def init_tparams(params):
    tparams = OrderedDict()
    for kk, pp in params.iteritems():
        tparams[kk] = theano.shared(params[kk], name=kk)
    return tparams

# load parameters
def load_params(path, params):
    pp = numpy.load(path)
    for kk, vv in params.iteritems():
        if kk not in pp:
            warnings.warn('%s is not in the archive'%kk)
            continue
        params[kk] = pp[kk]

    return params

# layers: 'name': ('parameter initializer', 'feedforward')
layers = {'ff': ('param_init_fflayer', 'fflayer'), 
          'gru': ('param_init_gru', 'gru_layer'),
          'gru_cond': ('param_init_gru_cond', 'gru_cond_layer'),
          }

def get_layer(name):
    fns = layers[name]
    return (eval(fns[0]), eval(fns[1]))

# some utilities
def ortho_weight(ndim):
    W = numpy.random.randn(ndim, ndim)
    u, s, v = numpy.linalg.svd(W)
    return u.astype('float32')

def norm_weight(nin,nout=None, scale=0.01, ortho=True):
    if nout == None:
        nout = nin
    if nout == nin and ortho:
        W = ortho_weight(nin)
    else:
        W = scale * numpy.random.randn(nin, nout)
    return W.astype('float32')

def tanh(x):
    return tensor.tanh(x)

def linear(x):
    return x

def concatenate(tensor_list, axis=0):
    """
    Alternative implementation of `theano.tensor.concatenate`.
    This function does exactly the same thing, but contrary to Theano's own
    implementation, the gradient is implemented on the GPU.
    Backpropagating through `theano.tensor.concatenate` yields slowdowns
    because the inverse operation (splitting) needs to be done on the CPU.
    This implementation does not have that problem.
    :usage:
        >>> x, y = theano.tensor.matrices('x', 'y')
        >>> c = concatenate([x, y], axis=1)
    :parameters:
        - tensor_list : list
            list of Theano tensor expressions that should be concatenated.
        - axis : int
            the tensors will be joined along this axis.
    :returns:
        - out : tensor
            the concatenated tensor expression.
    """
    concat_size = sum(tt.shape[axis] for tt in tensor_list)

    output_shape = ()
    for k in range(axis):
        output_shape += (tensor_list[0].shape[k],)
    output_shape += (concat_size,)
    for k in range(axis + 1, tensor_list[0].ndim):
        output_shape += (tensor_list[0].shape[k],)

    out = tensor.zeros(output_shape)
    offset = 0
    for tt in tensor_list:
        indices = ()
        for k in range(axis):
            indices += (slice(None),)
        indices += (slice(offset, offset + tt.shape[axis]),)
        for k in range(axis + 1, tensor_list[0].ndim):
            indices += (slice(None),)

        out = tensor.set_subtensor(out[indices], tt)
        offset += tt.shape[axis]

    return out

# feedforward layer: affine transformation + point-wise nonlinearity
def param_init_fflayer(options, params, prefix='ff', nin=None, nout=None, ortho=True):
    if nin == None:
        nin = options['dim_proj']
    if nout == None:
        nout = options['dim_proj']
    params[_p(prefix,'W')] = norm_weight(nin, nout, scale=0.01, ortho=ortho)
    params[_p(prefix,'b')] = numpy.zeros((nout,)).astype('float32')

    return params

def fflayer(tparams, state_below, options, prefix='rconv', activ='lambda x: tensor.tanh(x)', **kwargs):
    return eval(activ)(tensor.dot(state_below, tparams[_p(prefix,'W')])+tparams[_p(prefix,'b')])

# GRU layer
def param_init_gru(options, params, prefix='gru', nin=None, dim=None):
    if nin == None:
        nin = options['dim_proj']
    if dim == None:
        dim = options['dim_proj']
    W = numpy.concatenate([norm_weight(nin,dim),
                           norm_weight(nin,dim)], axis=1)
    params[_p(prefix,'W')] = W
    U = numpy.concatenate([ortho_weight(dim),
                           ortho_weight(dim)], axis=1)
    params[_p(prefix,'U')] = U
    params[_p(prefix,'b')] = numpy.zeros((2 * dim,)).astype('float32')

    Wx = norm_weight(nin, dim)
    params[_p(prefix,'Wx')] = Wx
    Ux = ortho_weight(dim)
    params[_p(prefix,'Ux')] = Ux
    params[_p(prefix,'bx')] = numpy.zeros((dim,)).astype('float32')

    return params

def gru_layer(tparams, state_below, options, prefix='gru', mask=None, **kwargs):
    nsteps = state_below.shape[0]
    if state_below.ndim == 3:
        n_samples = state_below.shape[1]
    else:
        n_samples = 1

    dim = tparams[_p(prefix,'U')].shape[0]

    if mask == None:
        mask = tensor.alloc(1., state_below.shape[0], 1)

    def _slice(_x, n, dim):
        if _x.ndim == 3:
            return _x[:, :, n*dim:(n+1)*dim]
        return _x[:, n*dim:(n+1)*dim]

    state_below_ = tensor.dot(state_below, tparams[_p(prefix, 'W')]) + tparams[_p(prefix, 'b')]
    state_belowx = tensor.dot(state_below, tparams[_p(prefix, 'Wx')]) + tparams[_p(prefix, 'bx')]

    def _step(m_, x_, xx_, h_):
        preact = tensor.dot(h_, tparams[_p(prefix, 'U')])
        preact += x_

        r = tensor.nnet.sigmoid(_slice(preact, 0, dim))
        u = tensor.nnet.sigmoid(_slice(preact, 1, dim))

        preactx = tensor.dot(h_, tparams[_p(prefix, 'Ux')])
        preactx = preactx * r
        preactx = preactx + xx_

        h = tensor.tanh(preactx)

        h = u * h_ + (1. - u) * h
        h = m_[:,None] * h + (1. - m_)[:,None] * h_

        return h, r, u, preact, preactx

    rval, updates = theano.scan(_step, 
                                sequences=[mask, state_below_, state_belowx],
                                outputs_info = [tensor.alloc(0., n_samples, dim),
                                                None, None, None, None],
                                name=_p(prefix, '_layers'),
                                n_steps=nsteps)
    return rval

# Conditional GRU layer with Attention
def param_init_gru_cond(options, params, prefix='gru_cond', nin=None, dim=None, dimctx=None, nin_below=None):
    if nin == None:
        nin = options['dim']
    if dim == None:
        dim = options['dim']
    if dimctx == None:
        dimctx = options['dim']
    if nin_below == None:
        nin_below = options['nfeature']

    params = param_init_gru(options, params, prefix, nin=nin, dim=dim)

    # context to GRU
    Wc = norm_weight(nin_below,dim*2)
    params[_p(prefix,'Wc')] = Wc

    Wcx = norm_weight(nin_below,dim)
    params[_p(prefix,'Wcx')] = Wcx

    # attention: context -> hidden (internal dimension of alignment network)
    Wc_att = norm_weight(nin_below, dimctx)
    params[_p(prefix,'Wc_att')] = Wc_att

    # attention: GRU-> hidden (internal dimension of alignment network)
    Wd_att = norm_weight(dim,dimctx)
    params[_p(prefix,'Wd_att')] = Wd_att

    # attention: hidden bias
    b_att = numpy.zeros((dimctx,)).astype('float32')
    params[_p(prefix,'b_att')] = b_att

    # attention: 
    V_att = norm_weight(dimctx,nin_below)
    params[_p(prefix,'V_att')] = V_att
    c_att = numpy.zeros((nin_below,)).astype('float32')
    params[_p(prefix, 'c_att')] = c_att

    return params

def gru_cond_layer(tparams, state_below, options, prefix='gru', 
                    mask=None, context=None, one_step=False, 
                    init_memory=None, init_state=None, 
                    context_mask=None,
                    **kwargs):

    assert context, 'Context must be provided'

    if one_step:
        assert init_state, 'previous state must be provided'

    nsteps = state_below.shape[0]
    if state_below.ndim == 3:
        n_samples = state_below.shape[1]
    else:
        n_samples = 1

    # mask
    if mask == None:
        mask = tensor.alloc(1., state_below.shape[0], 1)

    dim = tparams[_p(prefix, 'U')].shape[0]

    # initial/previous state
    if init_state == None:
        init_state = tensor.alloc(0., n_samples, dim)

    # projected context 
    assert context.ndim == 2, 'Context must be 2-d: #samples x # annotations'
    pctx_ = tensor.dot(context, tparams[_p(prefix,'Wc_att')]) + tparams[_p(prefix,'b_att')]
    
    # drop out Bernoulli vector
    d = tensor.switch(use_noise, trng.binomial(pctx_.shape, n=1, p=prob, dtype=pctx_.dtype), prob)

    # projected x
    state_below_ = tensor.dot(state_below, tparams[_p(prefix, 'W')]) + tparams[_p(prefix, 'b')]
    state_belowx = tensor.dot(state_below, tparams[_p(prefix, 'Wx')]) + tparams[_p(prefix, 'bx')]
        
    def _slice(_x, n, dim):
        if _x.ndim == 3:
            return _x[:, :, n*dim:(n+1)*dim]
        return _x[:, n*dim:(n+1)*dim]

    def _step(m_, x_, xx_, h_, ctx_, pctx_, d):

        # attention
        pstate_ = tensor.dot(h_, tparams[_p(prefix,'Wd_att')])
        pctx__ = pctx_ + pstate_ 
        pctx__ = tensor.tanh(pctx__)
        
        pctx__ = pctx__ * d # dropout here
        
        alpha = tensor.dot(pctx__, tparams[_p(prefix,'V_att')]) + tparams[_p(prefix, 'c_att')]
        
        # softmax like processing
        alpha = tensor.exp(alpha)
        alpha = alpha/alpha.sum(axis=1, keepdims=True) 
       
        if context_mask:
            alpha = alpha * context_mask        
           
       
        ctx_ = context * alpha # current context

        preact = tensor.dot(h_, tparams[_p(prefix, 'U')])
        preact += x_
        preact += tensor.dot(ctx_, tparams[_p(prefix, 'Wc')])

        r = tensor.nnet.sigmoid(_slice(preact, 0, dim))
        u = tensor.nnet.sigmoid(_slice(preact, 1, dim))

        preactx = tensor.dot(h_, tparams[_p(prefix, 'Ux')])
        preactx *= r
        preactx += xx_
        preactx += tensor.dot(ctx_, tparams[_p(prefix, 'Wcx')])

        h = tensor.tanh(preactx)

        h = u * h_ + (1. - u) * h
        h = m_[:,None] * h + (1. - m_)[:,None] * h_

        return h, ctx_, alpha, pstate_, preact, preactx, r, u

    if one_step:
        rval = _step(mask, state_below_, state_belowx, init_state, None, None, pctx_)
    else:
        rval, updates = theano.scan(_step, 
                                    sequences=[mask, state_below_, state_belowx],
                                    outputs_info = [init_state, 
                                                    tensor.alloc(0., n_samples, context.shape[1]),                                                    
                                                    None, None, None, None, None, None],
                                    non_sequences=[pctx_, d],
                                    name=_p(prefix, '_layers'),
                                    n_steps=nsteps)
    return rval

# initialize all parameters
def init_params(options):
    params = OrderedDict()
    
    # decoder
    params = get_layer(options['decoder'])[0](options, params, prefix='decoder', 
                                              nin=options['dim_word'], dim=options['dim'], 
                                              dimctx=options['dim_hid'], nin_below=options['nfeature'])
    # readout
    params = get_layer('ff')[0](options, params, prefix='ff_logit_lstm', nin=options['dim'], nout=options['dim_hid'], ortho=False)
    params = get_layer('ff')[0](options, params, prefix='ff_logit_ctx', nin=options['nfeature'], nout=options['dim_hid'], ortho=False)
    params = get_layer('ff')[0](options, params, prefix='ff_logit', nin=options['dim_hid'], nout=options['n_words'])

    return params

# build a training model
def build_model(corr_file,tparams, options):
      
    # description string: #words x #samples
    x = tensor.matrix('x', dtype='float32')
    x_mask = tensor.matrix('x_mask', dtype='float32')
    y = tensor.matrix('y', dtype='int64')
    y_mask = tensor.matrix('y_mask', dtype='float32')
    

    n_timesteps_trg = y.shape[0]
    n_samples = x.shape[1]

    init_memory = None
    
    # word embedding (target)
    import scipy.io
    from sklearn.decomposition import PCA
    whos_mat_data = scipy.io.whosmat(corr_file)
    matlab_data = scipy.io.loadmat(corr_file)
    correlations = matlab_data[whos_mat_data[0][0]]
    pca = PCA(n_components=options['dim_word'])
    pca.fit(correlations)
    correlations_reduced = pca.transform(correlations)
    n_clusters, dim_reduced = correlations_reduced.shape
    Wemb = numpy.zeros((n_clusters+1, dim_reduced), dtype=numpy.float32)
    Wemb[1:,:] = numpy.array(correlations_reduced, dtype=numpy.float32)
    Wemb_tensor = tensor.constant(Wemb, dtype=numpy.float32)
    
    emb = Wemb_tensor[y.flatten()].reshape([n_timesteps_trg, n_samples, options['dim_word']])
    emb_shifted = tensor.zeros_like(emb)
    emb_shifted = tensor.set_subtensor(emb_shifted[1:], emb[:-1])
    emb = emb_shifted
    
    # decoder
    proj = get_layer(options['decoder'])[1](tparams, emb, options, 
                                            prefix='decoder', 
                                            mask=y_mask, context=x.T, 
                                            context_mask=x_mask.T,
                                            one_step=False, 
                                            init_state=None,
                                            init_memory=init_memory)
    proj_h = proj[0]
    if options['decoder'].startswith('lstm'):
        ctxs = proj[2]
        alphas = proj[3]
    else:
        ctxs = proj[1]
        alphas = proj[2]
    
    proj_h = dropout_layer(proj_h) # Drop out here
    
    # compute word probabilities
    logit_lstm = get_layer('ff')[1](tparams, proj_h, options, prefix='ff_logit_lstm', activ='linear')
    logit_ctx = get_layer('ff')[1](tparams, ctxs, options, prefix='ff_logit_ctx', activ='linear')
    logit = tensor.tanh(logit_lstm + logit_ctx)
    logit = dropout_layer(logit) # Dropout here
    logit = get_layer('ff')[1](tparams, logit, options, prefix='ff_logit', activ='linear')
    logit_shp = logit.shape
    probs = tensor.nnet.softmax(logit.reshape([logit_shp[0]*logit_shp[1], logit_shp[2]]))
    
    # cost
    y_flat = y.flatten()
    y_flat_idx = tensor.arange(y_flat.shape[0]) * options['n_words'] + y_flat
    cost = -tensor.log(probs.flatten()[y_flat_idx]+1e-8)
    cost = cost.reshape([y.shape[0],y.shape[1]])
    cost = (cost * y_mask).sum(0)
    cost = cost.mean()   
    
    return x, x_mask, y, y_mask, alphas, cost


# optimizers
# name(hyperp, tparams, grads, inputs (list), cost) = f_grad_shared, f_update
def adam(lr, tparams, grads, inp, cost):
    gshared = [theano.shared(p.get_value() * 0., name='%s_grad'%k) for k, p in tparams.iteritems()]
    gsup = [(gs, g) for gs, g in zip(gshared, grads)]

    f_grad_shared = theano.function(inp, cost, updates=gsup)

    lr0 = 0.0002
    b1 = 0.1
    b2 = 0.001
    e = 1e-8

    updates = []

    i = theano.shared(numpy.float32(0.))
    i_t = i + 1.
    fix1 = 1. - b1**(i_t)
    fix2 = 1. - b2**(i_t)
    lr_t = lr0 * (tensor.sqrt(fix2) / fix1)

    for p, g in zip(tparams.values(), gshared):
        m = theano.shared(p.get_value() * 0.)
        v = theano.shared(p.get_value() * 0.)
        m_t = (b1 * g) + ((1. - b1) * m)
        v_t = (b2 * tensor.sqr(g)) + ((1. - b2) * v)
        g_t = m_t / (tensor.sqrt(v_t) + e)
        p_t = p - (lr_t * g_t)
        updates.append((m, m_t))
        updates.append((v, v_t))
        updates.append((p, p_t))
    updates.append((i, i_t))

    f_update = theano.function([lr], [], updates=updates, on_unused_input='ignore')

    return f_grad_shared, f_update

def adadelta(lr, tparams, grads, inp, cost):
    zipped_grads = [theano.shared(p.get_value() * numpy.float32(0.), name='%s_grad'%k) for k, p in tparams.iteritems()]
    running_up2 = [theano.shared(p.get_value() * numpy.float32(0.), name='%s_rup2'%k) for k, p in tparams.iteritems()]
    running_grads2 = [theano.shared(p.get_value() * numpy.float32(0.), name='%s_rgrad2'%k) for k, p in tparams.iteritems()]

    zgup = [(zg, g) for zg, g in zip(zipped_grads, grads)]
    rg2up = [(rg2, 0.95 * rg2 + 0.05 * (g ** 2)) for rg2, g in zip(running_grads2, grads)]

    f_grad_shared = theano.function(inp, cost, updates=zgup+rg2up)
    
    updir = [-tensor.sqrt(ru2 + 1e-6) / tensor.sqrt(rg2 + 1e-6) * zg for zg, ru2, rg2 in zip(zipped_grads, running_up2, running_grads2)]
    ru2up = [(ru2, 0.95 * ru2 + 0.05 * (ud ** 2)) for ru2, ud in zip(running_up2, updir)]
    param_up = [(p, p + ud) for p, ud in zip(itemlist(tparams), updir)]

    f_update = theano.function([lr], [], updates=ru2up+param_up, on_unused_input='ignore')

    return f_grad_shared, f_update

def rmsprop(lr, tparams, grads, inp, cost):
    zipped_grads = [theano.shared(p.get_value() * numpy.float32(0.), name='%s_grad'%k) for k, p in tparams.iteritems()]
    running_grads = [theano.shared(p.get_value() * numpy.float32(0.), name='%s_rgrad'%k) for k, p in tparams.iteritems()]
    running_grads2 = [theano.shared(p.get_value() * numpy.float32(0.), name='%s_rgrad2'%k) for k, p in tparams.iteritems()]

    zgup = [(zg, g) for zg, g in zip(zipped_grads, grads)]
    rgup = [(rg, 0.95 * rg + 0.05 * g) for rg, g in zip(running_grads, grads)]
    rg2up = [(rg2, 0.95 * rg2 + 0.05 * (g ** 2)) for rg2, g in zip(running_grads2, grads)]

    f_grad_shared = theano.function(inp, cost, updates=zgup+rgup+rg2up)

    updir = [theano.shared(p.get_value() * numpy.float32(0.), name='%s_updir'%k) for k, p in tparams.iteritems()]
    updir_new = [(ud, 0.9 * ud - 1e-4 * zg / tensor.sqrt(rg2 - rg ** 2 + 1e-4)) for ud, zg, rg, rg2 in zip(updir, zipped_grads, running_grads, running_grads2)]
    param_up = [(p, p + udn[1]) for p, udn in zip(itemlist(tparams), updir_new)]
    f_update = theano.function([lr], [], updates=updir_new+param_up, on_unused_input='ignore')

    return f_grad_shared, f_update

def sgd(lr, tparams, grads, x, mask, y, cost):
    gshared = [theano.shared(p.get_value() * 0., name='%s_grad'%k) for k, p in tparams.iteritems()]
    gsup = [(gs, g) for gs, g in zip(gshared, grads)]

    f_grad_shared = theano.function([x, mask, y], cost, updates=gsup)

    pup = [(p, p - lr * g) for p, g in zip(itemlist(tparams), gshared)]
    f_update = theano.function([lr], [], updates=pup)

    return f_grad_shared, f_update
def read_data(smri_fname,fmri_fname):
    fMRI_file = open(fmri_fname, 'r')
    sMRI_file = open(smri_fname, 'r')   
    
    data_fMRI = csv.reader(fMRI_file, delimiter=',')
    data_sMRI = csv.reader(sMRI_file, delimiter=',')
    table_fMRI = [row for row in data_fMRI]
    table_sMRI = [row for row in data_sMRI]
    
    fMRI_sentences = [map(int, s) for s in table_fMRI]
    sMRI_sentences = [map(float, s) for s in table_sMRI]
    fMRI_sentences = [numpy.asarray(arr, dtype=numpy.int64) for arr in fMRI_sentences]
    sMRI_sentences = [numpy.asarray(arr, dtype=numpy.float32) for arr in sMRI_sentences]
    
    assert len(fMRI_sentences)==len(sMRI_sentences), 'Number of subjects should be same in both modality'        
    return zip(sMRI_sentences, fMRI_sentences) # structural data as x and functional as y
 
def generate_batch(data, batch_size=1):
      
    sMRI_sentences, fMRI_sentences = zip(*data)
    
    sMRI_sentences = list(sMRI_sentences)
    fMRI_sentences = list(fMRI_sentences)
    
    nbatches = len(fMRI_sentences)//batch_size
    rem = len(fMRI_sentences)%batch_size
    pos = 0
    fMRI = []
    sMRI = []
    for _ in range(nbatches):
        fMRI.append(fMRI_sentences[pos:pos+batch_size])
        sMRI.append(sMRI_sentences[pos:pos+batch_size])
        pos = pos + batch_size
    
    if rem != 0:
        fMRI.append(fMRI_sentences[pos:])
        sMRI.append(sMRI_sentences[pos:])
    
    data = zip(sMRI, fMRI)    
          
    return data
    
#     return data
def prepare_data(seq_x, seq_y):
    # check for number of examples in input and output sequence
    assert len(seq_x) == len(seq_y), 'Number of samples in input and output should be same'
    n_samples = len(seq_x)
    
    lengths_x = [len(s) for s in seq_x]
    lengths_y = [len(s) for s in seq_y]    
    lengths_x = numpy.array(lengths_x).astype(numpy.int64)
    lengths_y = numpy.array(lengths_y).astype(numpy.int64)
    
    maxlen_x = lengths_x.max()
    maxlen_y = lengths_y.max()    
    
    assert all(lengths_x == maxlen_x), 'All sequences should have same length'
    assert all(lengths_y == maxlen_y), 'All sequences should have same length'
    
    x = numpy.zeros((maxlen_x, n_samples)).astype(numpy.float32)
    x_mask = numpy.ones_like(x, dtype=numpy.float32)
    y = numpy.zeros((maxlen_y, n_samples)).astype(numpy.int64)
    y_mask = numpy.ones_like(y, dtype=numpy.float32)
    
    # Fill matrix x and y, #n_timesteps x #n_samples
    for i in range(n_samples):
        x[:,i] = seq_x[i]
        y[:,i] = seq_y[i]
        
    return x, x_mask, y, y_mask      

def train(smri_file,fmri_file,corr_file_name,temp_dir="Results\\",dim_word=4, # word vector dimensionality
          dim=50, # the number of GRU units
          dim_hid=50,
          decoder='gru_cond',
          patience=10,
          max_epochs=100,
          dispFreq=202,
          lrate=0.01,         
          n_words=6,
          maxlen=150, # maximum length of the description
          optimizer='rmsprop', 
          batch_size = 16,
          valid_batch_size = 1,
          num_runs=100, reload_=False):

    # Model options
    model_options = locals().copy()

    print 'Loading data...'
    data = read_data(smri_file,fmri_file)
    print 'Done'
    model_options['nfeature'] = len(data[0][0])

    print 'Building model...'
    params = init_params(model_options)
    tparams = init_tparams(params)
    x, x_mask, y, y_mask, \
          alphas, \
          cost = \
          build_model(corr_file_name, tparams, model_options)
    print 'Done'
    
    inps = [x, x_mask, y, y_mask]

    # Regularizer
    cost = cost - tensor.alloc(0.5)*((tensor.sum(alphas**2, axis=2)).mean(0)).mean()
    
    print 'Building Alignment...',
    get_alignment = theano.function(inps, alphas)
    print 'Done'
    
    print 'Computing gradient...',
    grads = tensor.grad(cost, wrt=itemlist(tparams))
    print 'Done'
  
    lr = tensor.scalar(name='lr')
    print 'Building optimizers...',
    f_grad_shared, f_update = eval(optimizer)(lr, tparams, grads, inps, cost)
    print 'Done'
  
    train = generate_batch(data, batch_size)       
    
    for run_no in range(num_runs):
        print 'Processing %d-th run...' %(run_no + 1)
        params = init_params(model_options)
        zipp(params, tparams)
        uidx = 0
        for eidx in xrange(max_epochs):
            n_samples = 0
            train_cost = 0.        
            use_noise.set_value(1)
            for x, y in train:
                n_samples += len(x)
                uidx += 1       
                x, x_mask, y, y_mask = prepare_data(x, y)            
       
                cost1 = f_grad_shared(x, x_mask, y, y_mask)            
                f_update(lrate)
                train_cost += cost1      
            
            print 'Epoch ', eidx, 'Update ', uidx, 'Cost ', train_cost/len(train)
            
            
        if (not (os.path.exists(temp_dir))):
            os.mkdir(temp_dir)			
	saveto = os.path.join(temp_dir, 'model_run%d.txt' %(run_no + 1)) #'Results/model_run%d.txt' %(run_no + 1)    
        params = unzip(tparams)
        numpy.savez(saveto, **params)
        
        use_noise.set_value(0)
        a_file = open(os.path.join(temp_dir, 'alignment_run%d.txt' %(run_no + 1)), 'w')
        for x,y in data:
            x, x_mask, y, y_mask = prepare_data([x], [y])
            alignment = get_alignment(x, x_mask, y, y_mask)
            alignment = alignment.squeeze()
            a_file.write(str(alignment.shape)+'\n')
            numpy.savetxt(a_file, alignment, fmt = '%.4f')
            #a_file.write('\n')
            
        a_file.close()
        print 'Files written'

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Multi-modality application. Inputs are smri loadings and dfnc information (correlation and centroids information')
    parser.add_argument('--smri',type=str)
    parser.add_argument('--dfnc',type=str)
    parser.add_argument('--states',type=str)
    parser.add_argument('--outdir', type=str, default="results\\")
    parser.add_argument('--dim_word', type=int, default=4)
    parser.add_argument('--dim', type=int, default=50)
    parser.add_argument('--dim_hid', type=int, default=50)
    parser.add_argument('--patience', type=int, default=10)
    parser.add_argument('--dispFreq', type=int, default=202)
    parser.add_argument('--lrate', type=float, default=0.01)
    parser.add_argument('--n_words', type=int, default=6)
    parser.add_argument('--maxlen', type=int, default=150)
    parser.add_argument('--max_epochs', type=int, default=100)
    parser.add_argument('--num_runs', type=int, default=100)

    args = parser.parse_args()
    
    # call method
    train(args.smri, args.dfnc, args.states, temp_dir=args.outdir, dim_word=args.dim_word, dim=args.dim, dim_hid=args.dim_hid, patience=args.patience, dispFreq=args.dispFreq, \
          n_words=args.n_words, maxlen=args.maxlen, max_epochs=args.max_epochs, num_runs=args.num_runs)
