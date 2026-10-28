"""A fast training step of DeepGO-SE (model.DeepGOModel): the same loss, gradients and Adam update as the official
PyTorch step, computed by Triton kernels and replayed as a CUDA graph. Needs a CUDA GPU and Triton
(which needs the Python development headers, Python.h, on CPATH or the default include path).

One step: BCE on the batch plus all four EL losses over every normal form, with Adam (lr 5e-4, betas 0.9/0.999, eps 1e-8,
no weight decay), in FP32. Classifier and MLP gradients are written out by hand; the EL loss and its gradient through
the nine BatchNorm populations (one per class column of NF1-NF4, as the official el_loss calls go_norm nine times) are
computed by Triton kernels from an edge cache; the class table gets a fused Adam update, every other parameter torch
Adam(capturable=True). Full batches replay a CUDA graph (StreamingStep.capture); the caller recaptures when the learning
rate changes and runs partial batches eagerly (StreamingStep._step).

Checked against the official PyTorch step on five worlds trained both ways: development micro-AUPRC within -0.011 to
+0.013 per world and 0.002 on average (ordinary world-to-world variation); runs are byte-reproducible.
"""
import copy

import numpy as np
import torch
import triton as tr
import triton.language as tl


# ---- the constraint graph
@tr.jit
def _radius_backward(Active,Edges,Signs,Ptr,Rad,Grad,V:tl.constexpr,B:tl.constexpr=128):
    i=tl.program_id(0);start=tl.load(Ptr+i);end=tl.load(Ptr+i+1)
    acc=tl.full((B,),0.,tl.float32)
    for p in range(start,end,B):
        q=p+tl.arange(0,B);e=tl.load(Edges+q,q<end,0);s=tl.load(Signs+q,q<end,0)
        acc+=tl.load(Active+e)*s
    r=tl.load(Rad+i);sgn=tl.where(r>0,1.,tl.where(r<0,-1.,0.))
    tl.store(Grad+i,tl.sum(acc,0)*sgn)


def tensor(a,device,dtype=None):
    return torch.as_tensor(np.ascontiguousarray(a),device=device,dtype=dtype)


def routes(rows,edges,signs,n,device):
    order=np.argsort(rows,kind='stable');rows=rows[order]
    edges=edges[order];signs=signs[order]
    ptr=np.r_[0,np.cumsum(np.bincount(rows,minlength=n))]
    starts=[];ends=[];chunks=[0];dest=[];multirows=[];multiptr=[0]
    for row,(a,b) in enumerate(zip(ptr[:-1],ptr[1:])):
        k=(b-a+15)//16
        if k>1:multirows.append(row)
        for j in range(a,b,16):
            starts.append(j);ends.append(min(j+16,b))
            dest.append(multiptr[-1]+(j-a)//16 if k>1 else -row-1)
        if k>1:multiptr.append(multiptr[-1]+k)
        chunks.append(len(starts))
    return dict(edges=tensor(edges,device,torch.int32),signs=tensor(signs,device,torch.float32),
        starts=tensor(starts,device,torch.int32),ends=tensor(ends,device,torch.int32),
        chunks=tensor(chunks,device,torch.int32),ptr=tensor(ptr,device,torch.int32),n=n,
        dest=tensor(dest,device,torch.int32),multirows=tensor(multirows,device,torch.int32),
        multiptr=tensor(multiptr,device,torch.int32),npartial=multiptr[-1])


class GraphPlan:
    def __init__(self,model,forms):
        self.device=model.go_embed.weight.device;self.d=model.embed_dim
        self.v=model.go_embed.num_embeddings;self.nr=model.rel_embed.num_embeddings
        a,b,c,d=[np.asarray(x.cpu() if torch.is_tensor(x) else x,dtype=np.int64) for x in forms]
        assert all(len(x)>1 for x in (a,b,c,d))
        cols=[a[:,0],a[:,1],b[:,0],b[:,1],b[:,2],c[:,1],c[:,2],d[:,0],d[:,2]]
        ids=[];counts=[];inverse=[];self.offset=[0];role=[]
        for s,col in enumerate(cols):
            u,inv,ct=np.unique(col,return_inverse=True,return_counts=True)
            inverse.append(inv+self.offset[-1]);ids.extend(u);counts.extend(ct);role.extend([s]*len(u));self.offset.append(len(ids))
        self.u=len(ids);self.ids=tensor(ids,self.device,torch.int64)
        self.counts=tensor(counts,self.device,torch.float32);self.role=tensor(role,self.device,torch.int32)
        self.ns=tensor([len(x) for x in cols],self.device,torch.float32)
        classrows=np.full((self.v,9),-1,np.int32)
        for j,(i,s) in enumerate(zip(ids,role)):classrows[i,s]=j
        self.classrows=tensor(classrows,self.device)
        aa=[];bb=[];rr=[];rs=[];ra=[];rb=[];sa=[];sb=[];scale=[]
        def add(x,y,r,sgn,r1,r2,s1,s2,n):
            aa.extend(x);bb.extend(y);rr.extend(np.full(n,-1) if r is None else r);rs.extend([sgn]*n)
            ra.extend(r1);rb.extend(r2);sa.extend([s1]*n);sb.extend([s2]*n);scale.extend([1/n]*n)
        i=inverse
        add(i[0],i[1],None,0,a[:,0],a[:,1],1,-1,len(a))
        add(i[2],i[3],None,0,b[:,0],b[:,1],-1,-1,len(b))
        add(i[4],i[2],None,0,b[:,0],b[:,1],-1,0,len(b))
        add(i[4],i[3],None,0,b[:,0],b[:,1],0,-1,len(b))
        add(i[5],i[6],c[:,0],1,c[:,1],c[:,2],1,-1,len(c))
        # Preserve the author's relation-column radius lookup in NF4.
        add(i[7],i[8],d[:,1],-1,d[:,1],d[:,2],-1,-1,len(d))
        self.e=len(aa);ei=np.arange(self.e,dtype=np.int32)
        self.edge=[tensor(q,self.device,torch.int32) for q in (aa,bb,rr,rs,ra,rb)]
        self.edge += [tensor(q,self.device,torch.float32) for q in (sa,sb,scale)]
        self.zroutes=routes(np.r_[aa,bb],np.r_[ei,ei],np.r_[np.ones(self.e),-np.ones(self.e)],self.u,self.device)
        rr=np.array(rr);rs=np.array(rs);mask=rr>=0
        self.rroutes=routes(rr[mask],ei[mask],rs[mask],self.nr,self.device)
        self.radroutes=routes(np.r_[ra,rb],np.r_[ei,ei],np.r_[sa,sb],self.v,self.device)
        self.margin=model.margin;self.bn=model.go_norm
        assert self.bn.affine and self.bn.momentum is not None
        self.summary=dict(classes=self.v,dimensions=self.d,role_occurrences=sum(map(len,cols)),
            role_unique=self.u,distance_edges=self.e,gradient_chunks=len(self.zroutes['starts']))


# ---- step operations

@tr.jit
def _bn_running(Mu,Var,N,RunningMu,RunningVar,Tracked,D:tl.constexpr,MOMENTUM:tl.constexpr,B:tl.constexpr=256):
    d=tl.program_id(0)*B+tl.arange(0,B);mask=d<D
    m=tl.load(RunningMu+d,mask,0);v=tl.load(RunningVar+d,mask,0)
    for s in tl.static_range(9):
        n=tl.load(N+s)
        m=(1-MOMENTUM)*m+MOMENTUM*tl.load(Mu+s*D+d,mask,0)
        v=(1-MOMENTUM)*v+MOMENTUM*tl.load(Var+s*D+d,mask,0)*n/(n-1)
    tl.store(RunningMu+d,m,mask);tl.store(RunningVar+d,v,mask)
    if tl.program_id(0)==0:tl.store(Tracked,tl.load(Tracked)+9)


@tr.jit
def _bias_parameters(Step,Bias,LR):
    step=tl.load(Step)+1;tl.store(Step,step)
    n=step.to(tl.float32)
    tl.store(Bias,LR/(1-tl.exp(n*(-0.1053605156578263))))
    tl.store(Bias+1,tl.sqrt(1-tl.exp(n*(-0.0010005003335835344))))


@tr.jit
def _bce(Logits,Offset,Rad,Y,Grad,Part,BATCH:tl.constexpr,K:tl.constexpr,B:tl.constexpr=256):
    idx=tl.program_id(0)*B+tl.arange(0,B);valid=idx<BATCH*K;col=idx%K;row=idx//K
    z=tl.load(Logits+idx,valid,0)+tl.load(Offset+row,valid,0)+tl.abs(tl.load(Rad+col,valid,0))
    y=tl.load(Y+idx,valid,0);p=tl.sigmoid(z)
    loss=-y*tl.maximum(tl.log(p),-100.)-(1-y)*tl.maximum(tl.log(1-p),-100.)
    a=p*(1-p)
    grad=(p-y)*a/tl.maximum(a,1.e-12)/(BATCH*K)
    tl.store(Grad+idx,grad,valid)
    tl.store(Part+tl.program_id(0),tl.sum(tl.where(valid,loss,0.),0)/(BATCH*K))


@torch.no_grad()
def _classifier_gradients(h,w,rel,rad,y):
    logits=h@w.T;offset=h@rel.T;g=torch.empty_like(logits)
    part=torch.empty(tr.cdiv(logits.numel(),256),device=h.device,dtype=h.dtype)
    _bce[(len(part),)](logits,offset,rad,y,g,part,len(h),len(w))
    rows=g.sum(1,keepdim=True)
    return part.sum(),g@w+rows*rel,g.T@h,rows.T@h,g.sum(0).view(-1,1)*rad.sign()


def _block_forward(block,x):
    if not isinstance(block.activation,torch.nn.ReLU) or block.layer_norm is None or block.dropout is None:
        raise ValueError('Native ReLU / LayerNorm / dropout block required')
    linear=block.linear;norm=block.layer_norm
    pre=torch.addmm(linear.bias,x,linear.weight.T)
    rectified=torch.relu(pre)
    value,mu,inv=torch.ops.aten.native_layer_norm(rectified,list(norm.normalized_shape),norm.weight,norm.bias,norm.eps)
    out,mask=torch.ops.aten.native_dropout(value,block.dropout.p,True)
    return out,(x,pre,rectified,mu,inv,mask)


def _block_backward(block,cache,grad,input_gradient=True):
    x,pre,rectified,mu,inv,mask=cache;norm=block.layer_norm;linear=block.linear
    drop=torch.ops.aten.native_dropout_backward(grad,mask,1./(1.-block.dropout.p))
    g,norm.weight.grad,norm.bias.grad=torch.ops.aten.native_layer_norm_backward(
        drop,rectified,list(norm.normalized_shape),mu,inv,norm.weight,norm.bias,[True,True,True])
    g=torch.where(pre>0,g,0.)
    linear.weight.grad=g.T@x;linear.bias.grad=g.sum(0)
    return g@linear.weight if input_gradient else None


def capture_step(self,x,y):
    if self.graph is not None:raise RuntimeError('Step already captured')
    state={k:v.detach().clone() for k,v in self.model.state_dict().items()}
    m=self.m.clone();v=self.v.clone();counter=self.step_tensor.clone()
    other=copy.deepcopy(self.other.state_dict())
    cpu_rng=torch.get_rng_state();cuda_rng=torch.cuda.get_rng_state()
    self.graph_x=x.clone();self.graph_y=y.clone()
    stream=torch.cuda.Stream();stream.wait_stream(torch.cuda.current_stream())
    with torch.cuda.stream(stream):
        for _ in range(3):self._step(self.graph_x,self.graph_y)
    torch.cuda.current_stream().wait_stream(stream);torch.cuda.synchronize()
    self.graph=torch.cuda.CUDAGraph()
    with torch.cuda.graph(self.graph,stream=stream):self.graph_loss=self._step(self.graph_x,self.graph_y)
    self.model.load_state_dict(state);self.m.copy_(m);self.v.copy_(v);self.step_tensor.copy_(counter)
    # Preserve captured addresses when restoring Adam state after warmup.
    if other['state']:
        for param,key in zip(self.other.param_groups[0]['params'],other['param_groups'][0]['params']):
            old=other['state'][key]
            for key,value in old.items():self.other.state[param][key].copy_(value)
    else:
        for entry in self.other.state.values():
            for value in entry.values():
                if torch.is_tensor(value):value.zero_()
    self._initial_stats();torch.set_rng_state(cpu_rng);torch.cuda.set_rng_state(cuda_rng)

# ---- the streaming step

@tr.jit
def _shared_moments(W, C, Mu, Part, D: tl.constexpr, V: tl.constexpr,
                    NB: tl.constexpr, CENTER: tl.constexpr,
                    BR: tl.constexpr = 128, BD: tl.constexpr = 128):
    block = tl.program_id(0)
    rows = block * BR + tl.arange(0, BR)
    d = tl.program_id(1) * BD + tl.arange(0, BD)
    x = tl.load(W + rows[:, None] * D + d[None, :],
                (rows[:, None] < V) & (d[None, :] < D), 0)
    # One class-table read serves all nine occurrence populations.
    for s in tl.static_range(9):
        c = tl.load(C + rows * 9 + s, rows < V, 0)
        z = x
        if CENTER:
            delta = x - tl.load(Mu + s * D + d, d < D, 0)[None, :]
            z = delta * delta
        val = tl.sum(z * c[:, None], 0)
        tl.store(Part + (s * NB + block) * D + d, val, d < D)


@tr.jit
def _edge_values(W, R, Mu, Inv, Gamma, Beta, A, B, AS, BS, Rel, RS,
                 e, d, D: tl.constexpr):
    a = tl.load(A + e); b = tl.load(B + e)
    sa = tl.load(AS + e); sb = tl.load(BS + e)
    ha = (tl.load(W + a * D + d, d < D, 0)
          - tl.load(Mu + sa * D + d, d < D, 0)) * tl.load(Inv + sa * D + d, d < D, 0)
    hb = (tl.load(W + b * D + d, d < D, 0)
          - tl.load(Mu + sb * D + d, d < D, 0)) * tl.load(Inv + sb * D + d, d < D, 0)
    gamma = tl.load(Gamma + d, d < D, 0); beta = tl.load(Beta + d, d < D, 0)
    za = gamma * ha + beta; zb = gamma * hb + beta
    r = tl.load(Rel + e); rs = tl.load(RS + e)
    rel = tl.load(R + tl.maximum(r, 0) * D + d, (d < D) & (r >= 0), 0)
    # Preserve the native NF3 and NF4 parentheses, not an algebraic reassociation.
    v = tl.where(rs > 0, (za + rel) - zb, tl.where(rs < 0, za - (zb + rel), za - zb))
    return v, ha, hb


@tr.jit
def _edge_coefficients(W, R, Rad, Mu, Inv, Gamma, Beta, A, B, AS, BS, Rel, RS,
                 RA, RB, SA, SB, Scale, Loss, Active, Coef, D: tl.constexpr,
                 MARGIN: tl.constexpr, BD: tl.constexpr):
    e = tl.program_id(0); d = tl.arange(0, BD)
    v, ha, hb = _edge_values(W, R, Mu, Inv, Gamma, Beta, A, B, AS, BS, Rel, RS, e, d, D)
    norm = tl.sqrt(tl.sum(v * v, 0))
    ra = tl.load(RA + e); rb = tl.load(RB + e)
    hinge = norm + tl.load(SA + e) * tl.abs(tl.load(Rad + ra)) + tl.load(SB + e) * tl.abs(tl.load(Rad + rb)) - MARGIN
    scale = tl.load(Scale + e)
    active = tl.where(hinge > 0, scale, 0.)
    coef = tl.where(norm > 0, active / tl.maximum(norm, 1.e-30), 0.)
    tl.store(Loss + e, tl.maximum(hinge, 0.) * scale)
    tl.store(Active + e, active)
    tl.store(Coef+e,coef)


@tr.jit
def _constraint_summaries(W,R,Mu,Inv,Gamma,Beta,A,B,AS,BS,Rel,RS,Coef,EdgeG,
                          Edges,Starts,Ends,Summary,D:tl.constexpr,BD:tl.constexpr=128,BR:tl.constexpr=4):
    group = tl.program_id(0); d = tl.program_id(1) * BD + tl.arange(0, BD)
    start = tl.load(Starts + group); end = tl.load(Ends + group)
    first=tl.load(Edges+start);sa=tl.load(AS+first);sb=tl.load(BS+first)
    r=tl.load(Rel+first);rs=tl.load(RS+first)
    ma=tl.load(Mu+sa*D+d,d<D,0);ia=tl.load(Inv+sa*D+d,d<D,0)
    mb=tl.load(Mu+sb*D+d,d<D,0);ib=tl.load(Inv+sb*D+d,d<D,0)
    gamma=tl.load(Gamma+d,d<D,0);beta=tl.load(Beta+d,d<D,0)
    re=tl.load(R+tl.maximum(r,0)*D+d,(d<D)&(r>=0),0)
    ga=tl.full((BR,BD),0.,tl.float32);ha_sum=tl.full((BR,BD),0.,tl.float32);hb_sum=tl.full((BR,BD),0.,tl.float32)
    for begin in range(start,end,BR):
        pos=begin+tl.arange(0,BR);valid=pos<end;e=tl.load(Edges+pos,valid,0)
        a=tl.load(A+e,valid,0);b=tl.load(B+e,valid,0);mask=valid[:,None]&(d[None,:]<D)
        ha=(tl.load(W+a[:,None]*D+d[None,:],mask,0)-ma[None,:])*ia[None,:]
        hb=(tl.load(W+b[:,None]*D+d[None,:],mask,0)-mb[None,:])*ib[None,:]
        za=gamma[None,:]*ha+beta[None,:];zb=gamma[None,:]*hb+beta[None,:]
        v=tl.where(rs>0,(za+re[None,:])-zb,tl.where(rs<0,za-(zb+re[None,:]),za-zb))
        u=tl.load(Coef+e,valid,0)[:,None]*v
        tl.store(EdgeG+e[:,None]*D+d[None,:],u,mask)
        ga+=u;ha_sum+=u*ha;hb_sum-=u*hb
    g=tl.sum(ga,0)
    tl.store(Summary+(group*5+0)*D+d,g,d<D)
    tl.store(Summary+(group*5+1)*D+d,tl.sum(ha_sum,0),d<D)
    tl.store(Summary+(group*5+2)*D+d,-g,d<D)
    tl.store(Summary+(group*5+3)*D+d,tl.sum(hb_sum,0),d<D)
    tl.store(Summary+(group*5+4)*D+d,rs*g,d<D)


@tr.jit
def _summary_reduce(Summary, Slots, Ptr, Out, D: tl.constexpr,
                    BR: tl.constexpr = 32, BD: tl.constexpr = 128):
    row = tl.program_id(0); d = tl.program_id(1) * BD + tl.arange(0, BD)
    start = tl.load(Ptr + row); end = tl.load(Ptr + row + 1)
    acc = tl.full((BR, BD), 0., tl.float32)
    for p in range(start, end, BR):
        pos = p + tl.arange(0, BR)
        slot = tl.load(Slots + pos, pos < end, 0)
        acc += tl.load(Summary + slot[:, None] * D + d[None, :],
                       (pos[:, None] < end) & (d[None, :] < D), 0)
    tl.store(Out + row * D + d, tl.sum(acc, 0), d < D)


@tr.jit
def _incidences(EdgeG, Edges, Roles, Signs, K, start, end, d,
                D: tl.constexpr, BD: tl.constexpr):
    acc = tl.full((BD,), 0., tl.float32)
    for pos in range(start, end):
        e = tl.load(Edges + pos); role = tl.load(Roles + pos)
        u = tl.load(EdgeG + e * D + d, d < D, 0)
        acc += tl.load(Signs + pos) * u * tl.load(K + role * D + d, d < D, 0)
    return acc


@tr.jit
def _hub_partials(EdgeG, Edges, Roles, Signs, K, Starts, Ends, Part,
                  D: tl.constexpr, BD: tl.constexpr):
    chunk = tl.program_id(0); d = tl.program_id(1) * BD + tl.arange(0, BD)
    g = _incidences(EdgeG, Edges, Roles, Signs, K, tl.load(Starts + chunk),
                    tl.load(Ends + chunk), d, D, BD)
    tl.store(Part + chunk * D + d, g, d < D)


@tr.jit
def _stream_adam(W, New, EdgeG, Edges, Roles, Signs, Ptr, HubPtr, Part,
                 C, K, AA, BB, BC, M, V, STEP_SIZE, BC2,
                 D: tl.constexpr, NPRED: tl.constexpr, BD: tl.constexpr):
    i = tl.program_id(0); d = tl.program_id(1) * BD + tl.arange(0, BD); mask = d < D
    step_size=tl.load(STEP_SIZE);bc2=tl.load(BC2)
    lo = tl.load(HubPtr + i); hi = tl.load(HubPtr + i + 1)
    g = tl.full((BD,), 0., tl.float32)
    if hi > lo:
        for part in range(lo, hi):
            g += tl.load(Part + part * D + d, mask, 0)
    else:
        g = _incidences(EdgeG, Edges, Roles, Signs, K, tl.load(Ptr + i),
                        tl.load(Ptr + i + 1), d, D, BD)
    x = tl.load(W + i * D + d, mask, 0)
    g += tl.load(BC + i * D + d, mask & (i < NPRED), 0)
    for s in tl.static_range(9):
        c = tl.load(C + i * 9 + s)
        g += c * (x * tl.load(AA + s * D + d, mask, 0) + tl.load(BB + s * D + d, mask, 0))
    m = tl.load(M + i * D + d, mask, 0); v = tl.load(V + i * D + d, mask, 0)
    m += .1 * (g - m); v = .999 * v + .001 * g * g
    new = x - step_size * m / (tl.sqrt(v) / bc2 + 1.e-8)
    tl.store(M + i * D + d, m, mask); tl.store(V + i * D + d, v, mask)
    # All cross-class reads use the completed edge cache. In-place class writes
    # are safe: this kernel reads raw W only for its own row.
    tl.store(New + i * D + d, new, mask)


class StreamingPlan:
    """Compile fixed constraint routing once; no role-by-dimension tables."""

    def __init__(self, base):
        self.base = base
        device = base.device
        ids = base.ids.cpu().numpy(); role = base.role.cpu().numpy()
        ae, be, rel, rs = [x.cpu().numpy() for x in base.edge[:4]]
        a, b, sa, sb = ids[ae], ids[be], role[ae], role[be]
        self.endpoints = [tensor(x, device, torch.int32) for x in (a, b, sa, sb)]
        counts = np.zeros((base.v, 9), np.float32)
        counts[ids, role] = base.counts.cpu().numpy()
        self.counts = tensor(counts, device)

        order = np.lexsort((rs, rel, sb, sa)).astype(np.int32)
        keys = np.stack((sa[order], sb[order], rel[order], rs[order]), axis=1)
        boundaries = np.r_[0, np.flatnonzero(np.any(keys[1:] != keys[:-1], axis=1)) + 1, base.e]
        starts, ends, groups = [], [], []
        for lo, hi in zip(boundaries[:-1], boundaries[1:]):
            for pos in range(int(lo), int(hi), 64):
                starts.append(pos); ends.append(min(pos + 64, int(hi)))
                groups.append(keys[pos])
        self.groups = len(starts)
        self.group_edges, self.starts, self.ends = [tensor(x, device, torch.int32) for x in (order, starts, ends)]

        # Rows 0:9 = S1, 9:18 = S2, 18: = relation gradient.
        slots = [[] for _ in range(18 + base.nr)]
        for j, (s, t, r, _) in enumerate(groups):
            slots[s].append(j * 5); slots[s + 9].append(j * 5 + 1)
            slots[t].append(j * 5 + 2); slots[t + 9].append(j * 5 + 3)
            if r >= 0:
                slots[18 + r].append(j * 5 + 4)
        self.summary_slots = tensor([v for row in slots for v in row], device, torch.int32)
        self.summary_ptr = tensor(np.r_[0, np.cumsum([len(row) for row in slots])], device, torch.int32)

        # Keep two incidences even if a constraint joins a class to itself:
        # its two endpoints may belong to different normalization populations.
        rows = np.r_[a, b]; order = np.argsort(rows, kind='stable')
        ei = np.arange(base.e, dtype=np.int32)
        self.edges = tensor(np.r_[ei, ei][order], device)
        self.roles = tensor(np.r_[sa, sb][order], device, torch.int32)
        self.signs = tensor(np.r_[np.ones(base.e), -np.ones(base.e)][order], device, torch.float32)
        ptr = np.r_[0, np.cumsum(np.bincount(rows, minlength=base.v))]
        self.ptr = tensor(ptr, device, torch.int32)
        hub_starts, hub_ends, hub_ptr = [], [], [0]
        for lo, hi in zip(ptr[:-1], ptr[1:]):
            if hi - lo > 32:
                for pos in range(int(lo), int(hi), 32):
                    hub_starts.append(pos); hub_ends.append(min(pos + 32, int(hi)))
            hub_ptr.append(len(hub_starts))
        self.hubs = len(hub_starts)
        self.hub_starts, self.hub_ends, self.hub_ptr = [tensor(x, device, torch.int32) for x in (hub_starts, hub_ends, hub_ptr)]


class StreamingStep:
    """Full BCE + ontology + Adam update; exclusive owner of model optimization.

    Do not attach another optimizer to this model. The class Parameter object is
    retained and updated in place. Save model.state_dict()
    together with this object's state_dict() to resume training.
    """

    def __init__(self, model, plan, lr=5.e-4):
        w = model.go_embed.weight
        if not w.is_cuda or w.dtype != torch.float32 or not w.is_contiguous():
            raise ValueError('StreamingStep requires contiguous CUDA FP32 parameters')
        if any(p.dtype != torch.float32 or p.device != w.device for p in model.parameters()):
            raise ValueError('All parameters must share the class table device and FP32 dtype')
        if not plan.bn.track_running_stats or plan.bn.momentum is None:
            raise ValueError('Expected native BatchNorm running statistics and fixed momentum')
        self.model = model; self.plan = StreamingPlan(plan); self.lr = lr; self.graph=None
        self.m = torch.zeros_like(w); self.v = torch.zeros_like(w)
        self.step_tensor=torch.zeros((),device=w.device,dtype=torch.int64)
        self.bias=torch.empty(2,device=w.device,dtype=w.dtype)
        self.other = torch.optim.Adam([p for p in model.parameters() if p is not w], lr=lr, weight_decay=0,capturable=True)
        p = self.plan; d = plan.d; kw = dict(device=w.device, dtype=w.dtype)
        self.moment_parts = torch.empty((9, tr.cdiv(plan.v, 128), d), **kw)
        self.summaries = torch.empty((p.groups, 5, d), **kw)
        self.reduced = torch.empty((18 + plan.nr, d), **kw)
        self.hub_parts = torch.empty((p.hubs, d), **kw)
        self.loss = torch.empty(plan.e, **kw); self.active = torch.empty_like(self.loss)
        self.coef=torch.empty_like(self.loss)
        self.edge_grad = torch.empty((plan.e, d), **kw)
        self.radius_grad = torch.empty_like(model.go_rad.weight)
        self.mu=torch.empty((9,d),**kw);self.var=torch.empty_like(self.mu);self.inv=torch.empty_like(self.mu)

    @torch.no_grad()
    def _initial_stats(self):
        p=self.plan;b=p.base;d=b.d;nb=self.moment_parts.shape[1]
        _shared_moments[(nb,tr.cdiv(d,128))](self.model.go_embed.weight,p.counts,self.mu,self.moment_parts,d,b.v,nb,False)
        self.mu.copy_(self.moment_parts.sum(1)/b.ns[:,None])
        _shared_moments[(nb,tr.cdiv(d,128))](self.model.go_embed.weight,p.counts,self.mu,self.moment_parts,d,b.v,nb,True)
        self.var.copy_(self.moment_parts.sum(1)/b.ns[:,None]);self.inv.copy_(torch.rsqrt(self.var+b.bn.eps))

    @torch.no_grad()
    def _step(self, x, y):
        model = self.model; p = self.plan; b = p.base
        if not model.training or not model.go_norm.training:
            raise RuntimeError('StreamingStep is a training step, not an evaluation forward')
        if torch.is_autocast_enabled():
            raise RuntimeError('StreamingStep preserves the FP32 recipe; disable autocast')
        model.zero_grad(set_to_none=True)
        w = model.go_embed.weight; r = model.rel_embed.weight; rad = model.go_rad.weight
        gamma = model.go_norm.weight; beta = model.go_norm.bias
        first=model.net[0];second=model.net[1].fn
        hidden,cache0=_block_forward(first,x);residual,cache1=_block_forward(second,hidden)
        bce,gh,pred_grad,gr,grad_rad=_classifier_gradients(hidden+residual,w[:model.nb_gos],
            r[model.nb_rels:model.nb_rels+1],rad[:model.nb_gos],y)
        _block_backward(first,cache0,gh+_block_backward(second,cache1,gh),False)
        r.grad=torch.zeros_like(r);r.grad[model.nb_rels:model.nb_rels+1].copy_(gr)
        rad.grad=torch.zeros_like(rad);rad.grad[:model.nb_gos].copy_(grad_rad)
        with torch.no_grad():
            d = b.d; bd = tr.next_power_of_2(d); nb = self.moment_parts.shape[1]
            self._initial_stats();mu=self.mu;var=self.var;inv=self.inv
            _bn_running[(tr.cdiv(d,256),)](mu,var,b.ns,b.bn.running_mean,b.bn.running_var,b.bn.num_batches_tracked,d,b.bn.momentum)
            _edge_coefficients[(b.e,)](w, r, rad, mu, inv, gamma, beta, *p.endpoints,
                *b.edge[2:], self.loss, self.active, self.coef, d, b.margin, bd, num_warps=8)
            _constraint_summaries[(p.groups,tr.cdiv(d,128))](w,r,mu,inv,gamma,beta,*p.endpoints,
                b.edge[2],b.edge[3],self.coef,self.edge_grad,p.group_edges,p.starts,p.ends,self.summaries,d,num_warps=4)
            _summary_reduce[(18 + b.nr, tr.cdiv(d, 128))](self.summaries, p.summary_slots, p.summary_ptr, self.reduced, d)
            s1 = self.reduced[:9]; s2 = self.reduced[9:18]
            k = gamma[None, :] * inv
            aa = -k * inv * s2 / b.ns[:, None]
            bb = k * (mu * inv * s2 - s1) / b.ns[:, None]
            if p.hubs:
                _hub_partials[(p.hubs,tr.cdiv(d,256))](self.edge_grad, p.edges, p.roles, p.signs, k,
                    p.hub_starts, p.hub_ends, self.hub_parts, d, 256, num_warps=4)
            route = b.radroutes
            _radius_backward[(b.v,)](self.active, route['edges'], route['signs'], route['ptr'], rad, self.radius_grad, b.v)
            _bias_parameters[(1,)](self.step_tensor,self.bias,self.lr)
            _stream_adam[(b.v,tr.cdiv(d,256))](w, w, self.edge_grad, p.edges, p.roles, p.signs,
                p.ptr, p.hub_ptr, self.hub_parts, p.counts, k, aa, bb, pred_grad, self.m, self.v,
                self.bias[:1],self.bias[1:],
                d, model.nb_gos, 256, num_warps=4)
            for param, grad in ((r, self.reduced[18:]), (rad, self.radius_grad),
                                (gamma, s2.sum(0)), (beta, s1.sum(0))):
                if param.grad is None:
                    param.grad = grad
                else:
                    param.grad.add_(grad)
            self.other.step()
            self.last_bce = bce.detach()
            self.last_el = self.loss.sum()
            total = self.last_bce + self.last_el
        return total

    def __call__(self,x,y):
        if self.graph is None:return self._step(x,y)
        if x.shape!=self.graph_x.shape or y.shape!=self.graph_y.shape:raise ValueError('Captured batch shape differs')
        self.graph_x.copy_(x);self.graph_y.copy_(y);self.graph.replay()
        self.last_bce = self.graph_bce
        self.last_el = self.graph_el
        return self.graph_loss

    def capture(self,x,y):
        capture_step(self,x,y)
        self.graph_bce = self.last_bce
        self.graph_el = self.last_el

    def state_dict(self):
        return dict(engine='deepgose_streaming_v2', step=self.step_tensor, lr=self.lr,
                    m=self.m, v=self.v, other=self.other.state_dict())

    def load_state_dict(self, state):
        if self.graph is not None:raise RuntimeError('Restore optimizer before capturing')
        if state['engine'] != 'deepgose_streaming_v2':
            raise ValueError('Incompatible streaming optimizer state')
        if state['m'].shape != self.m.shape or state['v'].shape != self.v.shape:
            raise ValueError('Optimizer state does not match ontology table shape')
        self.m.copy_(state['m']); self.v.copy_(state['v'])
        self.step_tensor.copy_(state['step']);self.lr = float(state['lr'])
        self.other.load_state_dict(state['other'])
