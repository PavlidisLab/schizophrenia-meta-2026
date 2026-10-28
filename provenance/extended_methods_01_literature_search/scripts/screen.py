"""Two-pass batch screen. Usage: PYTHONPATH=<pylib> python3 screen.py pilot|full
Key is read from the ANTHROPIC_API_KEY environment variable (credential lookup removed in this copy)."""
import json, os, sys, time, random, hashlib, subprocess, datetime, csv
import prompts as P
HERE = os.path.dirname(os.path.abspath(__file__)) + '/'
SEED = 20261009; PILOT_N = 200; PILOT_CORPUS = 50
PROJ_LIMIT = 60.0; HARD_LIMIT = 75.0; CHUNK = 5000
def log(m):
    open(HERE+'run_log.txt','a').write(f"{datetime.datetime.now().isoformat(timespec='seconds')} {m}\n"); print(m, flush=True)
PSHA = hashlib.sha256(open(HERE+'prompts.py','rb').read()).hexdigest()
def get_key():
    # Repository copy: the credential lookup of the original was removed; set ANTHROPIC_API_KEY in the environment.
    if not os.environ.get('ANTHROPIC_API_KEY'): sys.exit('ERROR: set ANTHROPIC_API_KEY')
pool = [json.loads(l) for l in open(HERE+'pool.jsonl')]
byk = {r['key']: r for r in pool}
cid = {r['key']: f"r{i}" for i, r in enumerate(pool)}; kofc = {v:k for k,v in cid.items()}
def load(p):
    f = HERE+p
    return {json.loads(l)['key']: json.loads(l) for l in open(f)} if os.path.exists(f) else {}
RES = {1: load('results_pass1.jsonl'), 2: load('results_pass2.jsonl')}
def save_res(n, rec):
    RES[n][rec['key']] = rec
    open(HERE+f'results_pass{n}.jsonl','a').write(json.dumps(rec)+'\n')
def cost(model, u, batch=True):
    pi, po = P.PRICE[model]; f = 0.5 if batch else 1
    inp = u['input']*pi + u['cache_read']*pi*0.1 + u['cache_creation']*pi*1.25   # standard cache multipliers (assumption; not given in task)
    return (inp + u['output']*po) * f / 1e6
def usage_of(u):
    return dict(input=u.input_tokens or 0, output=u.output_tokens or 0, cache_read=getattr(u,'cache_read_input_tokens',0) or 0, cache_creation=getattr(u,'cache_creation_input_tokens',0) or 0)
def spent():
    f = HERE+'usage.tsv'
    if not os.path.exists(f): return 0.0
    return sum(float(r['cost_usd']) for r in csv.DictReader(open(f),delimiter='\t'))
def req(n, key):
    r = byk[key]; st = 'pass%d'%n
    sysp = P.PASS1_SYSTEM if n==1 else P.PASS2_SYSTEM
    user = (P.PASS1_USER if n==1 else P.PASS2_USER).format(title=r['title'] or '(no title)', abstract=r['abstract'] or '(no abstract available)')
    return dict(custom_id=cid[key], params=dict(model=P.MODELS[st], max_tokens=P.MAX_TOKENS[st],
        system=[{"type":"text","text":sysp,"cache_control":{"type":"ephemeral"}}],
        messages=[{"role":"user","content":user}],
        output_config={"effort":P.EFFORT[st],"format":{"type":"json_schema","schema":P.PASS1_SCHEMA if n==1 else P.PASS2_SCHEMA}}))
def run_batch(client, n, keys, stage):
    """Submit keys in chunks; resubmit only non-succeeded; returns after all done (max 3 rounds)."""
    todo = [k for k in keys if k not in RES[n]]
    for rnd in range(3):
        if not todo: break
        failed = []
        for i in range(0, len(todo), CHUNK):
            ch = todo[i:i+CHUNK]
            if spent() > HARD_LIMIT: sys.exit(f'ABORT: spent {spent():.2f} > {HARD_LIMIT}')
            b = client.messages.batches.create(requests=[req(n,k) for k in ch])
            log(f"pass{n} {stage} round{rnd} submitted batch {b.id} n={len(ch)}")
            while True:
                b = client.messages.batches.retrieve(b.id)
                if b.processing_status == 'ended': break
                time.sleep(30)
            tot = dict(input=0,output=0,cache_read=0,cache_creation=0); cnt = {}
            for res in client.messages.batches.results(b.id):
                k = kofc[res.custom_id]; cnt[res.result.type] = cnt.get(res.result.type,0)+1
                if res.result.type != 'succeeded': failed.append(k); continue
                m = res.result.message; u = usage_of(m.usage)
                for kk in tot: tot[kk] += u[kk]
                txt = ''.join(c.text for c in m.content if c.type=='text')
                try:
                    out = json.loads(txt)
                    assert out['decision'] in ('include','exclude','unsure')
                except Exception:
                    failed.append(k); cnt['unparsed'] = cnt.get('unparsed',0)+1; continue
                save_res(n, dict(key=k, batch_id=b.id, stop_reason=m.stop_reason, usage=u, stage=stage, **out))
            model = P.MODELS['pass%d'%n]; c = cost(model, tot)
            new = not os.path.exists(HERE+'usage.tsv')
            with open(HERE+'usage.tsv','a') as f:
                w = csv.writer(f, delimiter='\t')
                if new: w.writerow(['batch_id','stage','pass','model','round','requests','results','input_tokens','output_tokens','cache_read_tokens','cache_creation_tokens','cost_usd','prompt_sha256'])
                w.writerow([b.id,stage,n,model,rnd,len(ch),json.dumps(cnt),tot['input'],tot['output'],tot['cache_read'],tot['cache_creation'],f"{c:.6f}",PSHA])
            log(f"pass{n} {stage} batch {b.id} results {cnt} cost ${c:.4f}; cumulative ${spent():.4f}")
            if spent() > HARD_LIMIT: sys.exit(f'ABORT: spent {spent():.2f} > {HARD_LIMIT}')
        todo = failed
    if todo: log(f"pass{n}: {len(todo)} records still failed after 3 rounds")
def pilot_keys():
    rng = random.Random(SEED)
    c = sorted(k for k,r in byk.items() if r['in_corpus_536']=='yes'); o = sorted(k for k,r in byk.items() if r['in_corpus_536']=='no')
    return rng.sample(c, PILOT_CORPUS) + rng.sample(o, PILOT_N-PILOT_CORPUS)
def p2keys(keys): return [k for k in keys if k in RES[1] and RES[1][k]['decision'] in ('include','unsure')]
def tokens_per_req(n, keys):
    rs = [RES[n][k]['usage'] for k in keys if k in RES[n]]
    return {x: sum(u[x] for u in rs)/len(rs) for x in ('input','output','cache_read','cache_creation')} if rs else None
def project(pk):
    import anthropic
    pool_c = [k for k,r in byk.items() if r['in_corpus_536']=='yes']; pool_o = [k for k in byk if k not in set(pool_c)]
    out = {}
    tot = 0
    for n in (1,2):
        pass
    t1 = tokens_per_req(1, pk); c1 = cost(P.MODELS['pass1'], t1)   # cost per request
    rate = {}
    for lab, ks in (('corpus',[k for k in pk if byk[k]['in_corpus_536']=='yes']),('other',[k for k in pk if byk[k]['in_corpus_536']=='no'])):
        rate[lab] = len(p2keys(ks))/len(ks)
    n2 = rate['corpus']*len(pool_c) + rate['other']*len(pool_o)
    t2 = tokens_per_req(2, p2keys(pk)); c2 = cost(P.MODELS['pass2'], t2) if t2 else 0
    p1 = c1*len(byk); p2 = c2*n2
    return dict(pool=len(byk), pass1_cost_per_request=c1, pass1_tokens_per_request=t1, pass2_forward_rate=rate, pass2_expected_requests=n2,
                pass2_tokens_per_request=t2, pass2_cost_per_request=c2, projected_pass1=p1, projected_pass2=p2, projected_total=p1+p2)
if __name__ == '__main__':
    mode = sys.argv[1]; get_key()
    import anthropic; client = anthropic.Anthropic(default_headers={"anthropic-workspace-id": "wrkspc_01EoYR4Sn9RMKy8wrTLewG3U"})
    pk = pilot_keys(); json.dump(pk, open(HERE+'pilot_keys.json','w'))
    if mode == 'pilot':
        run_batch(client, 1, pk, 'pilot'); run_batch(client, 2, p2keys(pk), 'pilot')
        pr = project(pk); pr['pilot_spent_usd'] = spent(); pr['prompt_sha256'] = PSHA
        json.dump(pr, open(HERE+'pilot_projection.json','w'), indent=1); log('pilot projection '+json.dumps(pr))
        import pilot_report; pilot_report.write()
    elif mode == 'full':
        pr = json.load(open(HERE+'pilot_projection.json'))
        if pr['projected_total'] > PROJ_LIMIT: sys.exit(f"STOP: projected ${pr['projected_total']:.2f} > ${PROJ_LIMIT}")
        allk = [r['key'] for r in pool]
        run_batch(client, 1, allk, 'full'); run_batch(client, 2, p2keys(allk), 'full')
        import report; report.write()
