"""Recalculate added results independently from stored forecasts and raw targets."""
from pathlib import Path
from collections import defaultdict
from statistics import mean,median
import json,csv,math,sys
root=Path(__file__).resolve().parents[1]
def read(name):return json.loads((root/name).read_text())
with (root/'inputs/original/Blood_data_Tema.csv').open(encoding='utf-8-sig',newline='') as f:raw=list(csv.DictReader(f))
def verify(prefix,fields,expected):
 errors=read(f'results/{prefix}-errors.json');pred=read(f'results/{prefix}-predictions.json');summary=read(f'results/{prefix}-summary.json')
 assert len(summary)==expected and len(errors)==17*expected and len(pred)==170*expected
 groups=defaultdict(list)
 for r in pred:groups[tuple(r[k] for k in ['origin',*fields])].append(r)
 delta=0
 for e in errors:
  rows=groups[tuple(e[k] for k in ['origin',*fields])]
  assert [r['target'] for r in rows]==list(range(e['origin']+1,94))
  for r in rows:assert r['truth']==float(raw[r['target']-1]['QTY_DEMANDED'])
  d=[abs(r['prediction']-r['truth']) for r in rows]
  metrics={'MAPE':100*mean(v/r['truth'] for v,r in zip(d,rows)),'MAE':mean(d),'RMSE':math.sqrt(mean(v*v for v in d))}
  delta=max(delta,*[abs(metrics[k]-e[k]) for k in metrics])
 for s in summary:
  z=[r for r in errors if all(r[k]==s[k] for k in fields)];assert len(z)==17
  values={'median_MAPE':median(r['MAPE'] for r in z),'mean_MAPE':mean(r['MAPE'] for r in z),'median_MAE':median(r['MAE'] for r in z),'median_RMSE':median(r['RMSE'] for r in z)}
  delta=max(delta,*[abs(values[k]-s[k]) for k in values])
 assert delta<1e-9
 index={tuple(r[k] for k in ['origin',*fields]):r for r in errors}
 for r in errors:
  if r['selection']=='V':
   o=dict(r,selection='O');assert index[tuple(o[k] for k in ['origin',*fields])]['MAPE']<=r['MAPE']+1e-9
 return {'metric_rows':len(errors),'forecast_rows':len(pred),'summary_rows':len(summary),'max_metric_or_summary_delta':delta}
result={'window':verify('window',['method','windows','spacing','selection'],24)}
selections=read('results/window-selections.json');candidate_groups=defaultdict(set)
for r in selections:
 n=r['origin'];h=93-n;j=r['windows'];s=r['spacing']
 ends=[n-h-s*i for i in reversed(range(j))]
 assert list(map(int,r['inner_ends'].split(';')))==ends and min(ends)>=24 and max(ends)+h<=n
 ids=r['candidate_ids'].split(';');assert len(ids)==r['candidates'] and r['selected'] in ids and r['oracle'] in ids
 candidate_groups[(n,r['method'])].add(tuple(ids))
 if r['method']=='KNN':
  min_end=n-h-18
  expected=[f'{lag}:{k}' for lag in range(1,7) for k in range(1,min(52,min_end-lag-h+1)+1)]
  assert ids==expected
 else:assert len(ids)==9
assert all(len(v)==1 for v in candidate_groups.values())
assert len(read('results/window-checks.json'))==34
assert max(r['native_primary_delta'] for r in read('results/window-checks.json'))<1e-9
assert all(r['same_selection'] and r['max_forecast_delta']==0 for r in read('results/window-invariance.json'))
old=read('results/factorial-errors.json'); win=read('results/window-errors.json')
matched=[r for r in win if r['windows']==3 and r['spacing']==6 and r['selection']=='V']
primary_delta=max(abs(r['MAPE']-next(v for v in old if v['origin']==r['origin'] and v['method']==r['method'] and v['preprocessing']=='P' and v['selection']=='V' and v['rule']=='matched' and v['seed']==1)['MAPE']) for r in matched)
result['window'].update(common_candidates=True,boundaries_pass=True,native_primary_reproduction=True,common_grid_primary_max_MAPE_change=primary_delta)
if '--windows-only' not in sys.argv:
 result['grnn']=verify('grnn',['preprocessing','strategy','selection'],9)
 for r in read('results/grnn-selections.json'):
  n=r['origin'];h=93-n
  assert list(map(int,r['inner_ends'].split(';')))==[n-h-12,n-h-6,n-h]
  assert r['selected'] in range(1,5) and r['oracle'] in range(1,5)
 assert all(r['same_selection'] and r['max_forecast_delta']==0 for r in read('results/grnn-invariance.json'))
 hist=read('reference/grnn-published-errors.json');current=read('results/grnn-original-runtime.json')
 d=max(abs(a[k]-b[k]) for a,b in zip(hist,current) for k in ('RMSE','MAE','MAPE'))
 assert len(hist)==len(current)==17 and d<5.1e-6
 assert max(r['point_delta'] for r in read('results/grnn-control.json'))<1e-10
 result['grnn'].update(historical_rounded_max_difference=d,original_control_flow_pass=True)
 flow=read('results/information-flow-checks.json')
 assert flow['lstm']['unchanged_training'] and flow['lstm']['changed_input_steps']==list(range(2,19)) and flow['lstm']['changed_inversion_steps']==list(range(2,19))
 result['information_flow']='PASS'
 result['status']='PASS'
 (root/'results/extension-verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
