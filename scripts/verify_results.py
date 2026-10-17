"""Independent metric aggregation; no forecasting libraries or R outputs as inputs to metrics."""
from pathlib import Path
from collections import defaultdict
from statistics import mean, median
import json, math, csv

ROOT = Path(__file__).resolve().parents[1]
def read(name):
    return json.loads((ROOT / name).read_text(encoding='utf-8'))
def key(row):
    return tuple(row[k] for k in ['origin','method','preprocessing','selection','rule','seed'])
preds = read('results/factorial-predictions.json')
errors = read('results/factorial-errors.json')
summary = read('results/factorial-summary.json')
groups = defaultdict(list)
for row in preds:
    groups[key(row)].append(row)
rawpath = ROOT/'inputs/original/Blood_data_Tema.csv'
with rawpath.open(encoding='utf-8-sig', newline='') as f:
    raw = list(csv.DictReader(f))
delta = 0.0
for row in errors:
    z = groups[key(row)]
    assert len(z) == row['horizon'] == 93-row['origin']
    assert [v['target'] for v in z] == list(range(row['origin']+1,94))
    for v in z:
        assert v['truth'] == float(raw[v['target']-1]['QTY_DEMANDED'])
    diffs = [abs(v['prediction']-v['truth']) for v in z]
    calc = {'MAPE':100*mean(d/v['truth'] for d,v in zip(diffs,z)),
            'MAE':mean(diffs),'RMSE':math.sqrt(mean(d*d for d in diffs))}
    delta = max(delta, *(abs(calc[k]-row[k]) for k in calc))
assert delta < 1e-9, delta
summary_delta = 0.0
for row in summary:
    z = [e for e in errors if all(e[k]==row[k] for k in ['method','preprocessing','selection','rule','seed'])]
    assert len(z)==17
    calc = dict(median_MAPE=median(v['MAPE'] for v in z),mean_MAPE=mean(v['MAPE'] for v in z),
                median_MAE=median(v['MAE'] for v in z),median_RMSE=median(v['RMSE'] for v in z))
    summary_delta = max(summary_delta, *(abs(calc[k]-row[k]) for k in calc))
assert summary_delta < 1e-9
indexed = {key(e):e for e in errors}
for e in errors:
    if e['selection']=='V':
        o=dict(e,selection='O')
        assert indexed[key(o)]['MAPE'] <= e['MAPE']+1e-9
for row in read('results/factorial-selections.json'):
    ends = list(map(int,row['inner_ends'].split(';')))
    assert max(ends)+row['inner_H']<=row['origin']
    if row['method']=='KNN':
        for label in ['selected','oracle']:
            lag,k=map(int,row[label].split(':'))
            assert 1<=lag<=6 and 1<=k<=52
            assert k<=min(ends)-lag-row['inner_H']+1
            assert k<=row['origin']-lag-(93-row['origin'])+1
    else:
        for label in ['selected','oracle']:
            p,r=map(int,row[label].split(':')); assert p in (1,2,3) and r in (25,30,35)
pub = read('reference/nnar-published-errors.json')
nn = sorted([e for e in errors if e['method']=='NNAR' and e['preprocessing']=='G' and e['selection']=='O' and e['rule']=='matched' and e['seed']==1], key=lambda e:e['origin'])
nn_diff = [dict(origin=e['origin'], published_MAPE=p['MAPE'],current_MAPE=e['MAPE'],difference=e['MAPE']-p['MAPE']) for e,p in zip(nn,pub)]
assert len(nn_diff)==17
pilot = ROOT/'reference'
baseline = json.loads((pilot/'baseline-errors.json').read_text(encoding='utf-8'))
legacy = json.loads((pilot/'origin-errors.json').read_text(encoding='utf-8'))
cross_checks = {}
for method in ['KNN','ARIMA']:
    current = sorted([e for e in errors if e['method']==method and e['preprocessing']=='G' and e['seed']==1 and e['rule']==('matched' if method=='KNN' else 'none') and e['selection']==('O' if method=='KNN' else 'none')],key=lambda e:e['origin'])
    original = baseline[method]
    assert len(current)==len(original)==17
    diff=max(abs(a[k]-b[k]) for a,b in zip(current,original) for k in ['MAPE','MAE','RMSE'])
    assert diff<1e-9,(method,diff)
    cross_checks[method+'_original_reproduction_max_difference']=diff
current=sorted([e for e in errors if e['method']=='KNN' and e['preprocessing']=='P' and e['selection']=='V' and e['rule']=='matched'],key=lambda e:e['origin'])
prior=sorted([e for e in legacy if e['method']=='C_KNN' and e['rule']=='matched'],key=lambda e:e['origin'])
assert len(current)==len(prior)==17
diff=max(abs(a[k]-b[k]) for a,b in zip(current,prior) for k in ['MAPE','MAE','RMSE'])
assert diff<1e-9,diff
cross_checks['historical_KNN_previous_pilot_max_difference']=diff
output=dict(status='PASS',metric_rows=len(errors),forecast_rows=len(preds),summary_rows=len(summary),
    unique_outer_months=len({v['target'] for v in preds}),max_metric_delta=delta,max_summary_delta=summary_delta,
    identical_raw_targets=True,oracle_bounds_pass=True,selection_feasibility_pass=True,
    nnar_published_comparison=nn_diff,
    nnar_max_published_MAPE_difference=max(abs(v['difference']) for v in nn_diff),
    nnar_published_median_MAPE=median(v['MAPE'] for v in pub),
    cross_checks=cross_checks,
    scope='Recalculates metrics from stored forecasts; does not refit models.')
(ROOT/'results/python-verification.json').write_text(json.dumps(output,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in output.items() if k!='nnar_published_comparison'},indent=2))
